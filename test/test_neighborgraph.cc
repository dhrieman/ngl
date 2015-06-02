/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright 2015  David H. Rieman (david.h.rieman@gmail.com)
 * All rights reserved.
 *
 * THE BSD LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/
#include <assert.h>
#include <stdio.h>

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "neighbor-graph-impl.h"
#include "ngl.h"

using std::ifstream;
using std::runtime_error;
using std::string;
using std::vector;

void readPoints(const char *src, int dims, vector<double>* points) {
  ifstream inp(src);
  if (!inp.is_open()) {
    fprintf(stderr, "Cannot open file %s\n", src);
    throw runtime_error("Cannot open file");
  }
  string line;
  while (!inp.eof()) {
    getline(inp, line);
    if (line.length() == 0) continue;
    const char *sline = line.c_str();
    char *tok = strtok((char*) sline, " \t\r\n");
    if (tok && tok[0] == '#') continue;
    int ndims = 0;
    vector<double> readPts;
    while (tok && ndims < dims) {
      readPts.push_back(atof(tok));
      points->push_back(atof(tok));
      tok = strtok(NULL, " \t\r\n");
      ndims++;
    }
    if (ndims < dims) {
      fprintf(stderr, "Error reading input file: num.dims=%d, expecting %d\n",
          ndims, dims);
      throw runtime_error("Error reading input file");
    }
  }
  inp.close();
}

void readEdges(const char * src, vector<int>* edges, int *expectedEdges) {
  ifstream inp(src);
  if (!inp.is_open()) {
    fprintf(stderr, "Cannot open file %s\n", src);
    throw runtime_error("Cannot open file");
  }
  string line;
  while (!inp.eof()) {
    getline(inp, line);
    if (line.length() == 0) continue;
    const char *sline = line.c_str();
    char *tok = strtok((char*) sline, " \t\r\n");
    if (tok && tok[0] == '#') continue;
    assert(tok);
    int i1 = atoi(tok);
    tok = strtok(NULL, " \t\r\n");
    int i2 = atoi(tok);
    edges[i1].push_back(i2);
    edges[i2].push_back(i1);
    (*expectedEdges)++;
  }
  inp.close();
}

bool edgeExists(int i1, int i2, vector<int>* edges) {
  if (find(edges[i1].begin(), edges[i1].end(), i2) != edges[i1].end()) {
    return true;
  }
  if (find(edges[i2].begin(), edges[i2].end(), i1) != edges[i2].end()) {
    return true;
  }
  return false;
}

template<typename NeighborhoodBuilder>
bool testNeighborGraphBuilder(const char *src,
                              const char *expected,
                              string testName,
                              int dims, double param = 0.0) {
  vector<double> data;
  readPoints(src, dims, &data);
  int num_points = data.size() / dims;
  vector<double*> points;
  for (int i = 0; i < num_points; ++i) {
    double* p = new double[dims];
    for (int k = 0; k < dims; ++k) {
      p[k] = data[i * dims + k];
    }
    points.push_back(p);
  }
  fprintf(stdout, "Test method %s, dims %d numPoints %d\n",
      testName.c_str(), dims, num_points);

  ngl::DoubleVectorSpace s(dims);
  NeighborhoodBuilder rngb(s);
  rngb.addPoints(points);
  rngb.setParam(param);

  ngl::NeighborGraphImpl neighborGraph;
  rngb.computeNeighborGraph(&neighborGraph);

  vector<int>* edges = new vector<int>[num_points];
  int numExpectedEdges = 0;
  readEdges(expected, edges, &numExpectedEdges);

  int numEdges = neighborGraph.size();
  if (numExpectedEdges != numEdges) {
    fprintf(stderr, "Obtained %d edges. Expecting %d\n", numEdges,
        numExpectedEdges);
    return false;
  }
  bool passed = true;
  for (int i = 0; i < numEdges; ++i) {
    int src;
    int dst;
    neighborGraph.getEdge(i, &src, &dst);
    if (!edgeExists(src, dst, edges)) {
      passed = false;
      break;
    }
  }
  delete[] edges;
  for (int i = 0; i < num_points; ++i) {
    delete points[i];
  }

  return passed;
}

struct TestCase {
  string points;
  string expected_edges;
  int dims;
};

vector<TestCase> tests;

void addTest(const char* points, const char* expected_edges, int dims) {
  TestCase test;
  test.points = string(points);
  test.expected_edges = string(expected_edges);
  test.dims = dims;
  tests.push_back(test);
}

using ngl::BSkeletonBuilderD;
using ngl::GabrielGraphBuilderD;
using ngl::RelativeNeighborGraphBuilderD;

int main(int argc, char* argv[]) {
  addTest("../test_data/p2d", "../test_data/graphs/p2d.2.b", 2);
  addTest("../test_data/p3d", "../test_data/graphs/p3d.2.b", 3);
  addTest("../test_data/p4d", "../test_data/graphs/p4d.2.b", 4);
  addTest("../test_data/p5d", "../test_data/graphs/p5d.2.b", 5);
  for (int i = 0; i < tests.size(); ++i) {
    bool passed = testNeighborGraphBuilder<RelativeNeighborGraphBuilderD>(
        tests[i].points.c_str(), tests[i].expected_edges.c_str(),
        "Relative Neighbor", tests[i].dims);
    if (passed) {
      fprintf(stderr, "[\033[32m PASSED \033[0m]\n ");
    } else {
      fprintf(stderr, "[\033[31m FAILED \033[0m]\n ");
    }
  }
  tests.clear();
  addTest("../test_data/p2d", "../test_data/graphs/p2d.1.b", 2);
  addTest("../test_data/p3d", "../test_data/graphs/p3d.1.b", 3);
  addTest("../test_data/p4d", "../test_data/graphs/p4d.1.b", 4);
  addTest("../test_data/p5d", "../test_data/graphs/p5d.1.b", 5);
  for (int i = 0; i < tests.size(); ++i) {
    bool passed = testNeighborGraphBuilder<GabrielGraphBuilderD>(
        tests[i].points.c_str(), tests[i].expected_edges.c_str(),
        "Gabriel", tests[i].dims);
    if (passed) {
      fprintf(stderr, "[\033[32m PASSED \033[0m]\n ");
    } else {
      fprintf(stderr, "[\033[31m FAILED \033[0m]\n ");
    }
  }

  tests.clear();
  addTest("../test_data/p2d", "../test_data/graphs/p2d.1.4.b", 2);
  addTest("../test_data/p3d", "../test_data/graphs/p3d.1.4.b", 3);
  addTest("../test_data/p4d", "../test_data/graphs/p4d.1.4.b", 4);
  addTest("../test_data/p5d", "../test_data/graphs/p5d.1.4.b", 5);
  for (int i = 0; i < tests.size(); ++i) {
    bool passed = testNeighborGraphBuilder<BSkeletonBuilderD>(
        tests[i].points.c_str(), tests[i].expected_edges.c_str(),
        "BSkeleton", tests[i].dims, 1.4);
    if (passed) {
      fprintf(stderr, "[\033[32m PASSED \033[0m]\n ");
    } else {
      fprintf(stderr, "[\033[31m FAILED \033[0m]\n ");
    }
  }

  tests.clear();
  addTest("../test_data/p2d", "../test_data/graphs/p2d.0.7.b", 2);
  addTest("../test_data/p3d", "../test_data/graphs/p3d.0.7.b", 3);
  addTest("../test_data/p4d", "../test_data/graphs/p4d.0.7.b", 4);
  addTest("../test_data/p5d", "../test_data/graphs/p5d.0.7.b", 5);
  for (int i = 0; i < tests.size(); ++i) {
    bool passed = testNeighborGraphBuilder<BSkeletonBuilderD>(
        tests[i].points.c_str(), tests[i].expected_edges.c_str(),
        "BSkeleton", tests[i].dims, 0.7);
    if (passed) {
      fprintf(stderr, "[\033[32m PASSED \033[0m]\n ");
    } else {
      fprintf(stderr, "[\033[31m FAILED \033[0m]\n ");
    }
  }
}
