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
#include <stdlib.h>

#include <string>
#include <unordered_map>
#include <vector>

#include "io-utils.h"
#include "neighbor-graph-impl.h"
#include "ngl.h"
#include "stringprintf.h"

using std::string;
using std::unordered_map;
using std::vector;

using ngl::BSkeletonBuilderD;
using ngl::DoubleVectorSpace;
using ngl::GabrielGraphBuilderD;
using ngl::NeighborGraphImpl;
using ngl::RelativeNeighborGraphBuilderD;

void usage(char* argv0) {
  printf("\nUsage:\n%s source dims method param\n", argv0);
  printf("  source\tSource file containing input points.\n");
  printf("  dims  \tNumber of dimensions.\n");
  printf("  method\tNeighborhood method");
  printf("  (RelativeNeighbor | Gabriel | BSkeleton).\n");
  printf("  param  \tOptional parameter. Beta for beta skeleton.\n");
}

ngl::NeighborhoodBuilderD* createNeighborhoodBuilderD(
    const DoubleVectorSpace& space,
    const string& method, double beta) {
  ngl::NeighborhoodBuilderD* builder;
  if (method == "RelativeNeighbor") {
    builder = new RelativeNeighborGraphBuilderD(space);
  } else if (method == "Gabriel") {
    builder = new GabrielGraphBuilderD(space);
  } else if (method == "BSkeleton") {
    builder = new BSkeletonBuilderD(space, beta);
  } else {
    fprintf(stderr, "Method %s not found.\n", method.c_str());
    return nullptr;
  }
  return builder;
}

void estimateScale(const vector<double*>& points, int dims,
                   ngl::NeighborhoodBuilderD& builder,
                   vector<double>* scales) {
  DoubleVectorSpace space(dims);

  NeighborGraphImpl neighborGraph;
  builder.addPoints(points);
  builder.computeNeighborGraph(&neighborGraph);
  
  vector<int>* neighbors = new vector<int>[points.size()];
  scales->clear();
  for (int i = 0; i < neighborGraph.size(); ++i) {
    int src;
    int dst;
    neighborGraph.getEdge(i, &src, &dst);
    neighbors[src].push_back(dst);
    neighbors[dst].push_back(src);
  }

  // Define local scale of a point as half the median distance to its
  // neighbors
  for (int i = 0; i < points.size(); ++i) {
    vector<double> distances;
    for (int j = 0; j < neighbors[i].size(); ++j) {
      int neighbor = neighbors[i][j];
      double d = space.distance(points[i], points[neighbor]);
      distances.push_back(d);
    }
    if (distances.size() == 0) continue;
    int neighbor = distances.size() / 2;
    std::nth_element(distances.begin(), distances.begin() + neighbor,
                     distances.end());
    double scale = distances[neighbor];
    scales->push_back(scale * 0.5);
  }
}

void perturbPoints(const vector<double*>& points, int dims,
                   const vector<double>& scales,
                   vector<double*>* perturbed_points) {
  assert(perturbed_points);
  assert(perturbed_points->size() == points.size());
  assert(scales.size() == points.size());
  for (int i = 0; i < points.size(); ++i) {
    //printf("Scales[%d] = %f\n", i, scales[i]);
    for (int d = 0; d < dims; ++d) {
      double delta =
          2.0 * (double) rand() / (double) RAND_MAX - 1.0;
      (*perturbed_points)[i][d] = points[i][d] + delta * scales[i];
    }
  }
}

void computePerturbedNeighborGraph(
    const vector<double*>& points, int dims,
    const vector<double>& scales,
    const string& method,
    double beta = 0.0) {
  unordered_map<string, int> edgeMap;
  int numIterations = 100;

  vector<double*> perturbed_points;
  for (int i = 0; i < points.size(); ++i) {
    double* p = new double[dims];
    perturbed_points.push_back(p);
  }

  DoubleVectorSpace space(dims);
  for (int it = 0; it < numIterations; ++it) {
    perturbPoints(points, dims, scales, &perturbed_points);

    ngl::NeighborhoodBuilderD* builder =
      createNeighborhoodBuilderD(space, method, beta);
    builder->addPoints(perturbed_points);

    NeighborGraphImpl neighborGraph;
    builder->computeNeighborGraph(&neighborGraph);

    for (int i = 0; i < neighborGraph.size(); ++i) {
      int src;
      int dst;
      neighborGraph.getEdge(i, &src, &dst);
      string hash = StringPrintf("%d %d", src, dst);
      if (edgeMap.find(hash) == edgeMap.end()) {
        edgeMap[hash] = 0;
      }
      edgeMap[hash] = edgeMap[hash] + 1;
    }
  }

  for (int i = 0; i < perturbed_points.size(); ++i) {
    delete perturbed_points[i];
  }
  
  for (auto it = edgeMap.begin(); it != edgeMap.end(); ++it) {
    string edge = it->first;
    int occurrences = it->second;
    double weight = (double) occurrences / (double) numIterations;
    printf("%s %f\n", edge.c_str(), weight);
  }
}

int main(int argc, char* argv[]) {
  if (argc < 4) {
    usage(argv[0]);
    return 1;
  }
  int dims = atoi(argv[2]);
  string method = string(argv[3]);
  double beta = 0.0;
  if (argc > 4) {
    beta = atof(argv[4]);
  }

  vector<double> data;
  readPoints(argv[1], dims, &data);
  int num_points = data.size() / dims;
  vector<double*> points;
  for (int i = 0; i < num_points; ++i) {
    double* p = new double[dims];
    for (int k = 0; k < dims; ++k) {
      p[k] = data[i * dims + k];
    }
    points.push_back(p);
  }

  DoubleVectorSpace space(dims);
  ngl::NeighborhoodBuilderD* scaleBuilder =
      createNeighborhoodBuilderD(space, method, beta);
  vector<double> scales;
  estimateScale(points, dims, *scaleBuilder, &scales);
  delete scaleBuilder;

  computePerturbedNeighborGraph(points, dims, scales, method, beta);
}
