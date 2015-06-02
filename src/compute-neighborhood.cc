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

#include <string>
#include <vector>

#include "io-utils.h"
#include "neighbor-graph-impl.h"
#include "ngl.h"

using std::string;
using std::vector;

void usage(char* argv0) {
  printf("\nUsage:\n%s source dims method param\n", argv0);
  printf("  source\tSource file containing input points.\n");
  printf("  dims  \tNumber of dimensions.\n");
  printf("  method\tNeighborhood method");
  printf("  (RelativeNeighbor | Gabriel | BSkeleton).\n");
  printf("  param  \tOptional parameter. Beta for beta skeleton.\n");
}

template<typename NeighborhoodBuilder>
void computeNeighborGraph(const vector<double*>& points, int dims,
      double param = 0.0) {
  ngl::DoubleVectorSpace space(dims);

  NeighborhoodBuilder builder(space);
  builder.setParam(param);
  builder.addPoints(points);

  ngl::NeighborGraphImpl neighborGraph;
  builder.computeNeighborGraph(&neighborGraph);
  for (int i = 0; i < neighborGraph.size(); ++i) {
    int src;
    int dst;
    neighborGraph.getEdge(i, &src, &dst);
    printf("%d %d\n", src, dst);
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

  if (method == "RelativeNeighbor") {
    computeNeighborGraph<ngl::RelativeNeighborGraphBuilderD>(points, dims);
  } else if (method == "Gabriel") {
    computeNeighborGraph<ngl::GabrielGraphBuilderD>(points, dims);
  } else if (method == "BSkeleton") {
    computeNeighborGraph<ngl::BSkeletonBuilderD>(points, dims, beta);
  } else {
    fprintf(stderr, "Method %s not found.\n", method.c_str());
    return 1;
  }
}
