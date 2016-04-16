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

using ngl::BSkeleton;
using ngl::BSkeletonBuilderD;
using ngl::DoubleVectorSpace;
using ngl::NeighborhoodBuilder;

void usage(char* argv0) {
  printf("\nUsage:\n%s source dims method param\n", argv0);
  printf("  source\tSource file containing input points.\n");
  printf("  dims  \tNumber of dimensions.\n");
  printf("  method\tNeighborhood method");
  printf("  (RelativeNeighbor | Gabriel | BSkeleton).\n");
  printf("  beta0  \tLower bound for beta.\n");
  printf("  beta1  \tUpper bound for beta.\n");
}


typedef BSkeleton<double*, double> BSkeletonD;
typedef double* Point;
const double TOLERANCE = 1e-3;

double estimateBeta(const DoubleVectorSpace& space,
    const Point& p, const Point& q, const Point &r, double beta0,
    double beta1) {
  assert(beta0 < beta1);

  // Terminate early if r shadowed by beta0-edge (p, q)
  BSkeletonD bskeleton0(space, beta0);
  bskeleton0.set(p, q);
  if (bskeleton0.shadows(r)) {
    return beta0;
  }

  // Terminate early if r not shadowed by beta0-edge (p, q)
  BSkeletonD bskeleton1(space, beta1);
  bskeleton1.set(p, q);
  if (!bskeleton1.shadows(r)) {
    return beta1;
  }

  double beta_lo = beta0;
  double beta_hi = beta1;
  double beta;
  do {
    beta = 0.5 * (beta_hi + beta_lo);
    BSkeletonD bskeleton(space, beta);
    bskeleton.set(p, q);
    bool shadows = bskeleton.shadows(r);
    if (shadows) {
      beta_hi = beta;
    } else {
      beta_lo = beta;
    }
  } while (beta_hi - beta_lo > TOLERANCE);
  return beta;
}

void computeProbabilisticNeighborGraph(const vector<double*>& points, int dims,
      double beta0, double beta1) {
  assert(beta0 < beta1);
  DoubleVectorSpace space(dims);
  
  // Compute upper bound graph, i.e., beta0-skeleton.
  BSkeletonBuilderD builder(space, beta0);
  builder.addPoints(points);

  ngl::NeighborGraphImpl neighborGraph;
  builder.computeNeighborGraph(&neighborGraph);

  for (int i = 0; i < neighborGraph.size(); ++i) {
    int src;
    int dst;
    neighborGraph.getEdge(i, &src, &dst);
    double min_beta = beta1;
    // Compute probability of each edge.
    for (int j = 0; j < points.size(); ++j) {
      if (j == src || j == dst) continue;
      // Find beta via binary search
      double beta = estimateBeta(space, points[src], points[j], points[dst],
          beta0, beta1);
      assert(beta >= beta0);
      assert(beta <= beta1);
      min_beta = fmin(min_beta, beta);
    }
    double probability = (min_beta - beta0) / (beta1 - beta0);
    printf("%d %d %g\n", src, dst, probability);
  }
}

int main(int argc, char* argv[]) {
  if (argc < 4) {
    usage(argv[0]);
    return 1;
  }
  int dims = atoi(argv[2]);
  string method = string(argv[3]);
  double beta0 = 1.0;
  double beta1 = 2.0;
  if (argc > 4) {
    beta0 = atof(argv[4]);
  }
  if (argc > 5) {
    beta1 = atof(argv[5]);
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

  computeProbabilisticNeighborGraph(points, dims, beta0, beta1);
}
