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

#ifndef NGL_NEIGHBORHOOD_BUILDER_H
#define NGL_NEIGHBORHOOD_BUILDER_H

#include <assert.h>
#include <vector>

#include "empty-region.h"
#include "neighbor-graph.h"
#include "vectorspace.h"

using std::vector;

namespace ngl {

template<typename Point, typename Scalar>
class NeighborhoodBuilder {
 public:
  NeighborhoodBuilder(const VectorSpace<Point, Scalar>& space) :
      space_(space), relaxed_(false) {
  }
  
  virtual ~NeighborhoodBuilder() {}
  
  virtual EmptyRegion<Point, Scalar>& getEmptyRegion() = 0;

  void addPoints(const vector<Point>& points) {
    points_ = points;
  }

  void computeNeighborGraph(NeighborGraph *neighbor_graph) {
    valid_ = new bool[points_.size()];
    for (int i = 0; i < points_.size(); ++i) {
      valid_[i] = true;
    }
    assert(neighbor_graph);
    neighbor_graph->clear();
    resetVisitedNodes();
    vector<int> neighbors;
    for (int i = 0; i < points_.size(); ++i) {
      invalidate(i);
      getNeighbors(points_[i], &neighbors);
      makeValid(i);
      for (int k = 0; k < neighbors.size(); ++k) {
        neighbor_graph->addEdge(i, neighbors[k]);
      }
    }
    postProcessEdges(neighbor_graph);
    delete valid_;
  }

  void getNeighbors(const Point &p, vector<int>* neighbors) {
    assert(neighbors);
    neighbors->clear();
    EmptyRegion<Point, Scalar>& empty_region = getEmptyRegion();
    int num_points = points_.size();
    int* index_pool = new int[num_points];
    int pool_ptr = num_points - 1;

    int k = 0;
    for (int j = 0; j < num_points; ++j) {
      int idx = j;
      if (!isValid(idx)) {
        continue;
      }
      index_pool[k++] = idx;
      pool_ptr = k;
    }
    int numIterations = 0;

    while (pool_ptr > 0) {
      int nearest = -1;
      Scalar min_distance = 1.0e20;
      for (int j = 0; j < pool_ptr; ++j) {
        int idx = index_pool[j];
        assert(isValid(idx));
        Scalar d = space_.distanceSqr(p, points_[idx]);
        incrementVisitedNodes();
        if (d < min_distance) {
          min_distance = d;
          nearest = idx;
        }
      }
      empty_region.set(p, points_[nearest]);
      neighbors->push_back(nearest);
      // remove points pruned by hyperplane between (i , nearest)
      int n = pool_ptr - 1;
      for (int j = n ; j >= 0 ; --j) {
        int idx = index_pool[j];
        incrementVisitedNodes();
        if (empty_region.shadows(points_[idx])
            || idx == nearest) {
          // in place deletion
          deleteInPlace(index_pool, j, pool_ptr);
        }
      }
      numIterations++;
    }

    delete index_pool;
  }

  void invalidate(int index) {
    valid_[index] = false;
  }
  void makeValid(int index) {
    valid_[index] = true;
  }
 private:
  bool isValid(int index) {
    return valid_[index];
  }
  void resetVisitedNodes() {}
  void incrementVisitedNodes() {}
  inline void deleteInPlace(int *p, int a, int &b) {
    p[a] = p[b-1];
    b--;
  }

  virtual void postProcessEdges(NeighborGraph *neighbor_graph) {
    if (relaxed_) return;
    EmptyRegion<Point, Scalar>& empty_region = getEmptyRegion();
    int numEdges = neighbor_graph->size();
    for (int i = numEdges - 1; i >= 0; i--) {
      int i1, i2;
      neighbor_graph->getEdge(i, &i1, &i2);
      Point &p = points_[i1];
      Point &q = points_[i2];
      bool isEdge = true;
      for (unsigned int j = 0; j < points_.size() && isEdge; j++) {
        if (j == i1 || j == i2) continue;
        empty_region.set(p, points_[j]);
        if (empty_region.shadows(q)) {
          isEdge = false;
        }
      }
      if (!isEdge) {
        neighbor_graph->removeEdge(i);
      }
    }
  }

  const VectorSpace<Point, Scalar>& space_;
  bool relaxed_;
  vector<Point> points_;
  bool* valid_;
};

};  // namespace ngl

#endif  // NGL_NEIGHBORHOOD_BUILDER_H
