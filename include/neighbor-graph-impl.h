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

#ifndef NEIGHBOR_GRAPH_IMPL_H
#define NEIGHBOR_GRAPH_IMPL_H

#include <vector>

#include "neighbor-graph.h"

using std::vector;

namespace ngl {

// Implements a neighbor graph using a vector of Edges.
class NeighborGraphImpl: public NeighborGraph {
  struct Edge{
    int src;
    int dst;
    Edge(int s, int d) {
      src = s;
      dst = d;
    }
    bool operator == (const Edge& m) const {
      return m.src == src && m.dst == dst;
    }
  };

 public:
  virtual void clear() {
    edges_.clear();
  }
  virtual void addEdge(int src, int dst) {
    edges_.push_back(Edge(src, dst));
  }
  virtual void getEdge(int index, int* src, int* dst) {
    *src = edges_[index].src;
    *dst = edges_[index].dst;
  }
  virtual void removeEdge(int index) {
    assert(index >= 0);
    assert(index < edges_.size());
    edges_.erase(edges_.begin() + index);
  }
  virtual int size() {
    return edges_.size();
  }

 protected:
  vector<Edge> edges_;
};

};  // namespace ngl

#endif  // NEIGHBOR_GRAPH_IMPL_H

