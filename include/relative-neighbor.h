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

#ifndef NGL_RELATIVE_NEIGHBOR_H
#define NGL_RELATIVE_NEIGHBOR_H

#include "geometric-test.h"
#include "neighborhood-builder.h"
#include "vectorspace.h"

namespace ngl {

template<typename Point, typename Scalar>
class RelativeNeighborEdge {
 public:
  Point p;
  Point q;
  Scalar square_length;

  explicit RelativeNeighborEdge(const VectorSpace<Point, Scalar>& space) :
      space_(space),
      p(NULL), q(NULL) {
    space_.allocate(&p);
    space_.allocate(&q);
  }

  virtual ~RelativeNeighborEdge() {
    space_.deallocate(p);
    space_.deallocate(q);
  }

  void reset(const Point& pin, const Point& qin) {
    assert(p);
    space_.set(p, pin);
    space_.set(q, qin);
    square_length = space_.distanceSqr(p, q);
  }

  const VectorSpace<Point, Scalar>& space_;
};


template<typename Point, typename Scalar>
class RelativeNeighbor :
    public GeometricTest<Point, Scalar, RelativeNeighborEdge<Point, Scalar> > {
 public:
  explicit RelativeNeighbor(const VectorSpace<Point, Scalar>& space):
      GeometricTest<Point, Scalar, RelativeNeighborEdge<Point, Scalar> >(),
      space_(space), edge_(space) {
  }
  virtual ~RelativeNeighbor() {}

  virtual Scalar shadowing(const RelativeNeighborEdge<Point, Scalar>& edge,
                           const Point& r) {
    Scalar dp = space_.distanceSqr(r, edge.p);
    Scalar dq = space_.distanceSqr(r, edge.q);
    return fmax(dq - dp, edge.square_length - dp);
  }

  virtual void setActiveEdge(const Point& p, const Point& q) {
    edge_.reset(p, q);
  }
  virtual RelativeNeighborEdge<Point, Scalar>& getActiveEdge() {
    return edge_;
  }

 private:
  const VectorSpace<Point, Scalar>& space_;
  RelativeNeighborEdge<Point, Scalar> edge_;
};



template<typename Point, typename Scalar>
class RelativeNeighborGraphBuilder :
    public NeighborhoodBuilder<Point, Scalar,
        RelativeNeighbor<Point, Scalar> > {
 public:
  explicit RelativeNeighborGraphBuilder(const VectorSpace<Point, Scalar>& space)
      : NeighborhoodBuilder<Point, Scalar,
          RelativeNeighbor<Point, Scalar> >(space) {
  }
};

};  // namespace ngl

#endif  // NGL_RELATIVE_NEIGBOR_H
