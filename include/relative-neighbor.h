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

#include "empty-region.h"
#include "neighborhood-builder.h"
#include "vectorspace.h"

namespace ngl {

template<typename Point, typename Scalar>
class RelativeNeighbor :
    public EmptyRegion<Point, Scalar> {
 public:
  explicit RelativeNeighbor(const VectorSpace<Point, Scalar>& space):
      EmptyRegion<Point, Scalar>(),
      space_(space) {
    space_.allocate(&p);
    space_.allocate(&q);
  }
  virtual ~RelativeNeighbor() {
    space_.deallocate(p);
    space_.deallocate(q);
  }

  virtual void set(const Point& pin, const Point& qin) {
    space_.set(p, pin);
    space_.set(q, qin);
    square_length = space_.distanceSqr(p, q);
  }

  virtual Scalar shadowing(const Point& r) {
    Scalar dp = space_.distanceSqr(r, p);
    Scalar dq = space_.distanceSqr(r, q);
    return fmax(dq - dp, square_length - dp);
  }

 private:
  const VectorSpace<Point, Scalar>& space_;
  Point p;
  Point q;
  Scalar square_length;
};



template<typename Point, typename Scalar>
class RelativeNeighborGraphBuilder :
    public NeighborhoodBuilder<Point, Scalar> {
 public:
  explicit RelativeNeighborGraphBuilder(const VectorSpace<Point, Scalar>& space)
      : NeighborhoodBuilder<Point, Scalar>(space), region_(space) {}
 private:
  RelativeNeighbor<Point, Scalar>& getEmptyRegion() {
    return region_;
  }
  RelativeNeighbor<Point, Scalar> region_;
};

};  // namespace ngl

#endif  // NGL_RELATIVE_NEIGBOR_H
