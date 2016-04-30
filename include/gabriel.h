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

#ifndef NGL_GABRIEL_H
#define NGL_GABRIEL_H

#include "empty-region.h"
#include "neighborhood-builder.h"
#include "vectorspace.h"

namespace ngl {

template<typename Point, typename Scalar>
class Gabriel :
    public EmptyRegion<Point, Scalar> {
 public:
  explicit Gabriel(const VectorSpace<Point, Scalar>& space):
      EmptyRegion<Point, Scalar>(),
      space_(space) {
    space.allocate(&p);
    space.allocate(&q);
    space.allocate(&tmp);
    space.allocate(&midpoint);
  }
  virtual ~Gabriel() {
    space_.deallocate(p);
    space_.deallocate(q);
    space_.deallocate(tmp);
    space_.deallocate(midpoint);
  }
  
  void set(const Point& pin, const Point& qin) {
    space_.set(p, pin);
    space_.set(q, qin);
    space_.interpolate(pin, qin, 0.5, midpoint);
  }

  virtual Scalar shadowing(const Point& r) {
    space_.interpolate(p, q, 2.0, tmp);
    Scalar dp = space_.distanceSqr(r, p);
    Scalar dq = space_.distanceSqr(r, tmp);
    return dq - dp;
  }

 private:
  const VectorSpace<Point, Scalar>& space_;
  Point p;
  Point q;
  Point tmp;
  Point midpoint;
  Scalar square_length;
};


template<typename Point, typename Scalar>
class GabrielGraphBuilder :
    public NeighborhoodBuilder<Point, Scalar> {
 public:
  explicit GabrielGraphBuilder(const VectorSpace<Point, Scalar>& space) :
      NeighborhoodBuilder<Point, Scalar>(space), region_(space) {}
 private:
  Gabriel<Point, Scalar>& getEmptyRegion() {
    return region_;
  }
  Gabriel<Point, Scalar> region_;
};


};

#endif  // NGL_GABRIEL_H
