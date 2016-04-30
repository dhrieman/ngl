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

#ifndef NGL_BSKELETON_H
#define NGL_BSKELETON_H

#include "empty-region.h"
#include "neighborhood-builder.h"
#include "vectorspace.h"

namespace ngl {

template<typename Point, typename Scalar>
class BSkeleton :
    public EmptyRegion<Point, Scalar> {
 public:
  explicit BSkeleton(const VectorSpace<Point, Scalar>& space, Scalar beta):
      EmptyRegion<Point, Scalar>(),
      space_(space), beta_(beta) {
    space_.allocate(&p);
    space_.allocate(&q);
    space_.allocate(&pq);
    space_.allocate(&planeCrossing);
    space_.allocate(&centroid);
    space_.allocate(&rq);
  }
  virtual ~BSkeleton() {
    space_.deallocate(p);
    space_.deallocate(q);
    space_.deallocate(pq);
    space_.deallocate(planeCrossing);
    space_.deallocate(centroid);
    space_.deallocate(rq);
  }

  void set(const Point& pin, const Point& qin) {
    assert(p);
    space_.set(p, pin);
    space_.set(q, qin);
    square_length = space_.distanceSqr(p, q);
    space_.subtract(q, p, pq);
    if (beta_ > 1) {
      space_.interpolate(p, q, 1.0 / beta_, planeCrossing);
      Scalar delta = 0.5 * (1.0 + 1.0 / (beta_ - 1.0));
      space_.interpolate(p, q, 1.0 - delta, centroid);
      distanceCentroidSqr = space_.distanceSqr(q, centroid);
    }
  }

  virtual Scalar shadowing(const Point& r) {
    if (beta_ <= 1.0) {
      space_.subtract(r, q, rq);
      Scalar dot = space_.dot(pq, rq);
      Scalar mag2sqr = space_.dot(rq, rq);
      Scalar test = (square_length > 0 && mag2sqr > 0)
          ? fabs(dot) * dot / (square_length * mag2sqr)
          : 1;
      Scalar maxdot = 1.0 - beta_ * beta_;
      return maxdot - test;
    } else {
      space_.subtract(r, planeCrossing, rq);
      Scalar dot = space_.dot(pq, rq);
      if (beta_ > 1) {
        Scalar drc2 = space_.distanceSqr(r, centroid);
        if (drc2 < distanceCentroidSqr) {
          return distanceCentroidSqr - drc2;
        }
      }
      return -dot;
    }
  }

 private:
  const VectorSpace<Point, Scalar>& space_;
  Scalar beta_;

  Point p;
  Point q;
  // Vector pq = q - p
  Point pq;
  // Point where umbra plane crosses edge.
  Point planeCrossing;
  Point centroid;
  Scalar square_length;
  // Distance squared centroid - q
  Scalar distanceCentroidSqr;
  // Vector rq = r - q, to be set on shadow queries
  Point rq;
};



template<typename Point, typename Scalar>
class BSkeletonBuilder :
    public NeighborhoodBuilder<Point, Scalar> {
 public:
  BSkeletonBuilder(const VectorSpace<Point, Scalar>& space, Scalar beta)
      : NeighborhoodBuilder<Point, Scalar>(space), region_(space, beta) {}
 private:
  BSkeleton<Point, Scalar>& getEmptyRegion() {
    return region_;
  }
  BSkeleton<Point, Scalar> region_;
};

};  // namespace ngl

#endif  // NGL_BSKELETON_H
