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

#include "geometric-test.h"
#include "neighborhood-builder.h"
#include "vectorspace.h"

namespace ngl {

template<typename Point, typename Scalar>
class BSkeletonEdge {
 public:
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

  explicit BSkeletonEdge(const VectorSpace<Point, Scalar>& space) :
      space_(space),
      p(NULL), q(NULL), pq(NULL), planeCrossing(NULL), centroid(NULL) {
    space_.allocate(&p);
    space_.allocate(&q);
    space_.allocate(&pq);
    space_.allocate(&planeCrossing);
    space_.allocate(&centroid);
  }

  virtual ~BSkeletonEdge() {
    space_.deallocate(p);
    space_.deallocate(q);
    space_.deallocate(pq);
    space_.deallocate(planeCrossing);
    space_.deallocate(centroid);
  }

  void reset(const Point& pin, const Point& qin, Scalar beta) {
    assert(p);
    space_.set(p, pin);
    space_.set(q, qin);
    square_length = space_.distanceSqr(p, q);
    space_.subtract(q, p, pq);
    if (beta > 1) {
      space_.interpolate(p, q, 1.0 / beta, planeCrossing);
      Scalar delta = 0.5 * (1.0 + 1.0 / (beta - 1.0));
      space_.interpolate(p, q, 1.0 - delta, centroid);
      distanceCentroidSqr = space_.distanceSqr(q, centroid);
    }
  }

  const VectorSpace<Point, Scalar>& space_;
};


template<typename Point, typename Scalar>
class BSkeleton :
    public GeometricTest<Point, Scalar, BSkeletonEdge<Point, Scalar> > {
 public:
  explicit BSkeleton(const VectorSpace<Point, Scalar>& space):
      GeometricTest<Point, Scalar, BSkeletonEdge<Point, Scalar> >(),
      space_(space), edge_(space) {
    space_.allocate(&rq);
  }
  virtual ~BSkeleton() {
    space_.deallocate(rq);
  }

  virtual Scalar shadowing(const BSkeletonEdge<Point, Scalar>& edge,
                           const Point& r) {
    Scalar beta =
        GeometricTest<Point, Scalar, BSkeletonEdge<Point, Scalar> >::getParam();
    if (beta <= 1.0) {
      space_.subtract(r, edge.q, rq);
      Scalar dot = space_.dot(edge.pq, rq);
      Scalar mag2sqr = space_.dot(rq, rq);
      Scalar test = (edge.square_length > 0 && mag2sqr > 0)
          ? fabs(dot) * dot / (edge.square_length * mag2sqr)
          : 1;
      Scalar maxdot = 1.0 - beta * beta;
      return maxdot - test;
    } else {
      space_.subtract(r, edge.planeCrossing, rq);
      Scalar dot = space_.dot(edge.pq, rq);
      if (beta > 1) {
        Scalar drc2 = space_.distanceSqr(r, edge.centroid);
        if (drc2 < edge.distanceCentroidSqr) {
          return edge.distanceCentroidSqr - drc2;
        }
      }
      return -dot;
    }
  }

  virtual void setActiveEdge(const Point& p, const Point& q) {
    Scalar beta =
        GeometricTest<Point, Scalar, BSkeletonEdge<Point, Scalar> >::getParam();
    edge_.reset(p, q, beta);
  }
  virtual BSkeletonEdge<Point, Scalar>& getActiveEdge() {
    return edge_;
  }

 private:
  const VectorSpace<Point, Scalar>& space_;
  BSkeletonEdge<Point, Scalar> edge_;
  Point rq;
};



template<typename Point, typename Scalar>
class BSkeletonBuilder :
    public NeighborhoodBuilder<Point, Scalar,
        BSkeleton<Point, Scalar> > {
 public:
  explicit BSkeletonBuilder(const VectorSpace<Point, Scalar>& space)
      : NeighborhoodBuilder<Point, Scalar,
          BSkeleton<Point, Scalar> >(space) {}
};

};  // namespace ngl

#endif  // NGL_BSKELETON_H
