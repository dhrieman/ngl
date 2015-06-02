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

#include "geometric-test.h"
#include "relative-neighbor.h"

namespace ngl {

template<typename Point, typename Scalar>
class GabrielEdge : public RelativeNeighborEdge<Point, Scalar> {
 public:
  Point midpoint;
  explicit GabrielEdge(const VectorSpace<Point, Scalar>& space) :
      RelativeNeighborEdge<Point, Scalar>(space),
      midpoint(NULL) {
    space.allocate(&midpoint);
  }

  virtual ~GabrielEdge() {
    RelativeNeighborEdge<Point, Scalar>::space_.deallocate(midpoint);
  }

  void reset(const Point& pin, const Point& qin) {
    RelativeNeighborEdge<Point, Scalar>::reset(pin, qin);
    RelativeNeighborEdge<Point, Scalar>::space_.interpolate(pin, qin, 0.5,
                                                            midpoint);
  }
};


template<typename Point, typename Scalar>
class Gabriel :
    public GeometricTest<Point, Scalar, GabrielEdge<Point, Scalar> > {
 public:
  explicit Gabriel(const VectorSpace<Point, Scalar>& space):
      GeometricTest<Point, Scalar, GabrielEdge<Point, Scalar> >(),
      space_(space), edge_(space) {
    space.allocate(&tmp);
  }
  virtual ~Gabriel() {
    space_.deallocate(tmp);
  }

  virtual Scalar shadowing(const GabrielEdge<Point, Scalar>& edge,
                           const Point& r) {
    space_.interpolate(edge.p, edge.q, 2.0, tmp);
    Scalar dp = space_.distanceSqr(r, edge.p);
    Scalar dq = space_.distanceSqr(r, tmp);
    return dq - dp;
  }

  virtual void setActiveEdge(const Point& p, const Point& q) {
    edge_.reset(p, q);
  }
  virtual GabrielEdge<Point, Scalar>& getActiveEdge() {
    return edge_;
  }

 private:
  const VectorSpace<Point, Scalar>& space_;
  GabrielEdge<Point, Scalar> edge_;
  Point tmp;
};


template<typename Point, typename Scalar>
class GabrielGraphBuilder :
    public NeighborhoodBuilder<Point, Scalar, Gabriel<Point, Scalar> > {
 public:
  explicit GabrielGraphBuilder(const VectorSpace<Point, Scalar>& space) :
      NeighborhoodBuilder<Point, Scalar, Gabriel<Point, Scalar> >(space) {
  }
};


};

#endif  // NGL_GABRIEL_H
