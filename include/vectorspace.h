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
#ifndef NGL_VECTORSPACE_H
#define NGL_VECTORSPACE_H

#include <math.h>

namespace ngl {

// Abstract class for a vector space using Point and Scalar for vectors and
// scalars, respectively.
template<typename Point, typename Scalar>
class VectorSpace {
 public:
  explicit VectorSpace(int dimensions) : dimensions_(dimensions) {}

  // L2-Distance between two points.
  inline Scalar distance(const Point &a, const Point &b) const {
    return sqrt(distanceSqr(a, b));
  }

  // Distance squared between two points.
  Scalar distanceSqr(const Point &a, const Point &b) const {
    Scalar dist2 = 0;
    for (unsigned int i = 0; i < dimensions_; ++i) {
      dist2 += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return dist2;
  }

  // Dot or inner product.
  Scalar dot(const Point &a, const Point &b) const {
    Scalar res = 0;
    for (unsigned int  i = 0; i < dimensions_; ++i) {
      res += a[i] * b[i];
    }
    return res;
  }

  // Normalize a point
  Scalar normalize(Point &a) const {
    Scalar adota = dot(a, a);
    Scalar lena = sqrt(adota);
    if (lena > 0) {
      for (unsigned int i = 0; i < dimensions_; ++i) {
        a[i] /= lena;
      }
    }
  }

  // c = a + b
  void add(const Point &a, const Point &b, Point &c) const {
    for (int i = 0; i < dimensions_; ++i) {
      c[i] = a[i] + b[i];
    }
  }

  // c = a - b
  void subtract(const Point &a, const Point &b, Point &c) const {
    for (int i = 0; i < dimensions_; ++i) {
      c[i] = a[i] - b[i];
    }
  }

  // dst = src
  void set(Point &dst, const Point &src) const {
    for (int i = 0; i < dimensions_; ++i) {
      dst[i] = src[i];
    }
  }

  // c = a + b * t
  void muladd(const Point &a, const Point &b, Scalar t, Point &c) const {
    for (int i = 0; i < dimensions_; ++i) {
      c[i] = a[i] + b[i] * t;
    }
  }

  // c = (1 - t) * a + t * b
  void interpolate(const Point &a, const Point &b, double t, Point &c) const {
    for (int i = 0; i < dimensions_; ++i) {
      c[i] = (1 - t) * a[i] + t * b[i];
    }
  }

  // b = t * a
  void mul(const Point &a, double t, Point &b) const {
    for (int i = 0; i < dimensions_; ++i) {
      b[i] = a[i] * t;
    }
  }

  // Allocate a point
  virtual void allocate(Point* p) const = 0;

  // Deallocate a point
  virtual void deallocate(Point p) const = 0;

  int dimensions() {
    return dimensions_;
  }
  int dimensions_;
};

};  // namespace ngl

#endif  // NGL_VECTORSPACE_H
