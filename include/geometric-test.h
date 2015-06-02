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
#ifndef NGL_GEOMETRIC_TEST_H
#define NGL_GEOMETRIC_TEST_H

namespace ngl {

// Abstract class for a geometric test. It provides the method
// shadowing(edge, point), which determines when an edge shadows a point.
template<typename Point, typename Scalar, typename Edge>
class GeometricTest {
 public:
  GeometricTest() {}

  virtual ~GeometricTest() {}

  // Returns true if the umbra region defined by edge shadows point r
  inline bool shadows(const Edge& edge, const Point &r) {
    return shadowing(edge, r) < 0;
  }

  // Returns a real value denoting the shadowing of a point r by the
  // umbra region defined by edge. As a convention, a negative value
  // denotes a point at the interior of the umbra, a positive value at
  // the outside, and zero at the boundary.
  virtual Scalar shadowing(const Edge& edge, const Point& r) = 0;

  // Sets the active edge to pq.
  virtual void setActiveEdge(const Point& p, const Point& q) = 0;

  // Returns the active edge.
  virtual Edge& getActiveEdge() = 0;
  
  void setParam(Scalar param) {
    param_ = param;
  }

  Scalar getParam() {
    return param_;
  }
  
 protected:
  // A parameterization value for this geometric test,
  // e.g., beta in Beta Skeletons
  Scalar param_;
};

};  // namespace ngl

#endif  // NGL_GEOMETRIC_TEST_H
