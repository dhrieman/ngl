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
#include <stdio.h>
#include <string>

#include "typed-vectorspace.h"

using ngl::TypedVectorSpace;
using ngl::VectorSpace;
using std::string;

template<typename Scalar>
bool verifyScalar(const Scalar& expected, const Scalar& computed,
                  double tolerance, const string& operation) {
  if (expected - computed < tolerance
      && computed - expected < tolerance) {
    fprintf(stderr, "[\033[32m PASSED \033[0m] %s\n", operation.c_str());
    return true;
  } else {
    printf("Tolerance = %g\n", tolerance);
    fprintf(stderr, "[\033[31m FAILED \033[0m] %s. The expected and computed"
      " points differ more than %g.",
      operation.c_str(), tolerance);
    return false;
  }
}

template<typename Point, typename Scalar>
bool verifyResult(const VectorSpace<Point, Scalar>& space,
                  const Point& expected,
                  const Point& computed, double tolerance,
                  const string& operation) {
  double d2 = space.distanceSqr(expected, computed);
  return verifyScalar<Scalar>(d2, 0, tolerance, operation);
}

template<typename Type>
void RunTest(int dimensions, double tolerance) {
  Type *a = new Type[dimensions];
  Type *b = new Type[dimensions];
  Type *c = new Type[dimensions];
  Type *a_plus_b = new Type[dimensions];
  Type *a_minus_b = new Type[dimensions];
  Type *a_copy = new Type[dimensions];
  Type *a_plus_2_b = new Type[dimensions];
  Type *a_times_3 = new Type[dimensions];
  Type *midpoint_a_b = new Type[dimensions];
  Type *a_quarter_b = new Type[dimensions];
  Type expectedDot = 0;
  Type lengthA = 0;
  Type lengthB = 0;
  for (int i = 0; i < dimensions; ++i) {
    a[i] = (Type) i;
    b[i] = (Type) (i*i) / 2;
    a_plus_b[i] = (Type) (a[i] + b[i]);
    a_minus_b[i] = (Type) (a[i] - b[i]);
    a_copy[i] = (Type) (a[i]);
    a_plus_2_b[i] = (Type) (a[i] + 2 * b[i]);
    a_times_3[i] = (Type) (3 * a[i]);
    midpoint_a_b[i] = (Type) (a[i] + b[i]) / 2;
    a_quarter_b[i] = (Type) (a[i] * 0.75 + b[i] * 0.25);
    expectedDot += a[i] * b[i];
    lengthA += a[i] * a[i];
    lengthB += b[i] * b[i];
  }
  TypedVectorSpace<Type> space(dimensions);
  space.add(a, b, c);
  verifyResult(space, a_plus_b, c, tolerance, "add");

  space.subtract(a, b, c);
  verifyResult(space, a_minus_b, c, tolerance, "subtract");

  space.set(c, a);
  verifyResult(space, a_copy, c, tolerance, "set");

  space.muladd(a, b, 2, c);
  verifyResult(space, a_plus_2_b, c, tolerance, "muladd");

  space.mul(a, 3, c);
  verifyResult(space, a_times_3, c, tolerance, "mul");

  space.interpolate(a, b, 0.5, c);
  verifyResult(space, midpoint_a_b, c, tolerance, "interpolate 0.5");

  space.interpolate(a, b, 0.25, c);
  verifyResult(space, a_quarter_b, c, tolerance, "interpolate 0.25");

  Type dot = space.dot(a, b);
  verifyScalar(expectedDot, dot, tolerance, "dot");

  delete a;
  delete b;
  delete c;
  delete a_plus_b;
  delete a_minus_b;
  delete a_copy;
  delete a_plus_2_b;
  delete a_times_3;
  delete midpoint_a_b;
  delete a_quarter_b;
}

int main(int argc, char* argv[]) {
  int dims[] = { 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 100};
  for (int k = 0; k < 11; ++k) {
    fprintf(stderr, "Testing geometry dims: %d, type: double\n", dims[k]);
    RunTest<double>(dims[k], 1e-8);
  }
  for (int k = 0; k < 11; ++k) {
    fprintf(stderr, "Testing geometry dims: %d, type: float\n", dims[k]);
    RunTest<float>(dims[k], 1e-6);
  }
  for (int k = 0; k < 11; ++k) {
    fprintf(stderr, "Testing geometry dims: %d, type: int\n", dims[k]);
    RunTest<int>(dims[k], 1e-8);
  }
}
