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
#ifndef NGL_H
#define NGL_H

#include "bskeleton.h"
#include "gabriel.h"
#include "neighborhood-builder.h"
#include "relative-neighbor.h"
#include "typed-vectorspace.h"
#include "neighbor-graph-impl.h"

namespace ngl {

typedef TypedVectorSpace<double> DoubleVectorSpace;

typedef NeighborhoodBuilder<double*, double> NeighborhoodBuilderD;
typedef RelativeNeighborGraphBuilder<double*, double>
    RelativeNeighborGraphBuilderD;
typedef GabrielGraphBuilder<double*, double> GabrielGraphBuilderD;
typedef BSkeletonBuilder<double*, double> BSkeletonBuilderD;

typedef TypedVectorSpace<float> FloatVectorSpace;
typedef NeighborhoodBuilder<float*, float> NeighborhoodBuilderF;
typedef RelativeNeighborGraphBuilder<float*, float>
    RelativeNeighborGraphBuilderF;
typedef GabrielGraphBuilder<float*, float> GabrielGraphBuilderF;
typedef BSkeletonBuilder<float*, float> BSkeletonBuilderF;

};

#endif  // NGL_H
