#ifndef NCrystal_Convolve_hh
#define NCrystal_Convolve_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/utils/NCSpan.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace MathUtils {

    ////////////////////////////////////////////////////////////////////////////
    //
    // Convolution h(z) = Int f(x) g(z-x) dx of two functions given as
    // piecewise linear interpolations of (xf,f) and (xg,g) (strictly
    // increasing grids, at least two points each, zero outside the grids),
    // evaluated at z. Candidate schemes (benchmarking stage):
    //
    //   DirectTrapezoid : trapezoidal rule over the xf points, with g linearly
    //                     interpolated (as in the Liquids plugin).
    //   ExactLinear     : exact convolution of the two interpolants, as
    //                     area(g)*f(z-c) (c is the centroid of g) plus local
    //                     corrections from the double antiderivative of g
    //                     (the narrower function is used as g).
    //   FFT             : both resampled on a common uniform grid of the given
    //                     step (0: smallest grid spacing of the inputs),
    //                     discrete convolution via FastConvolve, and the
    //                     result linearly interpolated at z.
    //
    // NB: Resolving a narrow input on a wide grid is costly for FFT (and its
    // accuracy is limited by the resampling), it is the caller's job to split
    // such inputs into analytic and smooth parts where needed.

    enum class ConvolveScheme { DirectTrapezoid, ExactLinear, FFT };

    VectD convolve( Span<const double> xf, Span<const double> f,
                    Span<const double> xg, Span<const double> g,
                    Span<const double> z, ConvolveScheme,
                    double fft_step = 0.0 );

  }
}

#endif
