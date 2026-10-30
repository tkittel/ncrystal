#ifndef NCrystal_Lines_hh
#define NCrystal_Lines_hh

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

#include "NCrystal/internal/math_utils/NCConvolve.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace MathUtils {

    ////////////////////////////////////////////////////////////////////////////
    //
    // Discrete (line) spectra, e.g. from Einstein oscillators, which can not
    // be handled by Fourier integrals (their I(t) is periodic), but combine
    // with a continuous S(beta) by shifting and adding.
    //
    // shiftAndAdd: h(z) = sum_k w_k g(z - s_k), with g given by (xg,g)
    // (strictly increasing grid, zero outside) and interpolated as in
    // convolve(..) (cf. NCConvolve.hh).
    //
    // oscillatorLines: lines of an Einstein oscillator of energy betaE
    // (in kT) and dimensionless strength A = alpha*w/betaE (w the oscillator's
    // weight in the normalised VDOS, alpha in units where the recoil is
    // alpha*w). The lines are at s_k = -k*betaE for integer k, with weights
    //
    //   w_k = exp(-A coth(betaE/2)) I_k(y) e^{k betaE/2},  y = A/sinh(betaE/2)
    //
    // (I_k modified Bessel functions; computed as exp(-A tanh(betaE/4) + k
    // betaE/2) times e^{-y} I_k(y) from Miller's backward recurrence,
    // normalised exactly by I_0 + 2 sum_{k>0} I_k = e^y). They satisfy
    // detailed balance w_{-k} = e^{-k betaE} w_k, sum_k w_k = 1 and
    // sum_k w_k s_k = -A betaE. Lines with w_k < reltol*max(w) are omitted.
    // Results are ordered by increasing shift.

    VectD shiftAndAdd( Span<const double> xg, Span<const double> g,
                       Span<const double> shifts, Span<const double> weights,
                       Span<const double> z,
                       ConvolveInterp = ConvolveInterp::Linear );

    struct Lines {
      VectD shifts;
      VectD weights;
    };
    Lines oscillatorLines( double A, double betaE, double reltol = 1e-14 );

  }
}

#endif
