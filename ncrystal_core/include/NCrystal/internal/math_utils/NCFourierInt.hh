#ifndef NCrystal_FourierInt_hh
#define NCrystal_FourierInt_hh

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
    // One-sided finite Fourier integrals of a function f = fr + i*fi sampled
    // on a uniform grid t_n = n*h (n=0..N-1, T=(N-1)*h), evaluated at w:
    //
    //   F(w) = (1/pi) Int_0^T [ fr(t) cos(wt) + fi(t) sin(wt) ] dt
    //
    // For f(-t) = conj(f(t)) this is the full (1/2pi) Int_{-T}^T f e^{-iwt} dt,
    // i.e. S(beta) from an intermediate function I(t) (fi empty means fi=0).
    //
    // Candidate schemes (benchmarking stage, the losers will be removed):
    //
    //   Trapezoid       : trapezoidal rule on f(t)e^{iwt}.
    //   TrapezoidEC     : as Trapezoid, plus the h^2/12 Euler-Maclaurin endpoint
    //                     correction (fr'(0)=0, other derivatives from
    //                     one-sided differences). Only valid for w*h << 1.
    //   PiecewiseLinear : f linear between samples, integrated exactly against
    //                     e^{iwt} (Filon-type, as in GASKET).
    //   CubicSpline     : f as a cubic spline (fr clamped fr'(0)=0, fi natural
    //                     at t=0, both natural at t=T), integrated exactly.
    //
    // Evaluation: Direct sums over all samples for each w (O(N) per w, uses
    // std::cos/sin). FFT* evaluates the periodic sample sums on a grid of
    // spacing 2pi/(M*h) (M = smallest power of two >= oversample*N) and
    // interpolates them (linearly or with 4-point Lagrange cubics), while the
    // exact piece integrals are still applied at w itself, so the piecewise
    // schemes remain valid beyond the Nyquist frequency pi/h in all modes.

    //
    // Optional tail (Trapezoid scheme only): Geometric continues the samples
    // beyond T as q_{N-1+k} = q_{N-1} r^k, r = q_{N-1}/q_{N-2} (q = fr - i*fi,
    // |r|<1 required), i.e. an exponential (or damped oscillating) decay, and
    // adds the exact trapezoid sum of that continuation for each w.

    enum class FourierScheme { Trapezoid, TrapezoidEC,
                               PiecewiseLinear, CubicSpline };
    enum class FourierEval { Direct, FFTLinear, FFTCubic };
    enum class FourierTail { None, Geometric };

    VectD fourierIntegral( double h, Span<const double> fr,
                           Span<const double> fi, Span<const double> w,
                           FourierScheme, FourierEval,
                           unsigned oversample = 8,
                           FourierTail = FourierTail::None );

  }
}

#endif
