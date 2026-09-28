#ifndef NCrystal_Math_FMA_hh
#define NCrystal_Math_FMA_hh

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

#include "NCrystal/internal/utils/NCFMADispatch.hh"
#include "NCrystal/internal/utils/NCMath.hh"

//The ordinary and MSVC/x86-only Windows-fast variants of NCMath.cc's
//expm1_reducedarg_taylor14/detail_stable_expm1 -- see NCFMADispatch.hh for
//the macros used below and the full mechanism/rationale, and
//NCMath_WINFMA.cc for the other half.
//
//expm1_reducedarg_taylor14 and the Cody-Waite constants are used only from
//within the including TU (by detail_stable_expm1 here, and by
//NC::stable_sinh's large-x branch in NCMath.cc), so unlike
//detail_stable_expm1 they need no WINFMA_DECLARE/FORWARD pair: written
//once, unconditionally, as plain anonymous-namespace entities, simply
//compiled under whatever flags the *including* translation unit uses --
//once NCMath_WINFMA.cc dispatches into detail_stable_expm1 via its own
///arch:AVX2 build, the helper (also compiled into that same TU) already
//benefits, with no separate twin needed.
//
//detail_stable_expm1 uses NCRYSTAL_FMADISPATCH_DECLARATOR_C (not the plain
//DECLARATOR): it is already extern "C" + NCRYSTAL_APPLY_C_NAMESPACE-wrapped
//in its ordinary form too, for the unrelated (Apple/Mach-O) rule-5 reason
//-- see docs/devel_fma_attribute.md. (The public NC::stable_exp is a
//trivial inline wrapper around stable_expm1 in NCMath.hh and needs none of
//this machinery.)

namespace NCRYSTAL_NAMESPACE {
  namespace {

    //Cody-Waite argument reduction constants for splitting x = n*ln(2) + r
    //with |r|<=ln(2)/2: ln2_hi carries only the top ~32 bits of ln(2) (all
    //lower bits exactly 0), so n*ln2_hi is exact for any |n| relevant here,
    //and the two-step fma reduction recovers r to near full double
    //precision even though x and n*ln(2) can be far larger than r itself.
    //Standard constants (as used by e.g. fdlibm's exp()). Also used by the
    //large-x branch of NC::stable_sinh in NCMath.cc:
    constexpr double ln2_hi = 6.93147180369123816490e-01;
    constexpr double ln2_lo = 1.90821492927058770002e-10;
    constexpr double invln2 = 1.44269504088896338700e+00;

    //expm1(r)=e^r-1 via a 14th order Taylor expansion (terms r^1/1! to
    //r^14/14!), Horner's method with std::fma throughout for reproducibility
    //(same technique as safe_xcothx in NCVDOSEval.cc). Only accurate for
    //the reduced range |r|<=ln(2)/2~0.3466 produced by the Cody-Waite
    //reduction in detail_stable_expm1 below: truncation error there is
    //~4e-18 relative (first omitted term, r^15/15!), comfortably below
    //double precision. Coefficients are exact rationals rounded to the
    //nearest double, so identical on every platform by construction:
    NCRYSTAL_FMADISPATCH_ATTR
    double expm1_reducedarg_taylor14( double r )
    {
      constexpr double c1 = 1.0;
      constexpr double c2 = 1.0/2;
      constexpr double c3 = 1.0/6;
      constexpr double c4 = 1.0/24;
      constexpr double c5 = 1.0/120;
      constexpr double c6 = 1.0/720;
      constexpr double c7 = 1.0/5040;
      constexpr double c8 = 1.0/40320;
      constexpr double c9 = 1.0/362880;
      constexpr double c10 = 1.0/3628800;
      constexpr double c11 = 1.0/39916800;
      constexpr double c12 = 1.0/479001600;
      constexpr double c13 = 1.0/6227020800.0;
      constexpr double c14 = 1.0/87178291200.0;
      double p = c14;
      p = std::fma( r, p, c13 );
      p = std::fma( r, p, c12 );
      p = std::fma( r, p, c11 );
      p = std::fma( r, p, c10 );
      p = std::fma( r, p, c9 );
      p = std::fma( r, p, c8 );
      p = std::fma( r, p, c7 );
      p = std::fma( r, p, c6 );
      p = std::fma( r, p, c5 );
      p = std::fma( r, p, c4 );
      p = std::fma( r, p, c3 );
      p = std::fma( r, p, c2 );
      p = std::fma( r, p, c1 );
      return r*p;
    }

  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, stable_expm1,
    ( double x )
  )
  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_stable_expm1,stable_expm1)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(stable_expm1,(x));
    // Evaluating expm1(x) by first finding integer n so that
    // x=n*ln2+r and |r|<ln2/2. Then use:
    //
    // exp(x)-1 = exp(n*ln2+r)-1
    //          = exp(r)*exp(n*ln2)-1
    //          = exp(r)*(2**n)-1
    //          = expm1(r)*(2**n) + (2**n-1)
    //
    // Evaluate expm1(r) via Taylor, and 2**n = std::ldexp(1,n).
    //
    // And use special branches for large |x|.

    if ( x >= 709.0 ) {
      //Here n below would reach 1024 (for x >~ 709.44), overflowing the
      //pow2n of the main-path formula even though the result itself can
      //still be finite. Use the exact 2^n*e^r form directly instead; the
      //"-1" of expm1 is relatively ~1e-308 here, far below 1 ulp, so it is
      //simply dropped. The first double whose correctly rounded expm1
      //overflows is (empirically, via mpmath) 709.78271289338397, and
      //returning infinity explicitly there keeps FE_OVERFLOW from being
      //raised under FPE-trapping tests:
      if ( x >= 709.78271289338397 )
        return kInfinity;
      const double n = std::round( x * invln2 );
      double r = std::fma( -n, ln2_hi, x );
      r = std::fma( -n, ln2_lo, r );
      const double exp_r = expm1_reducedarg_taylor14( r ) + 1.0;
      return std::ldexp( exp_r, static_cast<int>(n) );
    }
    if ( x <= -40.0 )
      return -1.0;
    if ( ncisnan(x) )
      NCRYSTAL_THROW(BadInput,"stable_expm1 called with NaN");
    const double n = std::round( x * invln2 );
    double r = std::fma( -n, ln2_hi, x );
    r = std::fma( -n, ln2_lo, r );
    //|r| is now less than ln2/2 ~= 0.34657359028, so ok for 14th order taylor:
    const double expm1_of_r = expm1_reducedarg_taylor14( r );
    //static cast to int is safe since x in (-40,710):
    const double pow2n = std::ldexp( 1.0, static_cast<int>(n) );
    return std::fma( pow2n, expm1_of_r, pow2n - 1.0 );
  }

}

#endif
