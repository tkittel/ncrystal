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

#if defined(_MSC_VER) && !defined(__clang__)
#  define NCMATHFMA_ALWAYS_INLINE __forceinline
#else
#  define NCMATHFMA_ALWAYS_INLINE inline __attribute__((always_inline))
#endif

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
    //The Horner sum p with expm1(r)=r*p, factored out (without the
    //dispatch attribute, so it can also inline fully into the ncerfc
    //machinery further below; every operation is an explicit std::fma,
    //so inlining into an fma clone cannot introduce contraction
    //differences). expm1_reducedarg_taylor14 wraps it with the exact
    //same operation sequence as always:
    NCMATHFMA_ALWAYS_INLINE
    double expm1_taylor14_over_r( double r )
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
      return p;
    }

    NCRYSTAL_FMADISPATCH_ATTR
    double expm1_reducedarg_taylor14( double r )
    {
      return r * expm1_taylor14_over_r( r );
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

  //////////////////////////////////////////////////////////////////////
  // ncerf/ncerfc: portable, libm-free erf and erfc (see NCMath.hh).   //
  // Rational approximations and region layout are those of W. J.      //
  // Cody's CALERF routine (netlib SPECFUN; W. J. Cody, "Rational      //
  // Chebyshev approximation for the error function", Math. Comp. 23   //
  // (1969) 631), with two modernisations: the exp(-x^2) tail factor   //
  // is computed via the fully portable detail_stable_exp above        //
  // (instead of libm exp), and the rounding of x^2 is compensated     //
  // with an exact std::fma residual split (instead of CALERF's        //
  // truncate-to-1/16ths trick). The rational recursions are written   //
  // as explicit std::fma Horner chains: CALERF's add-then-multiply    //
  // form is NOT contraction-proof (gcc's default -ffp-contract=fast   //
  // fuses each (a+C)*y multiply into the *following* addition across  //
  // statements, seen as last-bit Release-vs-Debug differences), so    //
  // the maximally-fused sequence is spelled out instead, making it    //
  // bit-identical on every platform and optimisation level.           //
  //////////////////////////////////////////////////////////////////////

  extern "C" double NCRYSTAL_APPLY_C_NAMESPACE(detail_ncerfc)( double x );
  extern "C" double NCRYSTAL_APPLY_C_NAMESPACE(detail_ncerf)( double x );

  namespace {

    //NB: the helpers below deliberately carry NCMATHFMA_ALWAYS_INLINE
    //and NOT NCRYSTAL_FMADISPATCH_ATTR: they must inline fully into both
    //clones of the dispatched detail_ncerf/detail_ncerfc entry points (a
    //decorated helper costs an ifunc call, measured 2-4x slower; without
    //forced inlining gcc leaves them as baseline-ISA functions with
    //library std::fma calls). Safe only because every contractable shape
    //in their bodies is already explicit std::fma:

    NCMATHFMA_ALWAYS_INLINE
    double ncerf_region1( double x )
    {
      //erf(x) for |x| <= 0.46875 (CALERF first interval):
      nc_assert( ncabs(x) <= 0.46875 );
      constexpr double A1 = 3.16112374387056560e0;
      constexpr double A2 = 1.13864154151050156e2;
      constexpr double A3 = 3.77485237685302021e2;
      constexpr double A4 = 3.20937758913846947e3;
      constexpr double A5 = 1.85777706184603153e-1;
      constexpr double B1 = 2.36012909523441209e1;
      constexpr double B2 = 2.44024637934444173e2;
      constexpr double B3 = 1.28261652607737228e3;
      constexpr double B4 = 2.84423683343917062e3;
      const double y = x * x;
      double xnum = std::fma( A5, y, A1 );
      double xden = y + B1;
      xnum = std::fma( xnum, y, A2 );
      xden = std::fma( xden, y, B2 );
      xnum = std::fma( xnum, y, A3 );
      xden = std::fma( xden, y, B3 );
      xnum = std::fma( xnum, y, A4 );
      xden = std::fma( xden, y, B4 );
      return x * xnum / xden;
    }

    NCMATHFMA_ALWAYS_INLINE
    double ncerfc_expneg_core( double z )
    {
      //exp(-z) for 0<=z<=~708, to full stable_exp-style precision:
      //the same Cody-Waite reduction (file-scope ln2_hi/ln2_lo/invln2
      //constants above) and 14-term expm1 Taylor as detail_stable_expm1
      //above, but reducing the negated argument directly so no
      //reciprocal division is needed, and defined locally (instead of
      //calling the exported entry point) so it inlines -- cf. the note
      //above:
      nc_assert( z >= 0.0 && z < 708.0 );
      const double n = std::nearbyint( z * invln2 );
      double t = std::fma( -n, ln2_hi, z );
      t = std::fma( -n, ln2_lo, t );
      const double r = -t;//so exp(-z)=2^(-n)*exp(r), |r|<=0.347
      const double p = expm1_taylor14_over_r( r );
      return std::ldexp( std::fma( r, p, 1.0 ), -int(n) );
    }

    NCMATHFMA_ALWAYS_INLINE
    double ncerfc_expmxsq_times( double y, double r )
    {
      //Evaluates exp(-y^2)*r without the naive y*y rounding loss:
      //y2lo=fma(y,y,-y2hi) is the exact rounding residual of y*y, and
      //exp(-y2hi-y2lo) = exp(-y2hi)*(1-y2lo) to within ~1e-27 (|y2lo|
      //<= 0.5*ulp(y2hi) <= ~4e-14 for the y^2<=705 relevant here):
      const double y2hi = y * y;
      const double y2lo = std::fma( y, y, -y2hi );
      const double e = ncerfc_expneg_core( y2hi );
      return std::fma( -y2lo, e, e ) * r;
    }

    NCMATHFMA_ALWAYS_INLINE
    double ncerfc_ypositive( double y )
    {
      //erfc(y) for y > 0.46875.
      nc_assert( !(y<=0.46875) );
      if ( y <= 4.0 ) {
        //CALERF second interval:
        constexpr double C1 = 5.64188496988670089e-1;
        constexpr double C2 = 8.88314979438837594e0;
        constexpr double C3 = 6.61191906371416295e1;
        constexpr double C4 = 2.98635138197400131e2;
        constexpr double C5 = 8.81952221241769090e2;
        constexpr double C6 = 1.71204761263407058e3;
        constexpr double C7 = 2.05107837782607147e3;
        constexpr double C8 = 1.23033935479799725e3;
        constexpr double C9 = 2.15311535474403846e-8;
        constexpr double D1 = 1.57449261107098347e1;
        constexpr double D2 = 1.17693950891312499e2;
        constexpr double D3 = 5.37181101862009858e2;
        constexpr double D4 = 1.62138957456669019e3;
        constexpr double D5 = 3.29079923573345963e3;
        constexpr double D6 = 4.36261909014324716e3;
        constexpr double D7 = 3.43936767414372164e3;
        constexpr double D8 = 1.23033935480374942e3;
        double xnum = std::fma( C9, y, C1 );
        double xden = y + D1;
        xnum = std::fma( xnum, y, C2 );
        xden = std::fma( xden, y, D2 );
        xnum = std::fma( xnum, y, C3 );
        xden = std::fma( xden, y, D3 );
        xnum = std::fma( xnum, y, C4 );
        xden = std::fma( xden, y, D4 );
        xnum = std::fma( xnum, y, C5 );
        xden = std::fma( xden, y, D5 );
        xnum = std::fma( xnum, y, C6 );
        xden = std::fma( xden, y, D6 );
        xnum = std::fma( xnum, y, C7 );
        xden = std::fma( xden, y, D7 );
        xnum = std::fma( xnum, y, C8 );
        xden = std::fma( xden, y, D8 );
        const double r = xnum / xden;
        return ncerfc_expmxsq_times( y, r );
      }
      if ( !( y < 26.543 ) ) {
        //Result would be below ~1e-308 (CALERF's XBIG; same policy as
        //erfcdiff's cutoff in NCMath.cc). The negated comparison also
        //routes NaN here (all its ordinary comparisons are false),
        //since it must not reach ncerfc_expneg_core below: the int
        //cast in its Cody-Waite reconstruction is undefined for NaN
        //(and its input assert rejects it):
        return ncisnan( y ) ? y : 0.0;
      }
      //CALERF third interval (4 < y < 26.543):
      constexpr double P1 = 3.05326634961232344e-1;
      constexpr double P2 = 3.60344899949804439e-1;
      constexpr double P3 = 1.25781726111229246e-1;
      constexpr double P4 = 1.60837851487422766e-2;
      constexpr double P5 = 6.58749161529837803e-4;
      constexpr double P6 = 1.63153871373020978e-2;
      constexpr double Q1 = 2.56852019228982242e0;
      constexpr double Q2 = 1.87295284992346047e0;
      constexpr double Q3 = 5.27905102951428412e-1;
      constexpr double Q4 = 6.05183413124413191e-2;
      constexpr double Q5 = 2.33520497626869185e-3;
      const double t = 1.0 / ( y * y );
      double xnum = std::fma( P6, t, P1 );
      double xden = t + Q1;
      xnum = std::fma( xnum, t, P2 );
      xden = std::fma( xden, t, Q2 );
      xnum = std::fma( xnum, t, P3 );
      xden = std::fma( xden, t, Q3 );
      xnum = std::fma( xnum, t, P4 );
      xden = std::fma( xden, t, Q4 );
      xnum = std::fma( xnum, t, P5 );
      xden = std::fma( xden, t, Q5 );
      double r = t * xnum / xden;
      r = ( kInvSqrtPi - r ) / y;
      return ncerfc_expmxsq_times( y, r );
    }

  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, ncerfc,
    ( double x )
  )
  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, ncerf,
    ( double x )
  )

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_ncerfc,ncerfc)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(ncerfc,(x));
    if ( ncabs(x) <= 0.46875 )
      return 1.0 - ncerf_region1( x );
    if ( x > 0.0 )
      return ncerfc_ypositive( x );
    //NB: written via std::fma (an *exact* product for b=-1.0, so
    //value-identical to plain subtraction) to prevent the inlined
    //multiply tail of ncerfc_ypositive from being contracted into
    //this subtraction on fma-capable targets. Also correctly
    //propagates NaN:
    return std::fma( -1.0, ncerfc_ypositive( -x ), 2.0 );
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_ncerf,ncerf)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(ncerf,(x));
    if ( ncabs(x) <= 0.46875 )
      return ncerf_region1( x );
    //No cancellation concern: erfc < 0.51 for the |x|>0.46875 args
    //reaching this point (fma with -1.0 for the same contraction
    //reason as in detail_ncerfc, value-identical to subtraction):
    const double e = std::fma( -1.0, ncerfc_ypositive( ncabs(x) ), 1.0 );
    return x > 0.0 ? e : -e;
  }

  //////////////////////////////////////////////////////////////////////
  // Fast trig (see NCMath.hh/.cc): the Taylor kernels behind          //
  // sincos_mpi2pi2, cos_mpipi, cos_mpi2pi2 and sin_mpi2pi2, written   //
  // as explicit std::fma Horner chains. Their original plain          //
  // add-then-multiply form was contracted differently by different    //
  // compilers/targets (gcc -ffp-contract=fast even across             //
  // statements, clang within expressions, baseline x86-64 not at      //
  // all), observed as last-bit F2 differences flipping knife-edge     //
  // counts in app_fillhkl on the aarch64 CI legs. The explicit        //
  // maximally-fused sequence below is the single canonical            //
  // evaluation on every platform and optimisation level.              //
  //////////////////////////////////////////////////////////////////////

  namespace {

    NCMATHFMA_ALWAYS_INLINE
    double trig_cos22_from_mx2( double mx2 )
    {
      //cos(x) = 1 + mx2*Horner (Taylor to 22nd order, mx2=-x^2),
      //coefficients 1/(2n)! exactly as in the original NCMath.cc:
      double q =           8.89679139245057328674889744250246834331248809e-22;
      q = std::fma( mx2, q, 4.1103176233121648584779906184361403746103695e-19 );
      q = std::fma( mx2, q, 1.56192069685862264622163643500573334235194041e-16 );
      q = std::fma( mx2, q, 4.77947733238738529743820749111754402759693765e-14 );
      q = std::fma( mx2, q, 1.14707455977297247138516979786821056662326504e-11 );
      q = std::fma( mx2, q, 2.08767569878680989792100903212014323125434237e-9 );
      q = std::fma( mx2, q, 2.75573192239858906525573192239858906525573192e-7 );
      q = std::fma( mx2, q, 2.48015873015873015873015873015873015873015873e-5 );
      q = std::fma( mx2, q, 1.38888888888888888888888888888888888888888889e-3 );
      q = std::fma( mx2, q, 4.16666666666666666666666666666666666666666667e-2 );
      q = std::fma( mx2, q, 0.5 );
      return std::fma( mx2, q, 1.0 );
    }

    NCMATHFMA_ALWAYS_INLINE
    double trig_sin19_over_x_from_mx2( double mx2 )
    {
      //sin(x)/x = 1 + mx2*Horner (Taylor to 19th order, mx2=-x^2),
      //coefficients 1/(2n+1)! exactly as in the original NCMath.cc:
      double p =           8.22063524662432971695598123687228074922073899e-18;
      p = std::fma( mx2, p, 2.81145725434552076319894558301032001623349274e-15 );
      p = std::fma( mx2, p, 7.64716373181981647590113198578807044415510024e-13 );
      p = std::fma( mx2, p, 1.60590438368216145993923771701549479327257105e-10 );
      p = std::fma( mx2, p, 2.50521083854417187750521083854417187750521084e-8 );
      p = std::fma( mx2, p, 2.75573192239858906525573192239858906525573192e-6 );
      p = std::fma( mx2, p, 1.98412698412698412698412698412698412698412698e-4 );
      p = std::fma( mx2, p, 8.33333333333333333333333333333333333333333333e-3 );
      p = std::fma( mx2, p, 1.66666666666666666666666666666666666666666667e-1 );
      return std::fma( mx2, p, 1.0 );
    }

  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    void, sincos_mpi2pi2,
    ( double A, double* cosA, double* sinA )
  )
  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, cos_mpipi,
    ( double A )
  )
  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, cos_mpi2pi2,
    ( double x )
  )
  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, sin_mpi2pi2,
    ( double x )
  )

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(void,detail_sincos_mpi2pi2,sincos_mpi2pi2)
  ( double A, double* cosA, double* sinA )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(sincos_mpi2pi2,(A,cosA,sinA));
    //NB: no nc_assert here or in the other extern "C" functions below
    //(it can throw, and MSVC assumes extern "C" functions do not:
    //C4297) -- input range asserts live in the NCMath.cc wrappers.
    //Evaluate at A/2 via (shorter, 15th/16th order) Taylor expansions
    //and get final results via double-angle formulas (exactly the
    //original NCMath.cc algorithm, with the Horner chains and the
    //double-angle reconstruction spelled out as explicit fma):
    const double x = 0.5 * A;
    const double mx2 = -( x * x );
    double p =           7.64716373181981647590113198578807044415510024e-13;
    p = std::fma( mx2, p, 1.60590438368216145993923771701549479327257105e-10 );
    p = std::fma( mx2, p, 2.50521083854417187750521083854417187750521084e-8 );
    p = std::fma( mx2, p, 2.75573192239858906525573192239858906525573192e-6 );
    p = std::fma( mx2, p, 1.98412698412698412698412698412698412698412698e-4 );
    p = std::fma( mx2, p, 8.33333333333333333333333333333333333333333333e-3 );
    p = std::fma( mx2, p, 1.66666666666666666666666666666666666666666667e-1 );
    const double s2 = x * std::fma( mx2, p, 1.0 );
    double q =           4.77947733238738529743820749111754402759693765e-14;
    q = std::fma( mx2, q, 1.14707455977297247138516979786821056662326504e-11 );
    q = std::fma( mx2, q, 2.08767569878680989792100903212014323125434237e-9 );
    q = std::fma( mx2, q, 2.75573192239858906525573192239858906525573192e-7 );
    q = std::fma( mx2, q, 2.48015873015873015873015873015873015873015873e-5 );
    q = std::fma( mx2, q, 1.38888888888888888888888888888888888888888889e-3 );
    q = std::fma( mx2, q, 4.16666666666666666666666666666666666666666667e-2 );
    q = std::fma( mx2, q, 0.5 );
    const double c2m1 = mx2 * q;//cos(A/2)-1
    const double k = 2.0 * c2m1;
    //sin(A)=2*sin(A/2)*cos(A/2)=(k+2)*s2; cos(A)=1+k*(c2m1+2)/... as in
    //the original, with each mul+add shape made explicit:
    *sinA = std::fma( 2.0, c2m1, 2.0 ) * s2;
    *cosA = std::fma( k, std::fma( mx2, q, 2.0 ), 1.0 );
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_cos_mpi2pi2,cos_mpi2pi2)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(cos_mpi2pi2,(x));
    return trig_cos22_from_mx2( -( x * x ) );
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_cos_mpipi,cos_mpipi)
  ( double A )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(cos_mpipi,(A));
    //abs/min/copysign tricks reduce the evaluation to [-pi/2,pi/2]:
    const double Aabs = ncabs( A );
    const double x = ncmin( Aabs, kPi - Aabs );
    const double c = trig_cos22_from_mx2( -( x * x ) );
    return std::copysign( c, kPiHalf - Aabs );
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_sin_mpi2pi2,sin_mpi2pi2)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(sin_mpi2pi2,(x));
    return x * trig_sin19_over_x_from_mx2( -( x * x ) );
  }

}

#endif
