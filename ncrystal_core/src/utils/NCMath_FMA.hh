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

namespace NCRYSTAL_NAMESPACE {
  namespace {

    //Cody-Waite argument reduction constants for splitting x = n*ln(2) + r
    //with |r|<=ln(2)/2: ln2_hi carries only the top ~32 bits of ln(2) (all
    //lower bits exactly 0), so n*ln2_hi is exact for any |n| relevant here,
    //and the two-step fma reduction recovers r to near full double
    //precision even though x and n*ln(2) can be far larger than r itself.
    //Standard constants:
    constexpr double ln2_hi = 6.93147180369123816490e-01;
    constexpr double ln2_lo = 1.90821492927058770002e-10;
    constexpr double invln2 = 1.44269504088896338700e+00;


    //Taylor expansion of (expm1(x)-x)/x^2, with the terms of expm1_taylor14
    //below (which is expm1(x) = x*(1+x*expm1_taylor14_p2(x))):
    NCRYSTAL_FMADISPATCH_INLINE double expm1_taylor14_p2( double x )
    {
      constexpr double c2 = 1.0/2.0;
      constexpr double c3 = 1.0/6.0;
      constexpr double c4 = 1.0/24.0;
      constexpr double c5 = 1.0/120.0;
      constexpr double c6 = 1.0/720.0;
      constexpr double c7 = 1.0/5040.0;
      constexpr double c8 = 1.0/40320.0;
      constexpr double c9 = 1.0/362880.0;
      constexpr double c10 = 1.0/3628800.0;
      constexpr double c11 = 1.0/39916800.0;
      constexpr double c12 = 1.0/479001600.0;
      constexpr double c13 = 1.0/6227020800.0;
      constexpr double c14 = 1.0/87178291200.0;
      double tmp = c14;
      tmp = std::fma( x, tmp, c13 );
      tmp = std::fma( x, tmp, c12 );
      tmp = std::fma( x, tmp, c11 );
      tmp = std::fma( x, tmp, c10 );
      tmp = std::fma( x, tmp, c9 );
      tmp = std::fma( x, tmp, c8 );
      tmp = std::fma( x, tmp, c7 );
      tmp = std::fma( x, tmp, c6 );
      tmp = std::fma( x, tmp, c5 );
      tmp = std::fma( x, tmp, c4 );
      tmp = std::fma( x, tmp, c3 );
      return std::fma( x, tmp, c2 );
    }

    //Taylor expansion of expm1(x)=exp(x)-1. Enough terms to be fully accurate
    //for |r|<ln2/2 and using std::fma for max portability.
    NCRYSTAL_FMADISPATCH_ATTR
    double expm1_taylor14( double x )
    {
      constexpr double c1 = 1.0;
      const double tmp = std::fma( x, expm1_taylor14_p2( x ), c1 );
      return x*tmp;
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
      const double exp_r = expm1_taylor14( r ) + 1.0;
      return std::ldexp( exp_r, static_cast<int>(n) );
    }
    if ( x <= -40.0 )
      return -1.0;
    const double n = std::round( x * invln2 );
    double r = std::fma( -n, ln2_hi, x );
    r = std::fma( -n, ln2_lo, r );
    //|r| is now less than ln2/2 ~= 0.34657359028, so ok for 14th order taylor:
    const double expm1_of_r = expm1_taylor14( r );
    //static cast to int is safe since x in (-40,710):
    const double pow2n = std::ldexp( 1.0, static_cast<int>(n) );
    return std::fma( pow2n, expm1_of_r, pow2n - 1.0 );
  }

  //////////////////////////////////////////////////////////////////////////////
  // ncerf/ncerfc. Written from scratch for NCrystal: no third-party erf/erfc
  // implementation was used or consulted. All constants in the generated block
  // are derived with mpmath by the NCDev script gen_ncerf_approx (see there
  // for details; it also simulates the code below in Python). Design:
  //
  //  |x|<1/2     : erf(x) = x*(c0+t*A(t)), t=x^2, c0=2/sqrt(pi) as hi+lo.
  //  1/2<=x<x0   : erfc(x) = exp(-x^2)*c3/(x+c3)*(1+eps(x)), c3~1/sqrt(pi)
  //                (exact asymptote), eps=P9/Q10 a small (<0.19) correction.
  //                exp(-x^2) = 2^-k*exp(r)*exp(c), r exact with |r|<~ln2/2
  //                and c tiny (incl. the rounding error of x*x).
  //  x>=x1 (x0)  : erf=1 (erfc=0), x1/x0 being the first doubles where the
  //                correctly rounded results reach those values.
  //
  // erf(-x)=-erf(x) and erfc(-x)=2-erfc(x). Every multiply-add is an explicit
  // std::fma and no plain product feeds a plain sum, so results are
  // bit-identical on all platforms. Worst errors vs. correctly rounded values
  // (50M points incl. dense windows at all thresholds): 1ULP for erf, 2ULP for
  // erfc (normal results; subnormal erfc results within 1 step).
  //////////////////////////////////////////////////////////////////////////////

  //FIXME: On targets without hardware FMA (WebAssembly has no fused
  //multiply-add instruction at all; also pre-2013 x86 CPUs) each std::fma here
  //becomes a libc software fma call: measured ~550ns/call in wasm (vs ~8ns
  //with hardware FMA, ~25-35x slower than musl's erf). Possible remedy (also
  //for all other explicit-fma code): an internal correctly rounded fma for
  //round-to-nearest and our argument ranges, built from error-free
  //transformations (cf. Boldo-Melquiond), verified bit-identical to hardware.

  namespace {

    //Generated by the NCDev script gen_ncerf_approx (do not edit):
    constexpr double ncerf_xsplit = 0.5;
    constexpr double ncerf_x1 = 5.9215871957945074;
    constexpr double ncerfc_x0 = 27.226017111108366;
    constexpr double ncerf_c0hi = 1.1283791670955126;
    constexpr double ncerf_c0lo = 1.5335459613165881e-17;
    constexpr double ncerf_c3 = 0.56418958354775628;
    constexpr double ncerf_invc3 = 1.7724538509055161;
    constexpr double ncerf_a0 = -0.37612638903183748;
    constexpr double ncerf_a1 = 0.1128379167095156;
    constexpr double ncerf_a2 = -0.026866170642132821;
    constexpr double ncerf_a3 = 0.0052239775293231074;
    constexpr double ncerf_a4 = -0.00085483118743757275;
    constexpr double ncerf_a5 = 0.00012054033933386649;
    constexpr double ncerf_a6 = -1.4863710204655132e-05;
    constexpr double ncerf_a7 = 1.4910186224329937e-06;
    constexpr double ncerf_p0 = 1.3002143073230382e-10;
    constexpr double ncerf_p1 = 0.6440746816773828;
    constexpr double ncerf_p2 = 1.2281051331750616;
    constexpr double ncerf_p3 = 1.1546262418511806;
    constexpr double ncerf_p4 = 0.68645270162227634;
    constexpr double ncerf_p5 = 0.28031988273418385;
    constexpr double ncerf_p6 = 0.080480474354960585;
    constexpr double ncerf_p7 = 0.015956492710316711;
    constexpr double ncerf_p8 = 0.0020249349889788631;
    constexpr double ncerf_p9 = 0.00013074452334787569;
    constexpr double ncerf_q1 = 3.4593893731651959;
    constexpr double ncerf_q2 = 5.5798104049836343;
    constexpr double ncerf_q3 = 5.5433376021453693;
    constexpr double ncerf_q4 = 3.7707246908040766;
    constexpr double ncerf_q5 = 1.8437666953360679;
    constexpr double ncerf_q6 = 0.66028814813330872;
    constexpr double ncerf_q7 = 0.17238436781955047;
    constexpr double ncerf_q8 = 0.03176078379702342;
    constexpr double ncerf_q9 = 0.0037944768362466832;
    constexpr double ncerf_q10 = 0.00023173863389154285;
    //End of generated block.

    //2^m for -1022<=m<=1023, from its bit pattern:
    NCRYSTAL_FMADISPATCH_INLINE double ncerf_pow2( int m )
    {
      const std::uint64_t b = static_cast<std::uint64_t>( m + 1023 ) << 52;
      double r;
      std::memcpy( &r, &b, sizeof(r) );
      return r;
    }

    //erf(x) for |x|<1/2 (monotonic, single final rounding):
    NCRYSTAL_FMADISPATCH_INLINE double ncerf_small( double x )
    {
      const double t = x*x;
      double a = ncerf_a7;
      a = std::fma( t, a, ncerf_a6 );
      a = std::fma( t, a, ncerf_a5 );
      a = std::fma( t, a, ncerf_a4 );
      a = std::fma( t, a, ncerf_a3 );
      a = std::fma( t, a, ncerf_a2 );
      a = std::fma( t, a, ncerf_a1 );
      a = std::fma( t, a, ncerf_a0 );
      //x*(c0hi+c0lo+t*A(t)):
      return std::fma( x, ncerf_c0hi, x * std::fma( t, a, ncerf_c0lo ) );
    }

    //For 1/2<=x<x0 returns v and k with erfc(x)=v*2^-k, v in ~[0.01,1]:
    struct ncerf_tail_t { double v; int k; };
    NCRYSTAL_FMADISPATCH_INLINE ncerf_tail_t ncerf_tail( double x )
    {
      //exp(-x^2) = 2^-k*(1+T), using x^2=h+l, -x^2 = -k*ln2 + r + c:
      const double h = x*x;
      const double l = std::fma( x, x, -h );
      const int k = static_cast<int>( std::fma( h, invln2, 0.5 ) );
      const double kd = static_cast<double>( k );
      const double r = std::fma( kd, ln2_hi, -h );//exact (k*ln2_hi exact)
      const double c = std::fma( kd, ln2_lo, -l );//|c|<~2e-7
      const double cc = std::fma( 0.5*c, c, c );//expm1(c)
      const double rq = r * expm1_taylor14_p2( r );
      const double em1 = std::fma( r, rq, r );//expm1(r)
      const double T = em1 + std::fma( cc, em1, cc );
      //eps(x):
      double p = ncerf_p9;
      p = std::fma( x, p, ncerf_p8 );
      p = std::fma( x, p, ncerf_p7 );
      p = std::fma( x, p, ncerf_p6 );
      p = std::fma( x, p, ncerf_p5 );
      p = std::fma( x, p, ncerf_p4 );
      p = std::fma( x, p, ncerf_p3 );
      p = std::fma( x, p, ncerf_p2 );
      p = std::fma( x, p, ncerf_p1 );
      p = std::fma( x, p, ncerf_p0 );
      double q = ncerf_q10;
      q = std::fma( x, q, ncerf_q9 );
      q = std::fma( x, q, ncerf_q8 );
      q = std::fma( x, q, ncerf_q7 );
      q = std::fma( x, q, ncerf_q6 );
      q = std::fma( x, q, ncerf_q5 );
      q = std::fma( x, q, ncerf_q4 );
      q = std::fma( x, q, ncerf_q3 );
      q = std::fma( x, q, ncerf_q2 );
      q = std::fma( x, q, ncerf_q1 );
      q = std::fma( x, q, 1.0 );
      const double eps = p / q;
      const double W = std::fma( eps, T, eps + T );//(1+eps)(1+T)-1
      //Anchor c3/(x+c3) = y0*(1+d), via exact 2Sum s+slo=x+c3 and the exact
      //division remainder c3-s*y0:
      const double s = x + ncerf_c3;
      const double sb = s - x;
      const double slo = ( x - ( s - sb ) ) + ( ncerf_c3 - sb );
      const double y0 = ncerf_c3 / s;
      const double d = std::fma( -slo, y0, std::fma( -s, y0, ncerf_c3 ) )
                       * ncerf_invc3;
      const double S = W + std::fma( W, d, d );//(1+W)(1+d)-1
      return { std::fma( y0, S, y0 ), k };
    }
  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, ncerf,
    ( double x )
  )
  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_ncerf,ncerf)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(ncerf,(x));
    if ( ncisnan( x ) )
      return x;
    const double ax = ncabs( x );
    if ( ax < ncerf_xsplit )
      return ncerf_small( x );
    if ( ax >= ncerf_x1 )
      return std::copysign( 1.0, x );
    const ncerf_tail_t tl = ncerf_tail( ax );
    return std::copysign( std::fma( -tl.v, ncerf_pow2( -tl.k ), 1.0 ), x );
  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, ncerfc,
    ( double x )
  )
  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_ncerfc,ncerfc)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(ncerfc,(x));
    if ( ncisnan( x ) )
      return x;
    const double ax = ncabs( x );
    if ( ax < ncerf_xsplit )
      return 1.0 - ncerf_small( x );//erfc>0.47, so no loss of precision
    if ( x >= ncerfc_x0 )
      return 0.0;
    if ( x <= -ncerf_x1 )
      return 2.0;
    const ncerf_tail_t tl = ncerf_tail( ax );
    if ( x < 0.0 )
      return std::fma( -tl.v, ncerf_pow2( -tl.k ), 2.0 );
    //Scale in two steps, so only the last (into subnormals) can be inexact:
    const int k1 = tl.k / 2;
    return ( tl.v * ncerf_pow2( -k1 ) ) * ncerf_pow2( k1 - tl.k );
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

    NCRYSTAL_FMADISPATCH_INLINE
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

    NCRYSTAL_FMADISPATCH_INLINE
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
    nc_assert( ncabs(A) <= kPiHalf );
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
    nc_assert( ncabs(x) <= kPiHalf );
    return trig_cos22_from_mx2( -( x * x ) );
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_cos_mpipi,cos_mpipi)
  ( double A )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(cos_mpipi,(A));
    //abs/min/copysign tricks reduce the evaluation to [-pi/2,pi/2]:
    const double Aabs = ncabs( A );
    nc_assert( Aabs <= kPi );
    const double x = ncmin( Aabs, kPi - Aabs );
    const double c = trig_cos22_from_mx2( -( x * x ) );
    return std::copysign( c, kPiHalf - Aabs );
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_sin_mpi2pi2,sin_mpi2pi2)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(sin_mpi2pi2,(x));
    nc_assert( ncabs(x) <= kPiHalf );
    return x * trig_sin19_over_x_from_mx2( -( x * x ) );
  }


}

#endif
