#ifndef NCrystal_FMA_hh
#define NCrystal_FMA_hh

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

#include "NCrystal/core/NCDefs.hh"

//ncforceinline: force-inlining, e.g. for small helpers called from
//NCRYSTAL_FMADISPATCH_ATTR functions, so they are compiled into each dispatch
//variant of their caller (a non-inlined helper would only get the baseline ISA,
//and a separately dispatched one costs an ifunc call). In this case helpers
//must obey the same rules as their callers (e.g. explicit std::fma
//everywhere). Of course, force-inlining might have other uses outside the FMA
//dispatching:
#if defined(__GNUC__) || defined(__clang__)
#  define ncforceinline inline __attribute__((always_inline))
#elif defined(_MSC_VER)
#  define ncforceinline __forceinline
#else
#  define ncforceinline inline
#endif

namespace NCRYSTAL_NAMESPACE {

  //Fused multiply-add a*b+c with a single rounding, to be used instead of
  //std::fma. On targets known to lack hardware FMA (wasm) it is softFMA, since
  //std::fma there is a slow library routine. Elsewhere it is std::fma, so
  //compilers can still inline and vectorise it as before. Force-inlined, so
  //NCRYSTAL_FMADISPATCH_ATTR functions get hardware FMA in their FMA variants.
  ncforceinline double ncfma( double, double, double ) noexcept;

  //Correctly rounded a*b+c without hardware FMA, giving results identical to
  //std::fma. Faster than many software std::fma implementations (wasm: 5.6ns vs
  //14-19ns, glibc: 3.8ns vs 6-11ns), but slower than hardware FMA (1.8ns) and
  //hinders auto vectorisation.
  double softFMA( double a, double b, double c ) noexcept;

}


////////////////////////////
// Inline implementations //
////////////////////////////

#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("fp-contract=off")))
#endif
inline double NCrystal::softFMA( double a, double b, double c ) noexcept
{
  //Algorithm 5 of S. Boldo and G. Melquiond, "Emulation of a FMA and correctly
  //rounded sums: proved algorithms using rounding to odd", IEEE Trans. Comput.
  //57 (2008) 462: exact product (Veltkamp/Dekker), TwoSum, and an addition
  //rounded to odd (by adjusting the last bit). Falls back to std::fma where the
  //error-free transformations might not be exact (over/underflow, inf, nan).
  //Requires no FP contraction (disabled by an attribute for GCC and a pragma
  //for Clang).
#ifdef __clang__
#  pragma clang fp contract(off)
#endif
  const double aa = std::fabs(a), ab = std::fabs(b), ac = std::fabs(c);
  const double p = aa * ab;
  if ( !( p < 1e300 && ac < 1e300 && aa < 1e290 && ab < 1e290
          && p > 1e-280 && ( ac == 0.0 || ac > 1e-280 ) ) )
    return std::fma( a, b, c );
  //Exact product uh+ul=a*b:
  const double uh = a * b;
  constexpr double splitter = 134217729.0;//2^27+1
  const double ta = splitter * a, ah = ta - ( ta - a ), al = a - ah;
  const double tb = splitter * b, bh = tb - ( tb - b ), bl = b - bh;
  const double ul = ( ( ah * bh - uh ) + ah * bl + al * bh ) + al * bl;
  //TwoSum th+tl=c+uh:
  const double th = c + uh;
  const double tbb = th - c;
  const double tl = ( c - ( th - tbb ) ) + ( uh - tbb );
  //v = tl+ul rounded to odd:
  double v = tl + ul;
  const double vb = v - tl;
  const double ve = ( tl - ( v - vb ) ) + ( ul - vb );
  if ( ve != 0.0 ) {
    std::uint64_t bits;
    std::memcpy( &bits, &v, sizeof(bits) );
    if ( !( bits & 1u ) ) {
      //v is even, so move it one ulp towards tl+ul:
      bits += ( ( ve > 0.0 ) == ( v > 0.0 ) ? 1u : std::uint64_t(-1) );
      std::memcpy( &v, &bits, sizeof(bits) );
    }
  }
  return th + v;
}

ncforceinline double NCrystal::ncfma( double a, double b,
                                           double c ) noexcept
{
#ifdef __wasm__
  return softFMA( a, b, c );
#else
  return std::fma( a, b, c );
#endif
}

#endif
