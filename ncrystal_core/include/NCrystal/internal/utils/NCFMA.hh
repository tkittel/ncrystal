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
//variant of their caller. In this case helpers must obey the same rules as
//their callers (e.g. explicit ncfma everywhere).  Of course, force-inlining
//might have uses outside the FMA dispatching:
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
  //std::fma. Faster than many software std::fma implementations (wasm: 6.6ns vs
  //19ns, glibc: 4.8ns vs 6-11ns), but slower than hardware FMA (1.8ns) and
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
  std::uint64_t ba, bb, bc;
  std::memcpy( &ba, &a, sizeof(ba) );
  std::memcpy( &bb, &b, sizeof(bb) );
  std::memcpy( &bc, &c, sizeof(bc) );
  const unsigned ea = unsigned(ba>>52)&0x7ffu, eb = unsigned(bb>>52)&0x7ffu;
  const unsigned ec = unsigned(bc>>52)&0x7ffu;
  //Exponents of a, b and (non-zero) c must be within +-400, so nothing below
  //can overflow or underflow:
  constexpr unsigned lo = 1023u-400u, hi = 1023u+400u;
  const bool cok = ( ec - lo <= hi - lo ) || ( ( bc << 1 ) == 0 );
  if ( !( ( ea - lo <= hi - lo ) & ( eb - lo <= hi - lo ) & cok ) )
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
  //v = tl+ul rounded to odd (if inexact and even, move one ulp towards
  //tl+ul), without branches since that is often mispredicted:
  double v = tl + ul;
  const double vb = v - tl;
  const double ve = ( tl - ( v - vb ) ) + ( ul - vb );
  std::uint64_t bits, bve;
  std::memcpy( &bits, &v, sizeof(bits) );
  std::memcpy( &bve, &ve, sizeof(bve) );
  const std::uint64_t adj = std::uint64_t( ( ve != 0.0 ) & !( bits & 1u ) );
  const std::uint64_t dir = ( bve ^ bits ) >> 63;
  bits += adj - 2u*( adj & dir );
  std::memcpy( &v, &bits, sizeof(bits) );
  return th + v;
}

ncforceinline double NCrystal::ncfma( double a, double b,
                                           double c ) noexcept
{
#ifdef __wasm__
  return softFMA( a, b, c );
#elif defined(__GNUC__) || defined(__clang__)
  return __builtin_fma( a, b, c );//as std::fma, minus a call in debug builds
#else
  return std::fma( a, b, c );
#endif
}

#endif
