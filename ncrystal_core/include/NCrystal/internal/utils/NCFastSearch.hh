#ifndef NCrystal_FastSearch_hh
#define NCrystal_FastSearch_hh

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

//Hints that the following branch condition is data-dependent, encouraging
//a conditional move (cmov) instead of a mispredictable jump. Only Clang
//needs and exposes this (__builtin_unpredictable); GCC already emits cmov
//unaided, and other compilers get the plain condition. Verified via
//disassembly on both GCC and Clang.
#ifdef __clang__
#  define NCRYSTAL_UNPREDICTABLE(x) __builtin_unpredictable(x)
#else
#  define NCRYSTAL_UNPREDICTABLE(x) (x)
#endif

namespace NCRYSTAL_NAMESPACE {

  //Drop-in replacements for std::lower_bound/std::upper_bound on a sorted
  //array of double, returning the same index that
  //(std::lower_bound/upper_bound(first,first+n,val)-first) would, but
  //written so that the inner comparison compiles to a branchless
  //conditional move rather than an ordinary conditional branch.
  //
  //Same precondition and behaviour as the standard library functions: first
  //must be sorted in non-decreasing order (ties/duplicate values are fine,
  //e.g. a cumulative-sum array with some zero-weight entries), n>=1.
  //
  //Rationale: a binary search's real-hardware cost is dominated by branch
  //mispredictions (each comparison is ~50/50 for random queries), which a
  //conditional move avoids. Measured 3-5x faster than
  //std::lower_bound/upper_bound for cache-resident arrays, fading
  //gradually for very large arrays but never slower. Results verified
  //identical to the standard library's (see tests/src/app_fastsearch).
  std::size_t fastLowerBoundIdx( const double* first, std::size_t n, double val ) noexcept;
  std::size_t fastUpperBoundIdx( const double* first, std::size_t n, double val ) noexcept;
  inline std::size_t fastLowerBoundIdx( Span<const double> s, double val ) noexcept
  {
    return fastLowerBoundIdx( s.data(), s.size(), val );
  }
  inline std::size_t fastUpperBoundIdx( Span<const double> s, double val ) noexcept
  {
    return fastUpperBoundIdx( s.data(), s.size(), val );
  }

}

////////////////////////////
// Inline implementations //
////////////////////////////

inline std::size_t NCrystal::fastLowerBoundIdx( const double* a, std::size_t n,
                                                double val ) noexcept
{
  std::size_t lo = 0;
  std::size_t len = n;
  while ( len > 1 ) {
    std::size_t half = len/2;
    bool advance = ( a[lo+half-1] < val );
    lo += NCRYSTAL_UNPREDICTABLE(advance) ? half : 0;
    len -= half;
  }
  if ( len && NCRYSTAL_UNPREDICTABLE( a[lo] < val ) )
    ++lo;
  return lo;
}

inline std::size_t NCrystal::fastUpperBoundIdx( const double* a, std::size_t n,
                                                double val ) noexcept
{
  std::size_t lo = 0;
  std::size_t len = n;
  while ( len > 1 ) {
    std::size_t half = len/2;
    bool advance = !( val < a[lo+half-1] );//a[lo+half-1] <= val
    lo += NCRYSTAL_UNPREDICTABLE(advance) ? half : 0;
    len -= half;
  }
  if ( len && NCRYSTAL_UNPREDICTABLE( !( val < a[lo] ) ) )
    ++lo;
  return lo;
}

#undef NCRYSTAL_UNPREDICTABLE

#endif
