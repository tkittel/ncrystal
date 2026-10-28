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
#ifdef NCRYSTAL_UNPREDICTABLE
#  undef NCRYSTAL_UNPREDICTABLE
#endif
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
  //identical to the standard lib functions.
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

  //Guide table for searches in a fixed sorted array (see RandIdxByWeightPicker
  //in NCRandUtils.hh for the method and references), but with its steps taken
  //in the bit pattern of the values rather than in the values. For positive
  //doubles the bit pattern orders like the values and is roughly linear in
  //log2(x), so a subtraction and a shift select the step (no log or division),
  //giving roughly log-spaced steps. Useful for sorted arrays of positive finite
  //values roughly evenly spaced in log (e.g. energy grids spanning orders of
  //magnitude), for any searched values. Results are identical to
  //fastLowerBoundIdx/fastUpperBoundIdx, also outside the array range. The
  //array is passed to the search functions, and must be the one the table was
  //built from.
  class LogGuideTable final {
  public:
    LogGuideTable() = default;
    LogGuideTable( const double* x, std::size_t n );
    std::size_t lowerBoundIdx( const double* x, std::size_t n, double val ) const ncnoexceptndebug;
    std::size_t upperBoundIdx( const double* x, std::size_t n, double val ) const ncnoexceptndebug;
  private:
    std::vector<std::uint32_t> m_guide;
    std::uint64_t m_base = 0;
    unsigned m_shift = 0;
    static std::uint64_t bits( double v ) noexcept
    {
      std::uint64_t b;
      std::memcpy( &b, &v, sizeof(b) );
      return b;
    }
  };

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

inline NCrystal::LogGuideTable::LogGuideTable( const double* x,
                                                       std::size_t n )
{
  static_assert( sizeof(double) == sizeof(std::uint64_t), "" );
  static_assert( std::numeric_limits<double>::is_iec559, "" );
#ifndef NDEBUG
  nc_assert_always( n >= 1 && n < std::numeric_limits<std::uint32_t>::max() );
  nc_assert_always( x[0] > 0.0 && std::isfinite( x[n-1] ) );
  nc_assert_always( std::is_sorted( x, x + n ) );
#endif
  m_base = bits( x[0] );
  const std::uint64_t range = bits( x[n-1] ) - m_base;
  //Smallest power-of-two bucket width giving at most ~4 buckets per point:
  while ( ( range >> m_shift ) >= 4 * static_cast<std::uint64_t>( n ) )
    ++m_shift;
  const std::size_t nbuckets = static_cast<std::size_t>( range >> m_shift ) + 1;
  m_guide.reserve( nbuckets + 1 );
  std::size_t j = 0;
  for ( auto b : ncrange( nbuckets ) ) {
    //Index of first element with bit pattern >= the bucket start:
    const std::uint64_t start = m_base + ( static_cast<std::uint64_t>( b ) << m_shift );
    while ( j < n && bits( x[j] ) < start )
      ++j;
    m_guide.push_back( static_cast<std::uint32_t>( j ) );
  }
  m_guide.push_back( static_cast<std::uint32_t>( n ) );
}

inline std::size_t NCrystal::LogGuideTable::lowerBoundIdx( const double* x,
                                                               std::size_t n,
                                                               double val ) const ncnoexceptndebug
{
  nc_assert( !m_guide.empty() );
  if ( !( val > x[0] && val <= x[n-1] ) )
    return fastLowerBoundIdx( x, n, val );
  //Elements before m_guide[b] are below val, those from m_guide[b+1] above:
  const std::size_t b = static_cast<std::size_t>( ( bits( val ) - m_base ) >> m_shift );
  nc_assert( b + 1 < m_guide.size() );
  const std::size_t lo = m_guide[b];
  const std::size_t hi = m_guide[b+1];
  return lo + ( hi > lo ? fastLowerBoundIdx( x + lo, hi - lo, val ) : 0 );
}

inline std::size_t NCrystal::LogGuideTable::upperBoundIdx( const double* x,
                                                               std::size_t n,
                                                               double val ) const ncnoexceptndebug
{
  nc_assert( !m_guide.empty() );
  if ( !( val > x[0] && val < x[n-1] ) )
    return fastUpperBoundIdx( x, n, val );
  const std::size_t b = static_cast<std::size_t>( ( bits( val ) - m_base ) >> m_shift );
  nc_assert( b + 1 < m_guide.size() );
  const std::size_t lo = m_guide[b];
  const std::size_t hi = m_guide[b+1];
  return lo + ( hi > lo ? fastUpperBoundIdx( x + lo, hi - lo, val ) : 0 );
}

#undef NCRYSTAL_UNPREDICTABLE

#endif
