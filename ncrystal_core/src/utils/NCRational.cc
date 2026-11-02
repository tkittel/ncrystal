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

#include "NCrystal/internal/utils/NCRational.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {

    using int_t = Rational::int_t;
    constexpr int_t kMax = Rational::kMaxValue;

    [[noreturn]] void throwRange()
    {
      NCRYSTAL_THROW(CalcError,"Rational number exceeds the supported range"
                     " (magnitudes of numerator and denominator up to 2^53)");
    }

    int_t absv( int_t a ) noexcept { return a < 0 ? -a : a; }

    int_t gcdv( int_t a, int_t b ) noexcept
    {
      a = absv( a );
      b = absv( b );
      while ( b ) {
        const int_t t = a % b;
        a = b;
        b = t;
      }
      return a;
    }

    //Products and sums of values in [-kMax,kMax] (kMax=2^53) are exact in
    //int64 if the result is checked against kMax before use (sums never
    //overflow, products are checked by division first):
    int_t mulc( int_t a, int_t b )
    {
      if ( a != 0 && absv( b ) > kMax / absv( a ) )
        throwRange();
      return a * b;
    }

    int_t addc( int_t a, int_t b )
    {
      nc_assert( absv(a) <= kMax && absv(b) <= kMax );
      const int_t s = a + b;
      if ( absv( s ) > kMax )
        throwRange();
      return s;
    }

    //Normalise n/d (d!=0), with intermediate values in [-kMax,kMax]:
    void normalise( int_t& n, int_t& d )
    {
      nc_assert( d != 0 );
      if ( d < 0 ) {
        n = -n;
        d = -d;
      }
      const int_t g = gcdv( n, d );
      if ( g > 1 ) {
        n /= g;
        d /= g;
      }
      if ( n == 0 )
        d = 1;
    }
  }
}

NC::Rational::Rational( int_t n )
  : m_n( n ), m_d( 1 )
{
  if ( absv( n ) > kMax )
    throwRange();
}

NC::Rational::Rational( int_t n, int_t d )
{
  if ( d == 0 )
    NCRYSTAL_THROW(CalcError,"Rational number with zero denominator");
  if ( absv( n ) > kMax || absv( d ) > kMax )
    throwRange();
  normalise( n, d );
  m_n = n;
  m_d = d;
}

NC::Rational::int_t NC::Rational::floor() const noexcept
{
  //Integer division truncates towards zero:
  int_t q = m_n / m_d;
  if ( m_n % m_d != 0 && m_n < 0 )
    --q;
  return q;
}

NC::Rational NC::Rational::wrappedToUnitInterval() const
{
  const int_t f = floor();
  return Rational( NoCheck(), m_n - f * m_d, m_d );//(|f*d| <= |n|+d)
}

NC::Rational NC::Rational::operator+( const Rational& o ) const
{
  //a/b+c/d = (a*(d/g)+c*(b/g))/(b/g*d), with g=gcd(b,d):
  const int_t g = gcdv( m_d, o.m_d );
  int_t n = addc( mulc( m_n, o.m_d / g ), mulc( o.m_n, m_d / g ) );
  int_t d = mulc( m_d / g, o.m_d );
  normalise( n, d );
  return Rational( NoCheck(), n, d );
}

NC::Rational NC::Rational::operator-( const Rational& o ) const
{
  return *this + ( -o );
}

NC::Rational NC::Rational::operator*( const Rational& o ) const
{
  //Cross-reduce first, so the results are already normalised:
  const int_t g1 = gcdv( m_n, o.m_d );
  const int_t g2 = gcdv( o.m_n, m_d );
  if ( m_n == 0 || o.m_n == 0 )
    return Rational();
  int_t n = mulc( m_n / g1, o.m_n / g2 );
  int_t d = mulc( m_d / g2, o.m_d / g1 );
  normalise( n, d );
  return Rational( NoCheck(), n, d );
}

NC::Rational NC::Rational::operator/( const Rational& o ) const
{
  if ( o.m_n == 0 )
    NCRYSTAL_THROW(CalcError,"Rational division by zero");
  return *this * Rational( NoCheck(), ( o.m_n < 0 ? -o.m_d : o.m_d ),
                           absv( o.m_n ) );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Exact comparison a/b < c/d of fractions in [0,1) (0<=a<b, 0<=c<d),
    //without products which could overflow: compare the reciprocals' integer
    //parts, and recurse on their fractional parts (as in Euclid's algorithm,
    //so the number of steps is logarithmic):
    bool fracLess( int_t a, int_t b, int_t c, int_t d ) ncnoexceptndebug
    {
      bool flipped = false;
      while ( true ) {
        nc_assert( a >= 0 && a < b && c >= 0 && c < d );
        if ( a == 0 || c == 0 ) {
          const bool less = ( a == 0 && c != 0 );
          const bool greater = ( c == 0 && a != 0 );
          return flipped ? greater : less;
        }
        //a/b < c/d <=> b/a > d/c:
        const int_t qa = b / a;
        const int_t qc = d / c;
        if ( qa != qc )
          return flipped ? ( qa < qc ) : ( qa > qc );
        const int_t ra = b % a;
        const int_t rc = d % c;
        //Now compare ra/a and rc/c, with the sense flipped:
        b = a;
        a = ra;
        d = c;
        c = rc;
        flipped = !flipped;
      }
    }
  }
}

bool NC::Rational::operator<( const Rational& o ) const noexcept
{
  //Integer parts first, then the fractional parts (both in [0,1)):
  const int_t fa = floor();
  const int_t fb = o.floor();
  if ( fa != fb )
    return fa < fb;
  return fracLess( m_n - fa * m_d, m_d, o.m_n - fb * o.m_d, o.m_d );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Whether the correctly rounded double value of p/q (|p|,q <= 2^53, q>0)
    //is below (-1), equal to (0) or above (+1) x:
    int roundedCmp( int_t p, int_t q, double x ) noexcept
    {
      const double v = static_cast<double>( p ) / static_cast<double>( q );
      return v < x ? -1 : ( v > x ? 1 : 0 );
    }

    //The value (a+k*c)/(b+k*d) for k>=1 (a mediant-type step in the
    //Stern-Brocot tree), or false if its numerator or denominator would
    //exceed the limits:
    bool stepValue( int_t a, int_t b, int_t c, int_t d, int_t k, int_t maxden,
                    int_t& p, int_t& q ) ncnoexceptndebug
    {
      nc_assert( k >= 1 && b >= 1 && d >= 1 );
      if ( k > ( maxden - b ) / d )
        return false;
      q = b + k * d;//<= maxden
      const int_t ac = ( c < 0 ? -c : c );
      if ( ac != 0 && k > ( kMax - ( a < 0 ? -a : a ) ) / ac )
        return false;
      p = a + k * c;
      return true;
    }
  }
}

NC::Optional<NC::Rational>
NC::rationalFromDouble( double x, Rational::int_t maxDenominator )
{
  if ( !std::isfinite( x ) || !( ncabs( x ) <= double( kMax ) )
       || maxDenominator < 1 )
    return NullOpt;
  maxDenominator = ncmin( maxDenominator, kMax );
  const double fl = std::floor( x );
  if ( fl == x )
    return Rational( static_cast<int_t>( fl ) );

  //Walk the Stern-Brocot tree between the integers below and above x,
  //always keeping the bounds a/b and c/d on either side of the set of
  //rationals which round to x (an interval, since rounding is monotonic).
  //The first mediant in that interval is the simplest rational in it. Runs
  //of steps in the same direction are done at once, by finding the largest
  //k for which (a+k*c)/(b+k*d) (or (c+k*a)/(d+k*b)) stays on the same side
  //(exponential and then binary search), which keeps the number of
  //iterations logarithmic like for continued fractions:
  int_t a = static_cast<int_t>( fl ), b = 1;
  int_t c = a + 1, d = 1;
  while ( true ) {
    int_t p, q;
    if ( !stepValue( a, b, c, d, 1, maxDenominator, p, q ) )
      return NullOpt;
    const int s = roundedCmp( p, q, x );
    if ( s == 0 )
      return Rational( p, q );//(mediants are already reduced)
    //Moving the lower bound up (s<0), or the upper bound down (s>0):
    const int_t lo_n = ( s < 0 ? a : c ), lo_d = ( s < 0 ? b : d );
    const int_t hi_n = ( s < 0 ? c : a ), hi_d = ( s < 0 ? d : b );
    auto sameSide = [&]( int_t k, int_t& pk, int_t& qk )
    {
      return ( stepValue( lo_n, lo_d, hi_n, hi_d, k, maxDenominator, pk, qk )
               && roundedCmp( pk, qk, x ) == s );
    };
    int_t kgood = 1;//(k=1 is the mediant just tested)
    int_t kbad = 2;
    int_t pk, qk;
    while ( sameSide( kbad, pk, qk ) ) {
      kgood = kbad;
      if ( kbad > kMax / 2 )
        return NullOpt;
      kbad *= 2;
    }
    while ( kbad - kgood > 1 ) {
      const int_t kmid = kgood + ( kbad - kgood ) / 2;
      if ( sameSide( kmid, pk, qk ) )
        kgood = kmid;
      else
        kbad = kmid;
    }
    nc_assert_always( stepValue( lo_n, lo_d, hi_n, hi_d, kgood,
                                 maxDenominator, pk, qk ) );
    if ( s < 0 ) {
      a = pk;
      b = qk;
    } else {
      c = pk;
      d = qk;
    }
  }
}

std::ostream& NC::operator<<( std::ostream& os, const Rational& r )
{
  os << r.numerator();
  if ( !r.isInteger() )
    os << '/' << r.denominator();
  return os;
}
