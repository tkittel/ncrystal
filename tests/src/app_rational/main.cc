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
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/core/NCFmt.hh"
#include <iostream>
#include <sstream>

namespace NC = NCrystal;
using R = NC::Rational;
using int_t = R::int_t;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  int_t refgcd( int_t a, int_t b )
  {
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;
    while ( b ) {
      int_t t = a % b;
      a = b;
      b = t;
    }
    return a;
  }

  //Reference: a/b normalised with plain integer arithmetic (small values):
  bool isNormalised( const R& r, int_t n, int_t d )
  {
    if ( d < 0 ) {
      n = -n;
      d = -d;
    }
    const int_t g = refgcd( n, d );
    if ( n == 0 )
      return r.numerator() == 0 && r.denominator() == 1;
    return r.numerator() == n / g && r.denominator() == d / g;
  }

  template<class TFct>
  void requireThrows( const char * descr, TFct fct )
  {
    bool ok = false;
    try {
      fct();
    } catch ( NC::Error::CalcError& e ) {
      std::cout << descr << ": CalcError: " << e.what() << '\n';
      ok = true;
    }
    REQUIRE( ok );
  }

  std::string str( const R& r )
  {
    std::ostringstream ss;
    ss << r;
    return ss.str();
  }
}

int main()
{
  //Construction and normalisation:
  REQUIRE( str( R() ) == "0" );
  REQUIRE( str( R(7) ) == "7" );
  REQUIRE( str( R(-7) ) == "-7" );
  REQUIRE( str( R(2,4) ) == "1/2" );
  REQUIRE( str( R(-2,-4) ) == "1/2" );
  REQUIRE( str( R(2,-4) ) == "-1/2" );
  REQUIRE( str( R(0,-5) ) == "0" );
  REQUIRE( R(0,-5).denominator() == 1 );
  REQUIRE( R(6,3).isInteger() && R(6,3).numerator() == 2 );
  REQUIRE( R(R::kMaxValue).numerator() == R::kMaxValue );
  REQUIRE( R(-R::kMaxValue,R::kMaxValue) == R(-1) );
  requireThrows( "zero denominator", [](){ R(1,0); } );
  requireThrows( "too large numerator", [](){ R(R::kMaxValue+1); } );
  requireThrows( "too large denominator", [](){ R(1,R::kMaxValue+1); } );
  requireThrows( "division by zero", [](){ R(1,2)/R(); } );

  //Arithmetic on a grid of small fractions against an integer reference:
  std::vector<std::pair<int_t,int_t>> vals;
  for ( int_t d = 1; d <= 24; ++d )
    for ( int_t n = -30; n <= 30; ++n )
      vals.emplace_back( n, d );
  std::size_t nchecks = 0;
  for ( auto& x : vals ) {
    const R a( x.first, x.second );
    REQUIRE( isNormalised( a, x.first, x.second ) );
    REQUIRE( a.toDouble() == double(x.first)/double(x.second) );
    REQUIRE( -(-a) == a );
    REQUIRE( a + R() == a && a * R(1) == a );
    //floor and wrapping:
    const int_t fl = a.floor();
    REQUIRE( R(fl) <= a && a < R(fl+1) );
    const R w = a.wrappedToUnitInterval();
    REQUIRE( R() <= w && w < R(1) && w == a - R(fl) );
    for ( auto& y : vals ) {
      const R b( y.first, y.second );
      const int_t an = x.first, ad = x.second, bn = y.first, bd = y.second;
      REQUIRE( isNormalised( a + b, an*bd + bn*ad, ad*bd ) );
      REQUIRE( isNormalised( a - b, an*bd - bn*ad, ad*bd ) );
      REQUIRE( isNormalised( a * b, an*bn, ad*bd ) );
      if ( bn != 0 )
        REQUIRE( isNormalised( a / b, an*bd, ad*bn ) );
      REQUIRE( ( a < b ) == ( an*bd < bn*ad ) );
      REQUIRE( ( a == b ) == ( an*bd == bn*ad ) );
      REQUIRE( ( a <= b ) == ( an*bd <= bn*ad ) );
      REQUIRE( a + b == b + a && a * b == b * a );
      REQUIRE( ( a - b ) + b == a );
      ++nchecks;
    }
  }
  std::cout << "Arithmetic and comparison of " << vals.size()
            << " small fractions in all " << nchecks
            << " combinations: OK" << std::endl;

  //Associativity and distributivity on a few triplets:
  const R t1(3,8), t2(-5,12), t3(7,9);
  REQUIRE( ( t1 + t2 ) + t3 == t1 + ( t2 + t3 ) );
  REQUIRE( t1 * ( t2 + t3 ) == t1 * t2 + t1 * t3 );
  REQUIRE( ( t1 * t2 ) * t3 == t1 * ( t2 * t3 ) );

  //The motivation: symmetry images of rational coordinates are exact, so
  //they round to the same double as the exact value would:
  REQUIRE( ( R(1) - R(9,10) ).toDouble() == 0.1 );
  REQUIRE( !( 1.0 - 0.9 == 0.1 ) );//(unlike double arithmetic)
  REQUIRE( ( R(1) - R(37,40) ).toDouble() == 0.075 );
  REQUIRE( ( R(2,3) - R(1,3) ) == R(1,3) );
  REQUIRE( ( R(-1,10) ).wrappedToUnitInterval().toDouble() == 0.9 );
  std::cout << "1-9/10 = " << ( R(1) - R(9,10) ) << " = "
            << NC::fmt( ( R(1) - R(9,10) ).toDouble() )
            << " (double arithmetic: " << NC::fmt( 1.0 - 0.9 ) << ")" << '\n';

  //Comparisons near the limits, where cross multiplication would overflow:
  constexpr int_t k = R::kMaxValue;
  REQUIRE( R(k-2,k-1) < R(k-1,k) );
  REQUIRE( !( R(k-1,k) < R(k-2,k-1) ) );
  REQUIRE( R(1,k) < R(1,k-1) );
  REQUIRE( R(-1,k-1) < R(-1,k) );
  REQUIRE( R(k-1,k) < R(1) && R(-1) < R(-(k-1),k) );
  REQUIRE( R(k-3,k-2) > R(k-4,k-3) );
  std::cout << "Comparisons near 2^53 (no overflow): OK" << std::endl;

  //Results outside the supported range throw rather than being inexact:
  requireThrows( "sum too large", [](){ R(k) + R(1); } );
  requireThrows( "product too large", [](){ R(k/2) * R(3); } );
  requireThrows( "denominators too large",
                 [](){ R(1,k) + R(1,k-1); } );
  REQUIRE( R(k/2) * R(2) == R(k) );//(exactly at the limit)

  //////////////////////////////////////////////////////////////////////////
  // rationalFromDouble                                                   //
  //////////////////////////////////////////////////////////////////////////

  auto fromDbl = []( double x, int_t maxden )
  {
    return NC::rationalFromDouble( x, maxden );
  };

  //All fractions a/b with |a|<=b<=300 are recovered from their doubles:
  std::size_t nfrac = 0;
  for ( int_t den = 1; den <= 300; ++den ) {
    for ( int_t num = -den; num <= den; ++num ) {
      const R exact( num, den );
      const auto rec = fromDbl( exact.toDouble(), 300 );
      REQUIRE( rec.has_value() && rec.value() == exact );
      //Also with a larger denominator bound:
      const auto rec2 = fromDbl( exact.toDouble(), 10000000 );
      REQUIRE( rec2.has_value() && rec2.value() == exact );
      ++nfrac;
    }
  }
  std::cout << "rationalFromDouble recovers all " << nfrac
            << " fractions a/b with |a|<=b<=300: OK" << std::endl;

  //Decimals (as parsed from text) are recovered as the exact decimal values,
  //all with 6 digits and a sample with 7 digits:
  std::size_t ndec = 0;
  for ( int_t m = 0; m < 1000000; ++m ) {
    const R exact( m, 1000000 );
    const auto rec = fromDbl( exact.toDouble(), 1000000 );
    REQUIRE( rec.has_value() && rec.value() == exact );
    ++ndec;
  }
  for ( int_t m = 1; m < 10000000; m += 997 ) {
    std::ostringstream ss;
    ss << "0." << std::string( 7 - std::to_string(m).size(), '0' ) << m;
    const std::string str = ss.str();
    const double x = NC::str2dbl( str );
    const auto rec = fromDbl( x, 10000000 );
    REQUIRE( rec.has_value() && rec.value() == R( m, 10000000 ) );
    ++ndec;
  }
  REQUIRE( fromDbl( 0.30625, 1000 ).value() == R(49,160) );
  REQUIRE( fromDbl( 0.3333333333333333, 1000 ).value() == R(1,3) );
  REQUIRE( fromDbl( -0.75, 4 ).value() == R(-3,4) );
  REQUIRE( !fromDbl( 1.0-0.9, 1000000 ).has_value() );//not 0.1
  std::cout << "rationalFromDouble recovers " << ndec << " decimals: OK"
            << std::endl;

  //The result has the smallest denominator of all rationals rounding to x
  //(checked by brute force over the denominators), and all values it returns
  //round to x:
  NC::RandXRSRImpl rng( 12345 );
  std::size_t nbrute = 0, nfound = 0;
  for ( int i = 0; i < 3000; ++i ) {
    double x;
    if ( i % 3 == 0 ) {
      x = rng.generate() * 4.0 - 2.0;//random double, no simple value
    } else {
      //a fraction with a random denominator (<=3000), perturbed or not:
      const int_t den = 1 + static_cast<int_t>( rng.generate() * 3000 );
      const int_t num = static_cast<int_t>( ( rng.generate() * 4.0 - 2.0 )
                                            * den );
      x = R( num, den ).toDouble();
      if ( i % 3 == 2 )
        x = std::nextafter( x, 10.0 );
    }
    const int_t maxden = 3000;
    const auto rec = fromDbl( x, maxden );
    int_t brute_den = 0;
    for ( int_t q = 1; q <= maxden && !brute_den; ++q ) {
      const int_t p0 = static_cast<int_t>( std::floor( x * q ) );
      for ( int_t p = p0 - 1; p <= p0 + 2; ++p )
        if ( double(p) / double(q) == x )
          brute_den = q;
    }
    if ( brute_den ) {
      REQUIRE( rec.has_value() && rec.value().denominator() == brute_den );
      REQUIRE( rec.value().toDouble() == x );
      ++nfound;
    } else {
      REQUIRE( !rec.has_value() );
    }
    ++nbrute;
  }
  std::cout << "rationalFromDouble agrees with brute force search on " << nbrute
            << " values (" << nfound << " with a rational): OK" << std::endl;

  //Edge cases:
  REQUIRE( fromDbl( 0.0, 1 ).value() == R() );
  REQUIRE( fromDbl( -0.0, 1 ).value() == R() );
  REQUIRE( fromDbl( 17.0, 1 ).value() == R(17) );
  REQUIRE( fromDbl( -3.0, 5 ).value() == R(-3) );
  REQUIRE( !fromDbl( 0.5, 1 ).has_value() );//denominator 2 > 1
  REQUIRE( fromDbl( 0.5, 2 ).value() == R(1,2) );
  const double nan = std::numeric_limits<double>::quiet_NaN();
  REQUIRE( !fromDbl( nan, 10 ).has_value() );
  REQUIRE( !fromDbl( NC::kInfinity, 10 ).has_value() );
  REQUIRE( !fromDbl( -NC::kInfinity, 10 ).has_value() );
  REQUIRE( !fromDbl( 1e300, 10 ).has_value() );
  REQUIRE( !fromDbl( 1e-300, 1000000000 ).has_value() );//(fast)
  REQUIRE( !fromDbl( 0.5, 0 ).has_value() );
  REQUIRE( fromDbl( 1e15+0.5, 2 ).value() == R(2000000000000001,2) );
  REQUIRE( fromDbl( 1e-9, 1000000000 ).value() == R(1,1000000000) );
  //With a large maxDenominator, a simpler rational than the intended one can
  //round to the same double (near 1 their spacing is below the ulp):
  REQUIRE( fromDbl( 0.999999999, 1000000000 ).value()
           == R(999999972,999999973) );
  std::cout << "rationalFromDouble edge cases: OK" << std::endl;

  //Agreement with NC::detectSimpleRationalNumbers (denominators up to 64,
  //values in [0,1] and integers) on all its fractions, on random values, and
  //on the doubles next to its fractions. On the latter, it is somewhat
  //lenient (its lookup key is x*1e18 truncated to an integer), and sometimes
  //detects a fraction although it does not round to x. Return values of
  //compare: 0 (agree), 1 (only detectSimpleRationalNumbers finds the
  //fraction), 2 (disagreement):
  auto compare = [&fromDbl]( double x )
  {
    const auto det = NC::detectSimpleRationalNumbers( x );
    const auto rec = fromDbl( x, 64 );
    if ( !rec.has_value() )
      return det.second == 0 ? 0 : 1;
    return ( rec.value().numerator() == static_cast<int_t>( det.first )
             && rec.value().denominator() == static_cast<int_t>( det.second )
             ? 0 : 2 );
  };
  std::size_t nagree = 0, nlenient = 0, nneighbours = 0;
  for ( int_t den = 1; den <= 64; ++den ) {
    for ( int_t num = 0; num <= den; ++num ) {
      const double x = R( num, den ).toDouble();
      REQUIRE( compare( x ) == 0 );
      ++nagree;
      const double neighbours[2] = { std::nextafter( x, 2.0 ),
                                     std::nextafter( x, -1.0 ) };
      for ( double y : neighbours ) {
        const int c = compare( y );
        REQUIRE( c != 2 );
        ++nneighbours;
        if ( c == 1 )
          ++nlenient;
      }
    }
  }
  for ( int i = 0; i < 100000; ++i ) {
    REQUIRE( compare( rng.generate() ) == 0 );
    ++nagree;
  }
  std::cout << "rationalFromDouble(x,64) agrees with"
            << " detectSimpleRationalNumbers on " << nagree
            << " fractions and random values: OK" << std::endl;
  std::cout << "On " << nneighbours << " doubles next to those fractions,"
            << " only detectSimpleRationalNumbers finds a fraction for "
            << nlenient << " (and they never disagree): OK" << std::endl;

  std::cout << "All OK" << std::endl;
  return 0;
}
