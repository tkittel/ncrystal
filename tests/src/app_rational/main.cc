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


#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/core/NCFmt.hh"
#include <iostream>
#include <sstream>

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    //Private to the sgsym component (NCSGExactOrbit.hh):
    Optional<std::pair<std::int64_t,std::int64_t>>
    rationalFromDouble( double x, std::int64_t maxDenominator );
  }
}

namespace NC = NCrystal;
using int_t = std::int64_t;
using Frac = std::pair<int_t,int_t>;//numerator, denominator

#define REQUIRE(x) nc_assert_always(x)

namespace {

  //Reference: n/d in lowest terms (d>0), with plain integer arithmetic:
  Frac F( int_t n, int_t d = 1 )
  {
    nc_assert_always( d != 0 );
    if ( d < 0 ) {
      n = -n;
      d = -d;
    }
    int_t a = n < 0 ? -n : n, b = d;
    while ( b ) {
      const int_t t = a % b;
      a = b;
      b = t;
    }
    return { n / a, d / a };
  }

  double val( const Frac& f )
  {
    return static_cast<double>( f.first ) / static_cast<double>( f.second );
  }

}

int main()
{
  auto fromDbl = []( double x, int_t maxden )
  {
    return NC::detail::rationalFromDouble( x, maxden );
  };

  //All fractions a/b with |a|<=b<=300 are recovered from their doubles:
  std::size_t nfrac = 0;
  for ( int_t den = 1; den <= 300; ++den ) {
    for ( int_t num = -den; num <= den; ++num ) {
      const Frac exact = F( num, den );
      const auto rec = fromDbl( val( exact ), 300 );
      REQUIRE( rec.has_value() && rec.value() == exact );
      //Also with a larger denominator bound:
      const auto rec2 = fromDbl( val( exact ), 10000000 );
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
    const Frac exact = F( m, 1000000 );
    const auto rec = fromDbl( val( exact ), 1000000 );
    REQUIRE( rec.has_value() && rec.value() == exact );
    ++ndec;
  }
  for ( int_t m = 1; m < 10000000; m += 997 ) {
    std::ostringstream ss;
    ss << "0." << std::string( 7 - std::to_string(m).size(), '0' ) << m;
    const std::string str = ss.str();
    const double x = NC::str2dbl( str );
    const auto rec = fromDbl( x, 10000000 );
    REQUIRE( rec.has_value() && rec.value() == F( m, 10000000 ) );
    ++ndec;
  }
  REQUIRE( fromDbl( 0.30625, 1000 ).value() == F(49,160) );
  REQUIRE( fromDbl( 0.3333333333333333, 1000 ).value() == F(1,3) );
  REQUIRE( fromDbl( -0.75, 4 ).value() == F(-3,4) );
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
      x = val( F( num, den ) );
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
      REQUIRE( rec.has_value() && rec.value().second == brute_den );
      REQUIRE( val( rec.value() ) == x );
      ++nfound;
    } else {
      REQUIRE( !rec.has_value() );
    }
    ++nbrute;
  }
  std::cout << "rationalFromDouble agrees with brute force search on " << nbrute
            << " values (" << nfound << " with a rational): OK" << std::endl;

  //Edge cases:
  REQUIRE( fromDbl( 0.0, 1 ).value() == F( 0 ) );
  REQUIRE( fromDbl( -0.0, 1 ).value() == F( 0 ) );
  REQUIRE( fromDbl( 17.0, 1 ).value() == F(17) );
  REQUIRE( fromDbl( -3.0, 5 ).value() == F(-3) );
  REQUIRE( !fromDbl( 0.5, 1 ).has_value() );//denominator 2 > 1
  REQUIRE( fromDbl( 0.5, 2 ).value() == F(1,2) );
  const double nan = std::numeric_limits<double>::quiet_NaN();
  REQUIRE( !fromDbl( nan, 10 ).has_value() );
  REQUIRE( !fromDbl( NC::kInfinity, 10 ).has_value() );
  REQUIRE( !fromDbl( -NC::kInfinity, 10 ).has_value() );
  REQUIRE( !fromDbl( 1e300, 10 ).has_value() );
  REQUIRE( !fromDbl( 1e-300, 1000000000 ).has_value() );//(fast)
  REQUIRE( !fromDbl( 0.5, 0 ).has_value() );
  REQUIRE( fromDbl( 1e15+0.5, 2 ).value() == F(2000000000000001,2) );
  REQUIRE( fromDbl( 1e-9, 1000000000 ).value() == F(1,1000000000) );
  //With a large maxDenominator, a simpler rational than the intended one can
  //round to the same double (near 1 their spacing is below the ulp):
  REQUIRE( fromDbl( 0.999999999, 1000000000 ).value()
           == F(999999972,999999973) );
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
    return ( rec.value().first == static_cast<int_t>( det.first )
             && rec.value().second == static_cast<int_t>( det.second )
             ? 0 : 2 );
  };
  std::size_t nagree = 0, nlenient = 0, nneighbours = 0;
  for ( int_t den = 1; den <= 64; ++den ) {
    for ( int_t num = 0; num <= den; ++num ) {
      const double x = val( F( num, den ) );
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
