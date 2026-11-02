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
  //(R::kMaxValue rather than k in the lambdas, for MSVC):
  requireThrows( "sum too large", [](){ R(R::kMaxValue) + R(1); } );
  requireThrows( "product too large", [](){ R(R::kMaxValue/2) * R(3); } );
  requireThrows( "denominators too large",
                 [](){ R(1,R::kMaxValue) + R(1,R::kMaxValue-1); } );
  REQUIRE( R(k/2) * R(2) == R(k) );//(exactly at the limit)

  std::cout << "All OK" << std::endl;
  return 0;
}
