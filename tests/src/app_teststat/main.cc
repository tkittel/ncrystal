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

#include "NCTestUtils/NCTestStat.hh"
#include <iostream>

//Tests of the Kolmogorov-Smirnov utilities in NCTestStat.hh, with reference
//values from the Python version in NCTestUtils/stat.py.

#define REQUIRE(x) nc_assert_always(x)

namespace NC = NCrystal;

namespace {
  NC::VectD fracSeq( double step, unsigned n, double power = 1.0 )
  {
    NC::VectD v;
    for ( unsigned i = 1; i <= n; ++i ) {
      const double x = i * step;
      v.push_back( std::pow( x - std::floor( x ), power ) );
    }
    return v;
  }
  void check( const char * name, const NC::VectD& x, const NC::VectD& y,
              double ref_d, double ref_p )
  {
    const double d = nctest::ksDistance( x, y );
    const double p = nctest::ksPValue( x, y );
    std::cout << name << ": D=" << NC::fmt( d, "%.6g" )
              << " p=" << NC::fmt( p, "%.6g" ) << std::endl;
    REQUIRE( NC::floateq( d, ref_d, 1e-12, 1e-15 ) );
    REQUIRE( NC::floateq( p, ref_p, 1e-9, 1e-15 ) );
  }
}

int main()
{
  const NC::VectD a = fracSeq( 0.6180339887498949, 1000 );
  check( "identical samples", a, a, 0.0, 1.0 );
  check( "disjoint samples", { 1.0, 2.0, 3.0 }, { 4.0, 5.0 }, 1.0,
         nctest::ksPValue( { 1.0, 2.0, 3.0 }, { 4.0, 5.0 } ) );
  check( "same distribution", a, fracSeq( 0.4142135623730951, 800 ),
         0.00275000000000003, 0.9999999999809577 );
  check( "different distributions", a,
         fracSeq( 0.4142135623730951, 800, 1.2 ),
         0.0685, 0.03087694814896776 );
  NC::VectD t7, t5;
  for ( int i = 0; i < 50; ++i )
    t7.push_back( i % 7 );
  for ( int i = 0; i < 40; ++i )
    t5.push_back( i % 5 );
  check( "with ties", t7, t5, 0.28, 0.061339816229094256 );
  return 0;
}
