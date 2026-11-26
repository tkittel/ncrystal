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

////////////////////////////////////////////////////////////////////////////////
// Checks the precision of GaussOnSphere::circleIntegral and the methods it
// uses (table and numerical integration), against independent
// high precision reference values from devel/misc/gen_mosaicrefs.py. Also
// checks which table (if any) is set up for each configuration.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/phys_utils/NCGaussOnSphere.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCTestUtils/NCTestFindData.hh"
#include <fstream>
#include <iostream>
#include <map>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  struct RefCase { double mosdeg, ntrunc, prec, alpha, gamma, ref; };

  std::vector<RefCase> loadRefs()
  {
    std::vector<RefCase> res;
    for ( auto& v : nctest::load_test_data_numbers( "mosaicrefs",
                                                    "circleint.txt", 6 ) )
      res.push_back( { v[0], v[1], v[2], v[3], v[4], v[5] } );
    REQUIRE( res.size() >= 200 );
    return res;
  }

  //Expected table sizes (alpha points) for the configurations in the
  //reference file, by mosaicity and precision (zero if no table, e.g. since
  //the closed form is valid in places, the table range is empty, or the
  //density underflows within the truncation):
  unsigned expectedTableSize( double mosdeg, double prec )
  {
    using Key = std::pair<double,double>;
    static const std::map<Key,unsigned> expected = {
      { { 1.0, 1e-7 }, 128 }, { { 10.0, 1e-4 }, 64 }, { { 20.0, 1e-4 }, 64 },
      { { 20.0, 1e-7 }, 128 }, { { 40.0, 1e-3 }, 64 } };
    auto it = expected.find( { mosdeg, prec } );
    return it == expected.end() ? 0 : it->second;
  }
}

int main()
{
  std::map<std::tuple<double,double,double>,
           std::unique_ptr<NC::GaussOnSphere>> goss;
  unsigned nnumint(0), ntable(0), nother(0);
  for ( auto& rc : loadRefs() ) {
    auto& gosptr = goss[ std::make_tuple( rc.mosdeg, rc.ntrunc, rc.prec ) ];
    const double sigma = NC::MosaicityFWHM{ rc.mosdeg * NC::kDeg }.sigma().dbl();
    if ( !gosptr ) {
      gosptr = std::make_unique<NC::GaussOnSphere>( sigma, rc.ntrunc * sigma,
                                                    rc.prec );
      const unsigned na = expectedTableSize( rc.mosdeg, rc.prec );
      REQUIRE( gosptr->hasTable() == ( na > 0 ) );
      REQUIRE( gosptr->tableSize().first == na );
      REQUIRE( gosptr->tableSize().second == 4 * na );
    }
    const auto& gos = *gosptr;
    const double cg( std::cos(rc.gamma) ), sg( std::sin(rc.gamma) );
    const double ca( std::cos(rc.alpha) ), sa( std::sin(rc.alpha) );
    //Errors relative to the reference, or to values which are not negligible
    //(as in GaussOnSphere::setupTable):
    const double wfloor = 1e-3 * gos.getNormFactor() * sigma;
    auto relerr = [&rc,wfloor]( double val )
    {
      return NC::ncabs( val - rc.ref ) / NC::ncmax( rc.ref, wfloor );
    };
    //Numerical integration and table aim at 10 times below the accuracy:
    const double acc = NC::ncclamp( rc.prec * 0.1, 1e-7, 1e-4 );
    const double num = gos.circleIntegralNumerical( cg, sg, ca, sa );
    REQUIRE( relerr( num ) <= 0.1 * acc );
    const bool table_range = ( gos.hasTable()
                               && rc.alpha >= 1.2 * gos.getTruncangle() );
    if ( table_range )
      REQUIRE( relerr( gos.circleIntegralTable( cg, sg, ca, sa ) ) <= acc );
    //The final result, where not from the closed form:
    const double val = gos.circleIntegral( cg, sg, ca, sa );
    if ( val == num && !table_range ) {
      ++nnumint;
    } else if ( table_range ) {
      ++ntable;
      REQUIRE( relerr( val ) <= acc );
    } else {
      ++nother;
    }
  }
  std::cout << "Cases using numerical integration: " << nnumint << ", table: "
            << ntable << ", other: " << nother << std::endl;
  REQUIRE( nnumint >= 20 && ntable >= 20 && nother >= 20 );
  return 0;
}
