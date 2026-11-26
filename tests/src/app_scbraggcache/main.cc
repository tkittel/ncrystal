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

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include <iostream>

// Tests of the SCBragg cache: Results only depend on the (snapped) input
// direction, never on the call history, and an unscattered neutron keeps its
// exact input direction.

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  //Narrow mosaicity, so cross sections vary rapidly with the direction:
  const char * cfgstr = "Ge_sg227.ncmat;dcutoff=1;comp=bragg;mos=0.01deg"
                        ";dir1=@crys_hkl:0,0,1@lab:0,0,1"
                        ";dir2=@crys_hkl:1,0,0@lab:1,0,0";

  void testHistoryIndependence()
  {
    auto sc = NC::FactImpl::createScatter( NC::MatCfg( cfgstr ) );
    auto rng = NC::getRNG();
    //Find a neutron state on the flank of a reflection:
    const NC::NeutronEnergy ekin{ NC::NeutronWavelength{ 2.0 } };
    NC::Vector d0;
    double xs0 = 0.0;
    for ( unsigned i = 0; i < 1000000 && !( xs0 > 0.0 ); ++i ) {
      d0 = NC::randIsotropicDirection( *rng );
      NC::CachePtr cp;
      xs0 = sc->crossSection( cp, ekin, d0.as<NC::NeutronDirection>() ).dbl();
    }
    REQUIRE( xs0 > 0.0 );
    //Steps of directions across many grid cells (2^-47~7e-15) of the cache:
    const NC::Vector perp = d0.cross( NC::Vector( 0.3, -0.2, 1.0 ) ).unit();
    NC::CachePtr cp_shared;
    unsigned ndiffer = 0;
    double xs_prev = -1.0;
    for ( int i = -200; i <= 200; ++i ) {
      const NC::Vector d = d0 + perp * ( i * 1.3e-15 );
      const auto nd = d.as<NC::NeutronDirection>();
      NC::CachePtr cp_fresh;
      const double xs_fresh = sc->crossSection( cp_fresh, ekin, nd ).dbl();
      const double xs_shared = sc->crossSection( cp_shared, ekin, nd ).dbl();
      REQUIRE( xs_fresh == xs_shared );
      //Also when coming back from a direction just next to it:
      const auto nd2 = ( d + perp * 1e-15 ).as<NC::NeutronDirection>();
      (void)sc->crossSection( cp_shared, ekin, nd2 );
      REQUIRE( sc->crossSection( cp_shared, ekin, nd ).dbl() == xs_fresh );
      if ( xs_fresh != xs_prev )
        ++ndiffer;
      xs_prev = xs_fresh;
    }
    //The cache is recalculated when crossing grid cells, giving (slightly)
    //different results:
    std::cout << "Number of distinct consecutive values: " << ndiffer
              << " (out of 401)" << std::endl;
    REQUIRE( ndiffer > 50 );
  }

  void testUnscatteredDirection()
  {
    //Below the Bragg threshold or without contributing planes, the outcome
    //is exactly the input direction:
    auto sc = NC::FactImpl::createScatter( NC::MatCfg( cfgstr ) );
    auto rng = NC::getRNG();
    NC::CachePtr cp;
    const NC::NeutronDirection indir{ 0.1, 0.2, 0.3 };//not normalised
    for ( double wl : { 1.0, 2.0, 20.0 } ) {
      const NC::NeutronEnergy ekin{ NC::NeutronWavelength{ wl } };
      if ( sc->crossSection( cp, ekin, indir ).dbl() > 0.0 )
        continue;
      const auto out = sc->sampleScatter( cp, *rng, ekin, indir );
      REQUIRE( out.direction == indir );
      REQUIRE( out.ekin == ekin );
    }
  }
}

int main()
{
  testHistoryIndependence();
  testUnscatteredDirection();
  return 0;
}
