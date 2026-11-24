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
// Checks that a RotatedScatter wrapping a process set up in the crystal frame
// is equivalent to the same process set up directly in the lab frame.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/extd_utils/NCRotatedScatter.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/internal/extd_utils/NCABIUtils.hh"
#include "NCrystal/internal/scbragg/NCSCBragg.hh"
#include "NCrystal/internal/lcbragg/NCLCBragg.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/factories/NCFactImpl.hh"

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  using ProcPtr = NC::ProcImpl::ProcPtr;

  NC::SCOrientation mkOrient( const NC::OrientDir& dir1,
                              const NC::OrientDir& dir2,
                              double dirtol = 1e-4 )
  {
    NC::SCOrientation sco;
    sco.setPrimaryDirection( dir1 );
    sco.setSecondaryDirection( dir2, dirtol );
    return sco;
  }

  NC::SCOrientation crystalFrame()
  {
    return mkOrient( { NC::CrystalAxis{0.,0.,1.}, NC::LabAxis{0.,0.,1.} },
                     { NC::CrystalAxis{1.,0.,0.}, NC::LabAxis{1.,0.,0.} } );
  }

  std::vector<NC::SCOrientation> testOrientations()
  {
    using NC::CrystalAxis;
    using NC::HKLPoint;
    using NC::LabAxis;
    return {
      crystalFrame(),
      mkOrient( { CrystalAxis{0.,0.,1.}, LabAxis{1.,0.,0.} },
                { CrystalAxis{1.,0.,0.}, LabAxis{0.,1.,0.} } ),
      mkOrient( { HKLPoint{5.,1.,1.}, LabAxis{0.,0.,1.} },
                { HKLPoint{0.,-1.,1.}, LabAxis{0.,1.,0.} }, NC::kPi ),
      mkOrient( { HKLPoint{1.,1.,0.}, LabAxis{0.3,-1.2,0.7} },
                { CrystalAxis{0.2,0.1,1.}, LabAxis{1.,1.,1.} }, NC::kPi ),
      mkOrient( { CrystalAxis{-0.3,0.9,0.1}, LabAxis{-1.,0.2,0.1} },
                { HKLPoint{1.,0.,2.}, LabAxis{0.,0.,-1.} }, NC::kPi ),
    };
  }

  bool closeDirs( const NC::NeutronDirection& a,
                  const NC::NeutronDirection& b )
  {
    return ( a.as<NC::Vector>() - b.as<NC::Vector>() ).mag2() < 1e-18;
  }

  //Compares a process with a given orientation (ref) with a rotated one
  //(rot). Returns number of non-vanishing cross sections encountered.
  unsigned compare( const ProcPtr& ref, const ProcPtr& rot, bool oriented,
                    unsigned ndirs )
  {
    REQUIRE( rot->isOriented() );
    REQUIRE( ref->isOriented() == oriented );
    REQUIRE( ref->domain().elow == rot->domain().elow );
    REQUIRE( ref->domain().ehigh == rot->domain().ehigh );
    auto rng = NC::createBuiltinRNG( 117 );
    NC::CachePtr cp_ref, cp_rot;
    unsigned nxs(0);
    for ( auto wl : { 0.4, 0.9, 1.7, 2.6, 4.1, 6.5 } ) {
      const NC::NeutronEnergy ekin{ NC::wl2ekin(wl) };
      for ( unsigned i = 0; i < ndirs; ++i ) {
        const auto dir = NC::randIsotropicNeutronDirection(rng);
        const double xs_ref = ref->crossSection( cp_ref, ekin, dir ).dbl();
        const double xs_rot = rot->crossSection( cp_rot, ekin, dir ).dbl();
        REQUIRE( NC::floateq( xs_ref, xs_rot, 1e-9, 1e-14 ) );
        if ( !xs_ref )
          continue;
        ++nxs;
        //Scatterings, with identical random numbers:
        auto rng_ref = NC::createBuiltinRNG( 1000 + i );
        auto rng_rot = NC::createBuiltinRNG( 1000 + i );
        for ( unsigned j = 0; j < 3; ++j ) {
          auto o_ref = ref->sampleScatter( cp_ref, rng_ref, ekin, dir );
          auto o_rot = rot->sampleScatter( cp_rot, rng_rot, ekin, dir );
          REQUIRE( NC::floateq( o_ref.ekin.dbl(), o_rot.ekin.dbl(),
                                1e-12, 1e-99 ) );
          if ( oriented ) {
            REQUIRE( closeDirs( o_ref.direction, o_rot.direction ) );
          } else {
            //Only the scattering angle is comparable:
            auto mu = []( const NC::NeutronDirection& a,
                          const NC::NeutronDirection& b )
            {
              return a.as<NC::Vector>().dot( b.as<NC::Vector>() );
            };
            REQUIRE( NC::floateq( mu( dir, o_ref.direction ),
                                  mu( dir, o_rot.direction ),
                                  1e-9, 1e-9 ) );
          }
        }
      }
    }
    return nxs;
  }

  void testCrystal( const char * cfgstr, bool layered, int lcmode = 0 )
  {
    auto info = NC::FactImpl::createInfo( NC::MatCfg(cfgstr) );
    const auto reci = NC::getReciprocalLatticeRot( info->getStructureInfo() );
    const NC::MosaicityFWHM mos{ 1.5 * NC::kDeg };
    const NC::LCAxis lcaxis{ 0., 0., 1. };
    auto mk = [&]( const NC::SCOrientation& sco ) -> ProcPtr
    {
      if ( layered )
        return NC::makeSO<NC::LCBragg>( info, sco, mos, lcaxis, lcmode );
      return NC::makeSO<NC::SCBragg>( info, sco, mos );
    };
    const ProcPtr crystal = mk( crystalFrame() );
    for ( auto& sco : testOrientations() ) {
      auto rot = NC::makeSO<NC::RotatedScatter>
        ( crystal, NC::getCrystal2LabRot( sco, reci ) );
      REQUIRE( rot->process() == crystal );
      REQUIRE( NC::ProcImpl::NewABI::isPureElasticScatter( rot ) );
      const unsigned nxs = compare( mk( sco ), rot, true,
                                    lcmode ? 20 : 500 );
      REQUIRE( nxs >= ( layered ? 50 : 15 ) );
      //Unscattered neutrons must be exactly unchanged:
      NC::CachePtr cp;
      auto rng = NC::createBuiltinRNG( 5 );
      const NC::NeutronDirection dir{ 0.1, -0.5, 0.7 };
      const NC::NeutronEnergy ekin{ NC::wl2ekin(50.0) };
      REQUIRE( rot->crossSection( cp, ekin, dir ).dbl() == 0.0 );
      auto outcome = rot->sampleScatter( cp, rng, ekin, dir );
      REQUIRE( outcome.ekin == ekin );
      REQUIRE( outcome.direction == dir );
    }
    //Description must expose the wrapped process:
    auto descr = NC::makeSO<NC::RotatedScatter>
      ( crystal, NC::getCrystal2LabRot( crystalFrame(), reci ) )
      ->jsonDescription();
    REQUIRE( NC::contains( descr, "\"name\":\"RotatedScatter\"" ) );
    REQUIRE( NC::contains( descr, layered ? "\"name\":\"LCBragg\""
                                          : "\"name\":\"SCBragg\"" ) );
  }

  void testIsotropic()
  {
    //Rotating an isotropic (and inelastic) process must be harmless.
    auto proc = NC::FactImpl::createScatter( "Al_sg225.ncmat;comp=inelas" );
    auto info = NC::FactImpl::createInfo( NC::MatCfg("Al_sg225.ncmat") );
    const auto reci = NC::getReciprocalLatticeRot( info->getStructureInfo() );
    REQUIRE( !NC::ProcImpl::NewABI::isPureElasticScatter( proc ) );
    for ( auto& sco : testOrientations() ) {
      auto rot = NC::makeSO<NC::RotatedScatter>
        ( proc, NC::getCrystal2LabRot( sco, reci ) );
      REQUIRE( !NC::ProcImpl::NewABI::isPureElasticScatter( rot ) );
      REQUIRE( compare( proc, rot, false, 10 ) == 60 );
    }
  }
}

int main()
{
  testCrystal( "Ge_sg227.ncmat;dcutoff=0.7", false );
  testCrystal( "Al2O3_sg167_Corundum.ncmat;dcutoff=0.9", false );
  testCrystal( "C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.6", true );
  testCrystal( "C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.9", true, 20 );
  testIsotropic();
  return 0;
}
