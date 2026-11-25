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
// Checks SCBragg and LCBragg instances sharing a model (set up in the crystal
// frame) for different crystal orientations.
////////////////////////////////////////////////////////////////////////////////

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

  std::vector<NC::SCOrientation> testOrientations()
  {
    using NC::CrystalAxis;
    using NC::HKLPoint;
    using NC::LabAxis;
    return {
      mkOrient( { CrystalAxis{0.,0.,1.}, LabAxis{0.,0.,1.} },
                { CrystalAxis{1.,0.,0.}, LabAxis{1.,0.,0.} } ),
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

  const NC::RotMatrix& identity()
  {
    static const double data[9] = { 1., 0., 0., 0., 1., 0., 0., 0., 1. };
    static const NC::RotMatrix m( data );
    return m;
  }

  NC::NeutronDirection rotate( const NC::RotMatrix& m,
                               const NC::NeutronDirection& d )
  {
    return ( m * d.as<NC::Vector>() ).as<NC::NeutronDirection>();
  }

  bool closeDirs( const NC::NeutronDirection& a,
                  const NC::NeutronDirection& b )
  {
    return ( a.as<NC::Vector>() - b.as<NC::Vector>() ).mag2() < 1e-18;
  }

  //Compares a process in the crystal frame (crys) with one rotated by
  //crystal2lab (rot), for neutrons rotated accordingly. If same is true,
  //rot2 must give identical results to rot. Returns number of non-vanishing
  //cross sections encountered.
  unsigned compare( const ProcPtr& crys, const ProcPtr& rot,
                    const ProcPtr& rot2, const NC::RotMatrix& crystal2lab,
                    unsigned ndirs )
  {
    REQUIRE( rot->isOriented() && crys->isOriented() );
    REQUIRE( crys->domain().elow == rot->domain().elow );
    REQUIRE( crys->domain().ehigh == rot->domain().ehigh );
    auto rng = NC::createBuiltinRNG( 117 );
    NC::CachePtr cp_crys, cp_rot, cp_rot2;
    unsigned nxs(0);
    for ( auto wl : { 0.4, 0.9, 1.7, 2.6, 4.1, 6.5 } ) {
      const NC::NeutronEnergy ekin{ NC::wl2ekin(wl) };
      for ( unsigned i = 0; i < ndirs; ++i ) {
        const auto dir = NC::randIsotropicNeutronDirection(rng);
        const auto dir_lab = rotate( crystal2lab, dir );
        const double xs_crys = crys->crossSection( cp_crys, ekin, dir ).dbl();
        const double xs_rot = rot->crossSection( cp_rot, ekin, dir_lab ).dbl();
        const double xs_rot2 = rot2->crossSection( cp_rot2, ekin,
                                                   dir_lab ).dbl();
        REQUIRE( NC::floateq( xs_crys, xs_rot, 1e-9, 1e-14 ) );
        REQUIRE( xs_rot == xs_rot2 );
        if ( !xs_crys )
          continue;
        ++nxs;
        //Scatterings, with identical random numbers:
        auto rng_crys = NC::createBuiltinRNG( 1000 + i );
        auto rng_rot = NC::createBuiltinRNG( 1000 + i );
        auto rng_rot2 = NC::createBuiltinRNG( 1000 + i );
        for ( unsigned j = 0; j < 3; ++j ) {
          auto o_crys = crys->sampleScatter( cp_crys, rng_crys, ekin, dir );
          auto o_rot = rot->sampleScatter( cp_rot, rng_rot, ekin, dir_lab );
          auto o_rot2 = rot2->sampleScatter( cp_rot2, rng_rot2, ekin, dir_lab );
          REQUIRE( o_crys.ekin == o_rot.ekin );
          REQUIRE( closeDirs( rotate( crystal2lab, o_crys.direction ),
                              o_rot.direction ) );
          REQUIRE( o_rot.direction == o_rot2.direction );
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
    std::shared_ptr<const NC::SCBragg::Model> sc_model;
    std::shared_ptr<const NC::LCBragg::Model> lc_model;
    if ( layered )
      lc_model = NC::LCBragg::createModel( info, mos, lcaxis, lcmode );
    else
      sc_model = NC::SCBragg::createModel( info, mos );
    auto mk = [&]( const NC::RotMatrix& m ) -> ProcPtr
    {
      if ( layered ) {
        auto p = NC::makeSO<NC::LCBragg>( lc_model, m );
        REQUIRE( p->model().get() == lc_model.get() );
        return p;
      }
      auto p = NC::makeSO<NC::SCBragg>( sc_model, m );
      REQUIRE( p->model().get() == sc_model.get() );
      return p;
    };
    const ProcPtr crys = mk( identity() );
    for ( auto& sco : testOrientations() ) {
      const auto crystal2lab = NC::getCrystal2LabRot( sco, reci );
      auto rot = mk( crystal2lab );
      REQUIRE( NC::ProcImpl::NewABI::isPureElasticScatter( rot ) );
      //Without sharing, via the constructor taking an orientation:
      ProcPtr rot2 = ( layered
                       ? ProcPtr( NC::makeSO<NC::LCBragg>( info, sco, mos,
                                                           lcaxis, lcmode ) )
                       : ProcPtr( NC::makeSO<NC::SCBragg>( info, sco, mos ) ) );
      const unsigned nxs = compare( crys, rot, rot2, crystal2lab,
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
      //Names, and descriptions with a shared model id:
      REQUIRE( std::string( rot->name() ) == ( layered ? "LCBragg" : "SCBragg" ) );
      REQUIRE( NC::contains( rot->jsonDescription(), "\"model_uid\":" ) );
      REQUIRE( NC::contains( rot->jsonDescription(), "\"crystal2lab\":" ) );
    }
  }
}

int main()
{
  testCrystal( "Ge_sg227.ncmat;dcutoff=0.7", false );
  testCrystal( "Al2O3_sg167_Corundum.ncmat;dcutoff=0.9", false );
  testCrystal( "C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.6", true );
  testCrystal( "C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.9", true, 20 );
  return 0;
}
