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


#include "NCrystal/internal/sabscatter/NCSABScatter.hh"
#include "NCrystal/internal/sab/NCSABExtended.hh"
#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCS = NCrystal::SAB;

namespace {

  std::vector<const NC::SABScatter*> findSAB( const NC::ProcImpl::Process& p )
  {
    std::vector<const NC::SABScatter*> res;
    if ( auto s = dynamic_cast<const NC::SABScatter*>( &p ) )
      res.push_back( s );
    if ( auto c = dynamic_cast<const NC::ProcImpl::ProcComposition*>( &p ) )
      for ( auto& comp : c->components() )
        if ( auto s = dynamic_cast<const NC::SABScatter*>( &*comp.process ) )
          res.push_back( s );
    return res;
  }

  void testSAB( const NC::SABScatter* s )
  {
    //The tabulated SAB cross section against the direct evaluation it
    //replaces: SABProcessor::crossSectionUnitSigmaBound up to Emax, and
    //k/E+extender above (see NCSABExtended.cc). Also check that the scaled
    //description agrees with the cross section of the process:
    const auto& ext = s->sabExtended();
    const auto& p = ext.processor();
    const auto emaxInfo = p.getEMaxInfo();
    const double emax = emaxInfo.ekin.dbl();
    const double k = emax * ( emaxInfo.crossSectionUnitSigmaBound.dbl()
                              - ext.extender().crossSection( emaxInfo.ekin ).dbl() );
    auto direct = [&]( double e )
    {
      if ( e <= emax )
        return p.crossSectionUnitSigmaBound( NC::NeutronEnergy{ e } ).dbl();
      return k / e + ext.extender().crossSection( NC::NeutronEnergy{ e } ).dbl();
    };
    const auto d = s->createIsoXSTable();
    nc_assert_always( d.has_value() );
    const auto& t = d.value();
    const auto bp = s->createIsoXSBreakpoints();
    nc_assert_always( bp.size() == 2 && bp.back() == emax );
    NC::CachePtr cp0;
    const double scale = s->crossSectionIsotropic( cp0, NC::NeutronEnergy{ 1.0 } ).dbl()
      / ext.crossSectionUnitSigmaBound( NC::NeutronEnergy{ 1.0 } ).dbl();
    NC::VectD energies = NC::geomspace( 1e-12, 1e3, 200001 );
    for ( auto e : t.E ) {
      energies.push_back( e );
      energies.push_back( std::nextafter( e, 0.0 ) );
      energies.push_back( std::nextafter( e, NC::kInfinity ) );
    }
    double maxdev_grid = 0.0, maxdev_ext = 0.0;
    NC::CachePtr cp;
    for ( auto e : energies ) {
      const double xs = s->crossSectionIsotropic( cp, NC::NeutronEnergy{ e } ).dbl();
      nc_assert_always( NC::floateq( xs, NC::evalIsoXS( t, e ), 1e-14, 0.0 ) );
      const double xsu = ext.crossSectionUnitSigmaBound( NC::NeutronEnergy{ e } ).dbl();
      const double ref = direct( e );
      const double dev = std::fabs( xsu / ref - 1.0 );
      ( e <= emax ? maxdev_grid : maxdev_ext ) = std::max( e <= emax ? maxdev_grid : maxdev_ext, dev );
      nc_assert_always( NC::floateq( xs, scale * xsu, 1e-14, 0.0 ) );
    }
    std::cout << "  " << ( maxdev_grid < 1e-13 ? "up to Emax within 1e-13" : "up to Emax NOT within 1e-13" )
              << ", " << ( maxdev_ext < 1e-9 ? "above within 1e-9" : "above NOT within 1e-9" )
              << " (SAB element with " << p.getEDivKTGrid().size() << " egrid points)" << std::endl;
    nc_assert_always( maxdev_grid < 1e-13 && maxdev_ext < 1e-9 );
  }

  void testMaterial( const char * cfgstr )
  {
    auto sc = NC::FactImpl::createScatter( NC::MatCfg( cfgstr ) );
    auto sabs = findSAB( *sc );
    nc_assert_always( !sabs.empty() );
    std::cout << cfgstr << ":" << std::endl;
    for ( auto s : sabs )
      testSAB( s );
  }

  void testFreeGasExtender( const char * cfgstr )
  {
    //Like testMaterial, but with the free-gas extender instead of the default
    //SCT extender, set up like the factories do it:
    auto info = NC::FactImpl::createInfo( NC::MatCfg( cfgstr ) );
    std::cout << cfgstr << " (free-gas extender):" << std::endl;
    for ( auto& di : info->getDynamicInfoList() ) {
      auto di_knl = dynamic_cast<const NC::DI_ScatKnl*>( di.get() );
      nc_assert_always( di_knl != nullptr );
      auto sabdata = NC::DynInfoUtils::extractSABDataFromDynInfo( di_knl );
      auto ext = NCS::SABExtended::createWithFGExtender( NC::SABCfg::createConfig( 3 ),
                                                         sabdata,
                                                         di_knl->energyGrid() );
      nc_assert_always( ext->extender().shortName().value_or("") == "freegas" );
      NC::SABScatter s( std::move( ext ), sabdata->boundXS() );
      testSAB( &s );
    }
  }

  void testNoXSSupport()
  {
    //SupportXS::NO: same samplings as with cross section support, while the
    //cross section table is unavailable:
    auto info = NC::FactImpl::createInfo( NC::MatCfg( "Al_sg225.ncmat" ) );
    auto di_knl = dynamic_cast<const NC::DI_ScatKnl*>
      ( info->getDynamicInfoList().front().get() );
    nc_assert_always( di_knl != nullptr );
    auto sabdata = NC::DynInfoUtils::extractSABDataFromDynInfo( di_knl );
    const auto teff = NC::DynInfoUtils::extractTeffFromDynInfo( di_knl );
    using SupportXS = NCS::SABExtended::SupportXS;
    auto mk = [&sabdata,&teff]( SupportXS sx )
    {
      return NCS::createSABExtendedWithCache( 3, sabdata, teff, nullptr, sx );
    };
    auto ext_yes = mk( SupportXS::YES );
    auto ext_no = mk( SupportXS::NO );
    nc_assert_always( ext_yes->supportsXS() && !ext_no->supportsXS() );
    nc_assert_always( &*ext_no != &*ext_yes );//different cache entries
    nc_assert_always( &*ext_no == &*mk( SupportXS::NO ) );//cache hit
    auto expectBadInput = []( std::function<void()> f )
    {
      bool threw = false;
      try { f(); } catch ( NC::Error::BadInput& ) { threw = true; }
      nc_assert_always( threw );
    };
    expectBadInput( [&ext_no](){ ext_no->xsEvaluator(); } );
    expectBadInput( [&ext_no](){ ext_no->xsBreakpoints(); } );
    expectBadInput( [&ext_no](){ NC::SABScatter( ext_no, 1.0 ); } );
    for ( double e : { 1e-4, 0.025, 0.5, 20.0 } ) {
      auto rng1 = NC::createBuiltinRNG( 1234 );
      auto rng2 = NC::createBuiltinRNG( 1234 );
      for ( int i = 0; i < 1000; ++i ) {
        auto o1 = ext_yes->sampleScatter( *rng1, NC::NeutronEnergy{ e } );
        auto o2 = ext_no->sampleScatter( *rng2, NC::NeutronEnergy{ e } );
        nc_assert_always( o1.ekin.dbl() == o2.ekin.dbl()
                          && o1.mu.dbl() == o2.mu.dbl() );
      }
    }
    std::cout << "SupportXS::NO: identical samplings, no cross sections"
              << std::endl;
  }
}

int main()
{
  testMaterial( "Al_sg225.ncmat;knllux=3;comp=inelas" );
  testMaterial( "Polyethylene_CH2.ncmat;knllux=3;comp=inelas" );
  testMaterial( "Al_sg225.ncmat;temp=20K;knllux=3;comp=inelas" );
  testMaterial( "Al_sg225.ncmat;temp=3000K;knllux=3;comp=inelas" );
  testMaterial( "LiquidHeavyWaterD2O_T293.6K.ncmat;knllux=3;comp=inelas" );
  testFreeGasExtender( "Al_sg225.ncmat" );
  testNoXSSupport();
  return 0;
}
