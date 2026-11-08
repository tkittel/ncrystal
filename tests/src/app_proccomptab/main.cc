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


#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCPI = NCrystal::ProcImpl;

namespace {

  //Elastic process without a table (not tabulated, unlike all built-in
  //processes of standard materials with the default settings):
  class UntabulatedScat final : public NCPI::ScatterIsotropicMat {
  public:
    const char* name() const noexcept final { return "UntabulatedScat"; }
    NC::CrossSect crossSectionIsotropic( NC::CachePtr&,
                                         NC::NeutronEnergy ekin ) const final
    {
      return NC::CrossSect{ 1.0 + 0.1 * std::sqrt( ekin.dbl() ) };
    }
    NC::ScatterOutcomeIsotropic sampleScatterIsotropic( NC::CachePtr&,
                                                        NC::RNG&,
                                                        NC::NeutronEnergy ekin
                                                        ) const final
    {
      return { ekin, NC::CosineScatAngle{ 0.123 } };
    }
  };

  void testComposition( const std::string& label, const NCPI::Process& sc )
  {
    //The cross section of the composition (from its merged table where
    //available) agrees with the scaled sum of the cross sections of its
    //components, as does its table when it has one. Sampling
    //works as before:
    auto pc = dynamic_cast<const NCPI::ProcComposition*>( &sc );
    nc_assert_always( pc != nullptr );
    const auto d = pc->createIsoXSTable();
    unsigned ntab = 0;
    for ( auto& c : pc->components() ) {
      if ( c.process->createIsoXSTable().has_value() )
        ++ntab;
    }
    NC::CachePtr cp;
    std::vector<NC::CachePtr> ccp( pc->components().size() );
    auto rng = NC::getRNG();
    for ( auto e : NC::geomspace( 1e-6, 20.0, 5001 ) ) {
      const NC::NeutronEnergy ekin{ e };
      const double xs = pc->crossSectionIsotropic( cp, ekin ).dbl();
      double sum = 0.0;
      for ( auto i : NC::ncrange( pc->components().size() ) ) {
        auto& c = pc->components()[i];
        if ( c.process->domain().contains( ekin ) )
          sum += c.scale * c.process->crossSectionIsotropic( ccp[i], ekin ).dbl();
      }
      nc_assert_always( NC::floateq( xs, sum, 1e-13, 0.0 ) );
      if ( d.has_value() )
        nc_assert_always( NC::floateq( xs, NC::evalIsoXS( d.value(), e ), 1e-13, 0.0 ) );
      (void)pc->sampleScatterIsotropic( cp, *rng, ekin );
    }
    std::cout << label << ": " << pc->components().size() << " components, "
              << ntab << " with tables" << ( d.has_value() ? ", composition has a table" : "" )
              << ": OK" << std::endl;
  }

  void testMaterial( const char * cfgstr )
  {
    testComposition( cfgstr,
                     *NC::FactImpl::createScatter( NC::MatCfg( cfgstr ) ) );
  }

  void testWithUntabulated( const char * cfgstr )
  {
    //Components of a material plus a process without a table, so the
    //composition has no table and samples by evaluating its components:
    auto sc = NC::FactImpl::createScatter( NC::MatCfg( cfgstr ) );
    auto pc = dynamic_cast<const NCPI::ProcComposition*>( &*sc );
    nc_assert_always( pc != nullptr );
    NCPI::ProcComposition::ComponentList cl;
    for ( auto& c : pc->components() )
      cl.emplace_back( c.scale, c.process );
    cl.emplace_back( 0.5, NC::makeSO<UntabulatedScat>() );
    auto comb = NCPI::ProcComposition::consumeAndCombine( std::move( cl ) );
    nc_assert_always( !comb->createIsoXSTable().has_value() );
    testComposition( std::string(cfgstr) + " + UntabulatedScat", *comb );
    //The untabulated process is picked with the right frequency:
    NC::CachePtr cp;
    auto rng = NC::getRNG();
    const NC::NeutronEnergy ekin{ 0.025 };
    const double xs_tot = comb->crossSectionIsotropic( cp, ekin ).dbl();
    const double xs_u = 0.5 * ( 1.0 + 0.1 * std::sqrt( 0.025 ) );
    const unsigned n = 100000;
    unsigned nu = 0;
    for ( auto i : NC::ncrange( n ) ) {
      (void)i;
      if ( comb->sampleScatterIsotropic( cp, *rng, ekin ).mu.dbl() == 0.123 )
        ++nu;
    }
    const double expected = n * xs_u / xs_tot;
    nc_assert_always( NC::ncabs( nu - expected ) < 5.0 * std::sqrt( expected ) );
    std::cout << "  -> untabulated process sampled with the expected frequency"
              << std::endl;
  }

  void testUntabulatedOnly()
  {
    //No component has a table:
    NCPI::ProcComposition::ComponentList cl;
    cl.emplace_back( 0.5, NC::makeSO<UntabulatedScat>() );
    cl.emplace_back( 2.0, NC::makeSO<UntabulatedScat>() );
    auto comb = NCPI::ProcComposition::consumeAndCombine( std::move( cl ) );
    testComposition( "2 x UntabulatedScat", *comb );
  }
}

int main()
{
  testMaterial( "Al_sg225.ncmat;knllux=3" );
  testMaterial( "Polyethylene_CH2.ncmat;knllux=3" );
  testMaterial( "Al_sg225.ncmat" );
  testMaterial( "Y2O3_sg206_Yttrium_Oxide.ncmat;dcutoff=0.3" );
  testMaterial( "phases<0.3*Al_sg225.ncmat;knllux=3&0.7*Polyethylene_CH2.ncmat;knllux=3>" );
  testWithUntabulated( "Al_sg225.ncmat;knllux=3" );
  testUntabulatedOnly();
  return 0;
}
