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
  void testMaterial( const char * cfgstr )
  {
    //The cross section of the composition (from its merged table where
    //available) agrees with the scaled sum of the cross sections of its
    //components, as does its isoXSDescr() table when it has one. Sampling
    //works as before:
    auto sc = NC::FactImpl::createScatter( NC::MatCfg( cfgstr ) );
    auto pc = dynamic_cast<const NCPI::ProcComposition*>( &*sc );
    nc_assert_always( pc != nullptr );
    const auto d = pc->isoXSDescr();
    unsigned ntab = 0;
    for ( auto& c : pc->components() ) {
      auto p = dynamic_cast<const NCPI::ScatterIsotropicMat*>( &*c.process );
      if ( p && p->isoXSDescr().table.has_value() )
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
      if ( d.table.has_value() )
        nc_assert_always( NC::floateq( xs, NC::evalIsoXS( d.table.value(), e ), 1e-13, 0.0 ) );
      (void)pc->sampleScatterIsotropic( cp, *rng, ekin );
    }
    std::cout << cfgstr << ": " << pc->components().size() << " components, "
              << ntab << " with tables" << ( d.table.has_value() ? ", composition has a table" : "" )
              << ": OK" << std::endl;
  }
}

int main()
{
  testMaterial( "Al_sg225.ncmat;knllux=3" );
  testMaterial( "Polyethylene_CH2.ncmat;knllux=3" );
  testMaterial( "Al_sg225.ncmat" );
  testMaterial( "Y2O3_sg206_Yttrium_Oxide.ncmat;dcutoff=0.3" );
  testMaterial( "phases<0.3*Al_sg225.ncmat;knllux=3&0.7*Polyethylene_CH2.ncmat;knllux=3>" );
  return 0;
}
