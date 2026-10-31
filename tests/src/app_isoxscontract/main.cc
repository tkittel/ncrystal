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
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCPI = NCrystal::ProcImpl;

namespace {
  void testMaterial( const char * cfgstr )
  {
    //Every table of the leaf processes agrees with crossSectionIsotropic (to
    //within rounding), and is zero outside domain():
    auto sc = NC::FactImpl::createScatter( NC::MatCfg( cfgstr ) );
    std::vector<const NCPI::Process*> leaves;
    if ( auto c = dynamic_cast<const NCPI::ProcComposition*>( &*sc ) ) {
      for ( auto& comp : c->components() )
        leaves.push_back( &*comp.process );
    } else {
      leaves.push_back( &*sc );
    }
    std::cout << cfgstr << ":";
    for ( auto leaf : leaves ) {
      auto p = dynamic_cast<const NCPI::ScatterIsotropicMat*>( leaf );
      if ( !p )
        continue;
      const auto d = p->isoXSDescr();
      if ( !d.table.has_value() )
        continue;
      const auto& t = d.table.value();
      const auto dom = p->domain();
      NC::VectD energies = NC::geomspace( 1e-10, 1e3, 20001 );
      for ( double e : { 0.0, dom.elow.dbl(), dom.ehigh.dbl() } ) {
        if ( !std::isfinite( e ) )
          continue;
        energies.push_back( e );
        energies.push_back( std::nextafter( e, 0.0 ) );
        energies.push_back( std::nextafter( e, NC::kInfinity ) );
      }
      NC::CachePtr cp;
      for ( auto e : energies ) {
        if ( !( e >= 0.0 ) )
          continue;
        const double v = NC::evalIsoXS( t, e );
        if ( !dom.contains( NC::NeutronEnergy{ e } ) ) {
          nc_assert_always( v == 0.0 );
          continue;
        }
        const double xs = p->crossSectionIsotropic( cp, NC::NeutronEnergy{ e } ).dbl();
        nc_assert_always( NC::floateq( xs, v, 1e-14, 0.0 ) );
      }
      std::cout << " " << p->name() << "=OK";
    }
    std::cout << std::endl;
  }
}

int main()
{
  testMaterial( "Al_sg225.ncmat" );
  testMaterial( "Polyethylene_CH2.ncmat" );
  testMaterial( "Al_sg225.ncmat;knllux=3" );
  testMaterial( "Polyethylene_CH2.ncmat;knllux=3" );
  testMaterial( "He_Gas_STP.ncmat" );
  testMaterial( "Y2O3_sg206_Yttrium_Oxide.ncmat;dcutoff=0.3" );
  testMaterial( "Al_sg225.ncmat;temp=20K;comp=coh_elas" );
  return 0;
}
