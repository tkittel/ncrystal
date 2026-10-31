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


#include "NCrystal/internal/freegas/NCFreeGas.hh"
#include "NCrystal/internal/phys_utils/NCFreeGasUtils.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;

int main()
{
  //FreeGas cross sections (from tabulations) against the formula in
  //FreeGasXSProvider::crossSectionUntabulated, and as described by isoXSDescr() (to
  //within rounding):
  double worst = 0.0;
  for ( double mass : { 1.008, 2.014, 12.011, 207.2, 238.03 } ) {
    for ( double temp : { 0.1, 20.0, 293.15, 3000.0 } ) {
      NC::FreeGas fg( NC::Temperature{ temp }, NC::AtomMass{ mass }, NC::SigmaFree{ 3.5 } );
      NC::FreeGasXSProvider ref( NC::Temperature{ temp }, NC::AtomMass{ mass }, NC::SigmaFree{ 3.5 } );
      const auto d = fg.isoXSDescr();
      nc_assert_always( d.table.has_value() && d.breakpts.empty() );
      const auto& t = d.table.value();
      NC::VectD energies = NC::geomspace( 1e-14, 1e6, 50003 );
      for ( auto e : t.E ) {
        energies.push_back( e );
        energies.push_back( std::nextafter( e, 0.0 ) );
        energies.push_back( std::nextafter( e, NC::kInfinity ) );
      }
      NC::CachePtr cp;
      for ( auto e : energies ) {
        const double xs = fg.crossSectionIsotropic( cp, NC::NeutronEnergy{ e } ).dbl();
        nc_assert_always( NC::floateq( xs, NC::evalIsoXS( t, e ), 1e-14, 0.0 ) );
        const double xsref = ref.crossSectionUntabulated( NC::NeutronEnergy{ e } ).dbl();
        worst = std::max( worst, std::fabs( xs / xsref - 1.0 ) );
      }
      NC::CachePtr cp0;
      nc_assert_always( std::isfinite( fg.crossSectionIsotropic( cp0, NC::NeutronEnergy{ 0.0 } ).dbl() ) );
    }
  }
  nc_assert_always( worst < 1e-10 );
  {
    NC::FreeGas fg( NC::Temperature{ 293.15 }, NC::AtomMass{ 1.008 }, NC::SigmaFree{ 20.0 } );
    std::cout << "Number of points: " << fg.isoXSDescr().table.value().E.size() << std::endl;
  }
  std::cout << "FreeGas tabulations within 1e-10 of the formula, and as"
            " described: OK" << std::endl;
  return 0;
}
