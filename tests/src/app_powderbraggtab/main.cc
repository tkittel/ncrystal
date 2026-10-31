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


#include "NCrystal/internal/powderbragg/NCPowderBragg.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include <iostream>

namespace NC = NCrystal;

namespace {
  void check( const char * label, const NC::PowderBragg& pb,
              const std::function<double(double)>& ref )
  {
    //Cross sections as described by isoXSDescr() (to within rounding), zero
    //below the first edge, C_i/E between edges with C_i increasing, and (if
    //available) within 1e-13 of an independent reference:
    const auto d = pb.isoXSDescr();
    nc_assert_always( d.table.has_value() );
    const auto& t = d.table.value();
    const auto& edges = d.breakpts;
    nc_assert_always( !edges.empty() && edges.front() == pb.domain().elow.dbl() );
    NC::VectD energies = NC::geomspace( 1e-6, 10.0, 4001 );
    for ( auto e : edges ) {
      energies.push_back( e );
      energies.push_back( std::nextafter( e, 0.0 ) );
      energies.push_back( std::nextafter( e, NC::kInfinity ) );
    }
    energies.push_back( 0.0 );
    energies.push_back( NC::kInfinity );
    std::sort( energies.begin(), energies.end() );
    NC::CachePtr cp;
    double prev_c = 0.0;
    for ( auto e : energies ) {
      const double xs = pb.crossSectionIsotropic( cp, NC::NeutronEnergy{ e } ).dbl();
      nc_assert_always( NC::floateq( xs, NC::evalIsoXS( t, e ), 1e-15, 0.0 ) );
      if ( e < edges.front() || !std::isfinite( e ) ) {
        nc_assert_always( xs == 0.0 );
        continue;
      }
      const double c = xs * e;
      nc_assert_always( c >= prev_c * ( 1.0 - 1e-15 ) );
      prev_c = c;
      if ( ref ) {
        nc_assert_always( NC::floateq( xs, ref( e ), 1e-13, 0.0 ) );
      }
    }
    std::cout << label << ": " << edges.size() << " edges, cross sections as"
              " described" << ( ref ? ", and as from the HKL list" : "" )
              << ": OK" << std::endl;
  }

  void testMaterial( const char * cfgstr )
  {
    auto info = NC::FactImpl::createInfo( NC::MatCfg( cfgstr ) );
    NC::PowderBragg pb( *info );
    //Independent reference: xs(E) = sum over planes with E >= E(2d) of
    //fsq*mult*d*wl^2/(2*V*n_atoms):
    const auto& si = info->getStructureInfo();
    auto ref = [&info,&si]( double e )
    {
      NC::StableSum sum;
      for ( const auto& hkl : info->hklList() )
        if ( e >= NC::wl2ekin( 2.0 * hkl.dspacing ) )
          sum.add( hkl.fsquared * hkl.multiplicity * hkl.dspacing );
      return sum.sum() * NC::ekin2wlsq( e ) * 0.5 / ( si.volume * si.n_atoms );
    };
    check( cfgstr, pb, ref );
  }
}

int main()
{
  testMaterial( "Al_sg225.ncmat" );
  testMaterial( "Na4Si3Al3O12Cl_sg218_Sodalite.ncmat" );
  testMaterial( "Al2O3_sg167_Corundum.ncmat;dcutoff=0.3" );
  {
    //A single plane (one edge, so the table needs an extra point above it):
    NC::PowderBragg pb( 10.0, NC::PowderBragg::VectDFM{ { 2.0, 1.0 } } );
    check( "single plane", pb, nullptr );
  }
  {
    //Merged instances (createMerged):
    auto i1 = NC::FactImpl::createInfo( NC::MatCfg( "Al_sg225.ncmat" ) );
    auto i2 = NC::FactImpl::createInfo( NC::MatCfg( "Cu_sg225.ncmat" ) );
    NC::PowderBragg pb1( *i1 ), pb2( *i2 );
    auto merged = pb1.createMerged( pb2, 0.3, 0.7 );
    nc_assert_always( merged != nullptr );
    check( "Al+Cu merged", dynamic_cast<const NC::PowderBragg&>( *merged ), nullptr );
  }
  return 0;
}
