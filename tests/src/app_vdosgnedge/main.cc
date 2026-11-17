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


//Gn functions evaluated exactly at the edges of their energy ranges, or up to
//8 ulps outside them (rounding of energies meant to be on the edge), must give
//the edge values of the raw spectra, and 0 further out, with both VDOSGn::eval
//and VDOSGn::evalMany, for linearly and cubically interpolated orders. The
//upper edge once gave the second-last value.

#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCV = NCrystal::VDOS;

namespace {

  void checkEdges( const char* label, const NC::VDOSData& vdosdata,
                   NCV::VDOSGn::Cfg cfg, unsigned maxorder )
  {
    NC::VDOS::VDOSEval ve( vdosdata );
    NCV::VDOSGn gn( ve, cfg );
    gn.growMaxOrder( maxorder );
    NC::VectD out, workbuf;
    for ( auto n : NC::ncrange( 1u, maxorder + 1 ) ) {
      const NCV::VDOSGn::Order order( n );
      const auto& spec = gn.getRawSpectrum( order );
      const auto range = gn.eRange( order );
      double specmax = 0.0;
      for ( double v : spec )
        specmax = NC::ncmax( specmax, v );
      //Exact for linear interpolation, rounding of the cubic coefficients
      //otherwise:
      auto near = [specmax]( double a, double b )
      { return NC::ncabs( a - b ) <= 1e-12 * specmax; };
      //Moves x by k ulps towards to:
      auto ulps = []( double x, double to, unsigned k )
      {
        for ( unsigned i = 0; i < k; ++i )
          x = std::nextafter( x, to );
        return x;
      };
      NC::VectD egrid;
      for ( unsigned k : { 32u, 8u, 4u, 1u } )
        egrid.push_back( ulps( range.first, -NC::kInfinity, k ) );
      egrid.push_back( range.first );
      egrid.push_back( range.second );
      for ( unsigned k : { 1u, 4u, 8u, 32u } )
        egrid.push_back( ulps( range.second, NC::kInfinity, k ) );
      for ( auto i : NC::ncrange( egrid.size() ) ) {
        const double e = egrid.at(i);
        const double v = gn.eval( order, e );
        if ( i == 0 || i + 1 == egrid.size() )
          nc_assert_always( v == 0.0 );//32 ulps outside
        else
          nc_assert_always( near( v, i < 5 ? spec.front() : spec.back() ) );
      }
      gn.evalMany( order, egrid, out, workbuf );
      for ( auto i : NC::ncrange( egrid.size() ) )
        nc_assert_always( out.at(i) == gn.eval( order, egrid.at(i) ) );
    }
    std::cout << label << ": Gn at the edges of their ranges OK for orders 1.."
              << maxorder << std::endl;
  }

}

int main()
{
  //Debye VDOS, whose G1 is large right up to its upper edge:
  const auto debye
    = NC::DynInfoUtils::createVDOSDebye( NC::DebyeTemperature{ 410.4 },
                                         NC::Temperature{ 200.0 },
                                         NC::SigmaBound{ 1.5 },
                                         NC::AtomMass{ 26.98 } );
  //Enough orders to include thinned (cubically interpolated) ones:
  checkEdges( "Debye, default", debye, NCV::VDOSGn::Cfg::Default, 40 );
  checkEdges( "Debye, legacy", debye, NCV::VDOSGn::Cfg::Legacy, 40 );
  return 0;
}
