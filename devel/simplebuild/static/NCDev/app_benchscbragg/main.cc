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

// Benchmark of cross section evaluations and scatterings in single crystals,
// for neutrons of random directions (so the caches of the models are of no
// use). Since only Bragg diffraction is of interest, "comp=bragg" is added to
// the cfg-string, which must include the mosaicity and orientation.

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include <chrono>
#include <iostream>

namespace NC = NCrystal;

int main( int argc, char** argv )
{
  if ( argc < 2 ) {
    std::cout << "Please provide args: <cfgstr> [ndirs] [wl1 wl2 ...]\n"
      "Example cfgstr: \"Ge_sg227.ncmat;mos=0.5deg;dir1=@crys_hkl:5,1,1"
      "@lab:0,0,1;dir2=@crys_hkl:0,-1,1@lab:0,1,0\"" << std::endl;
    return 1;
  }
  auto sc = NC::FactImpl::createScatter( std::string(argv[1])+";comp=bragg" );
  std::size_t ndirs = 10000;
  if ( argc >= 3 )
    ndirs = static_cast<std::size_t>( NC::StrView(argv[2]).toUInt64().value() );
  NC::VectD wls;
  for ( int i = 3; i < argc; ++i )
    wls.push_back( NC::StrView(argv[i]).toDbl().value() );
  if ( wls.empty() )
    wls = { 0.25, 0.5, 1.0, 2.0, 4.0 };

  auto rng = NC::getRNG();
  std::vector<NC::NeutronDirection> dirs;
  dirs.reserve( ndirs );
  for ( auto i : NC::ncrange(ndirs) ) {
    (void)i;
    dirs.push_back( NC::randIsotropicNeutronDirection(rng) );
  }

  auto seconds_since = []( std::chrono::steady_clock::time_point t0 )
  {
    return std::chrono::duration<double>
      ( std::chrono::steady_clock::now() - t0 ).count();
  };

  NC::CachePtr cp;
  std::cout << "wavelength[Aa]  <xs>[barn]  xs[us/call]  xs+scatter[us/call]\n";
  for ( auto wl : wls ) {
    const NC::NeutronEnergy ekin = NC::NeutronWavelength{ wl };
    double sumxs(0.0), t_xs(NC::kInfinity), t_xsscat(NC::kInfinity);
    for ( int rep = 0; rep < 3; ++rep ) {
      sumxs = 0.0;
      auto t0 = std::chrono::steady_clock::now();
      for ( auto& d : dirs )
        sumxs += sc->crossSection( cp, ekin, d ).get();
      t_xs = NC::ncmin( t_xs, seconds_since(t0) );
      t0 = std::chrono::steady_clock::now();
      for ( auto& d : dirs ) {
        if ( sc->crossSection( cp, ekin, d ).get() > 0.0 )
          sc->sampleScatter( cp, *rng, ekin, d );
      }
      t_xsscat = NC::ncmin( t_xsscat, seconds_since(t0) );
    }
    std::cout << wl << "  " << sumxs/ndirs << "  " << t_xs*1e6/ndirs
              << "  " << t_xsscat*1e6/ndirs << std::endl;
  }
  return 0;
}
