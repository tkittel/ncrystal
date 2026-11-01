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

//Benchmark of FreeGasSampler (ns/sample) over a grid of target masses and
//E/kT, plus a hash of the sampled bits (to verify that optimisations leave
//results bitwise unchanged). Optional arg: number of samples per point.

#include "NCrystal/internal/phys_utils/NCFreeGasSampler.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <chrono>
#include <cstring>
#include <iostream>

namespace NC = NCrystal;

namespace {
  std::uint64_t mixbits( std::uint64_t h, double x )
  {
    std::uint64_t b;
    std::memcpy( &b, &x, sizeof(b) );
    return ( h ^ b ) * 0x100000001b3ull;
  }
}

int main( int argc, char** argv )
{
  const std::uint64_t n = ( argc > 1
                            ? NC::StrView(argv[1]).toUInt64().value()
                            : 2000000 );
  const NC::Temperature temp{ 293.15 };
  std::uint64_t htot = 0xcbf29ce484222325ull;
  double ttot_ab(0.0), ttot_em(0.0);
  for ( double mass : { 1.008, 4.0026, 39.95, 207.2 } ) {
    const NC::FreeGasSampler fgs( temp, NC::AtomMass{ mass } );
    for ( double c : { 1e-3, 0.1, 1.0, 10.0, 1000.0 } ) {
      const NC::NeutronEnergy ekin{ c * temp.kT() };
      auto rng = NC::createBuiltinRNG( 12345 );
      std::uint64_t h = 0xcbf29ce484222325ull;
      auto t0 = std::chrono::steady_clock::now();
      for ( std::uint64_t i = 0; i < n; ++i ) {
        auto ab = fgs.sampleAlphaBeta( *rng, ekin );
        h = mixbits( mixbits( h, ab.first ), ab.second );
      }
      auto t1 = std::chrono::steady_clock::now();
      for ( std::uint64_t i = 0; i < n; ++i ) {
        auto em = fgs.sampleDeltaEMu( *rng, ekin );
        h = mixbits( mixbits( h, em.first ), em.second );
      }
      auto t2 = std::chrono::steady_clock::now();
      const double tab = std::chrono::duration<double>(t1-t0).count();
      const double tem = std::chrono::duration<double>(t2-t1).count();
      ttot_ab += tab;
      ttot_em += tem;
      htot = mixbits( htot, double(h) );
      std::cout << "A=" << NC::fmt(mass,"%-7g") << " E/kT=" << NC::fmt(c,"%-6g")
                << " ns/sample: alphabeta=" << NC::fmt(1e9*tab/n,"%5.1f")
                << " dEmu=" << NC::fmt(1e9*tem/n,"%5.1f")
                << "  hash=" << std::hex << h << std::dec << std::endl;
    }
  }
  std::cout << "AVERAGE ns/sample: alphabeta="
            << NC::fmt(1e9*ttot_ab/(20*n),"%.1f")
            << " dEmu=" << NC::fmt(1e9*ttot_em/(20*n),"%.1f")
            << "  TOTALHASH=" << std::hex << htot << std::dec << std::endl;
  return 0;
}
