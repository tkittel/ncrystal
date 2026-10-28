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

// Benchmark of MC-like usage of scatter processes (one CachePtr per material,
// as for separate Scatter objects). Per event: a fresh or chained energy, then
// crossSectionIsotropic followed by sampleScatterIsotropic at that energy.
//
// Args: mode cfgs nevents [chunk] [elo] [ehi] [pchain]
//   mode:   mc (xs+sample), xs (xs only), samp (sample only), rngcost; append
//           "dir" for random directions (oriented materials), "count" to count
//           random numbers
//   cfgs:   cfg-strings separated by '|' (materials are interleaved)
//   chunk:  number of consecutive events in each material (default 1)
//   elo/ehi: log-uniform energy range in eV (default 1e-5..1)
//   pchain: probability to continue with the outgoing energy (default 0.5)
// Prints the time per event and a checksum of all outcomes.

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include <iostream>
#include <chrono>
#include <cstring>

namespace NC = NCrystal;

namespace {
  class CountingRNG final : public NC::RNG {
  public:
    CountingRNG( NC::RNG& r ) : m_r(r) {}
    std::uint64_t n = 0;
    bool coinflip() override { ++n; return m_r.coinflip(); }
    std::uint64_t generate64RndmBits() override { ++n; return m_r.generate64RndmBits(); }
    std::uint32_t generate32RndmBits() override { ++n; return m_r.generate32RndmBits(); }
  protected:
    double actualGenerate() override { ++n; return m_r.generate(); }
  private:
    NC::RNG& m_r;
  };

  std::uint64_t mixbits( std::uint64_t h, double x )
  {
    std::uint64_t b;
    std::memcpy( &b, &x, sizeof(b) );
    h ^= b + 0x9e3779b97f4a7c15ull + ( h << 6 ) + ( h >> 2 );
    return h;
  }
}

int main( int argc, char** argv )
{
  if ( argc < 4 ) {
    std::cout << "Usage: mode cfgs nevents [chunk] [elo] [ehi] [pchain]"
              << std::endl;
    return 1;
  }
  const std::string mode( argv[1] );
  const auto cfgs = NC::split2( std::string(argv[2]), 0, '|' );
  auto toU = []( const char * s ) { return NC::StrView(s).toUInt64().value(); };
  auto toD = []( const char * s ) { return NC::StrView(s).toDbl().value(); };
  const std::uint64_t nevents = toU( argv[3] );
  const std::uint64_t chunk = argc > 4 ? toU( argv[4] ) : 1;
  const double elo = argc > 5 ? toD( argv[5] ) : 1e-5;
  const double ehi = argc > 6 ? toD( argv[6] ) : 1.0;
  const double pchain = argc > 7 ? toD( argv[7] ) : 0.5;

  auto rngso = NC::getRNG();
  NC::RNG& rng0 = *rngso;

  if ( mode == "rngcost" ) {
    NC::RandXRSRImpl inl( 1234 );
    double s1(0), s2(0);
    auto t0 = std::chrono::steady_clock::now();
    for ( std::uint64_t i = 0; i < nevents; ++i )
      s1 += rng0.generate();
    auto t1 = std::chrono::steady_clock::now();
    for ( std::uint64_t i = 0; i < nevents; ++i )
      s2 += inl.generate();
    auto t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double,std::nano> d1(t1-t0), d2(t2-t1);
    std::cout << "rng via RNG&: " << d1.count()/nevents << " ns, inline: "
              << d2.count()/nevents << " ns (" << s1+s2 << ")" << std::endl;
    return 0;
  }

  struct Mat { NC::ProcImpl::ProcPtr proc; NC::CachePtr cache; };
  std::vector<Mat> mats;
  for ( auto& c : cfgs )
    mats.push_back( Mat{ NC::FactImpl::createScatter( c ), NC::CachePtr() } );

  CountingRNG crng( rng0 );
  const bool count = ( mode.size() > 5 && mode.substr(mode.size()-5) == "count" );
  NC::RNG& rng = count ? static_cast<NC::RNG&>(crng) : rng0;
  const std::string m = count ? mode.substr(0,mode.size()-5) : mode;
  const bool do_dir = ( m.size() > 3 && m.substr(m.size()-3) == "dir" );
  const std::string m2 = do_dir ? m.substr(0,m.size()-3) : m;
  const bool do_xs = ( m2 == "mc" || m2 == "xs" );
  const bool do_samp = ( m2 == "mc" || m2 == "samp" );

  //Pre-generate the fresh energies and chain decisions with a separate stream,
  //so the measured RNG usage is only that of the processes:
  NC::RandXRSRImpl erng( 987 );
  const double logr = std::log( ehi / elo );
  std::uint64_t h = 0;
  double xssum = 0.0;
  std::vector<double> ekin( mats.size(), elo );
  auto t0 = std::chrono::steady_clock::now();
  for ( std::uint64_t i = 0; i < nevents; ++i ) {
    const std::size_t imat = ( i / chunk ) % mats.size();
    auto& mat = mats[imat];
    double& e = ekin[imat];
    if ( !( e >= elo && e <= ehi ) || erng.generate() >= pchain )
      e = elo * std::exp( logr * erng.generate() );
    const NC::NeutronEnergy ne{ e };
    if ( do_dir ) {
      //Random directions (fresh each event), for oriented materials:
      const double cth = 2.0 * erng.generate() - 1.0;
      const double phi = NC::k2Pi * erng.generate();
      const double sth = std::sqrt( 1.0 - cth * cth );
      const NC::NeutronDirection dir{ sth * std::cos(phi), sth * std::sin(phi), cth };
      if ( do_xs )
        xssum += mat.proc->crossSection( mat.cache, ne, dir ).dbl();
      if ( do_samp ) {
        auto out = mat.proc->sampleScatter( mat.cache, rng, ne, dir );
        h = mixbits( mixbits( h, out.ekin.dbl() ), out.direction[2] );
        e = out.ekin.dbl();
      }
      continue;
    }
    if ( do_xs )
      xssum += mat.proc->crossSectionIsotropic( mat.cache, ne ).dbl();
    if ( do_samp ) {
      auto out = mat.proc->sampleScatterIsotropic( mat.cache, rng, ne );
      h = mixbits( mixbits( h, out.ekin.dbl() ), out.mu.dbl() );
      e = out.ekin.dbl();
    }
  }
  auto t1 = std::chrono::steady_clock::now();
  std::chrono::duration<double,std::nano> dt( t1 - t0 );
  std::cout << "ns/event " << dt.count() / nevents
            << " xssum " << NC::fmt(xssum) << " hash " << h;
  if ( count )
    std::cout << " rng/event " << double(crng.n) / nevents;
  std::cout << std::endl;
  return 0;
}
