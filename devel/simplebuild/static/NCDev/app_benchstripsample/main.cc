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


#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/sab/NCSABBetaStripSampler.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include <iostream>
#include <chrono>
#include <cstdio>

// Times (alpha,beta) sampling of a kernel with SABProcessor, with the switch
// for beta-strip sampling off and on, and with a SABBetaStripSampler alone.
// For each E/kT=x, energies are log-uniform in [x/sqrt(5),x*sqrt(5)] (the
// strip sampler is built for the upper end). Prints setup times, then ns and
// tries per sample for each x.

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

namespace {
  double now()
  {
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double>( clock::now().time_since_epoch() )
      .count();
  }

  //Calls f(rng,e) for e in evals until at least 0.2s have passed, returning
  //ns/call and tries/call:
  template<class TFct>
  NC::PairDD timeit( TFct f, const NC::VectD& evals )
  {
    auto rng = NC::createBuiltinRNG( 12345 );
    std::uint64_t n = 0, ntries = 0;
    double sum_ab = 0.0;
    const double t0 = now();
    double t1;
    do {
      for ( auto e : evals ) {
        auto r = f( *rng, e );
        ntries += r.ntries;
        sum_ab += r.alpha + r.beta;
      }
      n += evals.size();
      t1 = now();
    } while ( t1 - t0 < 0.2 );
    if ( sum_ab == 0.123456789 )
      std::cout << "unlikely" << std::endl;//prevent optimising away
    return { 1e9 * ( t1 - t0 ) / double(n), double(ntries) / double(n) };
  }
}

int main( int argc, char** argv )
{
  if ( argc < 4 ) {
    std::cout << "Please provide args: [cfgstr] [displaylabel] [E/kT] .."
              << " [E/kT]" << std::endl;
    std::cout << "NB: displaylabel can be empty (\"\") for monoatomic materials"
              << std::endl;
    return 1;
  }
  const std::string cfgstr(argv[1]);
  NC::MatCfg matcfg(cfgstr);
  const std::string displaylbl(argv[2]);
  NC::VectD xvals;
  for ( int i = 3; i < argc; ++i ) {
    auto v = NC::StrView(argv[i]).toDbl();
    nc_assert_always( v.has_value() && v.value() > 0.0 );
    xvals.push_back( v.value() );
  }
  NC::Optional<std::string> dlbl;
  if ( !displaylbl.empty() )
    dlbl = displaylbl;
  auto info = NC::FactImpl::createInfo(cfgstr);
  auto di_knl = dynamic_cast<const NC::DI_ScatKnl*>(
    NC::InfoUtils::findDynInfo( info, dlbl ) );
  if (!di_knl)
    NCRYSTAL_THROW(BadInput,"Selected component does not provide a SAB knl");
  auto sab = NC::extractSABDataFromDynInfo(
    di_knl, NC::VDOS::VDOSLux( matcfg.get_vdoslux() ) );
  const double kT = sab->temperature().kT();
  const int knllux = ( matcfg.get_knllux() < 0
                       ? NC::SABCfg::sablux_default_luxury
                       : matcfg.get_knllux() );

  auto cfg_off = NC::SABCfg::createConfig( knllux );
  cfg_off.betaStripSampling = false;
  auto cfg_on = NC::SABCfg::createConfig( knllux );
  cfg_on.betaStripSampling = true;
  double t0 = now();
  NCS::SABProcessor poff( cfg_off, sab, nullptr );
  const double tsetup_off = now() - t0;
  t0 = now();
  NCS::SABProcessor pon( cfg_on, sab, nullptr );
  const double tsetup_on = now() - t0;
  const double xmax = pon.getEMaxInfo().E_div_kT;
  std::printf( "SABProcessor setup [ms]: beta-strip sampling off %.1f, on"
               " %.1f\n", tsetup_off*1e3, tsetup_on*1e3 );
  std::printf( "     E/kT |  off: tries    ns |   on: tries    ns |"
               " strip: tries    ns setup[ms]\n" );
  auto procfct = [kT]( const NCS::SABProcessor& p )
  {
    return [&p,kT]( NC::RNG& rng, double e )
    { return p.sampleScatterAlphaBeta( rng, NC::NeutronEnergy{ e*kT } ); };
  };
  const double sqrtwin = std::sqrt( 5.0 );
  for ( auto x : xvals ) {
    const double xhi = x * sqrtwin;
    if ( xhi > xmax )
      continue;
    NC::VectD evals;
    {
      auto rng = NC::createBuiltinRNG( 54321 );
      for ( int i = 0; i < 4000; ++i )
        evals.push_back( x * std::pow( 5.0, rng->generate() - 0.5 ) );
    }
    auto roff = timeit( procfct( poff ), evals );
    auto ron = timeit( procfct( pon ), evals );
    t0 = now();
    NCS::SABBetaStripSampler st( sab, xhi );
    const double tst = now() - t0;
    if ( !st.hasStrips() ) {
      std::printf( "%9.2g | %11.3g %5.0f | %11.3g %5.0f |  (no strips)\n",
                   x, roff.second, roff.first, ron.second, ron.first );
      continue;
    }
    auto rst = timeit( [&st]( NC::RNG& rng, double e )
                       { return st.sample( rng, e ); }, evals );
    std::printf( "%9.2g | %11.3g %5.0f | %11.3g %5.0f | %12.3g %5.0f %9.2g\n",
                 x, roff.second, roff.first, ron.second, ron.first,
                 rst.second, rst.first, tst*1e3 );
  }
  return 0;
}
