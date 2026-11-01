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

//Validates FreeGasSampler against the S(alpha,beta)-based sampler of NCrystal
//4.4 and earlier (FreeGasLegacySampler), with two-sample chi2 tests of the
//(alpha,beta) and (dE,mu) distributions over a grid of target masses and
//energies. Pass an argument to print detailed statistics.

#include "NCrystal/internal/phys_utils/NCFreeGasSampler.hh"
#include "NCFreeGasLegacySampler.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include <iostream>

namespace NC = NCrystal;

namespace {

  //Two-sample chi2/ndf (equal sizes), with nb bins equiprobable in the
  //reference sample ref (sorted in place):
  double chi2ndf( NC::VectD& ref, const NC::VectD& other, unsigned nb )
  {
    std::sort( ref.begin(), ref.end() );
    NC::VectD edges;
    for ( auto i : NC::ncrange( 1u, nb ) )
      edges.push_back( ref.at( ( ref.size() * i ) / nb ) );
    edges.erase( std::unique( edges.begin(), edges.end() ), edges.end() );
    auto bin = [&edges]( double x )
    {
      return std::size_t( std::upper_bound( edges.begin(), edges.end(), x )
                          - edges.begin() );
    };
    NC::VectD h1( edges.size()+1, 0.0 ), h2( edges.size()+1, 0.0 );
    for ( auto x : ref )
      h1.at( bin(x) ) += 1.0;
    for ( auto x : other )
      h2.at( bin(x) ) += 1.0;
    double chi2(0.0);
    unsigned ndf(0);
    for ( auto i : NC::ncrange( h1.size() ) ) {
      const double s = h1.at(i) + h2.at(i);
      if ( s > 0.0 ) {
        chi2 += NC::ncsquare( h1.at(i) - h2.at(i) ) / s;
        ++ndf;
      }
    }
    return ndf > 1 ? chi2 / ( ndf - 1 ) : 0.0;
  }

  //2D (beta-quantile x mu) binning folded into a single 1D variable:
  NC::VectD fold2D( const NC::VectD& b, const NC::VectD& mu, NC::VectD bref )
  {
    std::sort( bref.begin(), bref.end() );
    constexpr unsigned nb = 10;
    NC::VectD edges;
    for ( auto i : NC::ncrange( 1u, nb ) )
      edges.push_back( bref.at( ( bref.size() * i ) / nb ) );
    NC::VectD res;
    res.reserve( b.size() );
    for ( auto i : NC::ncrange( b.size() ) ) {
      auto ib = std::upper_bound( edges.begin(), edges.end(), b.at(i) )
        - edges.begin();
      auto imu = NC::ncmin( 9, int( ( mu.at(i) + 1.0 ) * 5.0 ) );
      res.push_back( double( ib * 10 + imu ) + 0.5 );
    }
    return res;
  }

}

namespace {
  //Smoke test at extreme inputs: results must be finite and kinematically ok.
  void testExtremes( NC::RNG& rng )
  {
    const NC::Temperature temp{ 1.0 };
    for ( double mass : { 0.5, 1.008, 250.0 } ) {
      const NC::FreeGasSampler ng( temp, NC::AtomMass{ mass } );
      for ( double c : { 0.0, 1e-12, 1e-5, 1e10, 1e14 } ) {
        const NC::NeutronEnergy ekin{ c * temp.kT() };
        for ( auto i : NC::ncrange( 2000 ) ) {
          (void)i;
          auto ab = ng.sampleAlphaBeta( rng, ekin );
          nc_assert_always( std::isfinite(ab.first) );
          nc_assert_always( std::isfinite(ab.second) );
          nc_assert_always( ab.second >= -c );
          auto alim = NC::getAlphaLimits( NC::EDivKT{c}, ab.second );
          nc_assert_always( ab.first >= alim.first && ab.first <= alim.second );
          auto em = ng.sampleDeltaEMu( rng, ekin );
          nc_assert_always( std::isfinite(em.first) );
          nc_assert_always( em.first >= -ekin.dbl() );
          nc_assert_always( em.second >= -1.0 && em.second <= 1.0 );
        }
      }
    }
    std::cout << "Extreme-input smoke tests OK" << std::endl;
  }
}

int main( int argc, char** )
{
  const bool verbose = argc > 1;
  const NC::Temperature temp{ 293.15 };
  const unsigned nsample = 50000;
  auto rng = NC::createBuiltinRNG( 123456 );
  testExtremes( *rng );
  double worst(0.0);
  for ( double mass : { 1.008, 4.0026, 39.95, 207.2 } ) {
    const NC::AtomMass m{ mass };
    const NC::FreeGasSampler ng( temp, m );
    for ( double c : { 1e-3, 0.1, 1.0, 10.0, 1000.0, 1e5, 1e7 } ) {
      const NC::NeutronEnergy ekin{ c * temp.kT() };
      const NC::FreeGasLegacySampler ref( ekin, temp, m );
      NC::VectD a1, b1, a2, b2, e1, mu1, e2, mu2;
      for ( auto i : NC::ncrange( nsample ) ) {
        (void)i;
        auto x1 = ref.sampleAlphaBeta( *rng );
        auto x2 = ng.sampleAlphaBeta( *rng, ekin );
        a1.push_back( x1.first );
        b1.push_back( x1.second );
        a2.push_back( x2.first );
        b2.push_back( x2.second );
        auto y1 = ref.sampleDeltaEMu( *rng );
        auto y2 = ng.sampleDeltaEMu( *rng, ekin );
        e1.push_back( y1.first );
        mu1.push_back( y1.second );
        e2.push_back( y2.first );
        mu2.push_back( y2.second );
      }
      auto f1 = fold2D( e1, mu1, e1 );
      auto f2 = fold2D( e2, mu2, e1 );
      const double cb = chi2ndf( b1, b2, 50 );
      const double ca = chi2ndf( a1, a2, 50 );
      const double ce = chi2ndf( e1, e2, 50 );
      const double cm = chi2ndf( mu1, mu2, 50 );
      const double c2 = chi2ndf( f1, f2, 100 );
      const double w = NC::ncmax( NC::ncmax( cb, ca ),
                                  NC::ncmax( NC::ncmax( ce, cm ), c2 ) );
      worst = NC::ncmax( worst, w );
      if ( verbose )
        std::cout << "A=" << NC::fmt(mass) << " E/kT=" << NC::fmt(c)
                  << " chi2/ndf: beta=" << NC::fmt(cb,"%.2f")
                  << " alpha=" << NC::fmt(ca,"%.2f")
                  << " dE=" << NC::fmt(ce,"%.2f")
                  << " mu=" << NC::fmt(cm,"%.2f")
                  << " (dE,mu)=" << NC::fmt(c2,"%.2f") << std::endl;
    }
  }
  if ( verbose )
    std::cout << "worst chi2/ndf: " << NC::fmt(worst,"%.2f") << std::endl;
  nc_assert_always( worst < 2.0 );
  std::cout << "All free-gas sampler comparisons OK" << std::endl;
  return 0;
}
