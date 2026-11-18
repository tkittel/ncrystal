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


#include "NCrystal/internal/elincscatter/NCElIncScatter.hh"
#include "NCrystal/internal/phys_utils/NCElIncXS.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include <iostream>

namespace NC = NCrystal;

namespace {

  //Mono-atomic cross section, via stable_expm1 (relative precision ~1e-16):
  double refXS( double e, double msd, double sigma )
  {
    constexpr double kkk = 16.0 * NC::kPiSq * NC::ekin2wlsqinv(1.0);
    const double t = kkk * msd * e;
    return sigma * ( t < 1e-200 ? 1.0 : -NC::stable_expm1( -t ) / t );
  }

  NC::VectD testEnergies()
  {
    NC::VectD res = NC::geomspace( 1e-14, 1e6, 100003 );
    res.insert( res.begin(), { 0.0, 1e-300, 5e-324 } );
    return res;
  }

}

namespace {
  //Whether e*xs of a table of the form P_i(u)/e (as made by ElIncXS)
  //increases monotonically, checked at 50 points per piece (merged tables
  //and sampling need this, also where e*xs flattens at high energies):
  bool isMonotone( const NC::IsoXSTable& t )
  {
    for ( auto i : NC::ncrange( t.E.size() - 1 ) ) {
      const double* p = t.p.data() + 4 * i;
      double prev = p[0];
      for ( int k = 1; k <= 50; ++k ) {
        const double u = k / 50.0;
        const double v = p[0] + u * ( p[1] + u * ( p[2] + u * p[3] ) );
        if ( v < prev )
          return false;
        prev = v;
      }
    }
    return true;
  }
}

int main()
{
  const auto energies = testEnergies();

  //tabulateMonoAtomic within 1e-10 of the formula, bit-identical through
  //evalIsoXS and IsoXSEvaluator:
  double worst = 0.0;
  for ( double msd : { 0.0, 1e-12, 1e-6, 1e-4, 0.003, 0.0117, 0.05, 0.3, 3.0, 1e3 } ) {
    for ( double sigma : { 0.0, 0.07, 5.08, 82.02 } ) {
      const NC::IsoXSEvaluator ev( NC::ElIncXS::tabulateMonoAtomic( msd, NC::SigmaBound{ sigma } ) );
      const auto& t = ev.table();
      for ( auto e : energies ) {
        const double v = ev.eval( e );
        nc_assert_always( v == NC::evalIsoXS( t, e ) );
        const double ref = refXS( e, msd, sigma );
        if ( sigma == 0.0 ) {
          nc_assert_always( v == 0.0 );
          continue;
        }
        worst = std::max( worst, std::fabs( v / ref - 1.0 ) );
      }
    }
  }
  nc_assert_always( worst < 1e-10 );
  {
    auto t = NC::ElIncXS::tabulateMonoAtomic( 0.07, NC::SigmaBound{ 5.08 } );
    std::cout << "Number of points: " << t.E.size() << std::endl;
  }

  //Tiny msd values give exactly the constant table:
  {
    auto t0 = NC::ElIncXS::tabulateMonoAtomic( 0.0, NC::SigmaBound{ 1e11 } );
    for ( double msd : { 5e-324, 1e-310, 1e-300, 0.999999e-290 } ) {
      auto t = NC::ElIncXS::tabulateMonoAtomic( msd, NC::SigmaBound{ 1e11 } );
      nc_assert_always( t.E == t0.E && t.p == t0.p && t.q.empty() && t.r.empty() );
      for ( auto i : NC::ncrange( 5 ) )
        nc_assert_always( t.tail[i] == t0.tail[i] );
    }
  }
  //The threshold assumes |p| < 4*sigma/c (with sigma/c = tail[2]):
  {
    auto t = NC::ElIncXS::tabulateMonoAtomic( 0.01, NC::SigmaBound{ 5.08 } );
    for ( auto v : t.p )
      nc_assert_always( std::fabs( v ) < 4.0 * t.tail[2] );
  }
  //Just above the threshold, with the largest sigma, everything is finite:
  {
    auto t = NC::ElIncXS::tabulateMonoAtomic( 1.000001e-290, NC::SigmaBound{ 0.999999e12 } );
    for ( auto v : t.p )
      nc_assert_always( std::isfinite( v ) );
    nc_assert_always( std::isfinite( t.E.back() ) && std::isfinite( t.tail[2] ) );
    const NC::IsoXSEvaluator ev( std::move(t) );
    const auto& tt = ev.table();
    for ( double e : { 1.0, tt.E.front(), tt.E[tt.E.size()/2], tt.E.back(), 1e300 } )
      nc_assert_always( std::isfinite( ev.eval( e ) ) );
  }

  //ElIncScatter: cross sections are exactly those of its table, and within
  //1e-10 of the summed formula:
  const NC::VectD msd = { 0.0117, 0.0051, 0.0 }, bixs = { 5.08, 0.0, 2.1 },
    scale = { 0.25, 0.5, 0.25 };
  NC::ElIncScatter sc( msd, bixs, scale );
  const NC::ProcImpl::ScatterIsotropicMat& base = sc;
  auto d = base.createIsoXSTable();
  nc_assert_always( base.createIsoXSBreakpoints().empty() && d.has_value() );
  NC::CachePtr cp;
  NC::VectD many( energies.size() );
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
  sc.evalManyXSIsotropic( cp, energies.data(), energies.size(), many.data() );
#endif
  for ( auto i : NC::ncrange( energies.size() ) ) {
    const double e = energies[i];
    const double v = sc.crossSectionIsotropic( cp, NC::NeutronEnergy{ e } ).dbl();
    nc_assert_always( v == NC::evalIsoXS( d.value(), e ) );
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    nc_assert_always( v == many[i] );
#endif
    double ref = 0.0;
    for ( auto j : NC::ncrange( msd.size() ) )
      ref += scale[j] * refXS( e, msd[j], bixs[j] );
    nc_assert_always( NC::floateq( v, ref, 1e-10, 0.0 ) );
  }

  //Many elements (as from anisotropic sites, incl. one with zero msd): one
  //table, within 1e-10 of the summed formula:
  NC::VectD mm, mb, ms;
  for ( auto i : NC::ncrange( 4000 ) ) {
    const double x = std::fmod( i * 0.6180339887498949, 1.0 );
    mm.push_back( i ? 0.002 + 0.05 * x : 0.0 );
    mb.push_back( 1.0 + std::fmod( i * 0.4142135623730951, 1.0 ) );
    ms.push_back( 1.0 / 4000 );
  }
  {
    NC::ElIncXS ei( mm, mb, ms );
    const NC::IsoXSEvaluator ev( ei.tabulate() );
    double w = 0.0;
    for ( auto e : NC::geomspace( 1e-14, 1e6, 2003 ) )
      w = std::max( w, std::fabs( ev.eval( e ) / ei.evaluate( NC::NeutronEnergy{ e } ).dbl() - 1.0 ) );
    nc_assert_always( w < 1e-10 );
    nc_assert_always( isMonotone( ev.table() ) );
    std::cout << "Number of points (4000 elements): " << ev.table().E.size() << std::endl;
  }

  //Monotone pieces, also for few elements with different msd and after
  //merging tables with different energy ranges:
  {
    for ( double m : { 0.0, 1e-4, 0.0117, 0.5 } )
      nc_assert_always( isMonotone( NC::ElIncXS::tabulateMonoAtomic( m, NC::SigmaBound{ 5.08 } ) ) );
    auto t1 = NC::ElIncXS( { 0.03, 0.006, 0.004 }, { 80.3, 5.6, 4.2 }, { 0.5, 0.4, 0.1 } ).tabulate();
    auto t2 = NC::ElIncXS( { 0.008, 0.003 }, { 11.1, 1.0 }, { 0.5, 0.5 } ).tabulate();
    nc_assert_always( isMonotone( t1 ) && isMonotone( t2 ) && isMonotone( d.value() ) );
    const std::pair<double,const NC::IsoXSTable*> inputs[2] = { { 0.96, &t1 }, { 0.04, &t2 } };
    auto tm = NC::mergeIsoXSTables( inputs );
    nc_assert_always( isMonotone( tm ) );
    std::cout << "Monotone pieces: OK" << std::endl;
  }

  //Sampled mu follow the exact distribution, (summed) density
  //w*exp(-a*(1-mu)) with a=2k^2*msd:
  {
    NC::ElIncScatter scm( mm, mb, ms );
    auto rng = NC::getRNG();
    auto cdf = []( double e, const NC::VectD& m, const NC::VectD& b,
                    const NC::VectD& s, double mu )
    {
      constexpr double kk = 8.0 * NC::kPiSq * NC::ekin2wlsqinv(1.0);
      double num = 0.0, den = 0.0;
      for ( auto j : NC::ncrange( m.size() ) ) {
        const double a = kk * m[j] * e, wt = b[j] * s[j];
        auto F = [a]( double x ) { return a < 1e-200 ? 1.0 + x : std::exp( -2.0 * a ) * NC::stable_expm1( a * ( 1.0 + x ) ) / a; };
        num += wt * F( mu );
        den += wt * F( 1.0 );
      }
      return num / den;
    };
    const unsigned N = 200000;
    for ( int which : { 0, 1 } ) {
      const NC::ElIncScatter& s = which ? scm : sc;
      const NC::VectD& m = which ? mm : msd;
      const NC::VectD& b = which ? mb : bixs;
      const NC::VectD& f = which ? ms : scale;
      for ( double e : { 1e-5, 0.025, 0.3, 5.0 } ) {
        NC::VectD mus;
        mus.reserve( N );
        for ( unsigned i = 0; i < N; ++i )
          mus.push_back( s.sampleScatterIsotropic( cp, *rng, NC::NeutronEnergy{ e } ).mu.dbl() );
        std::sort( mus.begin(), mus.end() );
        double dmax = 0.0;
        for ( double mu : NC::linspace( -1.0, 1.0, 41 ) ) {
          const double emp = double( std::upper_bound( mus.begin(), mus.end(), mu ) - mus.begin() ) / N;
          dmax = std::max( dmax, std::fabs( emp - cdf( e, m, b, f, mu ) ) );
        }
        nc_assert_always( dmax < 0.006 );
      }
    }
  }
  std::cout << "All ok" << std::endl;
}
