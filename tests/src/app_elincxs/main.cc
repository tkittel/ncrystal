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
  std::cout << "All ok" << std::endl;
}
