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


//Tests of the search for the lowest energy grid point, Emin, of SAB kernels,
//using synthetic f(E)=Sint(E)/sqrt(E) curves and a simple kernel.

#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/sab/NCSABCfg.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <functional>
#include <iostream>
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace detail {
      //Defined in NCSABProcessor.cc:
      PairDD diagonalLowEnergyLimit( const SABData&, const VectD& logS );
      double determineEMinDivKT( const SABCfg::Cfg&, double emax,
                                 const std::function<double(double)>& f_of_e,
                                 const PairDD& diaglimit );
    }
  }
}

namespace {

  using FctE = std::function<double(double)>;
  const double kEMax = 10.0;

  double findEMin( const FctE& f, double f0 = 1.0 )
  {
    NC::SABCfg::Cfg cfg;
    return NC::SABUtils::detail::determineEMinDivKT( cfg, kEMax, f,
                                                     NC::PairDD( f0, 1.0 ) );
  }

  //Reference: dense scan for the first E where |log(f(E)/f(0))| reaches 1e-3
  //of |log(f(EMax)/f(0))|. Returns the log(E) interval containing it:
  NC::PairDD denseScan( const FctE& f )
  {
    const double logref = std::log( f( 1e-30 ) );
    const double target = 1e-3 * NC::ncabs( std::log( f(kEMax) ) - logref );
    const double lo = std::log( 1e-15 );
    const double hi = std::log( kEMax );
    const unsigned n = 200000;
    double prev = lo;
    for ( auto i : NC::ncrange( 1u, n + 1 ) ) {
      const double x = lo + ( hi - lo ) * i / n;
      if ( NC::ncabs( std::log( f( std::exp(x) ) ) - logref ) >= target )
        return { prev, x };
      prev = x;
    }
    return { hi, hi };
  }

  void testCurve( const char * name, const FctE& f )
  {
    const double emin = findEMin( f );
    const auto ref = denseScan( f );
    const double x = std::log( emin );
    nc_assert_always( x >= ref.first - 1e-6 && x <= ref.second + 1e-6 );
    std::cout << name << ": Emin/kT = " << NC::fmt( emin, "%.4g" )
              << " (dense scan agrees)" << std::endl;
  }

  //Gaussian bump in log(E) at e0:
  double bump( double e, double e0, double height )
  {
    const double d = ( std::log( e ) - std::log( e0 ) ) / 0.3;
    return height * std::exp( -0.5 * d * d );
  }

  double sqrtLinIntegral( double a, double b )
  {
    return ( 2.0 / 3.0 ) * ( b * std::sqrt(b) - a * std::sqrt(a) );
  }

  void testDiagonalLimit()
  {
    //S constant, so f0=4*S*int sqrt(x)dx and <sqrt(alpha)>=int x/int sqrt(x),
    //over the diagonal from max(0,alpha_min) to min(alpha_max,beta_max):
    const double s = 0.7;
    NC::VectD a = { 0.1, 0.5, 1.0, 2.0, 5.0 };
    NC::VectD b = { -6.0, -1.0, 0.0, 0.3, 2.5, 4.0 };
    NC::VectD sab( a.size() * b.size(), s );
    NC::SABData data( std::move(a), std::move(b), std::move(sab),
                      NC::Temperature{ 293.15 }, NC::SigmaBound{ 1.0 },
                      NC::AtomMass{ 1.0 } );
    const auto logS = NC::SABUtils::createLogSTable( data );
    namespace NSD = NC::SABUtils::detail;
    const auto res = NSD::diagonalLowEnergyLimit( data, *logS );
    const double expect_f0 = 4.0 * s * sqrtLinIntegral( 0.1, 4.0 );
    const double expect_msa = 0.5 * ( 16.0 - 0.01 ) / sqrtLinIntegral( 0.1, 4.0 );
    nc_assert_always( NC::floateq( res.first, expect_f0, 1e-12, 0.0 ) );
    nc_assert_always( NC::floateq( res.second, expect_msa, 1e-9, 0.0 ) );
    std::cout << "diagonalLowEnergyLimit OK" << std::endl;
  }

  void testDiagonalLimitExtremeS()
  {
    //Subnormal S next to S=1 (S ratio above the double range) and rows
    //differing by 1e300: finite, and as a dense numerical integration of
    //4*sqrt(x)*S(x,x) along the diagonal:
    NC::VectD a = { 0.1, 1.0, 2.0 };
    NC::VectD b = { -0.5, 1.5, 2.5 };
    NC::VectD sab = { 1e-320, 1.0, 0.5,
                      1e-300, 2e-300, 3.0,
                      0.2, 0.4, 0.1 };
    auto sOf = [&a,&b,&sab]( double x )
    {
      //S(x,x), log-linear in alpha (logs of each value, as no ratio is
      //representable), linear in beta:
      const std::size_t ia = ( x < a[1] ? 0 : 1 );
      const std::size_t ib = ( x < b[1] ? 0 : 1 );
      auto row = [&]( std::size_t jb ) {
        const double fa = sab[jb*3+ia];
        const double fb = sab[jb*3+ia+1];
        const double r = ( x - a[ia] ) / ( a[ia+1] - a[ia] );
        return std::exp( std::log( fa ) + r * ( std::log(fb) - std::log(fa) ) );
      };
      const double u = ( x - b[ib] ) / ( b[ib+1] - b[ib] );
      return row( ib ) * ( 1.0 - u ) + row( ib + 1 ) * u;
    };
    NC::StableSum ref;
    const unsigned n = 400000;
    const double lo = 0.1, hi = 2.0, dx = ( hi - lo ) / n;
    for ( auto i : NC::ncrange( n ) ) {
      const double x = lo + ( i + 0.5 ) * dx;
      ref.add( 4.0 * std::sqrt( x ) * sOf( x ) * dx );
    }
    NC::SABData data( std::move(a), std::move(b), std::move(sab),
                      NC::Temperature{ 293.15 }, NC::SigmaBound{ 1.0 },
                      NC::AtomMass{ 1.0 } );
    const auto logS = NC::SABUtils::createLogSTable( data );
    namespace NSD = NC::SABUtils::detail;
    const auto res = NSD::diagonalLowEnergyLimit( data, *logS );
    nc_assert_always( std::isfinite( res.first ) );
    nc_assert_always( std::isfinite( res.second ) );
    nc_assert_always( NC::floateq( res.first, ref.sum(), 1e-7, 0.0 ) );
    std::cout << "diagonalLowEnergyLimit with extreme S values OK" << std::endl;
  }

  NC::PairDD diagLimit( NC::VectD a, NC::VectD b, NC::VectD sab )
  {
    NC::SABData data( std::move(a), std::move(b), std::move(sab),
                      NC::Temperature{ 293.15 }, NC::SigmaBound{ 1.0 },
                      NC::AtomMass{ 1.0 } );
    const auto logS = NC::SABUtils::createLogSTable( data );
    return NC::SABUtils::detail::diagonalLowEnergyLimit( data, *logS );
  }

  //Kernel with the same S(alpha) in all beta rows:
  NC::PairDD diagLimitBetaConst( const NC::VectD& a, const NC::VectD& b,
                                 const std::function<double(double)>& s )
  {
    NC::VectD sab;
    for ( auto ib : NC::ncrange( b.size() ) ) {
      (void)ib;
      for ( double x : a )
        sab.push_back( s( x ) );
    }
    return diagLimit( a, b, sab );
  }

  void testDiagonalLimitAnalytic()
  {
    //Diagonal from alpha=0 (sqrt(x) singularity), constant S:
    //f0=4*S*(2/3)*h^1.5 and <sqrt(alpha)>=(h^2/2)/((2/3)h^1.5):
    {
      const double s = 0.37, h = 3.0;
      auto res = diagLimitBetaConst( { 0.0, 0.01, 1.0, 3.0 },
                                     { -1.0, 0.5, 2.0, 7.0 },
                                     [s]( double ) { return s; } );
      nc_assert_always( NC::floateq( res.first,
                                     4.0 * s * ( 2.0 / 3.0 ) * h * std::sqrt(h),
                                     1e-14, 0.0 ) );
      nc_assert_always( NC::floateq( res.second, 0.75 * std::sqrt( h ),
                                     1e-14, 0.0 ) );
    }
    //Linear S (a zero corner), from alpha=0: S=x on [0,1], so f0=4*2/5 and
    //<sqrt(alpha)>=(1/3)/(2/5):
    {
      auto res = diagLimitBetaConst( { 0.0, 1.0 }, { -1.0, 2.0 },
                                     []( double x ) { return x; } );
      nc_assert_always( NC::floateq( res.first, 1.6, 1e-14, 0.0 ) );
      nc_assert_always( NC::floateq( res.second, 5.0 / 6.0, 1e-14, 0.0 ) );
    }
    //S=exp(k*alpha), with log changes up to ~690 within a cell (many pieces)
    //and the diagonal from alpha=0. Closed forms over [0,h] for k<0:
    //int sqrt(x)e^(kx) = [sqrt(x)e^(kx)/k]
    //                    + sqrt(pi)*erf(sqrt(-kx))/(2(-k)^1.5)
    //int x e^(kx) = [e^(kx)*(x/k-1/k^2)]:
    for ( double k : { -0.3, -5.0, -60.0, -230.0 } ) {
      const double h = 3.0;
      auto expk = [k]( double x ) { return std::exp( k * x ); };
      auto res = diagLimitBetaConst( { 0.0, 0.5, 1.5, 3.0 },
                                     { -1.0, 1.0, 4.0 }, expk );
      const double mk = -k;
      const double i1 = std::sqrt(h) * std::exp( k * h ) / k
        + std::sqrt( NC::kPi ) * std::erf( std::sqrt( mk * h ) )
        / ( 2.0 * mk * std::sqrt( mk ) );
      const double i2 = std::exp( k * h ) * ( h / k - 1.0 / ( k * k ) )
        + 1.0 / ( k * k );
      nc_assert_always( NC::floateq( res.first, 4.0 * i1, 1e-12, 0.0 ) );
      nc_assert_always( NC::floateq( res.second, i2 / i1, 1e-12, 0.0 ) );
    }
    //All S zero, and no overlap of the grid with the diagonal:
    {
      auto r1 = diagLimitBetaConst( { 0.1, 1.0 }, { -1.0, 2.0 },
                                    []( double ) { return 0.0; } );
      nc_assert_always( r1.first == 0.0 && r1.second == 0.0 );
      auto r2 = diagLimitBetaConst( { 2.0, 3.0 }, { -1.0, 1.5 },
                                    []( double ) { return 1.0; } );
      nc_assert_always( r2.first == 0.0 && r2.second == 0.0 );
    }
    std::cout << "diagonalLowEnergyLimit analytic cases OK" << std::endl;
  }
}

int main()
{
  testDiagonalLimit();
  testDiagonalLimitExtremeS();
  testDiagonalLimitAnalytic();

  testCurve( "monotonic", []( double e ) { return 1.0 + e; } );
  testCurve( "dip", []( double e ) { return 1.0 + e - bump( e, 1e-6, 0.01 ); } );
  testCurve( "peak", []( double e ) { return 1.0 + e + bump( e, 1e-6, 0.01 ); } );

  //No low energy limit (f0=0):
  nc_assert_always( findEMin( []( double e ) { return 1.0 + e; }, 0.0 ) == kEMax );
  std::cout << "f0=0: returns Emax" << std::endl;

  //f already changes at the first scan point:
  bool threw = false;
  try {
    findEMin( []( double e ) { return std::pow( e, 0.1 ); } );
  } catch ( NC::Error::CalcError& ) {
    threw = true;
  }
  nc_assert_always( threw );
  std::cout << "f changing at the first scan point: CalcError" << std::endl;
}
