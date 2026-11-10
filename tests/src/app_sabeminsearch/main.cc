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
#include "NCrystal/internal/utils/NCMath.hh"
#include <functional>
#include <iostream>
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace detail {
      //Defined in NCSABProcessor.cc:
      PairDD diagonalLowEnergyLimit( const SABData& );
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
    const auto res = NC::SABUtils::detail::diagonalLowEnergyLimit( data );
    const double expect_f0 = 4.0 * s * sqrtLinIntegral( 0.1, 4.0 );
    const double expect_msa = 0.5 * ( 16.0 - 0.01 ) / sqrtLinIntegral( 0.1, 4.0 );
    nc_assert_always( NC::floateq( res.first, expect_f0, 1e-12, 0.0 ) );
    nc_assert_always( NC::floateq( res.second, expect_msa, 1e-9, 0.0 ) );
    std::cout << "diagonalLowEnergyLimit OK" << std::endl;
  }
}

int main()
{
  testDiagonalLimit();

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
