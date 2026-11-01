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

//Tests of the IsoMacroXSTable creation with synthetic cross section curves: the
//tolerance must be met for complicated shapes, and problems which could lead
//to wrong tables must result in exceptions.

#include "NCrystal/internal/xstable/NCIsoMacroXSTable.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCX = NCrystal::XSTable;

namespace {

  using Fct = std::function<double(double)>;

  //Worst relative error of the table vs. the function, at log-spaced points
  //(from NCX::denseWavelengthMin(), below which the tolerance is not
  //required) and at points close to the given discontinuities. Where the
  //function is 0, the table must be exactly 0:
  double worstRelErr( const NCX::IsoMacroXSTable& t, const Fct& f,
                      const NC::VectD& discont, double wlmax )
  {
    NC::VectD wls = NC::geomspace( NCX::denseWavelengthMin(), wlmax, 200000 );
    for ( auto e : discont )
      for ( auto r : { 1-1e-8, 1+1e-8, 1-1e-5, 1+1e-5 } )
        wls.push_back( e * r );
    double worst = 0.0;
    for ( auto wl : wls ) {
      const double exact = f( wl );
      const double tab = NCX::evalIsoMacroXSTable( t.wavelengths().data(),
                                                   t.macroXS().data(),
                                                   t.slopes().data(), t.n(),
                                                   wl );
      const double err = ( exact > 0.0
                           ? NC::ncabs( tab - exact ) / exact
                           : ( tab == exact ? 0.0 : NC::kInfinity ) );
      worst = std::max( worst, err );
    }
    return worst;
  }

  //Numerical derivative, for smooth functions:
  Fct numSlope( const Fct& f )
  {
    return [f]( double wl )
    {
      const double w1 = wl * ( 1.0 - 1e-6 );
      const double w2 = wl * ( 1.0 + 1e-6 );
      return ( f( w2 ) - f( w1 ) ) / ( w2 - w1 );
    };
  }

  Fct constSlope( double s )
  {
    return [s]( double ) { return s; };
  }

  void testOK( const char* descr, const Fct& f, const Fct& df,
               NC::VectD discont = {}, double tol = 1e-3, double wlmax = 500.0 )
  {
    NCX::IsoMacroXSTableParams params;
    params.tol = tol;
    params.wlmax = wlmax;
    auto t = NCX::createIsoMacroXSTable( f, df, discont, params );
    const double worst = worstRelErr( t, f, discont, wlmax );
    nc_assert_always( worst <= tol );
    std::cout << descr << ": npts=" << t.n()
              << " xs(0)=" << std::round( t.macroXS().front() * 1e6 ) / 1e6
              << " worst/tol<1: OK" << std::endl;
  }

  template <class TError>
  void testError( const char* descr, const Fct& f, const Fct& df,
                  NC::VectD discont = {}, NCX::IsoMacroXSTableParams params = {} )
  {
    try {
      NCX::createIsoMacroXSTable( f, df, discont, params );
    } catch ( TError& e ) {
      std::cout << descr << ": got expected exception: " << e.what() << std::endl;
      return;
    }
    std::cout << descr << ": did not get the expected exception!" << std::endl;
    nc_assert_always( false );
  }

  double step( double wl, double e, double below, double above )
  {
    return wl < e ? below : above;
  }

}

int main()
{
  //Shapes which must be handled (with exact derivatives where simple):
  testOK( "constant", []( double ) { return 5.0; }, constSlope( 0.0 ) );
  testOK( "zero", []( double ) { return 0.0; }, constSlope( 0.0 ) );
  testOK( "a+b*wl", []( double wl ) { return 3.0 + 0.5 * wl; }, constSlope( 0.5 ) );
  testOK( "a/wl^0.1 at long wavelengths",
          []( double wl ) { return 1.0 + wl / std::pow( 1.0 + wl, 0.1 ); },
          []( double wl )
          {
            return ( std::pow( 1.0 + wl, -0.1 )
                     - 0.1 * wl * std::pow( 1.0 + wl, -1.1 ) );
          } );
  testOK( "vanishing at wl=0",
          []( double wl ) { return wl * wl / ( 1.0 + wl ); },
          []( double wl ) { return wl * ( wl + 2.0 ) / ( ( 1.0 + wl ) * ( 1.0 + wl ) ); } );
  testOK( "vanishing at wl=0 with large curvature",
          []( double wl ) { return 10.0 * wl * wl / ( 1.0 + wl * wl ); },
          []( double wl )
          {
            const double d = 1.0 + wl * wl;
            return 20.0 * wl / ( d * d );
          } );
  //No absolute scale is involved:
  testOK( "vanishing at wl=0, tiny",
          []( double wl ) { return 1e-30 * wl * wl / ( 1.0 + wl ); },
          []( double wl ) { return 1e-30 * wl * ( wl + 2.0 ) / ( ( 1.0 + wl ) * ( 1.0 + wl ) ); } );
  testOK( "vanishing at wl=0 like wl^4",
          []( double wl )
          {
            const double q = wl * wl;
            return q * q / ( 1.0 + q * q );
          },
          []( double wl )
          {
            const double d = 1.0 + wl * wl * wl * wl;
            return 4.0 * wl * wl * wl / ( d * d );
          } );
  testOK( "oscillating",
          []( double wl ) { return 2.0 + std::sin( 20.0 * wl ); },
          []( double wl ) { return 20.0 * std::cos( 20.0 * wl ); },
          {}, 1e-3, 20.0 );
  //Resonance-like peaks (Lorentzians in energy), as long as they are wider
  //than the dense sampling:
  Fct peaks = []( double wl )
  {
    const double e = NC::wl2ekin( wl );
    double v = 1.0 + 0.3 * wl;
    for ( double e0 : { 0.5, 3.0, 20.0 } ) {
      const double g = 0.01 * e0;//width relative to the energy: 1%
      v += 50.0 * g * g / ( ( e - e0 ) * ( e - e0 ) + g * g );
    }
    return v;
  };
  testOK( "peaks", peaks, numSlope( peaks ) );
  //Declared discontinuities (up and down), also with kinks:
  testOK( "steps",
          []( double wl )
          {
            return ( 1.0 + 0.1 * wl + step( wl, 2.0, 0.0, 3.0 )
                     + step( wl, 4.0, 5.0, 0.0 ) + step( wl, 7.0, 0.0, 0.2 * ( wl - 7.0 ) ) );
          },
          []( double wl ) { return 0.1 + step( wl, 7.0, 0.0, 0.2 ); },
          { 2.0, 4.0, 7.0 } );
  testOK( "drop to zero",
          []( double wl ) { return step( wl, 5.0, 2.0 + wl, 0.0 ); },
          []( double wl ) { return step( wl, 5.0, 1.0, 0.0 ); },
          { 5.0 } );

  //Problems which must result in exceptions:
  testError<NC::Error::CalcError>( "unexpected step",
                                   []( double wl ) { return 1.0 + step( wl, 2.0, 0.0, 3.0 ); },
                                   constSlope( 0.0 ) );
  testError<NC::Error::CalcError>( "no limit for wl->0",
                                   []( double wl ) { return 1.0 / wl; },
                                   []( double wl ) { return -1.0 / ( wl * wl ); } );
  testError<NC::Error::CalcError>( "negative",
                                   []( double wl ) { return 1.0 - wl; },
                                   constSlope( -1.0 ) );
  testError<NC::Error::CalcError>( "NaN",
                                   []( double wl ) { return wl > 3.0 ? std::numeric_limits<double>::quiet_NaN() : 1.0; },
                                   constSlope( 0.0 ) );
  testError<NC::Error::CalcError>( "infinite",
                                   []( double wl ) { return wl > 3.0 ? NC::kInfinity : 1.0; },
                                   constSlope( 0.0 ) );
  testError<NC::Error::CalcError>( "infinite slope",
                                   []( double ) { return 1.0; },
                                   []( double wl ) { return wl > 3.0 ? NC::kInfinity : 0.0; } );
  {
    NCX::IsoMacroXSTableParams p;
    p.tol = 0.0;
    testError<NC::Error::BadInput>( "tol=0", []( double ) { return 1.0; },
                                    constSlope( 0.0 ), {}, p );
    p = NCX::IsoMacroXSTableParams();
    p.wlmax = 1e-9;
    testError<NC::Error::BadInput>( "tiny wlmax", []( double ) { return 1.0; },
                                    constSlope( 0.0 ), {}, p );
  }

  //evalIsoMacroXSTable: pairs at discontinuities, cubic interpolation, and
  //clamped extrapolation:
  {
    const NC::VectD wl = { 0.0, 1.0, 2.0, 2.0, 3.0 };
    const NC::VectD xs = { 1.0, 2.0, 3.0, 1.0, 0.5 };
    const NC::VectD sl = { 1.0, 1.0, 1.0, -0.5, -0.5 };
    auto ev = [&]( double w )
    {
      return NCX::evalIsoMacroXSTable( wl.data(), xs.data(), sl.data(),
                                       wl.size(), w );
    };
    nc_assert_always( ev( -1.0 ) == 1.0 );
    nc_assert_always( ev( 0.5 ) == 1.5 );//cubic reduces to linear
    nc_assert_always( ev( 1.999 ) > 2.99 );
    nc_assert_always( ev( 2.0 ) == 1.0 );
    nc_assert_always( ev( 2.5 ) == 0.75 );
    nc_assert_always( ev( 3.5 ) == 0.25 );
    nc_assert_always( ev( 4.0 ) == 0.0 );//extrapolated: -0.0, clamped
    const NC::VectD sl2 = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    auto ev2 = [&]( double w )
    {
      return NCX::evalIsoMacroXSTable( wl.data(), xs.data(), sl2.data(),
                                       wl.size(), w );
    };
    nc_assert_always( ev2( 0.5 ) == 1.5 );//symmetric smoothstep
    nc_assert_always( ev2( 0.25 ) == 1.0 + 0.15625 );
    nc_assert_always( ev2( NC::kInfinity ) == 0.5 );
    std::cout << "evalIsoMacroXSTable: OK" << std::endl;
  }
  //reduceHermitePoints: the interpolation is never negative, even for a sample
  //with zeroes and sharp peaks:
  {
    NC::VectD x, y, slopes;
    for ( unsigned i = 0; i <= 2000; ++i ) {
      const double xi = 0.01 * i;
      x.push_back( xi );
      const double d = std::fmod( xi, 2.0 ) - 1.0;
      const bool peak = ( d > 0.0 && d < 0.2 );
      y.push_back( peak ? std::sin( 5.0 * NC::kPi * d ) * 10.0 : 0.0 );
      slopes.push_back( peak ? 50.0 * NC::kPi * std::cos( 5.0 * NC::kPi * d ) : 0.0 );
    }
    auto r = NCX::reduceHermitePoints( x, y, slopes, 1e-3 );
    NC::VectD tx, ty;
    for ( auto i : r.idx ) {
      tx.push_back( x[i] );
      ty.push_back( y[i] );
    }
    double ymin = 0.0;
    for ( auto w : NC::linspace( 0.0, 20.0, 200001 ) )
      ymin = std::min( ymin, NCX::evalIsoMacroXSTable( tx.data(), ty.data(),
                                                       r.slope.data(),
                                                       tx.size(), w ) );
    nc_assert_always( ymin >= 0.0 );
    std::cout << "reduceHermitePoints: " << r.idx.size() << " of " << x.size()
              << " points, never negative: OK" << std::endl;
  }
  return 0;
}
