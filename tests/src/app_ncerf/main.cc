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

// Tests NC::ncerf/ncerfc (NCMath.hh): accuracy vs. mpmath reference values,
// loose agreement with std::erf/erfc, symmetry/monotonicity/special values,
// absence of FP exceptions, and a table of exact output bit patterns (the
// functions are meant to be bit-identical on all platforms).

#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>
#include <iomanip>
#include <cfenv>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  std::uint64_t bits( double x )
  {
    std::uint64_t b;
    std::memcpy( &b, &x, sizeof(b) );
    return b;
  }

  //Number of representable doubles between a and b (+0 and -0 are equal):
  std::uint64_t ulpdist( double a, double b )
  {
    auto key = []( double x )
    {
      const std::uint64_t bx = bits( x );
      const std::uint64_t sgn = std::uint64_t(1) << 63;
      return ( bx & sgn ) ? sgn - ( bx & ~sgn ) : sgn + bx;
    };
    const std::uint64_t ka( key(a) ), kb( key(b) );
    return ka > kb ? ka - kb : kb - ka;
  }

  struct RefVal { double x, erf, erfc; };

  //Generated with "sb_ncdev_gen_ncerf_approx --testref":
  const RefVal refvals[] = {
    //{ x, erf(x), erfc(x) }, mpmath values rounded to nearest:
    { 1e-300, 1.1283791670955126e-300, 1.0 },
    { 9.9999999999999995e-21, 1.1283791670955125e-20, 1.0 },
    { 1e-08, 1.1283791670955126e-08, 0.99999998871620832 },
    { 0.001, 0.0011283787909692365, 0.99887162120903072 },
    { 0.0625, 0.070431977722387074, 0.92956802227761293 },
    { 0.10000000000000001, 0.1124629160182849, 0.88753708398171505 },
    { 0.25, 0.27632639016823696, 0.7236736098317631 },
    { 0.375, 0.40411690943482231, 0.59588309056517774 },
    { 0.45000000000000001, 0.47548171978692366, 0.52451828021307634 },
    { 0.49999999999999994, 0.52049987781304652, 0.47950012218695354 },
    { 0.5, 0.52049987781304652, 0.47950012218695348 },
    { 0.50000000000000011, 0.52049987781304663, 0.47950012218695337 },
    { 0.59999999999999998, 0.60385609084792591, 0.39614390915207409 },
    { 0.75, 0.71115563365351508, 0.28884436634648486 },
    { 1.0, 0.84270079294971489, 0.15729920705028513 },
    { 1.25, 0.92290012825645829, 0.07709987174354177 },
    { 1.5, 0.96610514647531076, 0.033894853524689274 },
    { 2.0, 0.99532226501895271, 0.0046777349810472662 },
    { 2.5, 0.99959304798255499, 0.00040695201744495892 },
    { 3.0, 0.99997790950300136, 2.2090496998585441e-05 },
    { 3.5, 0.99999925690162761, 7.4309837234141278e-07 },
    { 4.0, 0.99999998458274209, 1.541725790028002e-08 },
    { 4.5, 0.99999999980338394, 1.9661604415428876e-10 },
    { 5.0, 0.99999999999846256, 1.5374597944280349e-12 },
    { 5.5, 0.99999999999999267, 7.3578479179743983e-15 },
    { 5.7999999999999998, 0.99999999999999978, 2.3555893751564417e-16 },
    { 5.9000000000000004, 0.99999999999999989, 7.1904097835504783e-17 },
    { 6.0, 1.0, 2.1519736712498913e-17 },
    { 7.0, 1.0, 4.1838256077794142e-23 },
    { 8.0, 1.0, 1.1224297172982926e-29 },
    { 10.0, 1.0, 2.0884875837625449e-45 },
    { 12.0, 1.0, 1.3562611692059042e-64 },
    { 15.0, 1.0, 7.2129941724512068e-100 },
    { 18.0, 1.0, 6.0823692318163992e-143 },
    { 20.0, 1.0, 5.3958656116079012e-176 },
    { 22.0, 1.0, 1.6219058609334726e-212 },
    { 24.0, 1.0, 1.6489825831519335e-252 },
    { 25.0, 1.0, 8.3001725711965228e-274 },
    { 26.0, 1.0, 5.6631924088561432e-296 },
    { 26.5, 1.0, 2.2109076642637343e-307 },
    { 26.550000000000001, 1.0, 1.5552026941135507e-308 },
    { 26.600000000000001, 1.0, 1.0885125885442269e-309 },
    { 27.0, 1.0, 5.2370464393526292e-319 },
    { 27.199999999999999, 1.0, 9.8813129168249309e-324 },
    { -0.001, -0.0011283787909692365, 1.0011283787909693 },
    { -0.0625, -0.070431977722387074, 1.0704319777223872 },
    { -0.10000000000000001, -0.1124629160182849, 1.1124629160182848 },
    { -0.25, -0.27632639016823696, 1.2763263901682369 },
    { -0.375, -0.40411690943482231, 1.4041169094348223 },
    { -0.45000000000000001, -0.47548171978692366, 1.4754817197869237 },
    { -0.49999999999999994, -0.52049987781304652, 1.5204998778130465 },
    { -0.5, -0.52049987781304652, 1.5204998778130465 },
    { -0.50000000000000011, -0.52049987781304663, 1.5204998778130467 },
    { -0.59999999999999998, -0.60385609084792591, 1.6038560908479258 },
    { -0.75, -0.71115563365351508, 1.7111556336535152 },
    { -1.0, -0.84270079294971489, 1.8427007929497148 },
    { -1.25, -0.92290012825645829, 1.9229001282564582 },
    { -1.5, -0.96610514647531076, 1.9661051464753108 },
    { -2.0, -0.99532226501895271, 1.9953222650189528 },
    { -2.5, -0.99959304798255499, 1.999593047982555 },
    { -3.0, -0.99997790950300136, 1.9999779095030015 },
    { -3.5, -0.99999925690162761, 1.9999992569016276 },
    { -4.0, -0.99999998458274209, 1.9999999845827421 },
    { -4.5, -0.99999999980338394, 1.9999999998033839 },
    { -5.0, -0.99999999999846256, 1.9999999999984626 },
    { -5.5, -0.99999999999999267, 1.9999999999999927 },
    { -5.7999999999999998, -0.99999999999999978, 1.9999999999999998 },
    { -5.9000000000000004, -0.99999999999999989, 2.0 },
    { -6.0, -1.0, 2.0 },
  };

  void testRefValues()
  {
    std::uint64_t worst_erf(0), worst_erfc(0);
    for ( auto& r : refvals ) {
      const std::uint64_t de = ulpdist( NC::ncerf( r.x ), r.erf );
      const std::uint64_t dc = ulpdist( NC::ncerfc( r.x ), r.erfc );
      worst_erf = std::max( worst_erf, de );
      worst_erfc = std::max( worst_erfc, dc );
    }
    std::cout << "Reference values: " << sizeof(refvals)/sizeof(RefVal)
              << " points, worst ULP: erf "
              << worst_erf << ", erfc " << worst_erfc << std::endl;
    REQUIRE( worst_erf <= 2 );
    REQUIRE( worst_erfc <= 2 );
  }

  void testVsStd()
  {
    //Loose comparison with the platform libm (whose accuracy varies):
    const std::uint64_t tol = 32;
    unsigned n(0);
    for ( double x = -7.0; x < 28.0; x += 0.000731 ) {
      const double e = NC::ncerf( x );
      const double ec = NC::ncerfc( x );
      const double se = std::erf( x );
      const double sec = std::erfc( x );
      if ( std::fabs( se ) > 1e-250 )
        REQUIRE( ulpdist( e, se ) <= tol );
      if ( std::fabs( sec ) > 1e-250 )
        REQUIRE( ulpdist( ec, sec ) <= tol );
      ++n;
    }
    for ( int i = -820; i < 0; ++i ) {//|results|>1e-250
      const double x = std::ldexp( 1.2345, i );
      REQUIRE( ulpdist( NC::ncerf( x ), std::erf( x ) ) <= tol );
      REQUIRE( ulpdist( NC::ncerf( -x ), std::erf( -x ) ) <= tol );
      REQUIRE( ulpdist( NC::ncerfc( x ), std::erfc( x ) ) <= tol );
      REQUIRE( ulpdist( NC::ncerfc( -x ), std::erfc( -x ) ) <= tol );
      n += 2;
    }
    std::cout << "Compared with std::erf/erfc at " << n << " points"
              << std::endl;
  }

  void checkMonotonic( double x, double& prev_erf, double& prev_erfc )
  {
    const double e = NC::ncerf( x );
    const double ec = NC::ncerfc( x );
    if ( !( e >= prev_erf ) || !( ec <= prev_erfc ) )
      std::cout << "Monotonicity violated at x="
                << NC::fmt( x, "%.17g" ) << std::endl;
    REQUIRE( e >= prev_erf );
    REQUIRE( ec <= prev_erfc );
    prev_erf = e;
    prev_erfc = ec;
  }

  void testProperties()
  {
    unsigned n(0);
    //Symmetries and erf+erfc=1 on a sweep:
    for ( double x = 0.0; x < 30.0; x += 0.000913 ) {
      const double e = NC::ncerf( x );
      const double ec = NC::ncerfc( x );
      REQUIRE( NC::ncerf( -x ) == -e );
      REQUIRE( e >= 0.0 && e <= 1.0 );
      REQUIRE( ec >= 0.0 && ec <= 1.0 );
      REQUIRE( std::fabs( ( e + ec ) - 1.0 ) <= 3e-16 );
      REQUIRE( std::fabs( NC::ncerfc( -x ) - ( 2.0 - ec ) ) <= 5e-16 );
      ++n;
    }
    //Monotonicity on a sweep, and over consecutive doubles around each
    //boundary between approximations (+-1/2, +-x1 and x0):
    double pe( -1.0 ), pec( 2.0 );
    for ( double x = -7.0; x < 30.0; x += 0.000377 ) {
      checkMonotonic( x, pe, pec );
      ++n;
    }
    //Region boundaries (+-1/2, +-x1, x0) and windows elsewhere:
    std::vector<double> starts = { -5.9215871957945074, -0.5, 0.0, 0.5,
                                   5.9215871957945074, 27.226017111108366 };
    for ( int i = 0; i < 100; ++i )
      starts.push_back( -6.1 + i * 0.3337 );
    for ( auto edge : starts ) {
      double x = edge;
      for ( int i = 0; i < 1000; ++i )
        x = std::nextafter( x, -NC::kInfinity );
      pe = -1.0;
      pec = 2.0;
      for ( int i = 0; i < 2000; ++i ) {
        checkMonotonic( x, pe, pec );
        x = std::nextafter( x, NC::kInfinity );
        ++n;
      }
    }
    std::cout << "Checked properties at " << n << " points" << std::endl;
  }

  void testSpecialValues()
  {
    const double inf = NC::kInfinity;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::feclearexcept( FE_ALL_EXCEPT );
    REQUIRE( bits( NC::ncerf( 0.0 ) ) == bits( 0.0 ) );
    REQUIRE( bits( NC::ncerf( -0.0 ) ) == bits( -0.0 ) );
    REQUIRE( NC::ncerfc( 0.0 ) == 1.0 );
    REQUIRE( NC::ncerfc( -0.0 ) == 1.0 );
    REQUIRE( NC::ncerf( inf ) == 1.0 );
    REQUIRE( NC::ncerf( -inf ) == -1.0 );
    REQUIRE( NC::ncerfc( inf ) == 0.0 );
    REQUIRE( NC::ncerfc( -inf ) == 2.0 );
    REQUIRE( NC::ncerf( 1e300 ) == 1.0 );
    REQUIRE( NC::ncerfc( 1e300 ) == 0.0 );
    REQUIRE( NC::ncerfc( -1e300 ) == 2.0 );
    REQUIRE( NC::ncerf( std::numeric_limits<double>::max() ) == 1.0 );
    REQUIRE( NC::ncerfc( -std::numeric_limits<double>::max() ) == 2.0 );
    REQUIRE( std::isnan( NC::ncerf( nan ) ) );
    REQUIRE( std::isnan( NC::ncerfc( nan ) ) );
    REQUIRE( std::isnan( NC::ncerf( -nan ) ) );
    REQUIRE( std::isnan( NC::ncerfc( -nan ) ) );
    const double dmin = std::numeric_limits<double>::denorm_min();
    REQUIRE( NC::ncerf( dmin ) == dmin );//c0*dmin rounds to dmin
    REQUIRE( NC::ncerf( -dmin ) == -dmin );
    REQUIRE( NC::ncerfc( dmin ) == 1.0 );
    //Underflow into subnormals and to 0 around x0:
    const double x0 = 27.226017111108366;
    REQUIRE( NC::ncerfc( x0 ) == 0.0 );
    REQUIRE( NC::ncerfc( std::nextafter( x0, 0.0 ) ) == dmin );
    REQUIRE( NC::ncerfc( 26.6 ) < std::numeric_limits<double>::min() );
    REQUIRE( NC::ncerfc( 26.6 ) > 0.0 );
    //No overflow/invalid/divbyzero raised anywhere (underflow/inexact are
    //expected and fine):
    for ( double x = -40.0; x < 40.0; x += 0.0123 ) {
      (void)NC::ncerf( x );
      (void)NC::ncerfc( x );
    }
    REQUIRE( !std::fetestexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW ) );
    std::cout << "Special values and FP exceptions OK" << std::endl;
  }

  void printBitPatterns()
  {
    //Exact output bits at points in all regions (catches any platform
    //dependency):
    const double xs[] = { 1e-300, -3e-9, 0.123, -0.25, 0.4999, 0.5,
                          -0.5, 0.7, 1.0, -1.3, 1.9, 2.6, -3.3, 4.4, 5.5,
                          -5.9, 6.5, 9.1, 13.7, 19.3, 24.4, 26.1, 26.9,
                          27.22 };
    std::cout << "Bit patterns of ncerf(x) and ncerfc(x):" << std::endl;
    for ( auto x : xs ) {
      std::cout << "  x=" << std::setw(7) << x << " : 0x"
                << std::hex << std::setfill('0')
                << std::setw(16) << bits( NC::ncerf( x ) ) << " 0x"
                << std::setw(16) << bits( NC::ncerfc( x ) )
                << std::dec << std::setfill(' ') << std::endl;
    }
  }

}

int main()
{
  testRefValues();
  testVsStd();
  testProperties();
  testSpecialValues();
  printBitPatterns();
  return 0;
}
