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

#include "NCrystal/internal/utils/NCNumFmt.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/query/NCQuery.hh"
#include <iostream>

//Tests of the number formatting for data files (NCNumFmt.hh): the text
//conventions, the digits for relative tolerances, idempotence, the
//tolerances, the fractions (also next to and near them), grids, and the
//"numfmt" query with its hex encoded values.

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace {

  //Value of a text, which may be a fraction "p/q":
  double parse( const std::string& t )
  {
    const auto pos = t.find( '/' );
    const NC::StrView sv( t );
    if ( pos == std::string::npos )
      return NC::str2dbl( sv );
    return ( NC::str2dbl( sv.substr( 0, pos ) )
             / NC::str2dbl( sv.substr( pos + 1 ) ) );
  }

  std::uint64_t s_lcg = 12345;
  double rnd()
  {
    s_lcg = s_lcg * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>( s_lcg >> 11 ) / 9007199254740992.0;
  }

  //Random values of various kinds (decimals with few digits, fractions with
  //and without noise, calculated values, and values over a wide range):
  double randomValue( unsigned kind )
  {
    switch ( kind % 6 ) {
    case 0:
      return rnd();
    case 1: {
      const double d = std::pow( 10.0, std::floor( rnd() * 8 ) + 1 );
      return std::floor( rnd() * d ) / d;
    }
    case 2: {
      const double q = std::floor( rnd() * 48 ) + 1;
      return std::floor( rnd() * q ) / q;
    }
    case 3: {
      const double q = std::floor( rnd() * 48 ) + 1;
      const double noise = std::pow( 10.0, -17.0 + 5.0 * rnd() );
      return ( std::floor( rnd() * q ) + 1.0 ) / q
        * ( 1.0 + ( rnd() < 0.5 ? -noise : noise ) );
    }
    case 4:
      return ( std::floor( rnd() * 1000 ) / 1000
               + std::floor( rnd() * 1000 ) / 1000 ) / 3.0;
    default:
      return ( rnd() < 0.5 ? -1.0 : 1.0 ) * rnd()
        * std::pow( 10.0, -20.0 + 40.0 * rnd() );
    }
  }

  void testConventions()
  {
    std::cout << "Conventions:" << std::endl;
    const double vals[] = { 0.25, 1.0/3, 0.3333, 300.0, 1e-5, 0.1 + 0.2,
                            123456.789, 2.5e17, -0.0, -0.125, 17.0/48,
                            1e15, 1e16, 0.0001, 4.0495, 90.0 };
    for ( double rtol : { 1e-14, 1e-6, 1e-3, 0.0 } ) {
      std::cout << "  rtol=" << NC::fmt( rtol ) << " ("
                << NC::numFmtDigits( rtol ) << " digits):";
      for ( double x : vals )
        std::cout << ' ' << NC::fmtRelTol( x, rtol, true, 40 );
      std::cout << std::endl;
    }
    REQUIRE( NC::fmtRelTol( 0.25, 1e-14, false ) == "0.25" );
    REQUIRE( NC::fmtRelTol( 0.25, 1e-14, true ) == ".25" );
    REQUIRE( NC::fmtRelTol( -0.25, 1e-14, true ) == "-.25" );
    REQUIRE( NC::fmtRelTol( 0.25, 1e-14, true, 40 ) == "1/4" );
    //(no fractions for values outside (0,1)):
    REQUIRE( NC::fmtRelTol( -0.25, 1e-14, true, 40 ) == "-.25" );
    REQUIRE( NC::fmtRelTol( 1.5, 1e-14, true, 40 ) == "1.5" );
    REQUIRE( NC::numFmtDigits( 1e-14 ) == 15 );
    REQUIRE( NC::numFmtDigits( 9.9e-15 ) == 17 );
    REQUIRE( NC::numFmtDigits( 1e-6 ) == 7 );
    REQUIRE( NC::numFmtDigits( 9.9e-7 ) == 8 );
    REQUIRE( NC::numFmtDigits( 1e-3 ) == 4 );
    REQUIRE( NC::numFmtDigits( 0.0 ) == 17 );
  }

  void testIdempotenceAndTolerance()
  {
    //Formatting the value read back from a text gives the same text, and the
    //value is within the tolerance 10^(1-ndigits) (cf. numFmtDigits; half of
    //it without fractions, and exact for 17 digits):
    std::size_t n = 0;
    for ( unsigned ndig : { 1u, 2u, 4u, 7u, 10u, 15u, 17u } ) {
      for ( unsigned maxden : { 0u, 40u } ) {
        const double rtol = ( ndig == 17 ? 0.0
                              : ( maxden ? 1.0 : 0.5 )
                              * std::pow( 10.0, 1.0 - ndig ) );
        for ( unsigned i = 0; i < 60000; ++i ) {
          const double x = randomValue( i );
          const std::string t = NC::fmtSigDigits( x, ndig, true, maxden );
          const double y = parse( t );
          REQUIRE( NC::fmtSigDigits( y, ndig, true, maxden ) == t );
          if ( ndig == 17 )
            REQUIRE( y == x || ( x == 0.0 && y == 0.0 ) );
          else
            REQUIRE( NC::ncabs( y - x )
                     <= rtol * NC::ncabs( x ) * ( 1.0 + 1e-12 ) );
          ++n;
        }
      }
    }
    std::cout << "Idempotence and tolerances: " << n << " values OK"
              << std::endl;
  }

  void testFractions()
  {
    //Every p/q (q<=40) is written as the reduced fraction, at any number of
    //digits (with 17 only for the double of p/q itself), and never for q>40
    //(17/48):
    std::size_t n = 0, nulp = 0;
    for ( unsigned q = 2; q <= 40; ++q ) {
      for ( unsigned p = 1; p < q; ++p ) {
        unsigned a = p, b = q;
        while ( b ) {
          const unsigned t = a % b;
          a = b;
          b = t;
        }
        if ( a != 1 )
          continue;
        const double x = static_cast<double>( p ) / q;
        std::ostringstream ss;
        ss << p << '/' << q;
        for ( unsigned ndig : { 4u, 7u, 15u, 17u } )
          REQUIRE( NC::fmtSigDigits( x, ndig, true, 40 ) == ss.str() );
        //One ulp away, the fraction is written exactly if the plain text is
        //that of p/q (most often, but not if p/q is near a rounding boundary,
        //e.g. 3/7=0.42857142857142855 and the double below it):
        const std::string plain = NC::fmtSigDigits( x, 15, true );
        for ( double nb : { std::nextafter( x, 1.0 ),
                            std::nextafter( x, 0.0 ) } ) {
          const bool same = ( NC::fmtSigDigits( nb, 15, true ) == plain );
          REQUIRE( ( NC::fmtSigDigits( nb, 15, true, 40 ) == ss.str() )
                   == same );
          if ( same )
            ++nulp;
          REQUIRE( NC::fmtSigDigits( nb, 17, true, 40 ) != ss.str() );
        }
        //Far outside the tolerance, no fraction:
        REQUIRE( NC::fmtSigDigits( x * ( 1.0 + 1e-12 ), 15, true, 40 )
                 != ss.str() );
        ++n;
      }
    }
    REQUIRE( NC::fmtSigDigits( 17.0 / 48, 15, true, 40 )
             == ".354166666666667" );
    REQUIRE( NC::fmtSigDigits( 1.0 / 3, 15, true, 0 ) == ".333333333333333" );
    REQUIRE( nulp > n );//(i.e. most neighbours)
    std::cout << "Fractions: all " << n << " reduced fractions p/q with q<=40"
              << " OK" << std::endl;
  }

  void testShortest()
  {
    //17 digits give the shortest text reading back exactly: it reads back
    //exactly and has as many significant digits as found by a brute force
    //search, for values with all exponents (also subnormal) and with few
    //digits:
    auto bruteForceDigits = []( double x )
    {
      char buf[64];
      for ( unsigned m = 1; m <= 17; ++m ) {
        std::snprintf( buf, sizeof(buf), "%.*e", static_cast<int>( m - 1 ), x );
        if ( std::strtod( buf, nullptr ) == x )
          return m;
      }
      return 17u;
    };
    auto sigDigits = []( const std::string& t )
    {
      std::string d;
      for ( char c : t.substr( 0, t.find( 'e' ) ) )
        if ( c >= '0' && c <= '9' )
          d += c;
      d.erase( 0, d.find_first_not_of( '0' ) );
      d.erase( d.find_last_not_of( '0' ) + 1 );
      return static_cast<unsigned>( d.size() );
    };
    std::size_t n = 0;
    for ( unsigned i = 0; i < 200000; ++i ) {
      std::uint64_t u;
      double x;
      do {
        s_lcg = s_lcg * 6364136223846793005ULL + 1442695040888963407ULL;
        u = s_lcg;
        if ( i % 4 == 1 )
          u &= 0x800fffffffffffffULL;//subnormal
        std::memcpy( &x, &u, sizeof(x) );
        if ( i % 4 == 2 ) {
          //few digits:
          x = std::round( x * 1e3 ) / 1e3;
        }
      } while ( !std::isfinite( x ) || x == 0.0 );
      const std::string t = NC::fmtSigDigits( x, 17, false );
      REQUIRE( NC::str2dbl( t ) == x );
      REQUIRE( sigDigits( t ) == bruteForceDigits( x ) );
      ++n;
    }
    REQUIRE( NC::fmtSigDigits( 5e-324, 17, false ) == "5e-324" );
    REQUIRE( NC::fmtSigDigits( -5e-324, 17, false ) == "-5e-324" );
    REQUIRE( NC::fmtSigDigits( 2.2250738585072014e-308, 17, false )
             == "2.2250738585072014e-308" );
    REQUIRE( NC::fmtSigDigits( 1e-320, 17, false ) == "1e-320" );
    std::cout << "Shortest texts: " << n << " values OK" << std::endl;
  }

  void testGrids()
  {
    //Vectors are formatted as single values, and grids keep the differences
    //between neighbours within rtol (also for neighbours much closer than
    //the tolerance), with values still within rtol:
    auto checkGrid = []( const NC::VectD& g, double rtol )
    {
      const auto res = NC::fmtRelTolVector( g, rtol, true, true );
      REQUIRE( res.size() == g.size() );
      NC::VectD y;
      for ( auto i : NC::ncrange( g.size() ) ) {
        y.push_back( parse( res.at(i) ) );
        REQUIRE( NC::ncabs( y.back() - g.at(i) )
                 <= rtol * NC::ncabs( g.at(i) ) );
      }
      for ( auto i : NC::ncrange( g.size() - 1 ) ) {
        const double dx = g.at(i+1) - g.at(i);
        REQUIRE( NC::ncabs( ( y.at(i+1) - y.at(i) ) - dx )
                 <= rtol * NC::ncabs( dx ) );
        REQUIRE( y.at(i+1) > y.at(i) );
      }
      return res;
    };
    const double rtols[] = { 0.0, 5e-14, 1e-9, 5e-6, 1e-3 };
    const NC::VectD g0 = { 1.0, 1.00001, 1.00002, 2.0, 3.0 };
    for ( double rtol : rtols ) {
      std::cout << "Grid with rtol=" << NC::fmt( rtol ) << ":";
      for ( auto& t : checkGrid( g0, rtol ) )
        std::cout << ' ' << t;
      std::cout << std::endl;
      //Not as grid: as single values:
      auto v = NC::fmtRelTolVector( g0, rtol, false, true );
      for ( auto i : NC::ncrange( g0.size() ) )
        REQUIRE( v.at(i) == NC::fmtRelTol( g0.at(i), rtol, true ) );
    }
    std::size_t ngrids = 0;
    for ( unsigned i = 0; i < 300; ++i ) {
      NC::VectD g;
      double x = rnd();
      const unsigned n = 2 + static_cast<unsigned>( rnd() * 99 );
      for ( unsigned j = 0; j < n; ++j ) {
        g.push_back( x );
        x += std::pow( 10.0, -12.0 * rnd() ) * x;
      }
      for ( double rtol : rtols ) {
        checkGrid( g, rtol );
        ++ngrids;
      }
    }
    std::cout << "Grids: " << ngrids << " random grids OK" << std::endl;
  }

  //Hex text of doubles as expected by the queries (8 bytes little-endian):
  std::string hexDoubles( const NC::VectD& v )
  {
    std::string res;
    for ( double x : v ) {
      std::uint64_t u;
      std::memcpy( &u, &x, sizeof(u) );
      for ( unsigned k = 0; k < 8; ++k ) {
        const unsigned b = static_cast<unsigned>( ( u >> ( 8 * k ) ) & 0xff );
        res += "0123456789abcdef"[b >> 4];
        res += "0123456789abcdef"[b & 0xf];
      }
    }
    return res;
  }

  std::string query( const std::string& hex, const char * rtol = "0",
                     const char * grid = "0" )
  {
    std::ostringstream os;
    NC::JSONQuery( os, { "util", "numfmt", rtol, grid, "1", "40", hex } );
    return os.str();
  }

  void testQuery()
  {
    //The hex encoding is exact, also for extreme values:
    const NC::VectD v = { -0.0, 0.1, -1.5, 1.0/3, 5e-324,
                          2.2250738585072014e-308, 1.7976931348623157e308,
                          -123456.789 };
    std::string expected = "\"";
    for ( auto x : v )
      expected += ( expected.size() > 1 ? " " : "" )
        + NC::fmtRelTol( x, 0.0, true, 40 );
    expected += '"';
    const std::string res = query( hexDoubles( v ) );
    std::cout << "Query: " << res << std::endl;
    REQUIRE( res == expected );
    //Upper case hex digits are fine:
    std::string hexu = hexDoubles( v );
    for ( auto& c : hexu )
      c = static_cast<char>( std::toupper( static_cast<unsigned char>( c ) ) );
    REQUIRE( query( hexu ) == expected );
    //Grids via the query:
    const NC::VectD g3 = { 1.0, 1.00001, 2.0 };
    std::string g = "\"";
    for ( auto& t : NC::fmtRelTolVector( g3, 1e-3, true, true, 40 ) )
      g += ( g.size() > 1 ? " " : "" ) + t;
    g += '"';
    REQUIRE( query( hexDoubles( g3 ), "1e-3", "1" ) == g );
    //Malformed input is rejected:
    auto bad = [&]( const std::string& hex )
    {
      try {
        query( hex );
      } catch ( NC::Error::BadInput& ) {
        return true;
      }
      return false;
    };
    const NC::VectD vok = { 0.5, 0.25 };
    const std::string h = hexDoubles( vok );
    REQUIRE( !bad( h ) );
    REQUIRE( bad( h.substr( 1 ) ) );//length not a multiple of 16
    REQUIRE( bad( h + "0" ) );
    REQUIRE( bad( "g" + h.substr( 1 ) ) );//not a hex digit
    REQUIRE( bad( h.substr( 0, 15 ) + " " + h.substr( 16 ) ) );
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const NC::VectD vinf = { 0.5, inf };
    const NC::VectD vnan = { nan };
    REQUIRE( bad( hexDoubles( vinf ) ) );
    REQUIRE( bad( hexDoubles( vnan ) ) );
    std::cout << "Query: malformed input rejected" << std::endl;
  }
}

int main()
{
  testConventions();
  testIdempotenceAndTolerance();
  testFractions();
  testShortest();
  testGrids();
  testQuery();
  return 0;
}
