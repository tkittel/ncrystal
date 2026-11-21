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
#include "NCrystal/internal/utils/NCMath.hh"
#include <cstdio>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {

    //10^(1-n) for n=1..15 (literals, so exactly the doubles of these decimal
    //values):
    constexpr double s_units[15] = { 1.0, 0.1, 1e-2, 1e-3, 1e-4, 1e-5,
                                     1e-6, 1e-7, 1e-8, 1e-9, 1e-10,
                                     1e-11, 1e-12, 1e-13, 1e-14 };

    //printf of a finite double, as text with the decimal point always '.',
    //and (if wantValue, otherwise NaN) the value of that text. The value is
    //from strtod on the raw printf output, which uses the same (current)
    //locale, and is much faster than locale independent stream based
    //conversions:
    struct PrintedDbl {
      std::string text;
      double value;
    };
    PrintedDbl printfDbl( const char * fmtstr, int prec, double x,
                          bool wantValue )
    {
      char buf[64];
      const int n = std::snprintf( buf, sizeof(buf), fmtstr, prec, x );
      nc_assert_always( n > 0 && n < static_cast<int>( sizeof(buf) ) );
      PrintedDbl res;
      res.value = ( wantValue ? std::strtod( buf, nullptr )
                    : std::numeric_limits<double>::quiet_NaN() );
      //Any character other than digits, signs and 'e' is the decimal point:
      for ( char * c = buf; c != buf + n; ++c )
        if ( !( ( *c >= '0' && *c <= '9' ) || *c == '-' || *c == '+'
                || *c == 'e' ) )
          *c = '.';
      res.text.assign( buf, static_cast<std::size_t>( n ) );
      return res;
    }

    std::string withoutDotZero( std::string s )
    {
      if ( s.size() >= 2 && s.compare( s.size() - 2, 2, ".0" ) == 0 )
        s.resize( s.size() - 2 );
      return s;
    }

    //Text like "-1.2345000e+02" (from printf's %e or %g) in Python repr
    //style: without trailing zeros of the digits, and in fixed notation for
    //decimal exponents -4..15 (otherwise like "1.5e+16"):
    std::string reprStyle( const std::string& e )
    {
      const auto epos = e.find( 'e' );
      nc_assert_always( epos != std::string::npos );
      const int e10 = std::atoi( e.c_str() + epos + 1 );
      const bool neg = ( e[0] == '-' );
      std::string digits;
      for ( auto i : ncrange( std::size_t( neg ? 1 : 0 ), epos ) )
        if ( e[i] != '.' )
          digits += e[i];
      while ( digits.size() > 1 && digits.back() == '0' )
        digits.pop_back();
      std::string r;
      if ( e10 < -4 || e10 >= 16 ) {
        r = digits.substr( 0, 1 );
        if ( digits.size() > 1 )
          r += "." + digits.substr( 1 );
        r += e.substr( epos );
      } else if ( e10 >= 0 ) {
        const std::size_t nint = static_cast<std::size_t>( e10 ) + 1;
        if ( digits.size() <= nint ) {
          r = digits + std::string( nint - digits.size(), '0' );
        } else {
          r = digits.substr( 0, nint ) + "." + digits.substr( nint );
        }
      } else {
        r = "0." + std::string( static_cast<std::size_t>( -e10 - 1 ), '0' )
          + digits;
      }
      return neg ? "-" + r : r;
    }

    //Like Python's repr: the shortest digits reading back exactly (cf.
    //reprStyle). For normal values, any decimal with at most 15 significant
    //digits survives a round trip through a double, so if the shortest has
    //at most 15 digits, it is the 15 digit text without trailing zeros, and
    //otherwise the 16 or 17 digit text. Subnormal values have less precision,
    //so all digits are tried:
    std::string reprDbl( double x )
    {
      nc_assert( std::isfinite( x ) && x != 0.0 );
      PrintedDbl e;
      const bool normal = ( ncabs( x ) >= std::numeric_limits<double>::min() );
      for ( int m = ( normal ? 15 : 1 ); m <= 17; ++m ) {
        e = printfDbl( "%.*e", m - 1, x, m < 17 );
        if ( e.value == x )
          break;
      }
      return reprStyle( e.text );
    }

    //The text without fractions, and (if wantValue) its value:
    PrintedDbl fmtPlain( double x, unsigned ndigits, bool wantValue )
    {
      if ( ndigits >= 16 )
        return { withoutDotZero( reprDbl( x ) ), x };
      PrintedDbl s = printfDbl( "%.*g", static_cast<int>( ndigits ), x,
                                wantValue );
      const auto epos = s.text.find( 'e' );
      if ( epos != std::string::npos ) {
        //Same value, sometimes shorter in fixed notation (e.g. "300"
        //instead of "3e+02"), which reprStyle uses for exponents 0..15:
        const int e10 = std::atoi( s.text.c_str() + epos + 1 );
        if ( e10 >= 0 && e10 < 16 ) {
          std::string alt = withoutDotZero( reprStyle( s.text ) );
          if ( alt.size() < s.text.size() )
            s.text.swap( alt );
        }
      }
      s.text = withoutDotZero( std::move( s.text ) );
      return s;
    }

    //fmtSigDigits, also giving the value of the text (if wantValue,
    //otherwise possibly NaN):
    PrintedDbl fmtSigDigitsImpl( double x, unsigned ndigits,
                                 bool stripLeadingZero,
                                 unsigned maxDenominator, bool wantValue )
    {
      nc_assert_always( ndigits >= 1 && ndigits <= 17 );
      if ( !std::isfinite( x ) )
        return { std::isnan( x ) ? "nan" : ( x > 0.0 ? "inf" : "-inf" ), x };
      if ( x == 0.0 )
        return { "0", 0.0 };
      const bool fractions = ( maxDenominator > 0 && x > 0.0 && x < 1.0 );
      PrintedDbl res = fmtPlain( x, ndigits, wantValue || fractions );
      if ( fractions ) {
        //Candidates from the value of the text, so values with the same text
        //get the same result. With n<16 digits, a candidate with the same
        //text is within half a unit of the last digit, and so within
        //xt*10^(1-n), of the value xt of the text (a cheap check first):
        const double xt = res.value;
        const double dmax = ( ndigits >= 16 ? 0.0 : xt * s_units[ndigits-1] );
        for ( auto q : ncrange( 2u, maxDenominator + 1 ) ) {
          const double p = std::round( xt * q );
          if ( !( p >= 1.0 && p < q ) )
            continue;
          const double v = p / q;
          if ( !( ncabs( v - xt ) <= dmax ) )
            continue;
          if ( ndigits >= 16
               || fmtPlain( v, ndigits, false ).text == res.text ) {
            std::ostringstream ss;
            ss << static_cast<unsigned>( p ) << '/' << q;
            return { ss.str(), v };
          }
        }
      }
      if ( stripLeadingZero ) {
        const std::string& t = res.text;
        if ( t.compare( 0, 2, "0." ) == 0 )
          res.text = t.substr( 1 );
        else if ( t.compare( 0, 3, "-0." ) == 0 )
          res.text = "-" + t.substr( 2 );
      }
      return res;
    }
  }
}

unsigned NC::numFmtDigits( double rtol )
{
  nc_assert_always( rtol >= 0.0 && std::isfinite( rtol ) );
  for ( unsigned n = 1; n <= 15; ++n )
    if ( s_units[n-1] <= rtol )
      return n;
  return 17;
}

std::string NC::fmtSigDigits( double x, unsigned ndigits,
                              bool stripLeadingZero, unsigned maxDenominator )
{
  return fmtSigDigitsImpl( x, ndigits, stripLeadingZero,
                           maxDenominator, false ).text;
}

std::string NC::fmtRelTol( double x, double rtol, bool stripLeadingZero,
                           unsigned maxDenominator )
{
  return fmtSigDigits( x, numFmtDigits( rtol ), stripLeadingZero,
                       maxDenominator );
}

NC::VectS NC::fmtRelTolVector( Span<const double> x, double rtol, bool isGrid,
                               bool stripLeadingZero,
                               unsigned maxDenominator )
{
  const unsigned nd0 = numFmtDigits( rtol );
  const std::size_t n = x.size();
  const bool grid = ( isGrid && rtol > 0.0 && n >= 2 );
  VectS res;
  res.reserve( n );
  VectD y;//values of the texts (for grids)
  y.reserve( grid ? n : 0 );
  for ( auto v : x ) {
    auto f = fmtSigDigitsImpl( v, nd0, stripLeadingZero, maxDenominator,
                               grid );
    res.push_back( std::move( f.text ) );
    if ( grid )
      y.push_back( f.value );
  }
  if ( !grid )
    return res;
  //More digits where neighbours need them, until no more can be added:
  std::vector<unsigned> nd( n, nd0 );
  std::vector<char> bad( n );
  while ( true ) {
    std::fill( bad.begin(), bad.end(), 0 );
    bool anybad = false;
    for ( auto i : ncrange( n - 1 ) ) {
      const double dx = x[i+1] - x[i];
      if ( ncabs( ( y[i+1] - y[i] ) - dx ) > rtol * ncabs( dx ) ) {
        for ( auto j : { i, i + 1 } ) {
          if ( nd[j] < 17 ) {
            bad[j] = 1;
            anybad = true;
          }
        }
      }
    }
    if ( !anybad )
      return res;
    for ( auto i : ncrange( n ) ) {
      if ( !bad[i] )
        continue;
      nd[i] = ( nd[i] >= 15 ? 17 : nd[i] + 1 );
      auto f = fmtSigDigitsImpl( x[i], nd[i], stripLeadingZero,
                                 maxDenominator, true );
      res[i] = std::move( f.text );
      y[i] = f.value;
    }
  }
}
