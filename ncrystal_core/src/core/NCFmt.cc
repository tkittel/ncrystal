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

#include "NCrystal/core/NCFmt.hh"
namespace NC = NCrystal;

#include <streambuf>
#include <istream>
#include <locale>
#include <cerrno>
#include <cstdlib>
#include <clocale>
#include <cfenv>
#include <cstdio>

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    // Inspired by:
    // https://stackoverflow.com/questions/13059091/creating-an-input-stream-from-constant-memory/13059195#13059195
    struct nc_membuf : std::streambuf {
      //NB: base is memory buffer of given size, allocation/deallocation of
      //which should be handled by calling code.
      nc_membuf(char const* base, size_t size)
      {
        char* p = const_cast<char*>(base);
        this->setg(p, p, p + size);
      }
    };
    struct nc_imemstream : virtual nc_membuf, std::istream {
      nc_imemstream(char const* base, size_t size)
        : nc_membuf(base, size),
          std::istream(static_cast<std::streambuf*>(this))
      {
      }
    };

    namespace {
      static constexpr std::size_t denorm_str_maxsize = 64;

      Optional<double> raw_str2dbl_subnormal( const char * s_data,
                                              std::size_t s_size )
      {
        //Fallback for subnormal values, which some stream implementations
        //(e.g. libc++) fail to read: std::strtod, which rounds them correctly,
        //with '.' replaced by the decimal point of the C locale (as in
        //dbl2shortstr). Only accepted if all is consumed and the result is
        //nonzero and subnormal (true underflow and other failures stay
        //failures):
        if ( s_size > denorm_str_maxsize )
          return NullOpt;
        //Only plain decimal numbers (strtod also accepts e.g. "nan", "inf" and
        //hexadecimal numbers, which str2dbl rejects):
        for ( auto i : ncrange( s_size ) ) {
          const char c = s_data[i];
          if ( !( ( c >= '0' && c <= '9' ) || c == '.' || c == 'e'
                  || c == 'E' || c == '+' || c == '-' ) )
            return NullOpt;
        }
        const char * dp = std::localeconv()->decimal_point;
        const std::size_t ndp = std::strlen( dp );
        if ( !ndp || ndp > 8 )
          return NullOpt;
        char buf[denorm_str_maxsize * 8 + 1];
        std::size_t n = 0;
        for ( auto i : ncrange( s_size ) ) {
          if ( s_data[i] == '.' ) {
            std::memcpy( buf + n, dp, ndp );
            n += ndp;
          } else {
            buf[n++] = s_data[i];
          }
        }
        buf[n] = '\0';
        //strtod sets ERANGE and floating point exception flags (underflow,
        //overflow, e.g. for "1e400"), which must neither leak nor trap
        //(e.g. on i386 with traps enabled):
        const int saved_errno = errno;
        std::fenv_t fenv;
        const bool held = ( std::feholdexcept( &fenv ) == 0 );
        char * end = nullptr;
        const double v = std::strtod( buf, &end );
        if ( held )
          std::fesetenv( &fenv );
        errno = saved_errno;
        if ( end != buf + n || !std::isfinite( v ) || v == 0.0
             || !( std::fabs( v ) < std::numeric_limits<double>::min() ) )
          return NullOpt;
        return v;
      }

#ifndef NDEBUG
      static int dummy = []()
      {
        //Subnormals read back exactly (also on platforms where streams fail
        //for them, where raw_str2dbl relies on this):
        auto test = []( std::string s, double vref )
        {
          auto v = raw_str2dbl_subnormal( s.data(), s.size() );
          nc_assert_always( v.has_value() );
          nc_assert_always( v.value() == vref );
        };
        test("4.94066e-324",4.94066e-324);
        test("5e-320",5e-320);
        test("5e-324",5e-324);
        test("1.123456789123456789e-315",1.123456789123456789e-315);
        test("2.2250738585072009e-308",2.2250738585072009e-308);
        for ( std::uint64_t m : { std::uint64_t(1), std::uint64_t(3),
                                  std::uint64_t(123456789),
                                  ( std::uint64_t(1) << 52 ) - 1,
                                  ( std::uint64_t(1) << 51 ) + 12345 } ) {
          const double v = std::ldexp( static_cast<double>( m ), -1074 );
          char b[64];
          const int nb = std::snprintf( b, sizeof(b), "%.17g", v );
          nc_assert_always( nb > 0 && nb < 64 );
          for ( char * c = b; *c; ++c )
            if ( !( ( *c >= '0' && *c <= '9' ) || *c == '-' || *c == 'e' ) )
              *c = '.';//(locale independent)
          test( std::string( b ), v );
          test( std::string( "-" ) + b, -v );
        }
        //Not subnormal, or no value, is not accepted:
        for ( const char * bad : { "1e-400", "1e-300", "1.5", "abc", "1e-310x",
                                   "nan", "inf", "1e400", "0x1p-1070" } )
          nc_assert_always( !raw_str2dbl_subnormal( bad, std::strlen(bad) ).has_value() );
        return 1;
      }();
#endif

      Optional<double> raw_str2dbl_impl( const char * s_data,
                                         std::size_t s_size ) {
        //Using streams so we can specify the locale (TODO in c++17 we can
        //possibly use std::from_chars instead!). Using custom stream buffers to
        //reduce need for allocations:
        Optional<double> result;
        detail::nc_imemstream ss(s_data,s_size);
        ss.std::istream::imbue(std::locale::classic());
        double val;
        ss >> val;
        if ( ss.fail() ) {
          result = raw_str2dbl_subnormal( s_data, s_size );
        } else if ( ss.eof() ) {
          result = val;
        }
        return result;
      }
    }

    Optional<double> raw_str2dbl( const char * s_data, std::size_t s_size ) {
      return raw_str2dbl_impl(s_data,s_size);
    }
    Optional<std::int64_t> raw_str2int64( const char * s_data,
                                          std::size_t s_size ) {
      //Using streams so we can specify the locale (TODO in c++17 we can possibly
      //use std::from_chars instead!). Using custom stream buffers to reduce need
      //for allocations:
      detail::nc_imemstream ss(s_data,s_size);
      ss.std::istream::imbue(std::locale::classic());
      std::int64_t val;
      ss >> val;
      if ( !ss.fail() && ss.eof() )
        return val;
      else
        return NullOpt;
    }

    //Replaces the C locale's decimal-point string (possibly several bytes) in
    //snprintf output of a finite number by '.', returning the new length. The
    //output has only a sign, digits, the decimal point and "e+-dd", so any
    //other run of bytes is the decimal point. Not in NCFmt.hh, since it is
    //also needed in NCString.cc, where it is fwd declared.
    std::size_t normaliseDecimalPoint( char * buf, std::size_t n )
    {
      std::size_t j = 0;
      bool insep = false;
      for ( std::size_t i = 0; i < n; ++i ) {
        const char c = buf[i];
        if ( ( c >= '0' && c <= '9' ) || c == '-' || c == '+' || c == 'e' ) {
          buf[j++] = c;
          insep = false;
        } else if ( !insep ) {
          buf[j++] = '.';
          insep = true;
        }
      }
      return j;
    }

    //This one is not actually added to NCFMT.hh, since it is just needed in
    //NCString.cc, where it is fwd declared.
    Optional<std::uint64_t> raw_str2uint64( const char * s_data,
                                            std::size_t s_size ) {
      detail::nc_imemstream ss(s_data,s_size);
      ss.std::istream::imbue(std::locale::classic());
      std::uint64_t val;
      ss >> val;
      if ( !ss.fail() && ss.eof() ) {
        if ( val && s_data[0]=='-' )
          return NullOpt;//prevent negative values via underflow (but -0 is ok)
        return val;
      } else {
        return NullOpt;
      }
    }
  }
}

NC::ShortStrDbl NC::dbl2shortstr( double value, const char * fmtstr )
{
  //Consistent handling of special values (also means we won't have to worry
  //about nan!=nan issues below):
  if ( std::isnan( value ) )
    return ShortStrDbl::make("nan");
  if ( std::isinf( value ) )
    return value > 0.0 ? ShortStrDbl::make("inf") : ShortStrDbl::make("-inf");
  if ( value == 0.0 && !fmtstr )
    return ShortStrDbl::make("0");

  //printf's "%g" is really hard to replicate with ostreams, and the following
  //is anyway more efficient in terms of malloc's:
  //
  //TODO: When compiler supports it, and we move to C++17, we should likely
  //use std::to_char instead of these locale-dependent methods!

  //Use the delayed_init_t constructor of ShortStrDbl, to avoid a buffer copy
  //when fmtstr is given (this is complicated, but worth it in this case):
  char * buf;
  ShortStrDbl::size_type* buf_strlen_ptr;
  ShortStrDbl result( ShortStrDbl::delayed_init_t{}, &buf, &buf_strlen_ptr );
  ShortStrDbl::size_type& buf_strlen = *buf_strlen_ptr;

  auto dosnprintf = [&buf,&buf_strlen,value]( const char * fmtstr_ )
  {
    auto nwritten_strlen = std::snprintf(buf,ShortStrDbl::bufsize,fmtstr_,value);
    nc_assert( nwritten_strlen > 0 && std::int64_t(nwritten_strlen+1) <= std::int64_t(ShortStrDbl::bufsize) );
    nc_assert( buf[nwritten_strlen] == '\0' );
    buf_strlen = nwritten_strlen;
    nc_assert( (int)std::strlen(buf) == (int)buf_strlen );
  };

  if ( fmtstr ) {
    //honour request for specific fmt string:
    dosnprintf(fmtstr);
    //Protect against non-C locale creating commas instead of dots for decimal
    //separators:
    auto pcomma = std::memchr(buf,',',buf_strlen);
    if ( pcomma )
      *(char*)pcomma = '.';
  } else {
    //No specific format string requested, aim for lossless and "nice".

    constexpr const char* fmt_patterns[] = {
      "%.1g",  "%.2g",  "%.3g",  "%.4g",  "%.5g",
      "%.6g",  "%.7g",  "%.8g",  "%.9g", "%.10g",
      "%.11g", "%.12g", "%.13g", "%.14g", "%.15g",
      "%.16g", "%.17g", "%.18g", "%.19g", "%.20g",
      "%.21g", "%.22g", "%.23g", "%.24g", "%.25g",
      "%.26g", "%.27g", "%.28g", "%.29g", "%.30g",
    };

    unsigned nstrip = 0;
    double testval = value;
    while ( std::fabs(testval) < std::numeric_limits<double>::min() ) {
      ++nstrip;
      testval *= 10;
    }
    if ( nstrip > 0 )
      --nstrip;

    constexpr int idxfmt_std = std::numeric_limits<double>::max_digits10-1;
    const int idxfmt = std::max<int>(1,idxfmt_std-nstrip);
    const int idxfmtm2 = std::max<int>(0,(idxfmt_std-2)-nstrip);

    const char * fmt_full = fmt_patterns[idxfmt];
    const char * fmt_full_m2 = fmt_patterns[idxfmtm2];

    //First try with a bit less than maxdigits10, which can prevent some
    //originally nice looking values being presented in a less nice looking
    //fashion (e.g. 0.1 -> 0.09999999999999999). Round-trip check with strtod
    //on the raw output, which uses the same (C) locale as snprintf; if it
    //fails, fmt_full is exact anyway. Formats into a larger buffer, since the
    //locale's decimal point can have several bytes until normalised:
    char raw[64];
    auto rawsnprintf = [&raw,value]( const char * f )
    {
      const int n = std::snprintf( raw, sizeof(raw), f, value );
      nc_assert_always( n > 0 && n < static_cast<int>( sizeof(raw) ) );
      return static_cast<std::size_t>( n );
    };
    std::size_t n = rawsnprintf( fmt_full_m2 );
    const int saved_errno = errno;//strtod sets ERANGE for denormals
    if ( std::strtod( raw, nullptr ) != value )
      n = rawsnprintf( fmt_full );
    errno = saved_errno;
    n = detail::normaliseDecimalPoint( raw, n );
    nc_assert_always( n <= ShortStrDbl::max_strlen );
    std::memcpy( buf, raw, n );
    buf[n] = '\0';
    buf_strlen = static_cast<ShortStrDbl::size_type>( n );
#ifndef NDEBUG
    {
      //loss less check:
      auto backconv2 = detail::raw_str2dbl( buf, buf_strlen );
      nc_assert( backconv2.has_value() && backconv2.value() == value );
    }
#endif
  }
  return result;
}

std::pair<unsigned,unsigned> NC::detectSimpleRationalNumbers(double value)
{
  if (value<=0.0)
    return { 0, ( value == 0.0 ? 1 : 0 ) };
  if (value>=1.0) {
    if (value==1.0)
      return { 1, 1 };
    double intpart;
    if ( std::modf(value, &intpart) != 0.0
         || intpart >= (double)std::numeric_limits<unsigned>::max() )
      return {0,0};
    return { static_cast<unsigned>(intpart), unsigned(1) };
  }
  //ok, value is in (0,1).
  static std::map<uint64_t,std::pair<unsigned,unsigned>> s_cache = []()
  {
    std::map<uint64_t,std::pair<unsigned,unsigned>> result;
    for (unsigned b = 2; b <= 64; ++b) {
      for (unsigned a = 1; a < b; ++a) {
        uint64_t key = static_cast<uint64_t>((double(a)/b)*1e18);
        //NB: key fits in uint64_t since value of a/b is in (0,1)
        auto it = result.find(key);
        if (it==result.end())
          result[key] = { a,b };//only insert if not there (because if a1/b1 =
                                //a2/b2, we prefer the smallest b-value).
      }
    }
    return result;
  }();
  uint64_t key = static_cast<uint64_t>(value*1e18);
  auto it = s_cache.find(key);
  if ( it == s_cache.end() )
    return {0,0};
  return it->second;
}

void NC::detail::streamdblfmtg( std::ostream& os, double val )
{
  os << fmtg(val);
}

namespace NCRYSTAL_NAMESPACE {
  std::ostream& operator<<( std::ostream& os , const detail_FmtDblFrac& fd ) {
    auto ab = detectSimpleRationalNumbers(fd.val);
    if ( ab.second == 1 )
      os << ab.first;//integer
    else if ( ab.second == 0 )
      os << dbl2shortstr( fd.val, fd.fmtstr );
    else
      os << ab.first<<"/"<<ab.second;//fraction detected
    return os;
  }
}
