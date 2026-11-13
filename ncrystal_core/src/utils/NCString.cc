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

#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/core/NCFmt.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <istream>
#include <iomanip>
#include <locale>
namespace NC = NCrystal;

std::string NC::displayCharSafeQuoted( char ch_raw, char quote_char )
{
  signed char ch = *reinterpret_cast<signed char*>(&ch_raw);
  std::ostringstream os;
  os << quote_char;
  if ( ch == quote_char )
    os << '\\';//e.g. """ -> "\""
  if ( ch >= ' ' && ch < '\x7f' ) {//' ' = 32, 0x7f=127
    os << ch;
  } else {
    os << "\\x";
    int ch_val(ch);
    if ( ch_val < 10 )
      os << '0';
    os << ch_val;
  }
  os << quote_char;
  return os.str();
}

void NC::trim( std::string& input )
{
  char * cB = &input[0];
  const char * cE = cB + input.size();
  const char * c = cB;
  while (c!=cE&&(*c==' '||*c=='\t'||*c=='\n'||*c=='\r'))
    ++c;
  if (c==cE) {
    input.clear();
    return;
  }
  const char * clast = cE-1;
  while (clast>c&&(*clast==' '||*clast=='\t'||*clast=='\n'||*clast=='\r'))
    --clast;
  nc_assert(clast>=c);
  std::size_t len = (clast + 1) - c;
  if (c!=cB) {
    for (std::size_t i = 0; i < len; ++i)
      *(cB + i) = *(c + i);
  }
  if (len<input.size())
    input.resize(len);
}

NC::VectS NC::split2( const std::string& input, std::size_t maxsplit,  char sep )
{
  VectS output;
  split(output,input,maxsplit,sep);
  return output;
}

void NC::split(NC::VectS& output, const std::string& input, std::size_t maxsplit, char sep )
{
  bool keep_empty = (sep!=0);
  //Split input on whitespace (" \t\n") with minimal memory allocations and
  //copying, discarding whitespace which is either leading and trailing or
  //repeating.
  nc_assert(bool(sep!='\0')==bool(sep));
  output.clear();

  if (input.empty()) {
    if (sep)
      output.emplace_back();
    return;
  }

  if (sep&&input[0]==sep)
    output.emplace_back();
  const char * c = &input[0];
  const char * cE = c + input.size();
  const char * partbegin = 0;
  while (true) {
    if ( maxsplit && output.size() == maxsplit ) {
      output.emplace_back(c);
      return;
    }
    if (c==cE||(sep?*c==sep:(*c==' '||*c=='\t'||*c=='\n'||*c=='\r'))) {
      if (partbegin) {
        if (keep_empty||c>partbegin) {
          output.emplace_back(partbegin,c-partbegin);
        }
        partbegin= keep_empty ? c+1 : 0;
      }
      if (c==cE)
        return;
    } else if (!partbegin) {
      partbegin = c;
    }
    ++c;
  }
}

bool NC::startswith(const std::string& str, const std::string& substr)
{
  return str.size()>=substr.size() && str.compare(0, substr.size(), substr) == 0;
}

bool NC::endswith(const std::string& str, const std::string& substr)
{
  return str.size()>=substr.size() && str.compare(str.size()-substr.size(), substr.size(), substr) == 0;
}

bool NC::contains_any(const std::string& haystack, const std::string& needles)
{
  const char * cB = &needles[0];
  const char * cE = cB + needles.size();
  for (const char * c = cB; c!=cE; ++c)
    if (haystack.find(*c) != std::string::npos)
      return true;
  return false;
}

double NC::str2dbl( StrView sv, const char * errmsg)
{
  double result;
  if ( !safe_str2dbl(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not a double")<<": \""<<sv<<"\"");
  return result;
}

int NC::str2int( StrView sv, const char * errmsg)
{
  int result;
  if ( !safe_str2int(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not an integer")<<": \""<<sv<<"\"");
  return result;
}

std::int32_t NC::str2int32( StrView sv, const char * errmsg)
{
  std::int32_t result;
  if ( !safe_str2int(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not a 32bit signed integer")<<": \""<<sv<<"\"");
  return result;
}

std::int64_t NC::str2int64( StrView sv, const char * errmsg)
{
  std::int64_t result;
  if ( !safe_str2int(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not a 64bit signed integer")<<": \""<<sv<<"\"");
  return result;
}

std::uint64_t NC::str2uint64( StrView sv, const char * errmsg)
{
  std::uint64_t result;
  if ( !safe_str2uint(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not a 64bit unsigned integer")<<": \""<<sv<<"\"");
  return result;
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Shared by safe_str2dbl and Str2DblConverter, with rawConv(data,size)
    //giving the same result as detail::raw_str2dbl:
    template<class TRawConv>
    bool safe_str2dbl_impl( StrView s, double& result, TRawConv&& rawConv )
    {
      const auto n_size = s.size();
      if ( n_size==0 || isWhiteSpace( s.front() ) || isWhiteSpace( s.back() ) )
        return false;

      {
        auto try_conv = rawConv( s.data_begin(), n_size );
        if ( try_conv.has_value() ) {
          result = try_conv.value();
          return true;
        }
      }

      //Did not work, check for "+-inf"/"+-INF"/"nan"/"NaN","NAN" for
      //cross-platform consistency:
      if ( n_size == 3 ) {
        if ( isOneOf(s,"inf","INF") ) {
          result = kInfinity;
          return true;
        }
        if ( isOneOf(s,"nan","NAN","NaN") ) {
          static_assert(std::numeric_limits<double>::has_quiet_NaN,"");
          result = std::numeric_limits<double>::quiet_NaN();
          return true;
        }
      }
      if ( n_size == 4 ) {
        if ( isOneOf(s,"+inf","+INF") ) {
          result = kInfinity;
          return true;
        }
        if ( isOneOf(s,"-inf","-INF") ) {
          result = -kInfinity;
          return true;
        }
      }

      //Syntax error:
      return false;
    }

    //Stream reading from a memory buffer which can be changed between reads:
    struct ResettableMemBuf : std::streambuf {
      void reset( const char * data, std::size_t size )
      {
        char* p = const_cast<char*>(data);
        this->setg( p, p, p + size );
      }
    };
  }
}

bool NC::safe_str2dbl( StrView s, double& result )
{
  return safe_str2dbl_impl( s, result,
                            []( const char * d, std::size_t n )
                            { return detail::raw_str2dbl( d, n ); } );
}

//FIXME: Most of the remaining cost is the stream extraction itself (copying
//characters to a temporary string before strtod). The platform's own
//locale-independent strtod_l (_strtod_l on MSVC) would avoid it, if it gives
//identical results on all platforms. Also consider using this class in other
//readers of large text files (e.g. .laz/.lau), and for str2int.
struct NC::Str2DblConverter::Impl {
  ResettableMemBuf buf;
  std::istream is;
  Impl() : is( &buf ) { is.imbue( std::locale::classic() ); }
  Optional<double> rawConvert( const char * d, std::size_t n )
  {
    //Plain values are read as in detail::raw_str2dbl, but with the stream
    //reused. Anything else (syntax errors, denormals, ...) is rare, and is left
    //to detail::raw_str2dbl itself, to keep its handling in one place:
    buf.reset( d, n );
    is.clear();
    double val;
    is >> val;
    if ( !is.fail() && is.eof() )
      return val;
    return detail::raw_str2dbl( d, n );
  }
};

NC::Str2DblConverter::Str2DblConverter() = default;
NC::Str2DblConverter::~Str2DblConverter() = default;
NC::Str2DblConverter::Str2DblConverter( Str2DblConverter&& ) noexcept = default;
NC::Str2DblConverter&
NC::Str2DblConverter::operator=( Str2DblConverter&& ) noexcept = default;

NC::Optional<double> NC::Str2DblConverter::tryConvert( StrView s )
{
  auto& impl = *m_impl;
  double result;
  if ( !safe_str2dbl_impl( s, result,
                           [&impl]( const char * d, std::size_t n )
                           { return impl.rawConvert( d, n ); } ) )
    return NullOpt;
  return result;
}

double NC::Str2DblConverter::convert( StrView sv, const char * errmsg )
{
  auto result = tryConvert( sv );
  if ( !result.has_value() )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not a double")<<": \""<<sv<<"\"");
  return result.value();
}

bool NC::safe_str2int( StrView s, int32_t& result )
{
  //For robustness/simplicity simply parse as 64 bit value and check range.
  std::int64_t val64;
  static constexpr std::int64_t range_min
    = static_cast<std::int64_t>( std::numeric_limits<std::int32_t>::lowest() );
  static constexpr std::int64_t range_max
    = static_cast<std::int64_t>( std::numeric_limits<std::int32_t>::max() );
  if ( !safe_str2int(s,val64) || val64 > range_max || val64 < range_min )
    return false;
  result = static_cast<std::int32_t>(val64);
  return true;
}

bool NC::safe_str2int( StrView s, std::int64_t& result )
{
  const auto n_size = s.size();
  if ( n_size==0 || isWhiteSpace( s.front() ) || isWhiteSpace( s.back() ) )
    return false;
  auto try_conv = detail::raw_str2int64( s.data_begin(), n_size );
  if ( try_conv.has_value() ) {
    result = try_conv.value();
    return true;
  } else {
    return false;
  }
}
namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    //Fwd declared from NCFmt.cc:
    Optional<std::uint64_t> raw_str2uint64( const char * s_data,
                                            std::size_t s_size );
    std::size_t normaliseDecimalPoint( char * buf, std::size_t n );
  }
}

bool NC::safe_str2uint( StrView s, std::uint64_t& result )
{
  const auto n_size = s.size();
  if ( n_size==0 || isWhiteSpace( s.front() ) || isWhiteSpace( s.back() ) )
    return false;
  auto try_conv = detail::raw_str2uint64( s.data_begin(), n_size );
  if ( try_conv.has_value() ) {
    result = try_conv.value();
    return true;
  } else {
    return false;
  }
}

bool NC::contains_only(const std::string& haystack, const std::string& needles)
{
  const char * c = &haystack[0];
  const char * cE = c + haystack.size();
  for (;c!=cE;++c) {
    if (!contains(needles,*c))
      return false;
  }
  return true;
}

void NC::strreplace(std::string& str, const std::string& oldtxt, const std::string& newtxt)
{
  if( oldtxt.empty() )
    return;
  std::size_t offset = 0;
  while( ( offset = str.find(oldtxt, offset) ) != std::string::npos ) {
    str.replace( offset, oldtxt.length(), newtxt );
    offset += newtxt.size();
  }
}

std::string NC::joinstr( const Span<const StrView>& parts, StrView sep)
{
  const std::size_t n = parts.size();
  if ( n < 2 )
    return n ? parts.at(0).to_string() : std::string();

  std::size_t newsize( (n-1)*sep.size() );
  for (const auto& p : parts)
    newsize += p.size();
  std::string tmp;
  tmp.reserve(newsize);
  tmp.append(parts.front().data_begin(),parts.front().size());
  for (std::size_t i = 1; i < n; ++i) {
    tmp.append( sep.data_begin(), sep.size() );
    tmp.append( parts[i].data_begin(), parts[i].size() );
  }
  nc_assert(tmp.size()==newsize);
  return tmp;
}

unsigned NC::countTrailingDigits( const std::string& ss )
{
  auto nn = ss.size();
  nc_assert_always(static_cast<std::uint64_t>(nn)
                   <static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  int n = static_cast<int>(nn);

  int nTrailingDigits(0);
  while ( nTrailingDigits < n ) {
    int i = n-(nTrailingDigits+1);
    nc_assert(i>=0);
    auto c = ss.at( i );
    if ( 'c' < '0' || c > '9' )
      break;
    ++nTrailingDigits;
  }
  return static_cast<unsigned>(nTrailingDigits);
}

NC::PairSS NC::decomposeStrWithTrailingDigits( const std::string& ss )
{
  unsigned nTrailingDigits = countTrailingDigits(ss);
  if (nTrailingDigits==0)
    return {ss,std::string()};
  auto nn = ss.size() - nTrailingDigits;
  return { ss.substr(0,nn), ss.substr(nn) };
}

//FIXME: Move getenv stuff to core, and encapsulate this windows include better.
#if ( defined(_WIN32) || defined(WIN32) )
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif

namespace NCRYSTAL_NAMESPACE {

#if ( defined(_WIN32) || defined(WIN32) )
  //Forward declarations of the UTF-8<->UTF-16 conversion helpers defined in
  //the NCFileUtilsWin.cc:
  namespace winstr_details {
    std::wstring winimpl_str2wstr( const std::string& );
    std::string winimpl_wstr2str( const std::wstring& );
  }
#endif

  namespace {

#if ( defined(_WIN32) || defined(WIN32) )
    //Deliberately using the raw Win32 GetEnvironmentVariableW (wide/UTF-16)
    //API here rather than std::getenv: the CRT's own environment table
    //(which std::getenv reads from) is populated once per CRT module
    //instance and only kept in sync by that *same* instance's own
    //putenv/_putenv_s calls - it does not pick up SetEnvironmentVariable
    //changes made by a different CRT instance (confirmed to matter in
    //practice: a Debug-built NCrystal.dll loaded into a Release-built
    //Python interpreter never saw an env var set via Python's os.environ,
    //since Python's setter and this DLL's std::getenv end up backed by
    //two different CRT modules there). GetEnvironmentVariable instead
    //always reads the single, true, OS-maintained per-process value,
    //independent of which CRT (if any) last wrote it.
    //
    //Using GetEnvironmentVariableW for unicode support:
    Optional<std::string> platform_getenv( const char* name )
    {
      const std::wstring name_w = winstr_details::winimpl_str2wstr( name );
      const DWORD needed = GetEnvironmentVariableW( name_w.c_str(), nullptr, 0 );
      if ( needed == 0 )
        return NullOpt;//unset (or, rarely, a genuine API error - treated the same)
      std::wstring buf( static_cast<std::size_t>(needed), L'\0' );
      const DWORD written = GetEnvironmentVariableW( name_w.c_str(), &buf[0],
                                                      static_cast<DWORD>(buf.size()) );
      if ( written == 0 )
        return NullOpt;
      buf.resize( written );
      return winstr_details::winimpl_wstr2str( buf );
    }
#else
    Optional<std::string> platform_getenv( const char* name )
    {
      const char * v = std::getenv( name );
      if ( !v )
        return NullOpt;
      return std::string(v);
    }
#endif

    struct GetEnvResult {
      std::string varname;
      Optional<std::string> val;
    };
    GetEnvResult raw_getenv( std::string& v )
    {
      GetEnvResult res;
      //NB: keep platform_getenv here (not std::getenv): on Windows it
      //reads via the wide Win32 API, so changes made from a different
      //CRT module in the same process are seen:
      res.varname = ncgetenv_varname( v );
      res.val = platform_getenv( res.varname.c_str() );
      return res;
    }
  }
}

std::string NC::ncgetenv_varname( std::string v )
{
  nc_assert(!v.empty());
  nc_assert(!startswith(v,"NCRYSTAL_"));//common mistake
  std::string res;
#if defined(NCRYSTAL_NAMESPACE_PROTECTION) && defined( NCRYSTAL_NAMESPACED_ENVVARS )
  res.reserve( 64 );
  res = "NCRYSTAL";
  res += upperCase(ncrystal_xstr(NCRYSTAL_NAMESPACE_PROTECTION));
  res += '_';
#else
  res = "NCRYSTAL_";
#endif
  res += v;
  return res;
}

//Common access to environment variables (unset and empty vars both return
//empty strings / 0.0 / 0, depending on which function is used):
std::string NC::ncgetenv(std::string v, std::string defval)
{
  auto res = raw_getenv(v);
  return res.val.has_value() ? res.val.value() : defval;
}

double NC::ncgetenv_dbl(std::string v, double defval )
{
  auto res = raw_getenv(v);
  if (!res.val.has_value())
    return defval;
  double result;
  if ( !safe_str2dbl(res.val.value(), result ) )
    NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "
                    <<res.varname
                    <<" (expected a floating point number but got \""
                    <<res.val.value()<<"\").");
  return result;
}

int NC::ncgetenv_int(std::string v, int defval )
{
  auto res = raw_getenv(v);
  if (!res.val.has_value())
    return defval;
  int result;
  if ( !safe_str2int(res.val.value(), result ) )
    NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<res.varname
                    <<" (expected an integral number but got \""<<res.val.value()<<"\").");
  return result;
}

std::int64_t NC::ncgetenv_int64(std::string v, std::int64_t defval )
{
  auto res = raw_getenv(v);
  if (!res.val.has_value())
    return defval;
  std::int64_t result;
  if ( !safe_str2int(res.val.value(), result ) )
    NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<res.varname
                    <<" (expected an integral number but got \""<<res.val.value()<<"\").");
  return result;
}


bool NC::ncgetenv_bool(std::string v)
{
  auto res = raw_getenv(v);
  if (!res.val.has_value())
    return false;
  std::string evs(res.val.value());
  if (evs.size()==1) {
    if (evs[0]=='0')
      return false;
    if (evs[0]=='1')
      return true;
  }
  NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<res.varname
                  <<" (expected a Boolean value, \"0\" or \"1\", but got \""<<evs<<"\").");
}

std::string NC::bytes2hexstr(const std::vector<uint8_t>& v) {
  const char hexchars[] = "0123456789abcdef";
  std::ostringstream ss;
  for (auto e : v)
    ss << hexchars[ e >> 4 ] << hexchars[ e & 0x0F ];
  return ss.str();
}

std::vector<uint8_t> NC::hexstr2bytes(const std::string& v) {
  std::vector<uint8_t> res;
  auto hex2val = [](unsigned c) {
    unsigned val(99);
    constexpr unsigned val_a = static_cast<unsigned>('a');
    constexpr unsigned val_A = static_cast<unsigned>('A');
    constexpr unsigned val_0 = static_cast<unsigned>('0');
    static_assert( val_0 < val_A && val_A < val_a,"" );
    if ( c >= val_a )
      val = 10 + ( c - val_a );
    else if ( c >= val_A)
      val = 10 + ( c - val_A );
    else if ( c >= val_0)
      val = ( c - val_0 );
    if ( val > 15 )
      NCRYSTAL_THROW2(BadInput, "Invalid character encountered in hex string: "<<c<<" (numeric value)");
    return val;
  };

  auto it = v.begin();
  auto itE = v.end();
  const auto vlen = v.size();
  if ( vlen%2 != 0 ) {
    //0 pad, e.g. interpret (0x)"A" as (0x)"0A".
    res.reserve((vlen+1)/2);
    res.push_back(hex2val(*it++));
  } else {
    res.reserve( vlen/2 );
  }
  while ( it != itE ) {
    auto c = *it++;
    res.push_back ( hex2val(*it++) + 16* hex2val(c) );
  }
  return res;
}

void NC::streamJSON( std::ostream& os, double val )
{
  if ( std::isnan( val ) )
    NCRYSTAL_THROW(CalcError,"Can not represent not-a-number (NaN) values in JSON format!");

  if ( std::isinf(val) ) {
    //infinity is not supported in json, but at least the python json modules
    //parses an out of bounds value as such:
    os << ( val > 0 ? "1.0e99999" : "-1.0e99999" );
    return;
  }
  if ( val == 0.0 ) {
    os << "0.0";//avoid -0, 0.0, etc.
    return;
  }
  //Encode as string, but make sure that we present floating point numbers in a
  //type-preserving way ("5.0", not "5") so we won't end up mapping C++ floating
  //point types to integer types once the json string is decoded in e.g. Python.
  auto sstr = dbl2shortstr(val);
  auto sstr_v = sstr.to_view();
  if ( !sstr_v.contains_any( ".eE" ) && sstr_v.toInt().has_value() )
    os << fmt( val, "%.1f" );
  else
    os << sstr;
}

void NC::streamJSONFastNotPretty( std::ostream& os, double val )
{
  if ( !std::isfinite( val ) || val == 0.0 ) {
    streamJSON( os, val );
    return;
  }
  char buf[64];
  const int n = std::snprintf( buf, sizeof(buf), "%.17g", val );
  nc_assert_always( n > 0 && n < static_cast<int>( sizeof(buf) ) );
  const auto len = detail::normaliseDecimalPoint( buf,
                                                 static_cast<std::size_t>(n) );
  os.write( buf, static_cast<std::streamsize>( len ) );
  //%.17g has no '.' and no exponent exactly for integers below 1e17:
  if ( ncabs( val ) < 1e17 && val == std::floor( val ) )
    os << ".0";
}

void NC::streamJSON( std::ostream& os, StrView sv )
{
  //According to answer at:
  //https://stackoverflow.com/questions/3020094/how-should-i-escape-strings-in-json:
  //  "JSON is pretty liberal: The only characters you must escape are \, ", and
  //   control codes (anything less than U+0020)."
  //In addition we use the special forms for "\n" "\t" "\r" for more
  //readable/short results.

  os << '"';
  for ( auto ch : sv ) {
    //Code below written to be valid and hopefully non-warning producing for
    //both signed and unsigned char.
    switch ( ch ) {
    //Special cases first:
    case '\x00' : os << '"'; return;//null-char in data, immediately end!
    case '"'    : os << "\\\""; continue;
    case '\\'   : os << "\\\\"; continue;
    case '\n'   : os << "\\n"; continue;
    case '\r'   : os << "\\r"; continue;
    case '\t'   : os << "\\t"; continue;
    default:
      if ( ch < '\x20' && ch >= '\x01' ) {
        //Control character. Must encode hex value using \u00XX format. We use
        //snprintf rather than iomanip to avoid altering the state of the
        //stream, and to hopefully keep malloc's low:
        char buf[8];
        std::snprintf(buf, sizeof(buf)-1, "%04x", (unsigned char)(ch));
        os << "\\u"<<buf;
      } else {
        //Not control char or special case => simply pass through.
        os << ch;
      }
    }
  }
  os << '"';
}

NC::Optional<std::string> NC::findForbiddenChar( const StrView& teststr,
                                                 const StrView& forbidden_chars,
                                                 ExtraForbidOpt extra )
{
  Optional<char> badchar;
  if ( extra == ExtraForbidOpt::RequireSimpleASCII ) {
    for ( auto ch : teststr ) {
      if (!isSimpleASCII(ch,AllowTabs::Yes,AllowNewLine::Yes)) {
        badchar = ch;
        break;
      }
    }
  }
  if ( !badchar.has_value() && forbidden_chars.has_value() ) {
    auto pos = teststr.find_first_of(forbidden_chars);
    if ( pos != std::string::npos )
      badchar = teststr.at(pos);
  }
  if ( !badchar.has_value() )
    return NullOpt;
  return displayCharSafeQuoted(badchar.value());
}

std::string NC::fmtUInt64AsNiceDbl( std::uint64_t n )
{
  std::string res;
  if ( n >= 1000 ) {
    std::uint64_t k = 3;
    std::uint64_t tentok = 1000;
    std::uint64_t best_k = 0;
    std::uint64_t best_tentok = 0;
    const std::uint64_t ndiv10 = n/10;
    while( true ) {
      if ( n % tentok == 0 ) {
        best_k = k;
        best_tentok = tentok;
      }
      if ( tentok > ndiv10 )
        break;
      k += 1;
      tentok *= 10;
    }
    if ( best_k ) {
      res = std::to_string(n/best_tentok);
      res += 'e';
      res += std::to_string(best_k);
    }
  }
  if (res.empty())
    res = std::to_string(n);
  return res;
}

namespace NCRYSTAL_NAMESPACE {
  namespace HugeVectDB {
    struct DB {
      std::mutex mtx;
      bool enabled = false;
      std::map<std::string,VectD> data;
      std::uint64_t count = 0;
      std::string registerData( VectD&& v )
      {
        std::string key = std::to_string(++count);
        data.emplace(key,std::move(v));
        return key;
      }
      VectD retrieveData( const std::string& key )
      {
        auto it = data.find(key);
        if (it == data.end())
          NCRYSTAL_THROW2(BadInput,"Invalid key for JSON Huge Vector: \""
                          <<key<<"\"");
        VectD res = std::move(it->second);
        data.erase(it);
        return res;
      }
    };
    DB& getDB()
    {
      static DB db;
      static const bool dummy = []()
      {
        registerCacheCleanupFunction([]()
        {
          auto& thedb = getDB();
          NCRYSTAL_LOCK_GUARD(thedb.mtx);
          thedb.data.clear();
        });
        return false;
      }();
      (void) dummy;
      return db;
    }

    void enable( bool f )
    {
      auto& db = getDB();
      NCRYSTAL_LOCK_GUARD(db.mtx);
      db.enabled = f;
      //Discard any data left over from previous usage (e.g. because the
      //caller encountered an error before retrieving it):
      if ( f )
        db.data.clear();
    }

    VectD retrieveData( const std::string& key ) {
      auto& db = getDB();
      NCRYSTAL_LOCK_GUARD(db.mtx);
      return db.retrieveData(key);
    }

  }
}

void NC::streamJSONHugeDblVect( std::ostream& os, VectD&& v )
{
  auto& db = HugeVectDB::getDB();
  std::unique_lock<std::mutex> lock(db.mtx);
  if (db.enabled) {
    const std::size_t n = v.size();
    std::string key = db.registerData( std::move(v) );
    os << "\"__ncrystal__dblarray::"<<key<<"::"<<n<<"\"";
    return;
  } else {
    lock.unlock();
  }
  os << '[';
  bool first(true);
  for ( auto val : v ) {
    if ( first )
      first = false;
    else
      os << ',';
    streamJSONFastNotPretty( os, val );
  }
  os<<']';
}
