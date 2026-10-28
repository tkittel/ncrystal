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
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include <iostream>
#include <sstream>
#include <limits>
#include <cstring>
#include <cstdio>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(__LINE__ && x)

void dump(const char* settings,
          const std::string& s,
          const std::vector<std::string>& result)
{
  std::cout << '"' << s << "\" splits to "<<result.size()<<" parts (settings: "<<settings<<"):"<<std::endl;
  for (std::size_t i = 0; i<result.size(); ++i)
    std::cout << "       -->  \"" << result.at(i) << '"' << std::endl;
}

void test_split()
{
  std::vector<std::string> parts;
  std::string s;

  s="";
  NC::split(parts,s);
  dump("no options",s,parts);

  s="";
  NC::split(parts,s,0,';');
  dump("sep=;",s,parts);

  s="   ";
  NC::split(parts,s,0,';');
  dump("sep=;",s,parts);

  s="   ";
  NC::split(parts,s);
  dump("no options",s,parts);

  s="  hej med dig \t lala\t";
  NC::split(parts,s);
  dump("no options",s,parts);

  s=";hej;med;;dig;";
  NC::split(parts,s,0,';');
  dump("sep=;",s,parts);

  s=";hej;med; dig;;;hejsa";
  NC::split(parts,s,0,';');
  dump("sep=;",s,parts);

  s="hej;med;dig;";
  NC::split(parts,s,1,';');
  dump("maxsplit 1;sep=;",s,parts);

  s="hej med dig";
  NC::split(parts,s,1);
  dump("maxsplit 1",s,parts);

  s=" hej med dig";
  NC::split(parts,s,1);
  dump("maxsplit 1",s,parts);

  s=" hej med dig ";
  NC::split(parts,s,1);
  dump("maxsplit 1",s,parts);

  s="temp=77.0";
  NC::split(parts,s,0,'=');
  dump("sep='='",s,parts);
}

void require_bad_dbl(const char * c) {
  bool bad = false;
  try {
    NC::str2dbl(c);
  } catch ( NC::Error::BadInput& e ) {
    bad = true;
  }
  REQUIRE(bad);
}
void require_bad_int(const char * c) {
  bool bad = false;
  try {
    NC::str2int(c);
  } catch ( NC::Error::BadInput& e ) {
    bad = true;
  }
  REQUIRE(bad);
}
void require_bad_int32(const char * c) {
  bool bad = false;
  try {
    NC::str2int32(c);
  } catch ( NC::Error::BadInput& e ) {
    bad = true;
  }
  REQUIRE(bad);
}
void require_bad_int64(const char * c) {
  bool bad = false;
  try {
    NC::str2int64(c);
  } catch ( NC::Error::BadInput& e ) {
    bad = true;
  }
  REQUIRE(bad);
}
void require_bad_uint64(const char * c) {
  bool bad = false;
  try {
    NC::str2uint64(c);
  } catch ( NC::Error::BadInput& e ) {
    bad = true;
  }
  REQUIRE(bad);
}


void test_misc() {
  REQUIRE(!NC::isSimpleASCII("something\tsdf",NC::AllowTabs::No,NC::AllowNewLine::No));
  REQUIRE(!NC::isSimpleASCII("something\nsdf",NC::AllowTabs::No,NC::AllowNewLine::No));
  REQUIRE(!NC::isSimpleASCII("something\rsdf",NC::AllowTabs::No,NC::AllowNewLine::No));
  REQUIRE( NC::isSimpleASCII("something\t sdf",NC::AllowTabs::Yes,NC::AllowNewLine::No));
  REQUIRE( NC::isSimpleASCII("something\n sdf",NC::AllowTabs::No,NC::AllowNewLine::Yes));
  REQUIRE( NC::isSimpleASCII("something\r sdf",NC::AllowTabs::No,NC::AllowNewLine::Yes));
  REQUIRE( NC::isSimpleASCII("some\tthing\n sdf",NC::AllowTabs::Yes,NC::AllowNewLine::Yes));

  std::string t;
  t = "bla";NC::trim(t);REQUIRE(t=="bla");
  t = "\rbla";NC::trim(t);REQUIRE(t=="bla");
  t = "bla\t\n";NC::trim(t);REQUIRE(t=="bla");
  t = "\t\nb la\r";NC::trim(t);REQUIRE(t=="b la");
  t = "\t\nb l\t\na\r";NC::trim(t);REQUIRE(t=="b l\t\na");
  t = "b  la";NC::trim(t);REQUIRE(t=="b  la");
  t = " b l a ";NC::trim(t);REQUIRE(t=="b l a");

  REQUIRE(NC::startswith("hejsa","hej"));
  REQUIRE(!NC::startswith(" hejsa","hej"));
  REQUIRE(NC::startswith("hejsa",""));
  REQUIRE(NC::startswith("",""));
  REQUIRE(!NC::startswith("\n","\r"));
  REQUIRE(!NC::startswith("hej","hejsa"));

  REQUIRE(NC::contains("hejsa",'h'));
  REQUIRE(NC::contains("hejsa","h"));
  REQUIRE(NC::contains("hejsa","hej"));
  REQUIRE(NC::contains("hejsa","ejsa"));
  REQUIRE(!NC::contains("hejsa\n sdfsdf ds f$%^&^  \n lala","t"));
  REQUIRE(!NC::contains("hejsa\n sdfsdf ds f$%^&^  \n lala",'t'));
  REQUIRE(NC::contains("hejsa\n sdfsdf ds f$%^&^  \n lala","$"));
  REQUIRE(NC::contains("hejsa\n sdfsdf ds f$%^&^  \n lala",'$'));
  REQUIRE(!NC::contains("hejsa",'\0'));
  REQUIRE(!NC::contains("",'\0'));
  REQUIRE(NC::contains("",""));
  REQUIRE(NC::contains("sdfsdf",""));

  //Check if any of the chars in "needles" is present in the string (haystack):
  REQUIRE(NC::contains_any("sdfsdf","456457d"));
  REQUIRE(!NC::contains_any("","456457d"));
  REQUIRE(!NC::contains_any("sdfsd",""));
  REQUIRE(!NC::contains_any("hejsa","HEJSA"));
  REQUIRE(NC::contains_any("hej med dig!","'{}![]"));

  std::cout<<"Testing str2dbl"<<std::endl;
  REQUIRE(NC::str2dbl("2.0")==2.0);
  REQUIRE(NC::str2dbl("2.0e-3")==2.0e-3);
  REQUIRE(NC::str2dbl("-2.0E-3")==-2.0e-3);
  REQUIRE(NC::str2dbl("-1E-3")==-1.0e-3);
  require_bad_dbl("2.0a");
  require_bad_dbl(" 2.0");
  require_bad_dbl("2.0 ");
  require_bad_dbl("2e.0");
  require_bad_dbl("e-3");

  //Not sure if these should actually work or not, but they seem not to with gcc
  //6.3.1 so we put a test here for now. If they fail on some platforms, we
  //might have to add code in NC::str2dbl looking for strings like infinity
  //strings (to get consistent behaviour everywhere):

  REQUIRE(NC::ncisinf(NC::str2dbl("inf")));
  require_bad_dbl("\tinf");
  require_bad_dbl("\tINF");
  require_bad_dbl("  \t  inf");
  require_bad_dbl(" \t INF");
  require_bad_dbl("inf ");
  require_bad_dbl("INF ");
  require_bad_dbl("infinity");
  require_bad_dbl("INFINITY");

  std::cout<<"Testing str2int"<<std::endl;
  REQUIRE(NC::str2int("-1")==-1);
  REQUIRE(NC::str2int("34545")==34545);
  REQUIRE(NC::str2int("000")==0);
  require_bad_int("2.0");
  require_bad_int("1e2");
  require_bad_int("0.1");
  require_bad_int("2 ");
  require_bad_int(" 2");

  std::cout<<"Testing str2int32"<<std::endl;

  REQUIRE(NC::str2int32("-1")==-1);
  REQUIRE(NC::str2int32("34545")==34545);
  REQUIRE(NC::str2int32("000")==0);
  require_bad_int32("2.0");
  require_bad_int32("1e2");
  require_bad_int32("0.1");
  require_bad_int32("2 ");
  require_bad_int32(" 2");
  std::cout<<"Testing str2int32 out of bound values"<<std::endl;
  require_bad_int32("3000000000");
  require_bad_int32("-3000000000");

  std::cout<<"Testing str2int64"<<std::endl;
  REQUIRE(NC::str2int64("-1")==-1);
  REQUIRE(NC::str2int64("34545")==34545);
  REQUIRE(NC::str2int64("000")==0);
  REQUIRE(NC::str2int64("3000000000")==3000000000LL);
  REQUIRE(NC::str2int64("-3000000000")==-3000000000LL);
  require_bad_int64("2.0");
  require_bad_int64("1e2");
  require_bad_int64("0.1");
  require_bad_int64("2 ");
  require_bad_int64(" 2");

  std::cout<<"Testing str2int64 out of bound values"<<std::endl;
  REQUIRE(NC::str2int64("3000000000000")==3000000000000ULL);
  REQUIRE(NC::str2int64("-3000000000000")==-3000000000000LL);
  REQUIRE(NC::str2int64("9223372036854775807")==9223372036854775807LL);
  require_bad_int64("9223372036854775808");
  REQUIRE(NC::str2int64("-9223372036854775807")==-9223372036854775807LL);
  constexpr std::int64_t minvalint64 = -9223372036854775807LL;
  REQUIRE(NC::str2int64("-9223372036854775808")==(minvalint64-1));
  require_bad_int64("-9223372036854775809");

  std::cout<<"Testing str2uint64"<<std::endl;
  REQUIRE(NC::str2uint64("34545")==34545);
  REQUIRE(NC::str2uint64("000")==0);
  REQUIRE(NC::str2uint64("3000000000")==3000000000ULL);
  REQUIRE(NC::str2uint64("3000000000000")==3000000000000ULL);
  require_bad_uint64("2.0");
  require_bad_uint64("1e2");
  require_bad_uint64("0.1");
  require_bad_uint64("2 ");
  require_bad_uint64(" 2");
  require_bad_uint64("-1");
  require_bad_uint64("-3000000000000");
  require_bad_uint64("-3000000000");
  std::cout<<"Testing str2uint64 out of bound values"<<std::endl;
  REQUIRE(NC::str2uint64("18446744073709551615")==18446744073709551615ULL);
  require_bad_uint64("18446744073709551616");
  REQUIRE(NC::str2uint64("9223372036854775808")==9223372036854775808ULL);
}

namespace {
  //Str2DblConverter must agree with safe_str2dbl and str2dbl on everything
  //(values, failures and exception messages), also when good and bad inputs
  //are interleaved:
  struct ConverterCheck {
    NC::Str2DblConverter conv;
    std::size_t nchecks = 0;
    void check( const std::string& s )
    {
      double v1(0.0);
      const bool ok1 = NC::safe_str2dbl( s, v1 );
      const auto r2 = conv.tryConvert( s );
      REQUIRE( ok1 == r2.has_value() );
      const double v2 = r2.value_or( 0.0 );
      if ( ok1 )
        REQUIRE( NC::ncisnan(v1) ? NC::ncisnan(v2)
                 : std::memcmp( &v1, &v2, sizeof(double) ) == 0 );
      std::string e1, e2;
      try {
        NC::str2dbl( s );
      } catch ( NC::Error::BadInput& e ) {
        e1 = e.what();
      }
      try {
        conv.convert( s );
      } catch ( NC::Error::BadInput& e ) {
        e2 = e.what();
      }
      REQUIRE( e1 == e2 );
      ++nchecks;
    }
  };
}

void test_str2dblconverter()
{
  ConverterCheck cc;
  const char * fixed[] = {
    "2.0", "2.0e-3", "-2.0E-3", "-1E-3", "2.0a", " 2.0", "2.0 ", "2e.0", "e-3",
    "inf", "-inf", "+inf", "INF", "+INF", "-INF", "nan", "NaN", "NAN", "-nan",
    "infinity", "INFINITY", "\tinf", "inf ", "", "-", "+", ".", "1e", "1e+",
    "0x1p3", "0x10", "1,5", "1.5.2", "--1", "+-1", "0", "-0", "-0.0", "+0.",
    ".5", "5.", "00012", "1_000", "\t1", "1\n", "1e308", "1e309", "-1e309",
    "1.7976931348623157e308", "1.7976931348623159e308",
    "2.2250738585072014e-308", "4.94066e-324", "5e-324", "1e-320", "1e-301",
    "1.123456789123456789e-315", "-5e-320", "1e-400", "1e-3001",
    "3.14159265358979323846264338327950288419716939937510582097494459",
  };
  for ( auto s : fixed )
    cc.check( s );

  //Random bit patterns (all exponents, denormals, inf, nan) in several formats,
  //and random strings of number characters (mostly invalid):
  NC::RandXRSRImpl rng( 12345 );
  const char * fmtstrs[] = { "%.17g", "%.6e", "%g", "%.3f", "%a" };
  for ( int i = 0; i < 20000; ++i ) {
    const std::uint64_t bits = rng.genUInt64();
    double x;
    std::memcpy( &x, &bits, sizeof(double) );
    std::ostringstream ss;
    ss << NC::fmt( x );
    cc.check( ss.str() );
    for ( auto fs : fmtstrs ) {
      char buf[512];//"%.3f" of 1e308 needs ~315 chars
      std::snprintf( buf, sizeof(buf), fs, x );
      cc.check( buf );
    }
    const char chars[] = "0123456789.eE+-";
    std::string s;
    const auto n = 1 + rng.genUInt64() % 8;
    for ( std::uint64_t j = 0; j < n; ++j )
      s += chars[ rng.genUInt64() % ( sizeof(chars) - 1 ) ];
    cc.check( s );
  }

  //Still agrees after being moved:
  NC::Str2DblConverter conv2( std::move( cc.conv ) );
  cc.conv = std::move( conv2 );
  for ( auto s : fixed )
    cc.check( s );
  std::cout << "Str2DblConverter agrees with str2dbl in "
            << cc.nchecks << " checks" << std::endl;
}

int main(int,char**) {
  std::cout<<"test_split begin"<<std::endl;
  test_split();
  std::cout<<"test_split end"<<std::endl;
  std::cout<<"test_misc begin"<<std::endl;
  test_misc();
  std::cout<<"test_misc emd"<<std::endl;
  test_str2dblconverter();
  std::cout<<"all done"<<std::endl;
  return 0;
}
