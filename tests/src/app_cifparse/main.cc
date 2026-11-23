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

#include "NCrystal/internal/cif/NCCIFParse.hh"
#include "NCTestUtils/NCTestExpect.hh"
#include <iostream>
#include <sstream>

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  using Lenient = NCC::CIFParseLenient;
  using Prefixes = NC::Optional<NC::Span<const NC::StrView>>;

  NC::shared_obj<const NC::TextData> mkData( std::string s,
                                             const char* name = "test" )
  {
    return NC::makeSO<const NC::TextData>( NC::RawStrData( std::move( s ) ),
                                           NC::TextData::DataType{ "cif" },
                                           NC::DataSourceName( name ) );
  }

  NCC::Document parse( const std::string& s, Lenient l = Lenient::YES,
                       Prefixes filterValuesByPrefixes = NC::NullOpt )
  {
    return NCC::Document::parse( mkData( s ), l, filterValuesByPrefixes );
  }

  //Printable version of text (long texts abbreviated):
  std::string esc( NC::StrView s, std::size_t maxlen = 80 )
  {
    if ( s.size() > maxlen )
      return esc( s.substr( 0, 10 ) ) + "...(" + std::to_string( s.size() )
        + " chars)";
    std::string res;
    for ( char c : s ) {
      const auto u = static_cast<unsigned char>( c );
      if ( c == '\n' )
        res += "\\n";
      else if ( c == '\r' )
        res += "\\r";
      else if ( c == '\t' )
        res += "\\t";
      else if ( u < 32 || u > 126 ) {
        constexpr char hex[] = "0123456789ABCDEF";
        res += "\\x";
        res += hex[u / 16];
        res += hex[u % 16];
      } else {
        res += c;
      }
    }
    return res;
  }

  const char* kindStr( NCC::ValueKind k )
  {
    switch ( k ) {
    case NCC::ValueKind::Unquoted: return "U";
    case NCC::ValueKind::SingleQuoted: return "S";
    case NCC::ValueKind::DoubleQuoted: return "D";
    case NCC::ValueKind::TextField: return "T";
    case NCC::ValueKind::Unknown: return "?";
    case NCC::ValueKind::Inapplicable: return ".";
    }
    return "!";
  }

  std::string valStr( const NCC::Value& v )
  {
    return std::string( kindStr( v.kind ) ) + "[" + esc( v.text ) + "]@"
      + std::to_string( v.line );
  }

  //Full dump of the document (tags in file order with their values):
  std::string dumpStr( const NCC::Document& doc )
  {
    std::ostringstream os;
    for ( auto ib : NC::ncrange( doc.nBlocks() ) ) {
      auto b = doc.block( ib );
      os << "  data_" << esc( b.name() ) << " @" << b.line() << "\n";
      for ( auto it : NC::ncrange( b.nTags() ) ) {
        auto tag = b.tag( it );
        REQUIRE( b.has( tag ) );
        os << "    " << esc( tag );
        try {
          auto loop = b.loopOf( tag );
          auto vals = b.values( tag );
          if ( loop.has_value() ) {
            REQUIRE( vals.inLoop() );
            os << " (loop col " << loop.value().column( tag ).value() << "/"
               << loop.value().nTags() << ")";
          } else {
            REQUIRE( !vals.inLoop() && vals.size() == 1 );
          }
          os << ":";
          for ( auto i : NC::ncrange( vals.size() ) )
            os << " " << valStr( vals[i] );
        } catch ( NC::Error::LogicError& ) {
          os << ": (not loaded)";
        }
        os << "\n";
      }
    }
    return os.str();
  }

  //Parsing s (with args forwarded to parse) must throw BadInput:
  template<class... TArgs>
  void ensureParseError( const char* title, const std::string& s,
                         TArgs... args )
  {
    std::cout << title << ": ";
    nctest::ensureError<NC::Error::BadInput>( [&s,args...]()
                                              { parse( s, args... ); } );
  }

  //Accepted in both modes, with identical results:
  void showStrict( const char* title, const std::string& s )
  {
    const auto r1 = dumpStr( parse( s, Lenient::YES ) );
    const auto r2 = dumpStr( parse( s, Lenient::NO ) );
    REQUIRE( r1 == r2 );
    std::cout << title << " (both modes):\n" << r1;
  }

  //Accepted in both modes, with different results:
  void showBoth( const char* title, const std::string& s )
  {
    std::cout << title << " (lenient):\n" << dumpStr( parse( s ) );
    std::cout << title << " (strict):\n" << dumpStr( parse( s, Lenient::NO ) );
  }

  //Accepted only in lenient mode:
  void showLenientOnly( const char* title, const std::string& s )
  {
    std::cout << title << " (lenient):\n" << dumpStr( parse( s ) );
    ensureParseError( "  strict", s, Lenient::NO );
  }

  void testSyntax()
  {
    std::cout << "==== Syntax ====\n";
    //CIF 1.1 syntax, paragraphs 4, 7, 63 (single items and loops):
    showStrict( "items and loop", "data_a\n_x 1\n_y 2.5(3)\nloop_\n_l1 _l2\n"
                "a b\nc d\n_z -1e-3\n" );
    //Paragraphs 14, 15 (quoted strings may contain their quote char if not
    //followed by white space):
    showStrict( "quoted", "data_q\n_s1 'a dog's life'\n_s2 \"say \"hi\"!\"\n"
                "_s3 'it''s'\n_s4 ''\n_s5 \"\"\n_s6 'tab after'\t\n"
                "_s7 'x\"y' _s8 \"x'y\"\n_s9 ' '\n_s10 'a'b'\n" );
    //Paragraph 16 (no backslash escapes):
    showStrict( "backslash", "data_b\nloop_ _x _y 'an embedded \\' quote'\n" );
    //Paragraphs 17, 18, 20 (text fields):
    showStrict( "text fields", "data_t\n_t1\n;foo\n;\n_t2\n; foo\n  bar\n;\n"
                "_t5\n; a;b\n x ; y\n ;\n;\n"
                "_t6\n;trailing  \n;\n_t7 \n;'q' \"d\" # not a comment\n;\n" );
    //Paragraph 57 (unquoted strings, incl. ";" not at the start of a line,
    //"?" and "."):
    showStrict( "unquoted", "data_u\n_u1 ;abc\n_u2 a;b\n_u3 a'b\"c\n"
                "_u4 ?\n_u5 .\n_u6 '?'\n_u7 '.'\n_u8 \"?\"\n_u9 ;\n"
                "_u10 ?a\n_u11 .5\n_u12 a#b\n_u13 a_b\n_u14 ..\n" );
    //Paragraphs 11, 57 (reserved words only as exact matches, or as prefix
    //for data_ and save_; quoted they are fine):
    showStrict( "reserved", "data_r\n_r1 loop_x\n_r2 global_x\n_r3 stop_x\n"
                "_r4 'loop_'\n_r5 \"data_x\"\n_r6 xdata_\n_r7 'save_'\n"
                "_r8 dat_\n_r9\n;data_\n;\n" );
    //Paragraph 26 (case-insensitive reserved words and data names, while
    //case of values is kept):
    showStrict( "case", "DATA_Up\n_Tag VaLuE\nLoOp_ _L1 _l2 A b\n"
                "Data_x _y Z\n" );
    //Paragraphs 47, 60 (block names and tags of any non-blank chars):
    showStrict( "names", "data_a'b\"c#d$e;f[g]h_\n_t'a\"g#$;[]_ 1\n" );
    //Paragraphs 21, 45, 46, 58 (comments and white space everywhere):
    showStrict( "comments", "#\\#CIF_1.1\n# leading comment\n\n  \t\n"
                "data_c#x #c\n#c\n_a#b 1#c\n_b 2 #c\n"
                "_c # c\n  # c\n 3\nloop_ # c\n_l # c\n 1 # c\n # c\n 2\n"
                "#c\n 3\n_d 'x'#y' #c\n#final comment without eol" );
    showStrict( "comment after header", "data_x #c\n_a 1" );
    //Paragraph 45 ("#" inside values is data, e.g. labels in real files):
    showStrict( "# in loop values", "data_x\nloop_ _label _type\n"
                "N1T# Cp#2Yb\nC#1 'a#b' #c\n" );
    //Paragraph 24 (tabs and spaces are interchangeable):
    showStrict( "tabs", "data_tab\t_a\t1\tloop_\t_b\t_c\tx\t'y'\t\n" );
    //Paragraph 42 (final line without line ending):
    showStrict( "no final eol", "data_e\n_a 1\n_b last" );
    showStrict( "no final eol quoted", "data_e\n_a 'last'" );
    showStrict( "no final eol textfield", "data_e\n_a\n;x\n;" );
    showStrict( "loop at end", "data_e\nloop_ _a _b 1 2 3 4" );
    showStrict( "loop at end 2", "data_e\nloop_ _a _b 1 2 3 4\n" );
    //Paragraphs 25, 42 (LF, CRLF, CR and mixed line endings):
    showStrict( "CRLF", "data_w\r\n_a 1\r\n_t\r\n;x\r\ny\r\n;\r\n_b 2\r\n" );
    showStrict( "CR", "data_m\r_a 1\r_t\r;x\ry\r;\r_b 2\r" );
    showStrict( "mixed", "data_x\n_a 1\r\n_t\r;x\r\ny\nz\r\r;\n_b 2\r_c 3" );
    showStrict( "CR blank lines", "data_x\r\r\r_a 1\r\n\r\n_b 2\n\n_c 3" );
    //Multiple blocks (paragraph 58):
    showStrict( "blocks", "data_1 _a 1 data_2 _a 2\ndata_3\ndata_4 _b x" );
  }

  void testErrors()
  {
    std::cout << "==== Errors ====\n";
    ensureParseError( "cif2", "#\\#CIF_2.0\ndata_x _a 1\n" );
    ensureParseError( "save frame", "data_x\nsave_foo\n_a 1\nsave_\n" );
    ensureParseError( "save end", "data_x\n_a 1\nSAVE_\n" );
    ensureParseError( "global", "global_\n_a 1\ndata_x\n" );
    ensureParseError( "global in block", "data_x\n_a 1\nGlobal_\n" );
    ensureParseError( "stop", "data_x\nloop_ _a 1 2 stop_\n" );
    ensureParseError( "empty", "" );
    ensureParseError( "only comments", "#\\#CIF_1.1\n# nothing\n\n" );
    ensureParseError( "only whitespace", "  \n\t\r\n" );
    ensureParseError( "unterminated single", "data_x\n_a 'abc\n_b 1\n" );
    ensureParseError( "unterminated single (embedded)", "data_x\n_a 'abc'd\n" );
    ensureParseError( "unterminated single (comment)", "data_x\n_a 'abc'#c\n" );
    ensureParseError( "unterminated double", "data_x\n_a \"abc" );
    ensureParseError( "unterminated text field", "data_x\n_a\n;abc\n def\n" );
    ensureParseError( "unterminated text field (eof)", "data_x\n_a\n;abc" );
    ensureParseError( "text field end", "data_x\n_a\n;abc\n;x\n" );
    ensureParseError( "value instead of tag", "data_x\n_a 1 2\n" );
    ensureParseError( "value instead of tag (quoted)", "data_x\n_a 1\n'2'\n" );
    ensureParseError( "value instead of tag (text field)",
             "data_x\n_a 1\n;some text\nline 2\n;\n" );
    ensureParseError( "tag without value", "data_x\n_a\n_b 1\n" );
    ensureParseError( "tag without value (eof)", "data_x\n_a 1\n_b\n" );
    ensureParseError( "tag without value (eof, comment)", "data_x\n_a #c" );
    ensureParseError( "tag before block", "_a 1\ndata_x\n" );
    ensureParseError( "value before block", "# c\n 1\ndata_x\n" );
    ensureParseError( "loop before block", "loop_ _a 1\ndata_x\n" );
    ensureParseError( "loop without tags", "data_x\nloop_ 1 2\n" );
    ensureParseError( "loop without tags (eof)", "data_x\nloop_\n" );
    ensureParseError( "loop without values", "data_x\nloop_ _a _b\ndata_y\n" );
    ensureParseError( "loop with tag after tags",
                      "data_x\nloop_ _a _b\n_c 1\n" );
    ensureParseError( "loop without values (eof)", "data_x\nloop_ _a _b\n" );
    ensureParseError( "loop without values (loop)",
                      "data_x\nloop_ _a loop_ _b 1\n" );
    ensureParseError( "loop value count",
                      "data_x\nloop_ _a _b\n1 2\n3\n_c 1\n" );
    ensureParseError( "duplicate tag", "data_x\n_a 1\n_b 2\n_A 3\n" );
    ensureParseError( "duplicate tag (loop)",
                      "data_x\n_a 1\nloop_ _b _A 1 2\n" );
    ensureParseError( "duplicate tag (in loop)",
                      "data_x\nloop_ _Ab _aB 1 2\n" );
    ensureParseError( "duplicate tag (first reported)",
             "data_x\n_b 1\n_a 1\n_c 1\n_A 1\n_B 1\n" );
    ensureParseError( "duplicate block",
                      "data_ab\n_a 1\ndata_c\ndata_AB\n_a 1\n" );
    ensureParseError( "reserved word loop_ as value", "data_x\n_a loop_\n" );
    ensureParseError( "reserved word data_ as value", "data_x\n_a data_y\n" );
    ensureParseError( "reserved word save_ as value", "data_x\n_a save_\n" );
    ensureParseError( "reserved word global_ as value",
                      "data_x\n_a GLOBAL_\n" );
    ensureParseError( "reserved word stop_ as value", "data_x\n_a stop_\n" );
    ensureParseError( "reserved char [", "data_x\n_a [1,2]\n" );
    ensureParseError( "reserved char ]", "data_x\nloop_ _a x ]\n" );
    ensureParseError( "reserved char $", "data_x\n_a $frame\n" );
    ensureParseError( "empty block name", "data_\n_a 1\n" );
    ensureParseError( "empty tag name", "data_x\n_ 1\n" );
    ensureParseError( "non-ASCII in tag", "data_x\n_caf\xc3\xa9 1\n" );
    ensureParseError( "non-ASCII in block name", "data_caf\xc3\xa9\n_a 1\n" );
    ensureParseError( "control char in tag", "data_x\n_a\x0b"  "b 1\n" );
    ensureParseError( "non-ASCII before block", "\xc3\xa9" "data_x\n_a 1\n" );
  }

  void testLenient()
  {
    std::cout << "==== Lenient ====\n";
    //CIF 1.1 syntax, paragraph 28:
    const std::string line2048 = "data_x\n_a " + std::string( 2048 - 3, 'v' )
      + "\n";
    const std::string line2049 = "data_x\n_a " + std::string( 2048 - 2, 'v' )
      + "\n";
    REQUIRE( parse( line2048, Lenient::NO ).block( 0 ).values( "_a" )[0]
             .text.size() == 2045 );
    std::cout << "line of 2048 chars accepted in both modes\n";
    showLenientOnly( "line of 2049 chars", line2049 );
    showLenientOnly( "line of 2049 chars (CRLF)", "data_x\r\n_a "
                     + std::string( 2046, 'v' ) + "\r\n" );
    showLenientOnly( "last line of 2049 chars",
                     "data_x\n_a " + std::string( 2046, 'v' ) );
    showLenientOnly( "comment line of 2049 chars",
                     "data_x\n_a 1\n#" + std::string( 2048, 'c' ) + "\n" );
    showLenientOnly( "text field line of 2049 chars",
                     "data_x\n_a\n;" + std::string( 2048, 't' ) + "\n;\n" );
    //Paragraphs 29, 30:
    const std::string tag75 = "_" + std::string( 74, 't' );
    const std::string name75( 75, 'n' );
    showStrict( "tag and block name of 75 chars",
                "data_" + name75 + "\n" + tag75 + " 1\n" );
    showLenientOnly( "tag of 76 chars", "data_x\n" + tag75 + "t 1\n" );
    showLenientOnly( "loop tag of 76 chars",
                     "data_x\nloop_ _a " + tag75 + "t 1 2\n" );
    showLenientOnly( "block name of 76 chars", "data_" + name75 + "n\n" );
    //Paragraph 22 (lenient: non-ASCII bytes in values and comments):
    showLenientOnly( "non-ASCII unquoted", "data_x\n_a caf\xc3\xa9\n" );
    showLenientOnly( "non-ASCII first char",
                     "data_x\n_a \xc3\xa9t\xc3\xa9\n" );
    showLenientOnly( "non-ASCII single-quoted",
                     "data_x\n_a 'caf\xc3\xa9 x'\n" );
    showLenientOnly( "non-ASCII double-quoted",
                     "data_x\n_a \"\xe2\x84\xab\"\n" );
    showLenientOnly( "non-ASCII text field",
                     "data_x\n_a\n;\xc3\x85ngstr\xc3\xb6m\n;\n" );
    showLenientOnly( "non-ASCII comment", "data_x\n_a 1 # \xc3\xa9\n" );
    showLenientOnly( "non-ASCII comment before block",
                     "# \xc3\xa9\ndata_x\n_a 1\n" );
    //Lenient: no line ending after an opening ";" alone on its line
    //(paragraph 17 strictly keeps it):
    showBoth( "text field, empty first line",
              "data_x\n_t1\n;\nfoo\n;\n_t2\n;\n;\n_t3\n;\n\n;\n"
              "_t4\n;\r\nfoo\r\nbar\r\n;\n_t5\n;\rfoo\r;\n" );
    showStrict( "text field, content after opening ;",
                "data_x\n_t1\n;foo\nbar\n;\n_t2\n; \nfoo\n;\n" );
    //Lenient: a leading UTF-8 byte-order mark is skipped:
    showLenientOnly( "BOM", "\xef\xbb\xbf" "data_x\n_a 1\n" );
    showLenientOnly( "BOM and comment", "\xef\xbb\xbf" "#\\#CIF_1.1\n"
                     "data_x\n_a 1\n" );
    showLenientOnly( "BOM and text field",
                     "\xef\xbb\xbf" "data_x\n_a\n;\nx\n;" );
    ensureParseError( "BOM in block", "data_x\n\xef\xbb\xbf" "_a 1\n" );
    ensureParseError( "non-ASCII in reserved word as value",
                      "data_x\n_a save_caf\xc3\xa9\n" );
    showStrict( "no BOM", "data_x\n_a 1\n" );
    showLenientOnly( "VT in value", "data_x\n_a a\x0b" "b\n" );
    showLenientOnly( "FF in comment", "data_x\n_a 1 #\x0c\n" );
    showLenientOnly( "DEL in value", "data_x\n_a '\x7f'\n" );
  }

  void testLookup()
  {
    std::cout << "==== Lookup ====\n";
    auto doc = parse( "data_One\n_cell_a 1.0\n_Cell_B 2.0\nloop_\n"
                      "_atom_site_label\n_atom_site_x\nAl1 0.1\nAl2 0.2\n"
                      "Al3 0.3\ndata_two\n_cell_a 5.0\n_x ?\n" );
    REQUIRE( doc.nBlocks() == 2 );
    auto b1 = doc.findBlock( "one" );
    auto b2 = doc.findBlock( "TWO" );
    REQUIRE( b1.has_value() && b2.has_value() );
    REQUIRE( !doc.findBlock( "three" ).has_value() );
    REQUIRE( !doc.findBlock( "" ).has_value() );
    REQUIRE( !doc.findBlock( "on" ).has_value() );
    REQUIRE( doc.findBlock( "One" ).value().name() == "One" );
    std::cout << "block names: " << b1.value().name() << " "
              << b2.value().name() << "\n";
    for ( auto& b : { b1.value(), b2.value() } ) {
      auto a = b.values( "_CELL_A" );
      REQUIRE( a.size() == 1 && !a.inLoop() && !a.empty() );
      std::cout << "data_" << b.name() << " _cell_a: " << valStr( a[0] )
                << "\n";
    }
    auto bl = b1.value();
    REQUIRE( bl.has( "_cell_b" ) && !bl.has( "_cell_c" ) );
    REQUIRE( bl.values( "_cell_c" ).empty() );
    REQUIRE( !bl.loopOf( "_cell_c" ).has_value() );
    REQUIRE( !bl.loopOf( "_cell_a" ).has_value() );
    REQUIRE( !bl.has( "_cell" ) && !bl.has( "_cell_a_" ) );
    REQUIRE( !b2.value().has( "_atom_site_x" ) );
    auto x = bl.values( "_atom_site_X" );
    REQUIRE( x.size() == 3 && x.inLoop() );
    std::cout << "_atom_site_x:";
    for ( auto i : NC::ncrange( x.size() ) )
      std::cout << " " << valStr( x[i] );
    std::cout << "\n";
    auto loop = bl.loopOf( "_ATOM_SITE_label" );
    REQUIRE( loop.has_value() );
    auto& lv = loop.value();
    REQUIRE( lv.nTags() == 2 && lv.nRows() == 3 );
    REQUIRE( lv.column( "_atom_site_x" ).value() == 1 );
    REQUIRE( lv.column( "_Atom_Site_Label" ).value() == 0 );
    REQUIRE( !lv.column( "_cell_a" ).has_value() );
    REQUIRE( lv.tag( 1 ) == "_atom_site_x" );
    for ( auto r : NC::ncrange( lv.nRows() ) )
      std::cout << "row " << r << ": " << valStr( lv.value( r, 0 ) ) << " "
                << valStr( lv.value( r, 1 ) ) << "\n";
    auto col = lv.values( 0 );
    REQUIRE( col.size() == 3 && col[2].text == "Al3" );
    REQUIRE( doc.dataSourceName().str() == "test" );
  }

  void testFilter()
  {
    std::cout << "==== filterValuesByPrefixes ====\n";
    const std::string data = "data_k\n_cell_a 1\n_symmetry_x 'abc'\n"
      "_other\n;text\r\nfield\r\n;\nloop_ _atom_site_label _atom_site_x\n"
      "A 1 B 2\nloop_ _geom_a _geom_b\n1 2\n3 4\n_CELL_B 3\n"
      "loop_ _foo _Cell_c 1 2\ndata_k2\n_other 1\n_cell_a 2\n";
    constexpr NC::StrView prefixes[] = { NC::StrView::make( "_cell_" ),
                                         NC::StrView::make( "_ATOM_" ) };
    const NC::Span<const NC::StrView> prefixSpan( prefixes );
    const NC::Span<const NC::StrView> emptySpan;
    auto dAll = parse( data );
    auto dFilt = parse( data, Lenient::YES, prefixSpan );
    auto dNone = parse( data, Lenient::YES, emptySpan );
    std::cout << "no filter:\n" << dumpStr( dAll );
    std::cout << "filter:\n" << dumpStr( dFilt );
    std::cout << "empty filter:\n" << dumpStr( dNone );
    for ( auto ib : NC::ncrange( dAll.nBlocks() ) ) {
      auto ba = dAll.block( ib );
      auto bk = dFilt.block( ib );
      auto bn = dNone.block( ib );
      REQUIRE( ba.nTags() == bk.nTags() && ba.nTags() == bn.nTags() );
      for ( auto it : NC::ncrange( ba.nTags() ) ) {
        auto tag = ba.tag( it );
        REQUIRE( bk.tag( it ) == tag && bn.tag( it ) == tag );
        REQUIRE( bk.has( tag ) && bn.has( tag ) );
        bool threw = false;
        try {
          bn.values( tag );
        } catch ( NC::Error::LogicError& ) {
          threw = true;
        }
        REQUIRE( threw );
        threw = false;
        try {
          bn.loopOf( tag );
        } catch ( NC::Error::LogicError& ) {
          threw = true;
        }
        REQUIRE( threw );
        threw = false;
        try {
          auto vk = bk.values( tag );
          auto va = ba.values( tag );
          REQUIRE( vk.size() == va.size() );
          for ( auto i : NC::ncrange( va.size() ) )
            REQUIRE( vk[i].text == va[i].text && vk[i].kind == va[i].kind
                     && vk[i].line == va[i].line );
        } catch ( NC::Error::LogicError& ) {
          threw = true;
        }
        std::cout << "  " << tag << ( threw ? " filtered" : " kept" ) << "\n";
      }
    }
    try {
      dFilt.block( 0 ).values( "_geom_a" );
    } catch ( NC::Error::LogicError& e ) {
      std::cout << "LogicError: " << e.what() << "\n";
    }
    //Syntax errors are also found where values are not stored:
    ensureParseError( "error in filtered loop",
                      "data_x\n_a 1\nloop_ _b _c 1 2 3\n", Lenient::YES,
                      Prefixes( emptySpan ) );
  }

  void testMoves()
  {
    std::cout << "==== Moves ====\n";
    const std::string data = "data_m\n_a 'x'\n_t\n;l1\r\nl2\r\n;\n"
      "loop_ _l _v a 1 b 2\n";
    auto d1 = parse( data );
    const auto uid = d1.getUniqueID();
    const std::string ref = dumpStr( d1 );
    const char* tptr = d1.block( 0 ).values( "_t" )[0].text.data_begin();
    const char* aptr = d1.block( 0 ).values( "_a" )[0].text.data_begin();
    NCC::Document d2( std::move( d1 ) );
    REQUIRE( d2.getUniqueID() == uid );
    REQUIRE( dumpStr( d2 ) == ref );
    REQUIRE( d2.block( 0 ).values( "_t" )[0].text.data_begin() == tptr );
    REQUIRE( d2.block( 0 ).values( "_a" )[0].text.data_begin() == aptr );
    auto d3 = parse( "data_other _x 1" );
    REQUIRE( d3.getUniqueID() != uid );
    auto col = d2.block( 0 ).values( "_v" );
    d3 = std::move( d2 );
    REQUIRE( d3.getUniqueID() == uid );
    REQUIRE( dumpStr( d3 ) == ref );
    REQUIRE( col.size() == 2 && col[1].text == "2" );
    REQUIRE( d3.block( 0 ).values( "_t" )[0].text == "l1\nl2" );
    REQUIRE( d3.block( 0 ).values( "_a" )[0].text.data_begin() == aptr );
    std::cout << "after moves:\n" << ref;
  }

  void testMPNote()
  {
    //NB: Fallback for the materialsproject.org note, cf. Document::parse:
    std::cout << "==== MP note fallback ====\n";
    const std::string note = "# Note from NCrystal.cifutils: Data from"
      " materialsproject.org / mp-";
    showStrict( "note (added to the last block)",
                "data_a _x 1\r\ndata_b _y 2\r\n\r\n" + note + "149\r\n" );
    showStrict( "note without final line ending",
                "data_a _x 1\n" + note + "55" );
    showStrict( "tag already present", "data_a\n"
                "_[LOCAL]_ncrystal_mp_database_code mp-7\n" + note + "8\n" );
    showStrict( "note too early", "data_a _x 1\n" + note + "9\n#"
                + std::string( 500, 'x' ) + "\n" );
    showStrict( "not an id", "data_a _x 1\n" + note + "\n" + note + "1x\n" );
    constexpr NC::StrView prefixes[] = { NC::StrView::make( "_x" ) };
    std::cout << "filtered out:\n"
              << dumpStr( parse( "data_a _x 1\n" + note + "3\n", Lenient::NO,
                                 NC::Span<const NC::StrView>( prefixes ) ) );
  }

  void testLarge()
  {
    std::cout << "==== Large ====\n";
    std::string s;
    s.reserve( 20000000 );
    s += "# generated\n";
    const unsigned nblocks = 2000;
    for ( auto i : NC::ncrange( nblocks ) ) {
      s += "data_block" + std::to_string( i ) + "\n_cell_length_a "
        + std::to_string( i ) + ".5\n_other_item 'v" + std::to_string( i )
        + "'\n";
    }
    s += "data_big\n_cell_length_a 1\nloop_\n_atom_site_label\n"
      "_atom_site_fract_x\n_atom_site_fract_y\n";
    const unsigned nrows = 400000;
    for ( auto i : NC::ncrange( nrows ) ) {
      const auto si = std::to_string( i );
      s += "X" + si + " 0." + si + " " + si + "e-6\n";
    }
    s += "loop_ _junk_a _junk_b\n";
    for ( auto i : NC::ncrange( nrows ) )
      s += std::to_string( i ) + " ?\n";
    constexpr NC::StrView prefixes[] = { NC::StrView::make( "_atom_site_" ) };
    const NC::Span<const NC::StrView> prefixSpan( prefixes );
    for ( int ifilt : { 0, 1 } ) {
      auto doc = ( ifilt ? parse( s, Lenient::NO, prefixSpan )
                   : parse( s, Lenient::NO ) );
      std::cout << ( ifilt ? "with filter" : "no filter" ) << ": "
                << doc.nBlocks() << " blocks\n";
      REQUIRE( doc.nBlocks() == nblocks + 1 );
      auto b = doc.findBlock( "BLOCK1234" ).value();
      REQUIRE( b.has( "_cell_length_a" ) && b.has( "_other_item" ) );
      if ( !ifilt ) {
        std::cout << "  block1234: "
                  << valStr( b.values( "_cell_length_a" )[0] ) << " "
                  << valStr( b.values( "_other_item" )[0] ) << "\n";
      }
      auto big = doc.findBlock( "big" ).value();
      auto loop = big.loopOf( "_atom_site_fract_y" ).value();
      REQUIRE( loop.nRows() == nrows && loop.nTags() == 3 );
      for ( unsigned r : { 0u, 1u, 123456u, nrows - 1 } )
        std::cout << "  row " << r << ": " << valStr( loop.value( r, 0 ) )
                  << " " << valStr( loop.value( r, 1 ) ) << " "
                  << valStr( loop.value( r, 2 ) ) << "\n";
      auto y = big.values( "_atom_site_fract_y" );
      REQUIRE( y.size() == nrows && y[nrows - 1].text == "399999e-6" );
      REQUIRE( big.has( "_junk_b" ) );
      if ( ifilt ) {
        bool threw = false;
        try {
          big.values( "_junk_b" );
        } catch ( NC::Error::LogicError& ) {
          threw = true;
        }
        REQUIRE( threw );
      } else {
        auto jb = big.values( "_junk_b" );
        REQUIRE( jb.size() == nrows && jb[7].kind == NCC::ValueKind::Unknown );
        std::cout << "  _junk_a last: "
                  << valStr( big.values( "_junk_a" )[nrows - 1] ) << "\n";
      }
    }
  }
}

int main()
{
  testSyntax();
  testErrors();
  testLenient();
  testLookup();
  testFilter();
  testMoves();
  testMPNote();
  testLarge();
  return 0;
}
