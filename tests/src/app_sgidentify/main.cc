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

// Tests of NCSGIdentify: space group settings from operations and symbols.

#include "NCrystal/internal/sgsym/NCSGIdentify.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <iostream>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace {

  std::vector<NC::SpaceGroup> allSettings()
  {
    std::vector<NC::SpaceGroup> res;
    for ( auto i : NC::ncrange( NC::SpaceGroupSettingIdx::max_value ) )
      res.emplace_back( NC::SpaceGroupSettingIdx{
          static_cast<std::uint16_t>( i + 1 ) } );
    return res;
  }

  bool contains( const std::vector<NC::SpaceGroup>& v, const NC::SpaceGroup& sg )
  {
    return std::find( v.begin(), v.end(), sg ) != v.end();
  }

  //Symbol with doubled spaces and padding, and upper case lattice letter:
  std::string mangled( std::string s )
  {
    std::string res = "  ";
    for ( auto c : s ) {
      res += ( c == ' ' ? "  " : std::string( 1, c ) );
    }
    const auto i = res.find_first_not_of( " -" );
    res[i] = static_cast<char>( std::toupper( static_cast<unsigned char>( res[i] ) ) );
    return res + " ";
  }

  std::string withoutSpaces( std::string s )
  {
    s.erase( std::remove( s.begin(), s.end(), ' ' ), s.end() );
    return s;
  }

  std::string str( const std::vector<NC::SpaceGroup>& v )
  {
    std::string res = "[";
    for ( auto& sg : v )
      res += ( res.size() > 1 ? ", " : "" ) + sg.toString();
    return res + "]";
  }

  void testOperations()
  {
    unsigned nshared = 0;
    for ( auto& sg : allSettings() ) {
      auto ops = NC::SGSymmetry::get( sg ).operations();
      const auto found = NC::settingsFromOperations( ops );
      REQUIRE( contains( found, sg ) );
      if ( found.size() != 1 ) {
        REQUIRE( found.size() == 2 && sg.number() == 68 );
        ++nshared;
      }
      //Order and duplicates do not matter:
      std::reverse( ops.begin(), ops.end() );
      ops.push_back( ops.front() );
      REQUIRE( NC::settingsFromOperations( ops ) == found );
      //Missing operations (except for P1, the identity is always needed):
      if ( ops.size() > 2 ) {
        ops.pop_back();
        ops.pop_back();
        REQUIRE( NC::settingsFromOperations( ops ).empty() );
      }
    }
    std::cout << "operations: all 530 settings found, " << nshared
              << " share their operations with another setting" << std::endl;
    //An origin shift is not a tabulated setting:
    std::vector<NC::SymOp> shifted = { NC::SymOp( "x,y,z" ),
                                       NC::SymOp( "-x+1/2,-y+1/2,-z+1/2" ) };
    std::cout << "operations of P-1 shifted by (1/4,1/4,1/4): "
              << str( NC::settingsFromOperations( shifted ) ) << std::endl;
    std::cout << "no operations: "
              << str( NC::settingsFromOperations( {} ) ) << std::endl;
  }

  void testTableSymbols()
  {
    for ( auto& sg : allSettings() ) {
      //Hall symbols (with upper case lattice letter and extra spaces) give
      //exactly the settings with the same operations:
      const std::string hall = mangled( sg.hallSymbol() );
      const auto fh = NC::settingsFromHallSymbol( hall );
      REQUIRE( contains( fh, sg ) );
      REQUIRE( fh == NC::settingsFromOperations( NC::SGSymmetry::get( sg ).operations() ) );
      //Hermann-Mauguin symbols, as in the table, and without spaces:
      const std::string hm = sg.hmSymbol();
      REQUIRE( contains( NC::settingsFromHMSymbol( hm ), sg ) );
      const std::string hm_nospace = withoutSpaces( hm ), hm_mangled = mangled( hm );
      REQUIRE( NC::settingsFromHMSymbol( hm ) == NC::settingsFromHMSymbol( hm_nospace ) );
      REQUIRE( NC::settingsFromHMSymbol( hm ) == NC::settingsFromHMSymbol( hm_mangled ) );
    }
    std::cout << "symbols: all 530 settings found by Hall and Hermann-Mauguin symbols" << std::endl;
  }

  void testSymbols()
  {
    for ( const char * s : { "P 1 21/c 1", "P21/c", "P2_1/c", "p 21/C", "P 21/n", "P 1 1 21/a",
                             "P 2", "C 2/c", "R -3 m", "R -3 m:H", "R -3 m :H", "R -3 m H",
                             "R-3m:R", "R -3 m r", "P n n n", "P n n n:2", "Fd-3m", "F d -3 m:2",
                             "P n m a", "Pnma", "C m c e", "C m c a", "A e a 2", "B b e 2",
                             "C c c e", "C c c e:2", "C m m e", "A e m 2",
                             "C -1", "I m m a (-a,c,b)", "F d d 2 S", "", "  ", "x" } )
      std::cout << "H-M \"" << s << "\": " << str( NC::settingsFromHMSymbol( s ) ) << std::endl;
    for ( const char * s : { "-P 2ybc", "-p 2ybc", " -P  2ybc ", "P 2", "P 2y", "-R 3 2\"",
                             "P 2yb (1/2*x,y,-1/2*x+z)", "", "x" } )
      std::cout << "Hall \"" << s << "\": " << str( NC::settingsFromHallSymbol( s ) ) << std::endl;
  }
}

int main()
{
  testOperations();
  testTableSymbols();
  testSymbols();
  return 0;
}
