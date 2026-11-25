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

#include "NCrystal/internal/cif/NCCIFSetting.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCTestUtils/NCTestExpect.hh"
#include <iostream>

//Tests of CIF::chooseSetting: the space group setting of CIF data, from the
//listed symmetry operations (cross-checked against the symbol), or else from
//the symbol (several candidates resolved by the cell), or from an explicit
//override, in default and permissive mode.

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace {

  using Permissive = NCC::CIFPermissive;
  using VectS = std::vector<std::string>;

  const char * cubic = "5.43 5.43 5.43 90 90 90";
  const char * hexagonal = "3.0 3.0 10.0 90 90 120";
  const char * rhombohedral = "4.0 4.0 4.0 70 70 70";
  const char * orthorhombic = "5 6 7 90 90 90";

  VectS opsOf( const char * setting )
  {
    VectS res;
    NC::SpaceGroup sg( setting, NC::SpaceGroupParseMode::Strict );
    for ( auto& o : NC::SGSymmetry::get( sg ).operations() )
      res.push_back( o.toString() );
    return res;
  }

  struct Content {
    Content( std::string c, VectS o = {}, std::string ha = {},
             std::string h = {}, std::string n = {}, std::string e = {} )
      : cell( std::move( c ) ), ops( std::move( o ) ), hall( std::move( ha ) ),
        hm( std::move( h ) ), itnumber( std::move( n ) ),
        extra( std::move( e ) ) {}
    std::string cell;
    VectS ops;
    std::string hall, hm, itnumber, extra;
  };

  //CIF data with the content (empty items are absent):
  std::string cifData( const Content& c )
  {
    std::string res = "data_test\n" + c.extra;
    const char * names[] = { "length_a", "length_b", "length_c",
                             "angle_alpha", "angle_beta", "angle_gamma" };
    auto cellvals = NC::StrView( c.cell ).split();
    for ( auto i : NC::ncrange( 6 ) )
      res += std::string( "_cell_" ) + names[i] + " "
        + NC::vectAt( cellvals, i ).to_string() + "\n";
    if ( !c.hall.empty() )
      res += "_space_group_name_Hall '" + c.hall + "'\n";
    if ( !c.hm.empty() )
      res += "_space_group_name_H-M_alt '" + c.hm + "'\n";
    if ( !c.itnumber.empty() )
      res += "_space_group_IT_number " + c.itnumber + "\n";
    if ( !c.ops.empty() ) {
      res += "loop_\n_space_group_symop_operation_xyz\n";
      for ( auto& o : c.ops )
        res += "'" + o + "'\n";
    }
    return res + ( "loop_\n_atom_site_label\n_atom_site_fract_x\n"
                   "_atom_site_fract_y\n_atom_site_fract_z\nAl1 0 0 0\n" );
  }

  void test( const char * descr, const Content& c,
             Permissive perm = Permissive::NO, const char * ovr = nullptr )
  {
    VectS warnings;
    const NCC::CIFWarnFct warnfct = [&warnings]( const std::string& w )
    {
      warnings.push_back( w );
    };
    std::cout << descr << ": ";
    try {
      auto doc = NCC::Document::parse(
        NC::makeSO<const NC::TextData>( NC::RawStrData( cifData( c ) ),
                                        NC::TextData::DataType{ "cif" },
                                        NC::DataSourceName( "test.cif" ) ) );
      auto items = NCC::findUniqueStructProvItems( doc, perm, warnfct );
      const auto sg = NCC::chooseSetting( doc, items.structureItems, perm,
                                          ovr ? NC::StrView( ovr )
                                          : NC::StrView(), warnfct );
      std::cout << sg << std::endl;
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "rejected: " << e.what() << std::endl;
    }
    for ( auto& w : warnings )
      std::cout << "   warning: " << w << std::endl;
  }

  VectS concat( VectS a, const VectS& b )
  {
    a.insert( a.end(), b.begin(), b.end() );
    return a;
  }
}

int main()
{
  const auto ops227 = opsOf( "227:2" );
  const VectS ops227m1( ops227.begin(), ops227.end() - 1 );
  VectS ops227rev( ops227.rbegin(), ops227.rend() );
  const auto P = Permissive::YES;

  //A: listed operations (any order, with duplicates), cross-checked against
  //the symbol:
  test( "listed operations", { cubic, ops227 } );
  test( "listed operations, reversed with duplicates",
        { cubic, concat( ops227rev, VectS( ops227.begin(),
                                           ops227.begin() + 3 ) ) } );
  test( "listed operations and matching symbol",
        { cubic, ops227, "", "F d -3 m" } );
  test( "listed operations and unrecognised symbol (ignored)",
        { cubic, ops227, "", "F d -3 m (bla)" } );
  test( "shared operations of number 68", { orthorhombic, opsOf( "68:1" ) } );

  //A.2: unusable operations are an error, in permissive mode the symbol is
  //used:
  const std::pair<const char *,VectS> badops[] = {
    { "missing operation", ops227m1 },
    { "unparsable operation", concat( ops227m1, { "x+y+z" } ) },
    { "origin shift", { "x,y,z", "-x+1/2,-y+1/2,-z+1/2" } } };
  for ( auto& e : badops )
    for ( auto perm : { Permissive::NO, Permissive::YES } )
      test( ( std::string( e.first ) + ", permissive="
              + ( perm == P ? "yes" : "no" ) ).c_str(),
            { cubic, e.second, "", "F d -3 m:2" }, perm );
  test( "missing operation, without symbol, permissive",
        { cubic, ops227m1 }, P );
  test( "missing operation, pymatgen file",
        { cubic, ops227m1, "", "", "", "# generated using PyMatGen\n" } );

  //A.3: operations and symbol of different space groups:
  for ( auto perm : { Permissive::NO, Permissive::YES } )
    test( perm == P ? "symbol mismatch, permissive=yes"
          : "symbol mismatch, permissive=no",
          { cubic, ops227, "", "F m -3 m" }, perm );
  test( "Hall symbol mismatch", { cubic, ops227, "-F 4 2 3", "F d -3 m" } );

  //B: no listed operations, the symbol (Hall, H-M or number):
  test( "Hall symbol", { cubic, {}, "-F 4vw 2vw 3" } );
  test( "H-M symbol with origin choice", { cubic, {}, "", "F d -3 m:1" } );
  test( "H-M symbol, one setting", { cubic, {}, "", "Fm-3m" } );
  test( "unrecognised Hall symbol, H-M used", { cubic, {}, "bla", "Fm-3m" } );
  test( "IT number, one setting", { cubic, {}, "", "", "225" } );
  test( "IT number, not a number", { cubic, {}, "", "", "22x" } );
  test( "IT number, out of range", { cubic, {}, "", "", "231" } );
  test( "R -3 m, hexagonal axes", { hexagonal, {}, "", "R -3 m" } );
  test( "R -3 m, rhombohedral axes", { rhombohedral, {}, "", "R -3 m" } );
  test( "R -3 m, unfitting cell",
        { "5.43 5.43 5.43 90 90 100", {}, "", "R -3 m" } );
  test( "ambiguous origin choice", { cubic, {}, "", "F d -3 m" } );
  test( "ambiguous IT number", { cubic, {}, "", "", "227" } );
  test( "no symmetry information", { cubic } );
  test( "unrecognised symbol", { cubic, {}, "", "bla" } );
  test( "null symbols skipped", { cubic, {}, "?", "." , "225" } );

  //C: explicit override (ignoring the symmetry information of the data):
  test( "override, H-M symbol", { cubic, ops227, "", "F d -3 m" },
        Permissive::NO, "F d -3 m:1" );
  test( "override, setting", { cubic, {}, "", "F d -3 m" },
        Permissive::NO, "227:1" );
  test( "override, ambiguous", { cubic }, Permissive::NO, "F d -3 m" );
  test( "override, unknown", { cubic }, Permissive::NO, "bla" );
  test( "override, no such setting", { cubic }, Permissive::NO, "227:3" );
  return 0;
}
