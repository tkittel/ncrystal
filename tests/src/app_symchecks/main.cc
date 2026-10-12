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

#include "NCrystal/NCrystal.hh"
#include "NCrystal/internal/extd_utils/NCFillHKL.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include <iostream>
#include <sstream>

//The library has two independent checks that a structure has the symmetry of
//its declared space group: the exact test of the atom positions under the
//symmetry operations (sgsym checkStructureSymmetry, used at Info creation)
//and the comparison of |F|^2 within families of symmetry-equivalent planes
//(the FillHKL self-check). The positions are loaded without a space group
//here, so that both checks can be applied with a claimed space group number
//and their verdicts compared. They must agree, except where the point
//symmetry holds but a translation does not (the |F|^2 check only needs the
//former).

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace {

  struct Verdicts {
    double posdev;//largest deviation from the position check
    bool fsq_ok;//did the |F|^2 self-check pass?
  };

  Verdicts check( const char * name, unsigned sg, const char * cell,
                  const std::string& positions )
  {
    //Load the positions without a space group (which does not trigger any
    //symmetry check), then apply both checks with the claimed space group:
    std::string ncmat = "NCMAT v7\n@CELL\n";
    ncmat += cell;
    ncmat += "\n@ATOMPOSITIONS\n" + positions
      + "@DYNINFO\n  element Al\n  fraction 1\n  type vdosdebye\n"
      "  debye_temp 400\n";
    const std::string fn = std::string( name ) + ".ncmat";
    NC::registerInMemoryFileData( fn, std::move( ncmat ) );
    auto info = NC::createInfo( fn + ";dcutoff=0.5;inelas=0" );
    NC::StructureInfo si = info->getStructureInfo();
    REQUIRE( si.spacegroup == 0 );
    si.spacegroup = sg;
    std::vector<std::vector<NC::Vector>> pos;
    for ( auto& ai : info->getAtomInfos() ) {
      pos.emplace_back();
      for ( auto& p : ai.unitCellPositions() )
        pos.back().push_back( p.as<NC::Vector>() );
    }
    Verdicts v;
    const bool rhombo = NC::usesRhombohedralAxes( sg, si.alpha );
    v.posdev = NC::checkStructureSymmetry( sg, rhombo, pos ).maxDeviation;
    NC::FillHKLCfg cfg;
    cfg.dcutoff = 0.5;
    try {
      NC::calculateHKLPlanes( si, info->getAtomInfos(), cfg );
      v.fsq_ok = true;
    } catch ( NC::Error::BadInput& e ) {
      REQUIRE( std::string( e.what() ).find( "not consistent with space group" )
               != std::string::npos );
      v.fsq_ok = false;
    }
    std::cout << "  " << name << " as SG-" << sg
              << ": position check deviation " << NC::fmt( v.posdev, "%.3g" )
              << ", |F|^2 self-check "
              << ( v.fsq_ok ? "passed" : "failed" ) << std::endl;
    return v;
  }

  const char * cubic = "  lengths 4.05 4.05 4.05\n  angles 90 90 90";
  const char * hexagonal = "  lengths 4 4 6\n  angles 90 90 120";

  std::string fcc( const char * moved = nullptr )
  {
    std::string s = "  Al 0 0 0\n  Al 0 1/2 1/2\n  Al 1/2 0 1/2\n";
    s += ( moved ? moved : "  Al 1/2 1/2 0\n" );
    return s;
  }

  std::string p312orbit()
  {
    //One general orbit of P312 (SG-149, Laue class -31m):
    const double x = 0.1, y = 0.3, z = 0.2;
    const double ops[6][3] = { {x,y,z}, {-y,x-y,z}, {-x+y,-x,z},
                               {-y,-x,-z}, {-x+y,y,-z}, {x,x-y,-z} };
    std::ostringstream ss;
    ss.precision( 12 );
    for ( auto& p : ops ) {
      ss << "  Al";
      for ( double c : p )
        ss << ' ' << c - std::floor( c );
      ss << '\n';
    }
    return ss.str();
  }
}

int main()
{
  auto v = check( "fcc", 225, cubic, fcc() );
  REQUIRE( v.posdev < 1e-12 && v.fsq_ok );

  v = check( "fcc_slightly_moved", 225, cubic, fcc( "  Al 0.5003 1/2 0\n" ) );
  REQUIRE( v.posdev > 5e-4 && v.posdev < 1e-2 && v.fsq_ok );

  v = check( "fcc_moved", 225, cubic, fcc( "  Al 0.53 1/2 0\n" ) );
  REQUIRE( v.posdev > 1e-2 && !v.fsq_ok );

  v = check( "p312", 149, hexagonal, p312orbit() );
  REQUIRE( v.posdev < 1e-12 && v.fsq_ok );

  v = check( "p312", 150, hexagonal, p312orbit() );
  REQUIRE( v.posdev > 1e-2 && !v.fsq_ok );

  //A single atom on the body-centre: invariant under the m-3m rotations, so
  //|F|^2 is consistent within the Fm-3m families, but not under the
  //F-centring translations (the orbit of the 4b site has four positions):
  v = check( "body_centre", 225, cubic, "  Al 1/2 1/2 1/2\n" );
  REQUIRE( v.posdev > 1e-2 && v.fsq_ok );

  //A P2 structure with its 2-fold axis along c: the position check finds
  //the setting 3:c, whose symmetry-equivalent reflections differ from the
  //unique-axis-b ones NCrystal assumes, and the |F|^2 check fails as well:
  v = check( "p2_axis_c", 3, "  lengths 5 6 7\n  angles 90 90 90",
             "  Al 0.1 0.2 0.3\n  Al 0.9 0.8 0.3\n" );
  REQUIRE( v.posdev < 1e-12 && !v.fsq_ok );
  const std::vector<std::vector<NC::Vector>> p2pos
    = { { NC::Vector( 0.1, 0.2, 0.3 ), NC::Vector( 0.9, 0.8, 0.3 ) } };
  const auto res = NC::checkStructureSymmetry( 3, false, p2pos );
  const NC::SpaceGroup setting{ res.settingIdx.value() };
  std::cout << "  p2_axis_c: position check finds setting " << setting
            << std::endl;
  REQUIRE( setting == NC::SpaceGroup( "3:c", NC::SpaceGroupParseMode::Strict ) );
  std::cout << "All checks agree as expected" << std::endl;
  return 0;
}
