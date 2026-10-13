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
#include "NCrystal/internal/infobld/NCInfoBuilder.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
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
  const NC::SpaceGroup setting{ NC::SpaceGroupSettingIdx{ res.settingIdx } };
  std::cout << "  p2_axis_c: position check finds setting " << setting
            << std::endl;
  REQUIRE( setting == NC::SpaceGroup( "3:c", NC::SpaceGroupParseMode::Strict ) );
  std::cout << "All checks agree as expected" << std::endl;

  //The setting found by the check is stored on the StructureInfo (and
  //without positions, as for LAZ/LAU data, the one implied by the number and
  //cell, if unique):
  {
    using SG = NC::SpaceGroup;
    struct { SG operator()( const char * spec ) const { return SG( spec, NC::SpaceGroupParseMode::Relaxed ); } } sgr;
    auto asSG = []( NC::SpaceGroupSettingIdx hn )
    {
      std::cout << "  stored setting: "
                << ( hn.get() ? SG( hn ).toString() : "unknown" ) << std::endl;
      return hn.get() ? SG( hn ) : SG( "1", NC::SpaceGroupParseMode::Relaxed );
    };
    auto settingOf = [&asSG]( const char * name, unsigned sg,
                              const std::string& positions )
    {
      std::string ncmat = "NCMAT v7\n@CELL\n  lengths 5.43 5.43 5.43\n"
        "  angles 90 90 90\n@SPACEGROUP\n  " + std::to_string( sg )
        + "\n@ATOMPOSITIONS\n" + positions
        + "@DYNINFO\n  element Si\n  fraction 1\n  type vdosdebye\n"
        "  debye_temp 400\n";
      const std::string fn = std::string( name ) + "_setting.ncmat";
      NC::registerInMemoryFileData( fn, std::move( ncmat ) );
      auto info = NC::createInfo( fn + ";inelas=0" );
      return asSG( info->getStructureInfo().spacegroup_setting );
    };
    auto settingNoPos = [&asSG]( unsigned sg, double a, double c,
                                 double alpha, double gamma )
    {
      //Built directly, since NCMAT files always have positions:
      NC::InfoBuilder::SinglePhaseBuilder bldr;
      bldr.dataSourceName = NC::DataSourceName( "nopos" );
      NC::InfoBuilder::UnitCell uc;
      uc.structinfo.spacegroup = sg;
      uc.structinfo.lattice_a = uc.structinfo.lattice_b = a;
      uc.structinfo.lattice_c = c;
      uc.structinfo.alpha = uc.structinfo.beta = alpha;
      uc.structinfo.gamma = gamma;
      uc.structinfo.n_atoms = 2;
      bldr.unitcell = std::move( uc );
      bldr.composition = NC::createInfo( "Si_sg227.ncmat" )->getComposition();
      bldr.hklPlanes.emplace();
      bldr.hklPlanes.value().dspacingRange = NC::PairDD( 1.0, 2.0 );
      bldr.hklPlanes.value().source = NC::HKLList();
      auto info = NC::InfoBuilder::buildInfoPtr( std::move( bldr ) );
      return asSG( info->getStructureInfo().spacegroup_setting );
    };
    const std::string diamond1 = "  Si 0 0 0\n  Si 0 1/2 1/2\n  Si 1/2 0 1/2\n"
      "  Si 1/2 1/2 0\n  Si 1/4 1/4 1/4\n  Si 1/4 3/4 3/4\n  Si 3/4 1/4 3/4\n"
      "  Si 3/4 3/4 1/4\n";
    const std::string diamond2 = "  Si 1/8 1/8 1/8\n  Si 1/8 5/8 5/8\n"
      "  Si 5/8 1/8 5/8\n  Si 5/8 5/8 1/8\n  Si 7/8 7/8 7/8\n  Si 7/8 3/8 3/8\n"
      "  Si 3/8 7/8 3/8\n  Si 3/8 3/8 7/8\n";
    REQUIRE( settingOf( "diamond1", 227, diamond1 ) == sgr( "227" ) );
    REQUIRE( settingOf( "diamond2", 227, diamond2 ) == sgr( "227:2" ) );
    const std::string fccSi = "  Si 0 0 0\n  Si 0 1/2 1/2\n  Si 1/2 0 1/2\n"
      "  Si 1/2 1/2 0\n";
    REQUIRE( settingOf( "fcc", 225, fccSi ) == sgr( "225" ) );
    //Without positions: unique settings are known, others not (the two
    //settings of an R number are told apart by the cell):
    REQUIRE( settingNoPos( 225, 4, 4, 90, 90 ) == sgr( "225" ) );
    REQUIRE( settingNoPos( 227, 4, 4, 90, 90 ) == sgr( "1" ) );
    REQUIRE( settingNoPos( 166, 4, 10, 90, 120 ) == sgr( "166" ) );
    REQUIRE( settingNoPos( 166, 4, 4, 60, 60 ) == sgr( "166:r" ) );
  }

  //Sites: derived from positions, or positions derived from sites (the
  //diamond structure as its single site in both origin choices):
  {
    auto info = NC::createInfo( "Al_sg225.ncmat;inelas=0" );
    const auto& ai = info->getAtomInfos().front();
    REQUIRE( ai.sites().size() == 1 && ai.sites().front().multiplicity == 4 );
    REQUIRE( ai.sites().front().position == NC::AtomInfo::Pos( 0.0, 0.0, 0.0 ) );
    std::cout << "Al: 1 site (0,0,0) of multiplicity 4" << std::endl;
    for ( const char * site : { "0 0 0", "1/8 1/8 1/8" } ) {
      NC::InfoBuilder::SinglePhaseBuilder bldr;
      bldr.dataSourceName = NC::DataSourceName( "diamond_from_site" );
      NC::InfoBuilder::UnitCell uc;
      uc.structinfo.spacegroup = 227;
      //(a site at 1/8,1/8,1/8 is the 8a position only in origin choice 2;
      //in the default setting it would be the 16c position):
      if ( NC::StrView( site ) != "0 0 0" )
        uc.structinfo.spacegroup_setting = NC::SpaceGroup( "227:2", NC::SpaceGroupParseMode::Strict ).settingIdx();
      uc.structinfo.lattice_a = uc.structinfo.lattice_b = uc.structinfo.lattice_c = 5.43;
      uc.structinfo.alpha = uc.structinfo.beta = uc.structinfo.gamma = 90.0;
      NC::AtomInfo::AtomSites sites( 1 );
      const auto parts = NC::StrView( site ).split();
      for ( unsigned i = 0; i < 3; ++i ) {
        const auto frac = parts.at(i).split('/');
        sites.front().position[i] = ( frac.size() == 2
                                      ? frac.at(0).toDbl().value() / frac.at(1).toDbl().value()
                                      : frac.at(0).toDbl().value() );
      }
      uc.atomlist.emplace();
      uc.atomlist.value().emplace_back( NC::createInfo( "Si_sg227.ncmat" )
                                        ->getAtomInfos().front().indexedAtomData(),
                                        std::move( sites ), NC::NullOpt, NC::NullOpt );
      bldr.unitcell = std::move( uc );
      bldr.hklPlanes.emplace();
      bldr.hklPlanes.value().dspacingRange = NC::PairDD( 1.0, 2.0 );
      bldr.hklPlanes.value().source = NC::HKLList();
      auto info2 = NC::InfoBuilder::buildInfoPtr( std::move( bldr ) );
      const auto& ai2 = info2->getAtomInfos().front();
      const NC::SpaceGroup setting2( info2->getStructureInfo().spacegroup_setting );
      std::cout << "Diamond from site (" << site << "): " << ai2.numberPerUnitCell()
                << " positions, " << ai2.sites().size() << " site of multiplicity "
                << ai2.sites().front().multiplicity << ", setting " << setting2
                << std::endl;
      REQUIRE( ai2.numberPerUnitCell() == 8 && ai2.sites().size() == 1
               && ai2.sites().front().multiplicity == 8
               && info2->getStructureInfo().n_atoms == 8 );
      REQUIRE( setting2 == NC::SpaceGroup( NC::StrView( site ) == "0 0 0"
                                           ? "227:1" : "227:2",
                                           NC::SpaceGroupParseMode::Strict ) );
    }
  }
  return 0;
}
