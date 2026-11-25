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

#include "NCrystal/internal/cif/NCCIFSites.hh"
#include "NCrystal/internal/cif/NCCIFSetting.hh"
#include <iostream>

//Tests of CIF::getAtoms: the sites of CIF data (with positive occupancy)
//expanded into their orbits, with their elements and a single magnitude of
//their displacements, merged when at the same positions.

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace {

  using Permissive = NCC::CIFPermissive;

  //CIF data in the space group sg and cell, with the site rows (label, type
  //symbol, x, y, z, occupancy, Uiso) and extra data:
  std::string cifData( const std::string& sites, const std::string& extra = "",
                       const char * sg = "P -1",
                       const char * cell = "4 4 4 90 90 90" )
  {
    std::string res = "data_test\n_space_group_name_H-M_alt '";
    res += sg;
    res += "'\n";
    const char * names[] = { "length_a", "length_b", "length_c",
                             "angle_alpha", "angle_beta", "angle_gamma" };
    auto cellvals = NC::StrView( cell ).split();
    for ( auto i : NC::ncrange( 6 ) )
      res += std::string( "_cell_" ) + names[i] + " "
        + NC::vectAt( cellvals, i ).to_string() + "\n";
    return res + ( "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
                   "_atom_site_fract_x\n_atom_site_fract_y\n"
                   "_atom_site_fract_z\n_atom_site_occupancy\n"
                   "_atom_site_U_iso_or_equiv\n" ) + sites + extra;
  }

  void test( const char * descr, const std::string& data,
             Permissive perm = Permissive::NO, double tolSymprec = 0.01 )
  {
    std::vector<std::string> warnings;
    const NCC::CIFWarnFct warnfct = [&warnings]( const std::string& w )
    {
      warnings.push_back( w );
    };
    std::cout << descr << ":";
    try {
      auto doc = NCC::Document::parse(
        NC::makeSO<const NC::TextData>( NC::RawStrData( std::string( data ) ),
                                        NC::TextData::DataType{ "cif" },
                                        NC::DataSourceName( "test.cif" ) ) );
      auto items = NCC::findUniqueStructProvItems( doc, perm, warnfct );
      const auto sg = NCC::chooseSetting( doc, items.structureItems, perm,
                                          NC::StrView(), warnfct );
      const auto res = NCC::getAtoms( doc, items.structureItems, sg, perm,
                                      1e-3, tolSymprec, warnfct );
      const auto& atoms = res.atoms;
      if ( res.setting != sg )
        std::cout << " setting " << res.setting.toString();
      std::cout << '\n';
      for ( auto& a : atoms ) {
        std::cout << "   ";
        for ( auto& l : a.labels )
          std::cout << ( &l == &a.labels.front() ? "" : "+" ) << l;
        std::cout << ":";
        for ( auto& c : a.composition )
          std::cout << ' ' << NC::fmt( c.first ) << ' ' << c.second;
        std::cout << ", occupancy " << NC::fmt( a.occupancy )
                  << ", " << a.positions.size() << " positions (first "
                  << NC::fmt( a.positions.front().x() ) << ' '
                  << NC::fmt( a.positions.front().y() ) << ' '
                  << NC::fmt( a.positions.front().z() ) << ")";
        if ( a.uiso.has_value() )
          std::cout << ", uiso " << NC::fmt( a.uiso.value(), "%.6g" );
        if ( a.uij.has_value() ) {
          std::cout << ", uij";
          for ( auto u : a.uij.value() )
            std::cout << ' ' << NC::fmt( u, "%.6g" );
          nc_assert_always( a.uijs.size() == a.positions.size() );
        }
        if ( a.splitPositions )
          std::cout << ", split over " << a.splitPositions << " positions "
                    << NC::fmt( a.splitDistance, "%.3g" ) << " Aa apart";
        std::cout << '\n';
      }
    } catch ( NC::Error::BadInput& e ) {
      std::cout << " rejected: " << e.what() << std::endl;
    }
    for ( auto& w : warnings )
      std::cout << "   warning: " << w << std::endl;
  }

  std::string aniso( const std::string& rows, char form = 'U' )
  {
    std::string res = "loop_\n_atom_site_aniso_label\n";
    for ( auto ij : { "11", "22", "33", "12", "13", "23" } )
      res += std::string( "_atom_site_aniso_" ) + form + "_" + ij + "\n";
    return res + rows;
  }
}

int main()
{
  const auto P = Permissive::YES;
  std::cout << "==== Sites ====\n";
  test( "separate sites (general positions of P-1)",
        cifData( "Na1 Na 0.1 0.2 0.3 1 0.01\nCl1 Cl 0.3 0.1 0.2 1 0.01\n" ) );
  test( "unoccupied site skipped, occupancy absent or ?",
        cifData( "Na1 Na 0.1 0.2 0.3 0 0.01\nCl1 Cl 0.3 0.1 0.2 ? 0.01\n" ) );
  test( "invalid occupancy",
        cifData( "Na1 Na 0.1 0.2 0.3 x 0.01\n" ) );
  test( "occupancy above 1",
        cifData( "Na1 Na 0.1 0.2 0.3 1.01 0.01\n" ) );
  test( "occupancy slightly above 1 (rounding)",
        cifData( "Na1 Na 0.1 0.2 0.3 1.0000004 0.01\n" ) );
  test( "invalid coordinates",
        cifData( "Na1 Na 0.1 y 0.3 1 0.01\n" ) );
  test( "no sites with positive occupancy",
        cifData( "Na1 Na 0.1 0.2 0.3 0 0.01\n" ) );
  {
    //(Biso column added to the site loop):
    std::string d = cifData( "Na1 Na 0.1 0.2 0.3 1 x 0.8\n"
                             "Cl1 Cl 0.3 0.1 0.2 1 0 1.6(2)\n"
                             "K1 K 0.2 0.3 0.1 1 0.02 1.6\n" );
    d.replace( d.find( "_atom_site_U_iso_or_equiv\n" ),
               std::strlen( "_atom_site_U_iso_or_equiv\n" ),
               "_atom_site_U_iso_or_equiv\n_atom_site_B_iso_or_equiv\n" );
    test( "Biso used where Uiso is invalid or zero (B/8pi^2)", d );
  }
  std::cout << "==== Elements ====\n";
  test( "from type symbols (also isotopes)",
        cifData( "Na1 Na+ 0.1 0.2 0.3 1 0.01\nH1 D 0.3 0.1 0.2 1 0.01\n" ) );
  test( "from labels (with the formula)",
        cifData( "CaM1 ? 0.1 0.2 0.3 1 0.01\nO1 . 0.3 0.1 0.2 1 0.01\n",
                 "_chemical_formula_sum 'Ca O'\n" ) );
  test( "from labels (without formula), several examples",
        cifData( "Na1 ? 0.1 0.2 0.3 1 0.01\nCl1 ? 0.3 0.1 0.2 1 0.01\n"
                 "O(1) ? 0.2 0.3 0.1 1 0.01\nFe2 ? 0.4 0.4 0.1 1 0.01\n"
                 "Al ? 0.1 0.4 0.4 1 0.01\nNa2 ? 0.2 0.1 0.4 1 0.01\n" ) );
  test( "label not matching the formula",
        cifData( "Xx1 ? 0.1 0.2 0.3 1 0.01\n",
                 "_chemical_formula_sum 'Ca O'\n" ) );
  test( "label without formula",
        cifData( "CaM1 ? 0.1 0.2 0.3 1 0.01\n" ) );
  std::cout << "==== Anisotropic displacements ====\n";
  const std::string twosites = "Na1 Na 0.1 0.2 0.3 1 ?\nCl1 Cl 0.3 0.1 0.2 1"
    " 0.0101\n";
  test( "tensors (uiso from the Ueq, or the Uiso if consistent)",
        cifData( twosites, aniso( "Na1 0.012 0.009 0.015 0.002 -0.003 0.001\n"
                                  "Cl1 0.01 0.0101 0.0102 0 0 0\n" ) ) );
  test( "Uiso inconsistent with the Ueq",
        cifData( twosites, aniso( "Cl1 0.02 0.02 0.02 0 0 0\n" ) ) );
  test( "Uiso inconsistent with the Ueq, permissive",
        cifData( twosites, aniso( "Cl1 0.02 0.02 0.02 0 0 0\n" ) ), P );
  test( "all-zero tensor (placeholder for none)",
        cifData( twosites, aniso( "Cl1 0 0 0 0 0 0\n" ) ) );
  test( "invalid tensor values",
        cifData( twosites, aniso( "Cl1 0.01 0.01 x 0 0 0\n" ) ) );
  test( "tensor breaking the site symmetry",
        cifData( "Na1 Na 0 0 0 1 ?\n",
                 aniso( "Na1 0.0100 0.0200 0.0300 0.0010 0 0\n" ),
                 "P m -3 m" ) );
  test( "tensors given as B (converted to U, B = 8 pi^2 U)",
        cifData( twosites, aniso( "Na1 1.6 0.8 0.8 0 0 0.1\n"
                                  "Cl1 0.8 0.8 0.8 0 0 0\n", 'B' ) ) );
  test( "rows matched by label in order (labels reused)",
        cifData( "O1 O 0.1 0.2 0.3 1 ?\nO1 O 0.4 0.2 0.3 1 ?\n"
                 "Si1 Si 0.6 0.2 0.3 1 ?\n",
                 aniso( "O1 0.011 0.011 0.011 0 0 0\n"
                        "Si1 0.021 0.021 0.021 0 0 0\n"
                        "O1 0.031 0.031 0.031 0 0 0\n" ) ) );
  std::cout << "==== Merging ====\n";
  test( "mixed site (other position of the same orbit)",
        cifData( "Al1 Al 0.1 0.2 0.3 0.75 0.01\n"
                 "Cr1 Cr 0.9 0.8 0.7 0.25 0.01\n" ) );
  test( "partially occupied mixed site",
        cifData( "O1 O 0.1 0.2 0.3 0.4 0.01\nF1 F 0.1 0.2 0.3 0.4 0.01\n" ) );
  test( "rounding of total occupancy slightly above 1",
        cifData( "A Al 0.1 0.2 0.3 0.5000004 0.01\n"
                 "B Cr 0.1 0.2 0.3 0.5000004 0.01\n" ) );
  test( "total occupancy above 1",
        cifData( "A Al 0.1 0.2 0.3 0.7 0.01\nB Cr 0.1 0.2 0.3 0.7 0.01\n" ) );
  test( "same site given twice (total occupancy 2)",
        cifData( "A Al 0.1 0.2 0.3 1 0.01\nB Al 0.9 0.8 0.7 1 0.01\n" ) );
  test( "different Uiso at the same positions",
        cifData( "A Al 0.1 0.2 0.3 0.5 0.01\nB Cr 0.1 0.2 0.3 0.5 0.02\n" ) );
  test( "different tensors at the same positions",
        cifData( "A Al 0.1 0.2 0.3 0.5 ?\nB Cr 0.1 0.2 0.3 0.5 ?\n",
                 aniso( "A 0.01 0.01 0.01 0 0 0\n"
                        "B 0.02 0.02 0.02 0 0 0\n" ) ) );
  test( "positions at the cell boundary (periodic)",
        cifData( "A Al 0.99999 0.2 0.3 0.5 0.01\n"
                 "B Cr 0.00001 0.2 0.3 0.5 0.01\n" ) );
  const std::string counted = "loop_ _atom_type_symbol"
    " _atom_type_number_in_cell\nNa {}\nCl 2\n";
  for ( auto n : { "2", "2.01", "3", "?" } ) {
    std::string extra = counted;
    extra.replace( extra.find( "{}" ), 2, n );
    test( ( std::string( "atom count check, Na " ) + n ).c_str(),
          cifData( "Na1 Na 0.1 0.2 0.3 1 0.01\nCl1 Cl 0.3 0.1 0.2 1 0.01\n",
                   extra ) );
  }
  std::cout << "==== Higher symmetry in the same coordinates ====\n";
  //(P1 data with the inversion symmetry of P-1):
  const std::string p1bar = "Na1 Na 0.1 0.2 0.3 1 {u1}\n"
    "Na2 Na 0.9 0.8 0.7 1 {u2}\nCl1 Cl 0.3 0.1 0.2 1 0.01\n"
    "Cl2 Cl 0.7 0.9 0.8 1 0.01\n";
  auto withU = [&p1bar]( const char * u1, const char * u2 )
  {
    std::string r = p1bar;
    r.replace( r.find( "{u1}" ), 4, u1 );
    r.replace( r.find( "{u2}" ), 4, u2 );
    return r;
  };
  test( "P1 data with the symmetry of P-1",
        cifData( withU( "0.01", "0.01" ), "", "P 1" ) );
  test( "same, Uiso within tol_uij (mean used)",
        cifData( withU( "0.01", "0.010005" ), "", "P 1" ) );
  test( "same, different Uiso",
        cifData( withU( "0.01", "0.02" ), "", "P 1" ) );
  test( "same, different Uiso, permissive",
        cifData( withU( "0.01", "0.02" ), "", "P 1" ), P );
  test( "same, Uiso of only one site",
        cifData( withU( "0.01", "?" ), "", "P 1" ) );
  test( "same, tensors which are symmetry images (inversion: equal)",
        cifData( withU( "?", "?" ), aniso( "Na1 0.01 0.012 0.014 0.001 0 0\n"
                                         "Na2 0.01 0.012 0.014 0.001 0 0\n" ),
                 "P 1" ) );
  test( "same, tensors which are not symmetry images",
        cifData( withU( "?", "?" ), aniso( "Na1 0.01 0.012 0.014 0.001 0 0\n"
                                         "Na2 0.014 0.012 0.01 0.001 0 0\n" ),
                 "P 1" ) );
  test( "inversion symmetry only within 0.04 Aa (beyond tol_symprec)",
        cifData( "Na1 Na 0.1 0.2 0.3 1 0.01\nNa2 Na 0.91 0.8 0.7 1 0.01\n",
                 "", "P 1" ) );
  test( "same, with tol_symprec 0.1 Aa",
        cifData( "Na1 Na 0.1 0.2 0.3 1 0.01\nNa2 Na 0.91 0.8 0.7 1 0.01\n",
                 "", "P 1" ), Permissive::NO, 0.1 );
  test( "P1 data with the symmetry of P-1, not looked for",
        cifData( withU( "0.01", "0.01" ), "", "P 1" ), Permissive::NO, 0.0 );
  std::cout << "==== Split sites ====\n";
  //(two positions 0.054 Aa apart in Pm):
  const char * pm = "P 1 m 1";
  const char * pmcell = "5 5.4 6 90 90 90";
  test( "split site", cifData( "A Al 0.1 0.005 0.2 0.5 0.01\n", "", pm,
                               pmcell ) );
  test( "split site, permissive",
        cifData( "A Al 0.1 0.005 0.2 0.5 0.01\n", "", pm, pmcell ), P );
  test( "split site, total occupancy 0.5 of two elements, permissive",
        cifData( "A Al 0.1 0.005 0.2 0.25 0.01\nB Cr 0.1 0.005 0.2 0.25 0.01\n",
                 "", pm, pmcell ), P );
  test( "split site, total occupancy 0.6 of two elements",
        cifData( "A Al 0.1 0.005 0.2 0.3 0.01\nB Cr 0.1 0.005 0.2 0.3 0.01\n",
                 "", pm, pmcell ), P );
  return 0;
}
