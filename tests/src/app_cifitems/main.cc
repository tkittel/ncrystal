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


#include "NCrystal/internal/cif/NCCIFItems.hh"
#include "NCTestUtils/NCTestExpect.hh"
#include "NCTestUtils/NCTestFindData.hh"
#include <fstream>
#include <iostream>
#include <sstream>

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  using Permissive = NCC::CIFPermissive;

  //Warnings are passed to this function (and so appear in the log):
  const NCC::CIFWarnFct warnfct = []( const std::string& w )
  {
    std::cout << "WARNING-CALLBACK: " << w << std::endl;
  };

  NCC::Document parse( const std::string& s,
                       const char* name = "test.cif" )
  {
    return NCC::Document::parse(
      NC::makeSO<const NC::TextData>( NC::RawStrData( std::string( s ) ),
                                      NC::TextData::DataType{ "cif" },
                                      NC::DataSourceName( name ) ) );
  }

  std::string readTestData( const char* fn )
  {
    std::ifstream f( nctest::find_test_data( fn ) );
    REQUIRE( f.good() );
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
  }

  std::string str( const NCC::Value& v )
  {
    return "[" + v.text.to_string() + "]@" + std::to_string( v.line );
  }
  std::string str( const NC::Optional<NCC::Value>& v )
  {
    return v.has_value() ? str( v.value() ) : std::string( "-" );
  }
  std::string str( const NCC::ValueColumn& c )
  {
    if ( c.empty() )
      return "-";
    std::string res;
    for ( auto i : NC::ncrange( c.size() ) )
      res += ( i ? " " : "" ) + str( c[i] );
    return res;
  }
  std::string str( const std::vector<NCC::Value>& v )
  {
    if ( v.empty() )
      return "-";
    std::string res;
    for ( auto i : NC::ncrange( v.size() ) )
      res += ( i ? " " : "" ) + str( v[i] );
    return res;
  }
  std::string str( const NC::Optional<NCC::TaggedValue>& v )
  {
    return ( v.has_value() ? v.value().tag.to_string() + "="
             + str( v.value().value ) : std::string( "-" ) );
  }
  template<class T>
  std::string strAlt( const T& v )
  {
    std::string res;
    for ( auto& e : v )
      res += ( res.empty() ? "" : " | " ) + e.tag.to_string() + "="
        + str( e.value );
    return res.empty() ? "-" : res;
  }
  std::string strOps( const NC::SmallVector<NCC::TaggedColumn,2>& v )
  {
    std::string res;
    for ( auto& e : v )
      res += ( res.empty() ? "" : " | " ) + e.tag.to_string() + "="
        + str( e.values );
    return res.empty() ? "-" : res;
  }

  void dump( const NCC::StructureItems& s )
  {
    std::cout << "  cell:";
    for ( auto& v : s.cell )
      std::cout << " " << str( v );
    std::cout << "\n  sites: label " << str( s.sites.label )
              << "\n    x " << str( s.sites.fractX )
              << "\n    y " << str( s.sites.fractY )
              << "\n    z " << str( s.sites.fractZ )
              << "\n    type " << str( s.sites.typeSymbol )
              << "\n    uiso " << str( s.sites.uIsoOrEquiv )
              << "\n    biso " << str( s.sites.bIsoOrEquiv )
              << "\n    occ " << str( s.sites.occupancy ) << "\n";
    if ( s.aniso.has_value() ) {
      const auto& a = s.aniso.value();
      std::cout << "  aniso ("
                << ( a.form == NCC::StructureItems::AnisoForm::U ? "U" : "B" )
                << "): label " << str( a.label ) << "\n";
      for ( auto& c : a.components )
        std::cout << "    " << str( c ) << "\n";
      std::cout << "    skipped rows:";
      for ( auto r : a.skippedRows )
        std::cout << " " << r;
      std::cout << "\n";
    }
    std::cout << "  operations: " << strOps( s.operations )
              << "\n  hm: " << strAlt( s.hmSymbol )
              << "\n  hall: " << strAlt( s.hallSymbol )
              << "\n  number: " << strAlt( s.itNumber )
              << "\n  formula: " << str( s.formulaSum )
              << "\n  atom type numbers: " << str( s.atomTypeNumberInCell )
              << "\n  temperatures: " << str( s.ambientTemperature ) << " "
              << str( s.measurementTemperature ) << "\n";
  }

  void dump( const NCC::ProvenanceItems& p, const NCC::Document& doc )
  {
    const auto ob = p.otherBlock( doc );
    std::cout << "  (block data_" << p.block( doc ).name() << ", other block "
              << ( ob.has_value() ? "data_" + ob.value().name().to_string()
                   : std::string( "-" ) ) << ")"
              << "\n  title: " << str( p.title )
              << "\n  authors: "
              << ( p.authors.has_value()
                   ? p.authors.value().tag.to_string() + "="
                   + str( p.authors.value().values ) : std::string( "-" ) )
              << "\n  journal: " << str( p.journal )
              << "\n  year: " << str( p.year )
              << "\n  doi: " << str( p.doi )
              << "\n  cod: " << str( p.codDatabaseCode )
              << "\n  mp: " << str( p.mpDatabaseCode )
              << "\n  audit: " << str( p.auditCreationMethod ) << " "
              << str( p.auditCreationDate ) << " " << str( p.auditAuthors )
              << "\n  original formula: " << str( p.codOriginalFormulaSum )
              << "\n";
  }

  //Items of the unique structure block:
  void show( const char* title, const std::string& s,
             Permissive perm = Permissive::NO, const char* name = "test.cif" )
  {
    std::cout << title << ":\n";
    auto doc = parse( s, name );
    auto r = NCC::findUniqueStructProvItems( doc, perm, warnfct );
    std::cout << "  (block data_" << r.structureItems.block( doc ).name()
              << ")\n";
    dump( r.structureItems );
    if ( r.provenanceItems.has_value() ) {
      std::cout << " provenance:\n";
      dump( r.provenanceItems.value(), doc );
    } else {
      std::cout << " provenance: none\n";
    }
  }

  void ensureItemsError( const char* title, const std::string& s,
                         Permissive perm = Permissive::NO,
                         const char* name = "test.cif" )
  {
    std::cout << title << ": ";
    nctest::ensureError<NC::Error::BadInput>( [&s,perm,name]()
    {
      auto doc = parse( s, name );
      NCC::findUniqueStructProvItems( doc, perm, warnfct );
    } );
  }

  const std::string cell = "_cell_length_a 5.43(1)\n_cell_length_b 5.43\n"
    "_cell_length_c 6\n_cell_angle_alpha 90\n_cell_angle_beta 90.0\n"
    "_cell_angle_gamma 120\n";
  const std::string site = "loop_\n_atom_site_label\n_atom_site_fract_x\n"
    "_atom_site_fract_y\n_atom_site_fract_z\nAl1 0 0 0\nO1 0.25 0.25 0.25\n";

  void testStructure()
  {
    std::cout << "==== Structure ====\n";
    show( "minimal (absent angles)", "data_min\n_cell_length_a 4\n"
          "_cell_length_b 4\n_cell_length_c 4\n" + site );
    show( "full", "data_full\n" + cell
          + "_space_group_name_H-M_alt 'P 63/m m c'\n"
          "_symmetry_space_group_name_H-M 'P 63/m m c'\n"
          "_space_group_IT_number 194\n_symmetry_Int_Tables_number 149\n"
          "_symmetry_space_group_name_Hall '-P 6c 2c'\n"
          "loop_ _space_group_symop_operation_xyz x,y,z -x,-y,-z\n"
          "loop_ _symmetry_equiv_pos_as_xyz 'x, y, z' '-x, -y, -z'\n"
          "_chemical_formula_sum 'Al O'\n_diffrn_ambient_temperature 293(2)\n"
          "_cell_measurement_temperature ?\n"
          "loop_ _atom_type_symbol _atom_type_number_in_cell Al 2 O 2\n"
          "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
          "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
          "_atom_site_U_iso_or_equiv\n_atom_site_occupancy\n"
          "Al1 Al3+ 0 0 0 0.01(1) 1\nO1 O2- 0.25 0.25 0.25 ? .\n"
          "_publ_section_title\n;\nA title\n;\n"
          "loop_ _publ_author_name 'Doe, J.' 'Roe, R.'\n"
          "_journal_name_full 'J. Test'\n_journal_year 2001\n"
          "_journal_paper_doi 10.1000/x\n_cod_database_code 1234567\n"
          "_audit_creation_method 'by hand'\n_audit_creation_date 2020-01-01\n"
          "_cod_original_formula_sum 'Al2 O2'\n"
          "# Note from NCrystal.cifutils: Data from materialsproject.org"
          " / mp-149\n" );
    show( "single site (no loop), same symmetry operations in both tags",
          "data_one\n" + cell + "_atom_site_label Fe\n_atom_site_fract_x 0\n"
          "_atom_site_fract_y 0\n_atom_site_fract_z 1/2\n"
          "_space_group_symop_operation_xyz x,y,z\n"
          "_symmetry_equiv_pos_as_xyz x,y,z\n" );
    show( "only deprecated symmetry tags", "data_dep\n" + cell + site
          + "_symmetry_space_group_name_H-M 'P 1'\n_symmetry_Int_Tables_number"
          " 1\nloop_ _symmetry_equiv_pos_as_xyz x,y,z\n" );
  }

  void testAniso()
  {
    std::cout << "==== Aniso ====\n";
    const std::string u = "loop_\n_atom_site_aniso_label\n"
      "_atom_site_aniso_U_11\n_atom_site_aniso_U_22\n_atom_site_aniso_U_33\n"
      "_atom_site_aniso_U_12\n_atom_site_aniso_U_13\n_atom_site_aniso_U_23\n";
    show( "U", "data_u\n" + cell + site + u + "O1 1 2 3 4 5 6\n" );
    show( "B in the loop of the sites", "data_b\n" + cell
          + "loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n"
          "_atom_site_fract_z\n_atom_site_aniso_B_11\n_atom_site_aniso_B_22\n"
          "_atom_site_aniso_B_33\n_atom_site_aniso_B_12\n"
          "_atom_site_aniso_B_13\n_atom_site_aniso_B_23\n"
          "Al1 0 0 0 1 1 1 0 0 0\nO1 0.5 0.5 0.5 ? ? ? ? ? ?\n",
          Permissive::YES );
    const std::string allnull = "data_n\n" + cell + site + u
      + "Al1 1 2 3 4 5 6\nO1 ? . ? ? ? ?\n";
    ensureItemsError( "all unknown (strict)", allnull );
    show( "all unknown (permissive)", allnull, Permissive::YES );
    const std::string unknown = "data_x\n" + cell + site + u
      + "Al2 1 2 3 4 5 6\nAl1 1 2 3 4 5 6\n";
    ensureItemsError( "unknown label (strict)", unknown );
    show( "unknown label (permissive)", unknown, Permissive::YES );
    const std::string partial = "data_x\n" + cell + site
      + "loop_\n_atom_site_aniso_label\n_atom_site_aniso_U_33\nAl1 1\n";
    ensureItemsError( "missing component (strict)", partial );
    show( "missing component (permissive)", partial, Permissive::YES );
    ensureItemsError( "beta", "data_x\n" + cell + site
                      + "loop_ _atom_site_aniso_label _atom_site_aniso_beta_11"
                      " Al1 1\n" );
    ensureItemsError( "U and B", "data_x\n" + cell + site + u
                      + "Al1 1 2 3 4 5 6\n_atom_site_aniso_B_11 1\n" );
    const std::string nolabel = "data_x\n" + cell + site
      + "loop_\n_atom_site_aniso_type_symbol\n_atom_site_aniso_U_11\n"
      "_atom_site_aniso_U_22\n_atom_site_aniso_U_33\n_atom_site_aniso_U_12\n"
      "_atom_site_aniso_U_13\n_atom_site_aniso_U_23\nAl 1 2 3 4 5 6\n";
    ensureItemsError( "no aniso label (strict)", nolabel );
    show( "no aniso label (permissive)", nolabel, Permissive::YES );
  }

  void testErrors()
  {
    std::cout << "==== Errors ====\n";
    ensureItemsError( "no structure", "data_a _cell_length_a 4\n" );
    ensureItemsError( "two structures", "data_a\n" + cell + site
                      + "data_b\n" + cell + site );
    ensureItemsError( "broken second structure", "data_a\n" + cell + site
                      + "data_b\n_atom_site_label X\n" );
    ensureItemsError( "missing cell length", "data_a\n_cell_length_a 4\n"
                      "_cell_length_b 4\n" + site );
    ensureItemsError( "unknown cell length", "data_a\n_cell_length_a ?\n"
                      "_cell_length_b 4\n_cell_length_c 4\n"
                      "_atom_sites_solution_primary ?\n" );
    ensureItemsError( "unknown cell angle", "data_a\n_cell_length_a 4\n"
                      "_cell_length_b 4\n_cell_length_c 4\n"
                      "_cell_angle_beta .\n" + site );
    ensureItemsError( "cell length in a loop", "data_a\n_cell_length_b 4\n"
                      "_cell_length_c 4\nloop_ _cell_length_a 4 5\n" + site );
    ensureItemsError( "missing coordinate", "data_a\n" + cell
                      + "loop_ _atom_site_label _atom_site_fract_x "
                      "_atom_site_fract_y Al 0 0\n" );
    ensureItemsError( "unknown coordinate", "data_a\n" + cell
                      + "loop_ _atom_site_label _atom_site_fract_x "
                      "_atom_site_fract_y _atom_site_fract_z Al 0 ? 0\n" );
    ensureItemsError( "column in another loop", "data_a\n" + cell + site
                      + "loop_ _atom_site_occupancy 1 1\n" );
    ensureItemsError( "symmetry symbol in a loop", "data_a\n" + cell + site
                      + "loop_ _space_group_name_Hall 'P 1' 'P 2'\n" );
  }

  void testProvenance()
  {
    std::cout << "==== Provenance ====\n";
    //(values are case-sensitive, so "PRIMARY" is not the primary citation):
    const std::string cit = "loop_\n_citation_id\n_citation_title\n"
      "_citation_journal_full\n_citation_year\n_citation_doi\n"
      "other 'Other title' 'Other J.' 1990 10.1/other\n"
      "PRIMARY 'Upper title' 'Upper J.' 1998 10.1/upper\n"
      "primary 'Primary title' ? 1999 10.1/primary\n"
      "loop_ _citation_author_citation_id _citation_author_name\n"
      "other 'Other, A.' PRIMARY 'Upper, A.' primary 'Prim, A.'"
      " primary 'Prim, B.'\n";
    show( "primary citation fallback", "data_c\n" + cell + site + cit
          + "_journal_name_full 'J. Publ.'\n" );
    show( "no primary citation", "data_c\n" + cell + site
          + "loop_ _citation_id _citation_title a 'T'\n" );
    show( "no provenance", "data_n\n" + cell + site );
    show( "several values", "data_c\n" + cell + site
          + "loop_ _journal_year 2000 2001\n_journal_name_full 'J. X'\n" );
    auto doc = parse( "data_p _journal_year 2000\ndata_q _journal_year ?\n" );
    REQUIRE( !NCC::getStructureItems( doc.block( 0 ), Permissive::NO,
                                      warnfct ).has_value() );
    REQUIRE( !NCC::getProvenanceItems( doc.block( 1 ),
                                       warnfct ).has_value() );
    std::cout << "provenance only:\n";
    dump( NCC::getProvenanceItems( doc.block( 0 ), warnfct ).value(), doc );
  }

  void testOtherBlock()
  {
    std::cout << "==== Provenance from another block ====\n";
    //Synthetic files with the layout of older journal supplementary files,
    //where the publication is only in a data_global block:
    const auto one = readTestData( "cif_globalblock.cif" );
    show( "cif_globalblock.cif", one, Permissive::NO, "cif_globalblock.cif" );
    ensureItemsError( "cif_globalblock_two.cif",
                      readTestData( "cif_globalblock_two.cif" ),
                      Permissive::NO, "cif_globalblock_two.cif" );
    const std::string s = "data_s\n" + cell + site
      + "_audit_creation_method 'by hand'\n";
    show( "structure block items first", s + "data_g\n"
          "_audit_creation_method 'other'\n_journal_year 2000\n" );
    show( "no structure block provenance", "data_s\n" + cell + site
          + "data_g\n_journal_year 2000\n" );
    show( "two other blocks (not used)", s + "data_g\n_journal_year 2000\n"
          "data_h\n_journal_name_full 'J. H'\n" );
    show( "other block without usable items is ignored", s
          + "data_g\n_journal_year 2000\ndata_h\n_journal_name_full ?\n" );
    show( "other block items all present already", s
          + "data_g\n_audit_creation_method 'other'\n" );
  }

  void testBlockAccess()
  {
    std::cout << "==== Block access ====\n";
    auto doc = parse( "data_other _x 1\ndata_s\n" + cell + site );
    auto r = NCC::findUniqueStructProvItems( doc, Permissive::NO, warnfct );
    REQUIRE( r.structureItems.block( doc ).name() == "s" );
    REQUIRE( !r.provenanceItems.has_value() );
    NCC::Document moved( std::move( doc ) );
    REQUIRE( r.structureItems.block( moved ).index() == 1 );
    auto other = parse( "data_s\n" + cell + site );
    nctest::ensureError<NC::Error::LogicError>( [&r,&other]()
    {
      r.structureItems.block( other );
    } );
  }
}

int main()
{
  testStructure();
  testAniso();
  testErrors();
  testProvenance();
  testOtherBlock();
  testBlockAccess();
  return 0;
}
