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

#include "NCrystal/internal/cif/NCCIFDescribe.hh"
#include <iostream>

//Tests of CIF::getDescription: description lines, COD and MP ids and formula
//of CIF data, from its provenance items.

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace {

  void test( const char * descr, const std::string& items,
             std::uint64_t reqCod = 0, std::uint64_t reqMP = 0 )
  {
    std::vector<std::string> warnings;
    const NCC::CIFWarnFct warnfct = [&warnings]( const std::string& w )
    {
      warnings.push_back( w );
    };
    const std::string data = ( "data_test\n_cell_length_a 4\n_cell_length_b 4\n"
                               "_cell_length_c 4\n_atom_site_label Al1\n"
                               "_atom_site_fract_x 0\n_atom_site_fract_y 0\n"
                               "_atom_site_fract_z 0\n" ) + items;
    auto doc = NCC::Document::parse(
      NC::makeSO<const NC::TextData>( NC::RawStrData( std::string( data ) ),
                                      NC::TextData::DataType{ "cif" },
                                      NC::DataSourceName( "test.cif" ) ) );
    auto r = NCC::findUniqueStructProvItems( doc, NCC::CIFPermissive::NO,
                                             warnfct );
    const auto d = NCC::getDescription( r.structureItems, r.provenanceItems,
                                      reqCod, reqMP, warnfct );
    std::cout << descr << ": cod_id " << d.codID << ", mp_id " << d.mpID
              << ", formula \"" << d.formula << "\"\n";
    for ( auto& l : d.lines )
      std::cout << "   | " << l << '\n';
    for ( auto& w : warnings )
      std::cout << "   warning: " << w << '\n';
  }
}

int main()
{
  test( "nothing", "" );
  test( "one author, no journal",
        "_publ_section_title 'A title'\n_publ_author_name 'Doe, J.'\n" );
  test( "two authors, journal and year",
        "loop_ _publ_author_name 'Doe, J.' 'Roe, R.'\n"
        "_journal_name_full 'J. Test'\n_journal_year 2001\n" );
  test( "three authors, journal only",
        "loop_ _publ_author_name A B C\n_journal_name_full 'J. Test'\n" );
  test( "several authors without journal",
        "loop_ _publ_author_name A B\n" );
  test( "journal and year without authors",
        "_journal_name_full 'J. Test'\n_journal_year 2001\n" );
  test( "multi-line title (trailing spaces removed)",
        "_publ_section_title\n;\nFirst line   \nSecond line\n;\n" );
  test( "DOI and audit items",
        "_journal_paper_doi 10.1000/xyz\n_audit_creation_method 'by hand'\n"
        "_audit_creation_date 2001-01-01\n"
        "loop_ _audit_author_name 'A. B.' 'C. D.'\n" );
  test( "DOI of the primary citation",
        "loop_ _citation_id _citation_doi _citation_title primary"
        " 10.1000/abc 'Cited title'\n" );
  test( "COD and MP ids embedded",
        "_cod_database_code 1234567\n"
        "_[local]_ncrystal_mp_database_code mp-55\n" );
  test( "requested ids (the same, and others)",
        "_cod_database_code 1234567\n", 1234567, 66 );
  test( "requested COD id differing from the embedded one",
        "_cod_database_code 1234567\n", 7654321 );
  test( "requested ids without provenance items", "", 1000001, 0 );
  test( "formula, and the COD original formula as fallback",
        "_cod_original_formula_sum 'Al2 O3'\n" );
  test( "formula given", "_chemical_formula_sum 'Al'\n"
        "_cod_original_formula_sum 'Al2 O3'\n" );
  return 0;
}
