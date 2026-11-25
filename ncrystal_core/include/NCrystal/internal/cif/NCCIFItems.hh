#ifndef NCrystal_CIFItems_hh
#define NCrystal_CIFItems_hh

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
#include "NCrystal/core/NCSmallVector.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace CIF {

    //The data items of a block about its crystal structure and its
    //provenance, as stated in the data: validated for presence and
    //consistency of the tags as described below, but not interpreted (values
    //stay text, e.g. numbers are not parsed). Values "?" and "." of optional
    //single items are treated as absent, and so are provenance items with
    //several values (with a warning). All values refer into the Document,
    //which must outlive the items (moving it is fine).

    //Warnings are passed to a function supplied by the caller (never to
    //NCRYSTAL_WARN), so callers can record them (e.g. for NCMAT comments):
    using CIFWarnFct = std::function<void(const std::string&)>;

    //Whether getStructureItems accepts some flaws in the data (with
    //warnings), as described below. Errors which permissive mode would avoid
    //end with a hint to use it:
    enum class CIFPermissive { YES, NO };

    //A value with the tag it was found under (for items with alternative
    //tags):
    struct ItemsFactory;//(internal)

    struct TaggedValue {
      StrView tag;
      Value value;
    };
    struct TaggedColumn {
      StrView tag;
      ValueColumn values;
    };
    struct TaggedValues {
      StrView tag;
      std::vector<Value> values;
    };

    class StructureItems final : private MoveOnly {
    public:
      //The block the items were extracted from (pass the Document they were
      //extracted from):
      BlockView block( const Document& ) const;

      //Cell lengths (_cell_length_a/b/c) and angles (_cell_angle_alpha/beta/
      //gamma), required and not "?" or ".". An absent angle is
      //s_defaultAngle:
      std::array<Value,6> cell;
      //Dictionary default for an absent cell angle:
      static constexpr Value s_defaultAngle{ StrView::make("90"),
                                             ValueKind::Unquoted, 0 };

      //Sites, the _atom_site_ columns (all of the same size, in one loop
      //unless there is a single site). Label and fractional coordinates are
      //required and may not be "?" or ".". Absent optional columns are empty:
      struct Sites {
        ValueColumn label, fractX, fractY, fractZ;
        ValueColumn typeSymbol, uIsoOrEquiv, bIsoOrEquiv, occupancy;
      };
      Sites sites;

      //Anisotropic displacements, the _atom_site_aniso_ columns as U or B
      //(B = 8 pi^2 U) tensor components 11, 22, 33, 12, 13 and 23, with the
      //label column (_atom_site_aniso_label, or _atom_site_label if the
      //components are in the loop of the sites). All six components of the
      //form, and the label, are required (if missing: an error, unless in
      //permissive mode, where all anisotropic displacements are ignored with
      //a warning). Rows with a label which is not one of a site, or with
      //components "?" or ".", are an error, unless in permissive mode, where
      //they are ignored (with a warning) and listed in skippedRows.
      //Components given as beta, or as both U and B, are an error:
      enum class AnisoForm { U, B };
      struct Aniso {
        AnisoForm form;
        ValueColumn label;
        std::array<ValueColumn,6> components;
        std::vector<std::uint32_t> skippedRows;
      };
      Optional<Aniso> aniso;

      //Symmetry: the operations (_space_group_symop_operation_xyz or the
      //deprecated _symmetry_equiv_pos_as_xyz), the Hermann-Mauguin symbol
      //(_space_group_name_H-M_alt or _symmetry_space_group_name_H-M), the
      //Hall symbol (_space_group_name_Hall or _symmetry_space_group_name_Hall)
      //and the International Tables number (_space_group_IT_number or
      //_symmetry_Int_Tables_number). Each has an entry per tag present,
      //current tag first, but the deprecated one only if its text differs
      //(whether the two are consistent is not checked here):
      SmallVector<TaggedColumn,2> operations;
      SmallVector<TaggedValue,2> hmSymbol, hallSymbol, itNumber;

      Optional<Value> formulaSum;//_chemical_formula_sum
      ValueColumn atomTypeNumberInCell;//_atom_type_number_in_cell
      Optional<Value> ambientTemperature;//_diffrn_ambient_temperature
      Optional<Value> measurementTemperature;//_cell_measurement_temperature

      StructureItems( StructureItems&& ) = default;
      StructureItems& operator=( StructureItems&& ) = default;
    private:
      friend struct ItemsFactory;
      StructureItems( BlockView );
      UniqueIDValue m_docID;
      std::uint32_t m_block;
    };

    class ProvenanceItems final {
    public:
      //The block the items were extracted from (pass the Document they were
      //extracted from), and the other block from which absent items were
      //taken, if any (cf. findUniqueStructProvItems):
      BlockView block( const Document& ) const;
      Optional<BlockView> otherBlock( const Document& ) const;

      //The publication of the data (_publ_section_title, _publ_author_name,
      //_journal_name_full, _journal_year, _journal_paper_doi). Each absent
      //one is taken from the primary citation if any, i.e. the
      //_citation_ row with _citation_id "primary" (_citation_title,
      //_citation_author_name of the rows with _citation_author_citation_id
      //"primary", _citation_journal_full, _citation_year, _citation_doi):
      Optional<TaggedValue> title;
      Optional<TaggedValues> authors;
      Optional<TaggedValue> journal, year, doi;

      Optional<Value> codDatabaseCode;//_cod_database_code
      //_[local]_ncrystal_mp_database_code (cf. Document::parse):
      Optional<Value> mpDatabaseCode;
      Optional<Value> auditCreationMethod;//_audit_creation_method
      Optional<Value> auditCreationDate;//_audit_creation_date
      std::vector<Value> auditAuthors;//_audit_author_name
      Optional<Value> codOriginalFormulaSum;//_cod_original_formula_sum

    private:
      friend struct ItemsFactory;
      ProvenanceItems( BlockView );
      UniqueIDValue m_docID;
      std::uint32_t m_block;
      Optional<std::uint32_t> m_otherBlock;
    };

    //The structure items of the block, or NullOpt if it has no crystal
    //structure (no tag starting with _atom_site_ or _atom_sites_). Throws
    //BadInput, with the data source name, block and line, if the items are
    //not as described above:
    Optional<StructureItems> getStructureItems( BlockView, CIFPermissive,
                                                const CIFWarnFct& );

    //The provenance items of the block, or NullOpt if it has none (i.e. no
    //usable value of any of the provenance tags above):
    Optional<ProvenanceItems> getProvenanceItems( BlockView,
                                                  const CIFWarnFct& );

    //The structure items of the only block of the document with a crystal
    //structure (all blocks with one are validated), and its provenance items.
    //Items absent there are taken from the only block without a structure
    //which has provenance items, if there is exactly one such block (e.g. a
    //data_global block with the publication). Throws BadInput if there is
    //not exactly one block with a crystal structure:
    struct StructProvItems {
      StructureItems structureItems;
      Optional<ProvenanceItems> provenanceItems;
    };
    StructProvItems findUniqueStructProvItems( const Document&,
                                               CIFPermissive,
                                               const CIFWarnFct& );
  }
}

#endif
