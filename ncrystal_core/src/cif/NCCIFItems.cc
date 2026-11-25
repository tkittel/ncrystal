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

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    namespace {

      constexpr auto s_atomSiteLabel = StrView::make("_atom_site_label");

      //Location for messages, e.g. 'CIF data in "x.cif" (data_a, line 7)'
      //(line 0 for none):
      std::string where( BlockView b, std::uint32_t line )
      {
        return "CIF data in " + describeLocation( b, line );
      }

      [[noreturn]] void throwItemError( BlockView b, std::uint32_t line,
                                        const std::string& reason )
      {
        NCRYSTAL_THROW2( BadInput, "Invalid " << where( b, line ) << ": "
                         << reason );
      }

      //The value of a single item, NullOpt if absent or "?" or ".". Several
      //values (a loop) are an error, or ignored with a warning if a warning
      //function is given:
      Optional<Value> singleImpl( BlockView b, StrView tag,
                                  const CIFWarnFct* warnfct )
      {
        const auto v = b.values( tag );
        if ( v.empty() )
          return NullOpt;
        if ( v.size() != 1 ) {
          if ( !warnfct )
            throwItemError( b, v[0].line, "tag " + tag.to_string()
                            + " must have a single value" );
          (*warnfct)( "Ignoring " + tag.to_string() + " in "
                      + where( b, v[0].line ) + ", which has several values" );
          return NullOpt;
        }
        if ( isNull( v[0] ) )
          return NullOpt;
        return v[0];
      }

      Optional<Value> single( BlockView b, StrView tag )
      {
        return singleImpl( b, tag, nullptr );
      }

      Optional<Value> singleOrWarn( BlockView b, StrView tag,
                                    const CIFWarnFct& warnfct )
      {
        return singleImpl( b, tag, &warnfct );
      }

      //Columns: same size and in the same loop (or all single values) as
      //ref. Empty columns (absent tags) are ignored:
      void checkSameLoop( BlockView b, StrView reftag, const ValueColumn& ref,
                          StrView tag, const ValueColumn& col )
      {
        if ( col.empty() )
          return;
        bool ok = ( col.size() == ref.size() && col.inLoop() == ref.inLoop() );
        if ( ok && ref.inLoop() )
          ok = b.loopOf( reftag ).value().column( tag ).has_value();
        if ( !ok )
          throwItemError( b, col[0].line, "tag " + tag.to_string()
                          + " must be in the same loop as "
                          + reftag.to_string() );
      }

      //Entries for a current and a deprecated tag (the latter only if its
      //text differs):
      void addAlternatives( BlockView b, StrView cur, StrView dep,
                            SmallVector<TaggedValue,2>& out )
      {
        for ( auto tag : { cur, dep } ) {
          auto v = single( b, tag );
          if ( v.has_value()
               && ( out.empty() || out[0].value.text != v.value().text ) )
            out.push_back( TaggedValue{ tag, v.value() } );
        }
      }
      void addAlternatives( BlockView b, StrView cur, StrView dep,
                            SmallVector<TaggedColumn,2>& out )
      {
        for ( auto tag : { cur, dep } ) {
          auto v = b.values( tag );
          if ( v.empty() )
            continue;
          auto sameAsFirst = [&out,&v]()
          {
            if ( out.empty() || out[0].values.size() != v.size() )
              return false;
            for ( auto i : ncrange( v.size() ) )
              if ( out[0].values[i].text != v[i].text )
                return false;
            return true;
          };
          if ( !sameAsFirst() )
            out.push_back( TaggedColumn{ tag, v } );
        }
      }

      constexpr std::array<StrView,6> s_cellTags = {
        StrView::make("_cell_length_a"), StrView::make("_cell_length_b"),
        StrView::make("_cell_length_c"), StrView::make("_cell_angle_alpha"),
        StrView::make("_cell_angle_beta"), StrView::make("_cell_angle_gamma")
      };

      constexpr std::array<StrView,6> s_anisoU = {
        StrView::make("_atom_site_aniso_U_11"),
        StrView::make("_atom_site_aniso_U_22"),
        StrView::make("_atom_site_aniso_U_33"),
        StrView::make("_atom_site_aniso_U_12"),
        StrView::make("_atom_site_aniso_U_13"),
        StrView::make("_atom_site_aniso_U_23")
      };
      constexpr std::array<StrView,6> s_anisoB = {
        StrView::make("_atom_site_aniso_B_11"),
        StrView::make("_atom_site_aniso_B_22"),
        StrView::make("_atom_site_aniso_B_33"),
        StrView::make("_atom_site_aniso_B_12"),
        StrView::make("_atom_site_aniso_B_13"),
        StrView::make("_atom_site_aniso_B_23")
      };

      //Anisotropic displacements of block b (U or B form, not both; beta is
      //an error), checked against the site labels. NullOpt if absent, or if
      //permissive and tags are missing. Rows with an unknown label or "?" or
      //"." components are errors, or skipped (with a warning) if permissive:
      Optional<StructureItems::Aniso> getAniso( BlockView b,
                                                const ValueColumn& siteLabels,
                                                CIFPermissive permissive,
                                                const CIFWarnFct& warnfct )
      {
        using Aniso = StructureItems::Aniso;
        for ( auto i : { 11, 22, 33, 12, 13, 23 } ) {
          const std::string tag = "_atom_site_aniso_beta_"
            + std::to_string( i );
          if ( b.has( tag ) )
            throwItemError( b, b.values( tag ).empty() ? b.line()
                            : b.values( tag )[0].line, "anisotropic"
                            " displacements given as beta (tag " + tag
                            + ") are not supported" );
        }
        auto anyOf = [&b]( const std::array<StrView,6>& tags ) {
          for ( auto& t : tags )
            if ( b.has( t ) )
              return true;
          return false;
        };
        const bool hasU = anyOf( s_anisoU );
        const bool hasB = anyOf( s_anisoB );
        if ( !hasU && !hasB )
          return NullOpt;
        if ( hasU && hasB )
          throwItemError( b, b.values( s_anisoB[0] ).empty() ? b.line()
                          : b.values( s_anisoB[0] )[0].line, "anisotropic"
                          " displacements given both as U and as B" );
        const auto& tags = ( hasU ? s_anisoU : s_anisoB );
        //Missing tags: an error, or ignore all aniso data if permissive:
        auto missing = [&b,permissive,&warnfct]( StrView tag ) -> bool
        {
          const std::string what = ( "anisotropic displacements lack the tag "
                                     + tag.to_string() );
          if ( permissive == CIFPermissive::NO )
            throwItemError( b, 0, what + ". Use permissive mode to ignore the"
                            " anisotropic displacements." );
          warnfct( "Ignoring the anisotropic displacements in "
                   + where( b, 0 ) + ", which lack the tag "
                   + tag.to_string() );
          return true;
        };
        Aniso res;
        res.form = ( hasU ? StructureItems::AnisoForm::U
                     : StructureItems::AnisoForm::B );
        for ( auto i : ncrange( 6 ) )
          if ( ( res.components[i] = b.values( tags[i] ) ).empty()
               && missing( tags[i] ) )
            return NullOpt;
        for ( auto i : ncrange( 6 ) )
          checkSameLoop( b, tags[0], res.components[0], tags[i],
                         res.components[i] );
        constexpr auto anisoLabel = StrView::make("_atom_site_aniso_label");
        StrView labelTag = anisoLabel;
        res.label = b.values( anisoLabel );
        if ( res.label.empty() ) {
          const auto& c0 = res.components[0];
          if ( c0.inLoop()
               && b.loopOf( tags[0] ).value().column( s_atomSiteLabel )
               .has_value() ) {
            res.label = siteLabels;
            labelTag = s_atomSiteLabel;
          } else if ( missing( anisoLabel ) ) {
            return NullOpt;
          }
        }
        checkSameLoop( b, tags[0], res.components[0], labelTag, res.label );
        std::vector<StrView> known;
        known.reserve( siteLabels.size() );
        for ( auto i : ncrange( siteLabels.size() ) )
          known.push_back( siteLabels[i].text );
        std::sort( known.begin(), known.end() );
        for ( auto r : ncrange( res.label.size() ) ) {
          const auto& lbl = res.label[r];
          if ( !std::binary_search( known.begin(), known.end(), lbl.text ) ) {
            if ( permissive == CIFPermissive::NO )
              throwItemError( b, lbl.line, "the label \""
                              + lbl.text.to_string() + "\" of anisotropic"
                              " displacements is not the label of a site. Use"
                              " permissive mode to ignore such entries." );
            warnfct( "Ignoring anisotropic displacements of \""
                     + lbl.text.to_string() + "\" in " + where( b, lbl.line )
                     + ", which is not the label of a site" );
            res.skippedRows.push_back( static_cast<std::uint32_t>( r ) );
            continue;
          }
          bool anyNull = false;
          for ( auto& c : res.components )
            anyNull = anyNull || isNull( c[r] );
          if ( !anyNull )
            continue;
          if ( permissive == CIFPermissive::NO )
            throwItemError( b, lbl.line, "anisotropic displacements of \""
                            + lbl.text.to_string() + "\" have components"
                            " \"?\" or \".\". Use permissive mode to ignore"
                            " such entries." );
          warnfct( "Ignoring anisotropic displacements of \""
                   + lbl.text.to_string() + "\" in " + where( b, lbl.line )
                   + ", which have components \"?\" or \".\"" );
          res.skippedRows.push_back( static_cast<std::uint32_t>( r ) );
        }
        return res;
      }
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    struct ItemsFactory {
      static StructureItems structure( BlockView b ) { return { b }; }
      static ProvenanceItems provenance( BlockView b ) { return { b }; }

      static bool empty( const ProvenanceItems& p )
      {
        return !( p.title.has_value() || p.authors.has_value()
                  || p.journal.has_value() || p.year.has_value()
                  || p.doi.has_value() || p.codDatabaseCode.has_value()
                  || p.mpDatabaseCode.has_value()
                  || p.auditCreationMethod.has_value()
                  || p.auditCreationDate.has_value()
                  || !p.auditAuthors.empty()
                  || p.codOriginalFormulaSum.has_value() );
      }

      //Fill the items absent in dst from src (another block), recording
      //src's block in dst if any was used:
      static void fill( ProvenanceItems& dst, const ProvenanceItems& src )
      {
        bool used = false;
        auto f = [&used]( auto& d, const auto& s )
        {
          if ( !d.has_value() && s.has_value() ) {
            d = s;
            used = true;
          }
        };
        f( dst.title, src.title );
        f( dst.authors, src.authors );
        f( dst.journal, src.journal );
        f( dst.year, src.year );
        f( dst.doi, src.doi );
        f( dst.codDatabaseCode, src.codDatabaseCode );
        f( dst.mpDatabaseCode, src.mpDatabaseCode );
        f( dst.auditCreationMethod, src.auditCreationMethod );
        f( dst.auditCreationDate, src.auditCreationDate );
        f( dst.codOriginalFormulaSum, src.codOriginalFormulaSum );
        if ( dst.auditAuthors.empty() && !src.auditAuthors.empty() ) {
          dst.auditAuthors = src.auditAuthors;
          used = true;
        }
        if ( used )
          dst.m_otherBlock = src.m_block;
      }
    };
  }
}

NCC::StructureItems::StructureItems( BlockView b )
  : m_docID( b.document().getUniqueID() ),
    m_block( static_cast<std::uint32_t>( b.index() ) )
{
}

NCC::ProvenanceItems::ProvenanceItems( BlockView b )
  : m_docID( b.document().getUniqueID() ),
    m_block( static_cast<std::uint32_t>( b.index() ) )
{
}

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    namespace {
      BlockView blockFromID( const Document& doc, UniqueIDValue id,
                             std::uint32_t iblock )
      {
        if ( doc.getUniqueID() != id )
          NCRYSTAL_THROW( LogicError, "CIF items used with a Document other"
                          " than the one they were extracted from" );
        return doc.block( iblock );
      }
    }
  }
}

NCC::BlockView NCC::StructureItems::block( const Document& doc ) const
{
  return blockFromID( doc, m_docID, m_block );
}

NCC::BlockView NCC::ProvenanceItems::block( const Document& doc ) const
{
  return blockFromID( doc, m_docID, m_block );
}

NC::Optional<NCC::BlockView>
NCC::ProvenanceItems::otherBlock( const Document& doc ) const
{
  if ( !m_otherBlock.has_value() )
    return NullOpt;
  return blockFromID( doc, m_docID, m_otherBlock.value() );
}

NC::Optional<NCC::StructureItems>
NCC::getStructureItems( BlockView b, CIFPermissive permissive,
                        const CIFWarnFct& warnfct )
{
#ifndef NDEBUG
  nc_assert_always( bool(warnfct) );
#endif
  bool hasStructure = false;
  for ( auto i : ncrange( b.nTags() ) ) {
    const auto t = b.tag( i );
    if ( t.startswithUncased( "_atom_site_" )
         || t.startswithUncased( "_atom_sites_" ) ) {
      hasStructure = true;
      break;
    }
  }
  if ( !hasStructure )
    return NullOpt;
  auto res = ItemsFactory::structure( b );

  for ( auto i : ncrange( 6 ) ) {
    const auto tag = s_cellTags[i];
    const auto v = b.values( tag );
    if ( v.empty() ) {
      if ( i < 3 )
        throwItemError( b, 0, "missing tag " + tag.to_string() );
      res.cell[i] = StructureItems::s_defaultAngle;
      continue;
    }
    if ( v.size() != 1 || isNull( v[0] ) )
      throwItemError( b, v[0].line, "tag " + tag.to_string()
                      + " must have a single value which is not \"?\" or"
                      " \".\"" );
    res.cell[i] = v[0];
  }

  auto& s = res.sites;
  constexpr auto t_x = StrView::make("_atom_site_fract_x");
  constexpr auto t_y = StrView::make("_atom_site_fract_y");
  constexpr auto t_z = StrView::make("_atom_site_fract_z");
  constexpr auto t_type = StrView::make("_atom_site_type_symbol");
  constexpr auto t_uiso = StrView::make("_atom_site_U_iso_or_equiv");
  constexpr auto t_biso = StrView::make("_atom_site_B_iso_or_equiv");
  constexpr auto t_occ = StrView::make("_atom_site_occupancy");
  struct Col { StrView tag; ValueColumn* col; bool required; };
  const Col cols[] = { { s_atomSiteLabel, &s.label, true },
                       { t_x, &s.fractX, true }, { t_y, &s.fractY, true },
                       { t_z, &s.fractZ, true },
                       { t_type, &s.typeSymbol, false },
                       { t_uiso, &s.uIsoOrEquiv, false },
                       { t_biso, &s.bIsoOrEquiv, false },
                       { t_occ, &s.occupancy, false } };
  for ( auto& c : cols ) {
    *c.col = b.values( c.tag );
    if ( c.required && c.col->empty() )
      throwItemError( b, 0, "missing tag " + c.tag.to_string() );
    checkSameLoop( b, s_atomSiteLabel, s.label, c.tag, *c.col );
    if ( c.required )
      for ( auto i : ncrange( c.col->size() ) )
        if ( isNull( (*c.col)[i] ) )
          throwItemError( b, (*c.col)[i].line, "tag " + c.tag.to_string()
                          + " has a value \"?\" or \".\"" );
  }

  res.aniso = getAniso( b, s.label, permissive, warnfct );

  addAlternatives( b, "_space_group_symop_operation_xyz",
                   "_symmetry_equiv_pos_as_xyz", res.operations );
  addAlternatives( b, "_space_group_name_H-M_alt",
                   "_symmetry_space_group_name_H-M", res.hmSymbol );
  addAlternatives( b, "_space_group_name_Hall",
                   "_symmetry_space_group_name_Hall", res.hallSymbol );
  addAlternatives( b, "_space_group_IT_number",
                   "_symmetry_Int_Tables_number", res.itNumber );

  res.formulaSum = single( b, "_chemical_formula_sum" );
  res.atomTypeNumberInCell = b.values( "_atom_type_number_in_cell" );
  res.ambientTemperature = single( b, "_diffrn_ambient_temperature" );
  res.measurementTemperature = single( b, "_cell_measurement_temperature" );
  return Optional<StructureItems>( std::move( res ) );
}

NC::Optional<NCC::ProvenanceItems>
NCC::getProvenanceItems( BlockView b, const CIFWarnFct& warnfct )
{
#ifndef NDEBUG
  nc_assert_always( bool(warnfct) );
#endif
  auto res = ItemsFactory::provenance( b );
  auto single1 = [&b,&warnfct]( StrView tag )
  {
    return singleOrWarn( b, tag, warnfct );
  };

  //Row of the primary citation:
  constexpr auto primary = StrView::make("primary");
  Optional<std::size_t> citRow;
  {
    const auto ids = b.values( "_citation_id" );
    for ( auto i : ncrange( ids.size() ) )
      if ( !isNull( ids[i] ) && ids[i].text == primary ) {
        citRow = i;
        break;
      }
  }
  auto citation = [&b,&citRow]( StrView tag ) -> Optional<TaggedValue>
  {
    if ( !citRow.has_value() )
      return NullOpt;
    const auto v = b.values( tag );
    if ( v.size() != b.values( "_citation_id" ).size()
         || isNull( v[citRow.value()] ) )
      return NullOpt;
    return TaggedValue{ tag, v[citRow.value()] };
  };
  auto tagged = [&single1,&citation]( StrView tag, StrView citTag )
    -> Optional<TaggedValue>
  {
    auto v = single1( tag );
    if ( v.has_value() )
      return TaggedValue{ tag, v.value() };
    return citation( citTag );
  };
  auto nonNull = [&b]( StrView tag ) {
    std::vector<Value> out;
    const auto v = b.values( tag );
    for ( auto i : ncrange( v.size() ) )
      if ( !isNull( v[i] ) )
        out.push_back( v[i] );
    return out;
  };

  res.title = tagged( "_publ_section_title", "_citation_title" );
  {
    constexpr auto t_publ = StrView::make("_publ_author_name");
    constexpr auto t_cit = StrView::make("_citation_author_name");
    auto a = nonNull( t_publ );
    if ( !a.empty() ) {
      res.authors = TaggedValues{ t_publ, std::move( a ) };
    } else if ( citRow.has_value() ) {
      const auto ids = b.values( "_citation_author_citation_id" );
      const auto names = b.values( t_cit );
      std::vector<Value> v;
      if ( ids.size() == names.size() )
        for ( auto i : ncrange( ids.size() ) )
          if ( !isNull( ids[i] ) && !isNull( names[i] )
               && ids[i].text == primary )
            v.push_back( names[i] );
      if ( !v.empty() )
        res.authors = TaggedValues{ t_cit, std::move( v ) };
    }
  }
  res.journal = tagged( "_journal_name_full", "_citation_journal_full" );
  res.year = tagged( "_journal_year", "_citation_year" );
  res.doi = tagged( "_journal_paper_doi", "_citation_doi" );

  res.codDatabaseCode = single1( "_cod_database_code" );
  res.mpDatabaseCode = single1( "_[local]_ncrystal_mp_database_code" );
  res.auditCreationMethod = single1( "_audit_creation_method" );
  res.auditCreationDate = single1( "_audit_creation_date" );
  res.auditAuthors = nonNull( "_audit_author_name" );
  res.codOriginalFormulaSum = single1( "_cod_original_formula_sum" );
  if ( ItemsFactory::empty( res ) )
    return NullOpt;
  return res;
}

NCC::StructProvItems
NCC::findUniqueStructProvItems( const Document& doc, CIFPermissive permissive,
                                const CIFWarnFct& warnfct )
{
#ifndef NDEBUG
  nc_assert_always( bool(warnfct) );
#endif
  Optional<StructureItems> found;
  std::string names;
  std::size_t n = 0;
  std::vector<std::uint32_t> others;//blocks without a structure
  for ( auto i : ncrange( doc.nBlocks() ) ) {
    const auto b = doc.block( i );
    auto s = getStructureItems( b, permissive, warnfct );
    if ( !s.has_value() ) {
      others.push_back( static_cast<std::uint32_t>( i ) );
      continue;
    }
    names += ( n++ ? ", data_" : "data_" ) + b.name().to_string();
    if ( !found.has_value() )
      found = std::move( s );
  }
  if ( n != 1 ) {
    const std::string why
      = ( n ? "several data blocks have a crystal structure (" + names + ")"
          : "no data block has a crystal structure (tags _atom_site_...)" );
    NCRYSTAL_THROW2( BadInput, "Invalid CIF data in \""
                     << doc.dataSourceName() << "\": " << why );
  }
  const auto b = found.value().block( doc );
  auto prov = getProvenanceItems( b, warnfct );
  //Missing items from the only block without a structure with provenance:
  Optional<ProvenanceItems> other;
  for ( auto i : others ) {
    auto p = getProvenanceItems( doc.block( i ), warnfct );
    if ( !p.has_value() )
      continue;
    if ( other.has_value() ) {
      other.reset();
      break;
    }
    other = std::move( p );
  }
  if ( other.has_value() ) {
    if ( !prov.has_value() )
      prov = ItemsFactory::provenance( b );
    ItemsFactory::fill( prov.value(), other.value() );
  }
  return StructProvItems{ std::move( found.value() ), std::move( prov ) };
}

std::string NCC::describeLocation( BlockView b, std::uint32_t line )
{
  std::ostringstream ss;
  ss << '"' << b.document().dataSourceName() << "\" (data_" << b.name();
  if ( line )
    ss << ", line " << line;
  ss << ')';
  return ss.str();
}
