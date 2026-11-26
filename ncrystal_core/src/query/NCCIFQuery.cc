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

#include "NCCIFQuery.hh"
#include "NCrystal/internal/cif/NCCIFUtils.hh"
#include "NCrystal/internal/cif/NCCIFItems.hh"
#include "NCrystal/internal/cif/NCCIFSetting.hh"
#include "NCrystal/internal/cif/NCCIFSites.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    namespace NCC = CIF;

    //A value as [text,line], or null for "?" and ".":
    void streamValue( std::ostream& os, const NCC::Value& v )
    {
      if ( NCC::isNull( v ) ) {
        os << "null";
        return;
      }
      os << '[';
      streamJSON( os, v.text );
      os << ',' << v.line << ']';
    }

    //Values (all, or the rows not in skip, which is sorted) as a list:
    template<class TValues>
    void streamValues( std::ostream& os, const TValues& v,
                       const std::vector<std::uint32_t>& skip = {} )
    {
      os << '[';
      bool first = true;
      auto itSkip = skip.begin();
      for ( auto i : ncrange( v.size() ) ) {
        if ( itSkip != skip.end() && *itSkip == i ) {
          ++itSkip;
          continue;
        }
        if ( !first )
          os << ',';
        first = false;
        streamValue( os, v[i] );
      }
      os << ']';
    }

    //Entry "tag":value of a dict (with a leading comma unless first):
    void streamKey( std::ostream& os, StrView key, bool& first )
    {
      if ( !first )
        os << ',';
      first = false;
      streamJSON( os, key );
      os << ':';
    }

    template<class T>
    void streamAlternatives( std::ostream& os, const T& alts )
    {
      os << '[';
      bool first = true;
      for ( auto& a : alts ) {
        if ( !first )
          os << ',';
        first = false;
        os << "{\"tag\":";
        streamJSON( os, a.tag );
        os << ",\"value\":";
        streamValue( os, a.value );
        os << '}';
      }
      os << ']';
    }

    //Items of a structure block under their CIF tags: single items as
    //"_tag":value (cell angles absent from the data have line 0), the site
    //and aniso columns as "_tag":[values] in "atom_site" and "atom_site_aniso"
    //(with "form" "U" or "B", without skipped rows), and "symmetry" with an
    //entry per tag present (current tag first). Absent items are omitted:
    void streamStructure( std::ostream& os, const NCC::StructureItems& s )
    {
      constexpr const char* cellTags[] = {
        "_cell_length_a", "_cell_length_b", "_cell_length_c",
        "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma" };
      os << '{';
      bool first = true;
      for ( auto i : ncrange( 6 ) ) {
        streamKey( os, cellTags[i], first );
        streamValue( os, s.cell[i] );
      }

      streamKey( os, "atom_site", first );
      {
        const auto& t = s.sites;
        const std::pair<const char*,const NCC::ValueColumn*> cols[] = {
          { "_atom_site_label", &t.label },
          { "_atom_site_fract_x", &t.fractX },
          { "_atom_site_fract_y", &t.fractY },
          { "_atom_site_fract_z", &t.fractZ },
          { "_atom_site_type_symbol", &t.typeSymbol },
          { "_atom_site_U_iso_or_equiv", &t.uIsoOrEquiv },
          { "_atom_site_B_iso_or_equiv", &t.bIsoOrEquiv },
          { "_atom_site_occupancy", &t.occupancy } };
        os << '{';
        bool f2 = true;
        for ( auto& c : cols ) {
          if ( c.second->empty() )
            continue;
          streamKey( os, c.first, f2 );
          streamValues( os, *c.second );
        }
        os << '}';
      }

      if ( s.aniso.has_value() ) {
        const auto& a = s.aniso.value();
        const bool isU = ( a.form == NCC::StructureItems::AnisoForm::U );
        constexpr const char* ij[] = { "11", "22", "33", "12", "13", "23" };
        streamKey( os, "atom_site_aniso", first );
        os << "{\"form\":" << ( isU ? "\"U\"" : "\"B\"" )
           << ",\"_atom_site_aniso_label\":";
        streamValues( os, a.label, a.skippedRows );
        for ( auto i : ncrange( 6 ) ) {
          os << ",\"_atom_site_aniso_" << ( isU ? 'U' : 'B' ) << '_' << ij[i]
             << "\":";
          streamValues( os, a.components[i], a.skippedRows );
        }
        os << '}';
      }

      streamKey( os, "symmetry", first );
      os << "{\"operations\":[";
      {
        bool f2 = true;
        for ( auto& o : s.operations ) {
          if ( !f2 )
            os << ',';
          f2 = false;
          os << "{\"tag\":";
          streamJSON( os, o.tag );
          os << ",\"values\":";
          streamValues( os, o.values );
          os << '}';
        }
      }
      os << "],\"hm\":";
      streamAlternatives( os, s.hmSymbol );
      os << ",\"hall\":";
      streamAlternatives( os, s.hallSymbol );
      os << ",\"number\":";
      streamAlternatives( os, s.itNumber );
      os << '}';

      auto optional = [&os,&first]( const char* key,
                                    const Optional<NCC::Value>& v )
      {
        if ( !v.has_value() )
          return;
        streamKey( os, key, first );
        streamValue( os, v.value() );
      };
      optional( "_chemical_formula_sum", s.formulaSum );
      if ( !s.atomTypeNumberInCell.empty() ) {
        streamKey( os, "_atom_type_number_in_cell", first );
        streamValues( os, s.atomTypeNumberInCell );
      }
      optional( "_diffrn_ambient_temperature", s.ambientTemperature );
      optional( "_cell_measurement_temperature", s.measurementTemperature );
      os << '}';
    }

    void streamDoubles( std::ostream& os, const double* v, std::size_t n )
    {
      os << '[';
      for ( auto i : ncrange( n ) ) {
        if ( i )
          os << ',';
        streamJSON( os, v[i] );
      }
      os << ']';
    }

    //Atoms (cf. CIFAtom) as [{"labels","composition":[[fraction,marker],..],
    //"occupancy","uiso","uij","positions","uijs","split"},..], with null for
    //absent displacements, and "split" null or {"positions","distance"}:
    void streamAtoms( std::ostream& os, const std::vector<NCC::CIFAtom>& atoms )
    {
      os << '[';
      for ( auto ia : ncrange( atoms.size() ) ) {
        const auto& a = atoms[ia];
        if ( ia )
          os << ',';
        os << "{\"labels\":";
        streamJSON( os, a.labels );
        os << ",\"composition\":[";
        for ( auto i : ncrange( a.composition.size() ) ) {
          os << ( i ? ",[" : "[" );
          streamJSON( os, a.composition[i].first );
          os << ',';
          streamJSON( os, a.composition[i].second );
          os << ']';
        }
        os << "],\"occupancy\":";
        streamJSON( os, a.occupancy );
        os << ",\"uiso\":";
        if ( a.uiso.has_value() )
          streamJSON( os, a.uiso.value() );
        else
          os << "null";
        os << ",\"uij\":";
        if ( a.uij.has_value() )
          streamDoubles( os, a.uij.value().data(), 6 );
        else
          os << "null";
        os << ",\"positions\":[";
        for ( auto i : ncrange( a.positions.size() ) ) {
          if ( i )
            os << ',';
          const double p[3] = { a.positions[i].x(), a.positions[i].y(),
                                a.positions[i].z() };
          streamDoubles( os, p, 3 );
        }
        os << "],\"uijs\":";
        if ( a.uijs.empty() ) {
          os << "null";
        } else {
          os << '[';
          for ( auto i : ncrange( a.uijs.size() ) ) {
            if ( i )
              os << ',';
            streamDoubles( os, a.uijs[i].data(), 6 );
          }
          os << ']';
        }
        os << ",\"split\":";
        if ( a.splitPositions )
        {
          os << "{\"positions\":" << a.splitPositions << ",\"distance\":";
          streamJSON( os, a.splitDistance );
          os << '}';
        }
        else
          os << "null";
        os << '}';
      }
      os << ']';
    }

    //Provenance items under their CIF tags (as for streamStructure):
    void streamProvenance( std::ostream& os, const NCC::ProvenanceItems& p )
    {
      os << '{';
      bool first = true;
      auto tagged = [&os,&first]( const Optional<NCC::TaggedValue>& v )
      {
        if ( !v.has_value() )
          return;
        streamKey( os, v.value().tag, first );
        streamValue( os, v.value().value );
      };
      auto optional = [&os,&first]( const char* key,
                                    const Optional<NCC::Value>& v )
      {
        if ( !v.has_value() )
          return;
        streamKey( os, key, first );
        streamValue( os, v.value() );
      };
      tagged( p.title );
      if ( p.authors.has_value() ) {
        streamKey( os, p.authors.value().tag, first );
        streamValues( os, p.authors.value().values );
      }
      tagged( p.journal );
      tagged( p.year );
      tagged( p.doi );
      optional( "_cod_database_code", p.codDatabaseCode );
      optional( "_[local]_ncrystal_mp_database_code", p.mpDatabaseCode );
      optional( "_audit_creation_method", p.auditCreationMethod );
      optional( "_audit_creation_date", p.auditCreationDate );
      if ( !p.auditAuthors.empty() ) {
        streamKey( os, "_audit_author_name", first );
        streamValues( os, p.auditAuthors );
      }
      optional( "_cod_original_formula_sum", p.codOriginalFormulaSum );
      os << '}';
    }
  }
}

void NC::CIFQuery::JSONQuery( std::ostream& os, const Query& query )
{
  auto invalid = [&query](const char * reason){
    std::ostringstream ss;
    ss << "Invalid cif query: ";
    streamJSON( ss, query );//trick: use streamJSON for easy format.
    if ( reason )
      ss<< " ("<<reason<<')';
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  constexpr auto sv_cif = StrView::make("cif");
  if ( query.size() < 2 || query.front() != sv_cif ) {
    invalid(nullptr);
    return;
  }
  const auto& key = query.at(1).trimmed();
  auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
  const std::size_t nargs = static_cast<std::size_t>(query.size()-2);
  //Items of the last argument LIST, given as "RAWDATA",LIST from position i,
  //with LIST a separator character followed by the items separated by it
  //(e.g. "@0.1(2)@-3"). Empty if absent:
  constexpr auto sv_rawdata = StrView::make("RAWDATA");
  auto rawList = [&arg,&invalid,nargs,sv_rawdata]( std::size_t i,
                                                    const char * usage )
  {
    if ( nargs == i )
      return StrView::SplitVect();
    if ( nargs != i + 2 || arg(i) != sv_rawdata || arg(i+1).empty() )
      invalid( usage );
    const StrView l = arg(i+1);
    return l.substr(1).split( l.front() );
  };

  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_numbers = StrView::make("numbers");
  constexpr auto sv_elements = StrView::make("elements");
  constexpr auto sv_structure = StrView::make("structure");
  if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"cif\",\"list\"]");
    os << "[\"numbers\", \"elements\", \"structure\"]";
  } else if ( key == sv_numbers ) {
    //["cif","numbers","RAWDATA",LIST]: printed CIF numbers (e.g.
    //"0.0283(4)") as [[value,rounding error],...], with null for strings
    //which are not numbers (cf. parseCIFNumber):
    const auto items = rawList( 0, "correct usage: [\"cif\",\"numbers\","
                                "\"RAWDATA\",LIST]" );
    os << '[';
    for ( auto i : ncrange( items.size() ) ) {
      if ( i )
        os << ',';
      const auto x = parseCIFNumber( items[i] );
      if ( !x.has_value() ) {
        os << "null";
        continue;
      }
      os << '[';
      streamJSON( os, x.value().value );
      os << ',';
      streamJSON( os, x.value().roundingError );
      os << ']';
    }
    os << ']';
  } else if ( key == sv_elements ) {
    //["cif","elements","typesymbol","RAWDATA",LIST] and ["cif","elements",
    //"label","RAWDATA",FORMULA,"RAWDATA",LIST]: elements (or isotopes) of CIF
    //atom type symbols, or of site labels checked against a chemical formula,
    //as [[Z,A],...] (A=0 for elements), with null where not recognised (cf.
    //elementFromCIFTypeSymbol and elementFromCIFLabel):
    const auto usage = "correct usage: [\"cif\",\"elements\",\"typesymbol\","
      "\"RAWDATA\",LIST] or [\"cif\",\"elements\",\"label\",\"RAWDATA\","
      "FORMULA,\"RAWDATA\",LIST]";
    const bool fromlabel = ( nargs >= 3 && arg(0) == "label"
                             && arg(1) == sv_rawdata );
    if ( !fromlabel && !( nargs >= 1 && arg(0) == "typesymbol" ) )
      invalid( usage );
    const auto items = rawList( fromlabel ? 3 : 1, usage );
    os << '[';
    for ( auto i : ncrange( items.size() ) ) {
      if ( i )
        os << ',';
      const auto e = ( fromlabel ? elementFromCIFLabel( items[i], arg(2) )
                       : elementFromCIFTypeSymbol( items[i] ) );
      if ( e.has_value() )
        os << '[' << e.value().Z << ',' << e.value().A << ']';
      else
        os << "null";
    }
    os << ']';
  } else if ( key == sv_structure ) {
    //["cif","structure",PERMISSIVE,OPTION...,SRC]: the items of the only
    //block with a crystal structure (cf. findUniqueStructProvItems), its
    //space group setting (cf. chooseSetting) and atoms (cf. getAtoms), parsed
    //in lenient mode (cf. CIFParseLenient), with PERMISSIVE "strict" or
    //"permissive" (cf. CIFPermissive), OPTIONs KEY=VALUE (each optionally
    //preceded by "RAWDATA", for values with any characters, e.g. Hall symbols
    //like "P 62 2 (0 0 4)"; override_spacegroup=SG, tol_uij=X with default 1e-3, tol_symprec=X in
    //Aa with default 0.01, 0 to not look for a higher symmetry of the atoms,
    //and allow_cell_misfit=1 to accept a cell which
    //does not fit the setting, e.g. for a later refinement of the
    //structure), and SRC
    //either the name of a TextData or "RAWDATA",NAME,"RAWDATA",TEXT (the CIF
    //data and its name, which can be empty). Result: {"block":name,
    //"block_line":n,"structure":{...},"provenance":{...}|null,
    //"provenance_block":name|null,"spacegroup_setting":setting (a higher
    //symmetry of the atoms is used unless override_spacegroup is given, cf.
    //getAtoms),"cell":
    //{"a":..,...,"gamma":..}|null (snapped to the setting, cf. snapCell,
    //null if it does not fit and that is allowed),"atoms":
    //[{"labels":[..],"composition":[[fraction,"Al"],..],"occupancy":x,
    //"uiso":x|null,"uij":[6]|null,"positions":[[x,y,z],..],"uijs":[[6],..]
    //|null,"split":null|{"positions":k,"distance":d}},..],"warnings":[..]},
    //with values as [text,line] (line 0 for an absent cell
    //angle) or null for "?" and ".", and items under their CIF tags: the cell
    //and other single items "_tag":value, the site and aniso columns
    //"_tag":[values] in "atom_site" and "atom_site_aniso" (with "form" "U"
    //or "B", without rows ignored in permissive mode), and "symmetry":
    //{"operations":[{"tag":..,"values":[..]}],"hm":[{"tag":..,"value":..}],
    //"hall":[..],"number":[..]} (current tag first). "provenance_block" is
    //the other block from which absent provenance items were taken, if any,
    //"spacegroup_setting" e.g. "227:2", and "warnings" the list of warnings
    //(which are not emitted). Absent optional items are omitted:
    const auto usage = "correct usage: [\"cif\",\"structure\",PERMISSIVE,"
      "OPTION...,NAME] or [\"cif\",\"structure\",PERMISSIVE,OPTION...,"
      "\"RAWDATA\",NAME,\"RAWDATA\",TEXT], with PERMISSIVE \"strict\" or"
      " \"permissive\" and OPTIONs override_spacegroup=SG, tol_uij=X,"
      " tol_symprec=X or allow_cell_misfit=1";
    if ( nargs < 2 || ( arg(0) != "strict" && arg(0) != "permissive" ) )
      invalid( usage );
    const bool raw = ( nargs >= 5 && arg(nargs-4) == sv_rawdata
                       && arg(nargs-2) == sv_rawdata );
    const std::size_t nopts = nargs - ( raw ? 5 : 2 );
    if ( !raw && arg(nargs-1) == sv_rawdata )
      invalid( usage );
    constexpr auto sv_override = StrView::make("override_spacegroup=");
    constexpr auto sv_toluij = StrView::make("tol_uij=");
    constexpr auto sv_tolsymprec = StrView::make("tol_symprec=");
    constexpr auto sv_misfit = StrView::make("allow_cell_misfit=1");
    StrView overrideSG;
    Optional<double> tolUij, tolSymprec;
    bool allowMisfit = false;
    for ( std::size_t i = 0; i < nopts; ++i ) {
      if ( arg( 1 + i ) == sv_rawdata && i + 1 < nopts )
        ++i;
      const auto o = arg( 1 + i );
      if ( o.startswith( sv_override ) && !overrideSG.has_value()
           && o.size() > sv_override.size() ) {
        overrideSG = o.substr( sv_override.size() );
      } else if ( o == sv_misfit && !allowMisfit ) {
        allowMisfit = true;
      } else if ( o.startswith( sv_toluij ) && !tolUij.has_value() ) {
        tolUij = o.substr( sv_toluij.size() ).toDbl();
        if ( !tolUij.has_value() || !( tolUij.value() > 0.0 )
             || !std::isfinite( tolUij.value() ) )
          invalid( usage );
      } else if ( o.startswith( sv_tolsymprec ) && !tolSymprec.has_value() ) {
        tolSymprec = o.substr( sv_tolsymprec.size() ).toDbl();
        if ( !tolSymprec.has_value() || !( tolSymprec.value() >= 0.0 )
             || !std::isfinite( tolSymprec.value() ) )
          invalid( usage );
      } else {
        invalid( usage );
      }
    }
    const auto src = [&arg,nopts]( std::size_t i ) { return arg(1+nopts+i); };
    const bool permissive = ( arg(0) == "permissive" );
    const auto perm = ( permissive ? NCC::CIFPermissive::YES
                        : NCC::CIFPermissive::NO );
    const std::string rawname = ( raw && !src(1).empty() ? src(1).to_string()
                                  : std::string( "raw CIF data" ) );
    auto td = ( raw
                ? makeSO<const TextData>( RawStrData( src(3).to_string() ),
                                          TextData::DataType{ "cif" },
                                          DataSourceName( rawname ) )
                : FactImpl::createTextData(
                  TextDataPath( src(0).to_string() ) ) );
    const auto doc = NCC::Document::parse( std::move( td ),
                                           NCC::CIFParseLenient::YES );
    std::vector<std::string> warnings;
    const NCC::CIFWarnFct warnfct = [&warnings]( const std::string& w )
    {
      warnings.push_back( w );
    };
    const auto r = NCC::findUniqueStructProvItems( doc, perm, warnfct );
    const auto setting = NCC::chooseSetting( doc, r.structureItems, perm,
                                             overrideSG, warnfct );
    const auto snapped = NCC::snapCell( doc, r.structureItems, setting );
    if ( !snapped.cell.has_value() && !allowMisfit )
      NCRYSTAL_THROW2( BadInput, "Unit cell not compatible with the space"
                       " group of the CIF data: "
                       << snapped.errorMsg.value() );
    Optional<CellParameters> cell = snapped.cell;
    //(no higher symmetry is looked for if the setting was given):
    const auto res = NCC::getAtoms( doc, r.structureItems, setting, perm,
                                    tolUij.has_value() ? tolUij.value() : 1e-3,
                                    ( overrideSG.has_value() ? 0.0
                                      : ( tolSymprec.has_value()
                                          ? tolSymprec.value() : 0.01 ) ),
                                    warnfct );
    const auto& atoms = res.atoms;
    if ( res.setting != setting )
      cell = NCC::snapCell( doc, r.structureItems, res.setting ).cell;
    const auto b = r.structureItems.block( doc );
    os << "{\"block\":";
    streamJSON( os, b.name() );
    os << ",\"block_line\":" << b.line() << ",\"structure\":";
    streamStructure( os, r.structureItems );
    os << ",\"provenance\":";
    Optional<NCC::BlockView> other;
    if ( r.provenanceItems.has_value() ) {
      streamProvenance( os, r.provenanceItems.value() );
      other = r.provenanceItems.value().otherBlock( doc );
    } else {
      os << "null";
    }
    os << ",\"provenance_block\":";
    if ( other.has_value() )
      streamJSON( os, other.value().name() );
    else
      os << "null";
    os << ",\"spacegroup_setting\":";
    streamJSON( os, res.setting.toString() );
    os << ",\"cell\":";
    if ( cell.has_value() ) {
      const auto& c = cell.value();
      streamJSONDictEntry( os, "a", c.a, JSONDictPos::FIRST );
      streamJSONDictEntry( os, "b", c.b );
      streamJSONDictEntry( os, "c", c.c );
      streamJSONDictEntry( os, "alpha", c.alpha );
      streamJSONDictEntry( os, "beta", c.beta );
      streamJSONDictEntry( os, "gamma", c.gamma, JSONDictPos::LAST );
    } else {
      os << "null";
    }
    os << ",\"atoms\":";
    streamAtoms( os, atoms );
    os << ",\"warnings\":";
    streamJSON( os, warnings );
    os << '}';
  } else {
    invalid("unknown key");
  }
}
