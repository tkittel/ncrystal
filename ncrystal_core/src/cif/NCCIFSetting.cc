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
#include "NCrystal/internal/cif/NCCIFUtils.hh"
#include "NCrystal/internal/sgsym/NCSGIdentify.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    namespace {

      //The first value which is not "?" or ".", or a StrView without value:
      StrView firstValue( const SmallVector<TaggedValue,2>& alternatives )
      {
        for ( auto& a : alternatives )
          if ( !isNull( a.value ) )
            return a.value.text;
        return StrView();
      }

      std::string listSettings( const std::vector<SpaceGroup>& v )
      {
        std::string res;
        for ( auto& e : v )
          res += ( res.empty() ? "" : ", " ) + e.toString();
        return res;
      }

      //Whether the cell fits the setting (lenient tolerances, with the
      //rounding errors of the printed values):
      bool cellFits( const SpaceGroup& sg, const StructureItems& s )
      {
        double v[6], e[6];
        for ( auto i : ncrange( 6 ) ) {
          const auto x = parseCIFNumber( s.cell[i].text );
          if ( !x.has_value() )
            return false;
          v[i] = x.value().value;
          e[i] = x.value().roundingError;
        }
        using Tol = SGCellConstraints::Tolerances;
        return SGSymmetry::get( sg ).cellConstraints()
          .trySnap( { v[0], v[1], v[2], v[3], v[4], v[5] }, Tol::lenient(),
                    { e[0], e[1], e[2], e[3], e[4], e[5] } ).cell.has_value();
      }

      //The one of the candidate settings which the cell fits (what describes
      //their source, empty if none):
      SpaceGroup resolve( const std::vector<SpaceGroup>& settings,
                          const StructureItems& s, const std::string& what )
      {
        if ( settings.empty() )
          NCRYSTAL_THROW2( BadInput, "Could not determine space group from CIF"
                           << ( what.empty() ? std::string( " data" )
                                : " data (" + what + " not recognised)" ) );
        if ( settings.size() == 1 )
          return settings.front();
        std::vector<SpaceGroup> fit;
        for ( auto& e : settings )
          if ( cellFits( e, s ) )
            fit.push_back( e );
        if ( fit.size() == 1 )
          return fit.front();
        if ( fit.empty() )
          NCRYSTAL_THROW2( BadInput, "The unit cell does not fit any of the"
                           " settings of the " << what << " ("
                           << listSettings( settings ) << ")" );
        const auto f0 = fit.front().toString();
        NCRYSTAL_THROW2( BadInput, "The " << what << " has several settings"
                         " fitting the unit cell (" << listSettings( fit )
                         << "). Please specify which one to use (e.g."
                         " override_spacegroup=\"" << f0 << "\", or"
                         " --spacegroup=" << f0 << " for ncrystal_cif2ncmat)" );
      }

      //Settings from the space group symbol (Hall, else Hermann-Mauguin, else
      //number), and a description of the symbol used (or of the last one not
      //recognised, empty if there is none):
      struct SymbolSettings {
        std::vector<SpaceGroup> settings;
        std::string descr;
      };
      SymbolSettings symbolSettings( const StructureItems& s )
      {
        SymbolSettings res;
        const auto hall = firstValue( s.hallSymbol );
        if ( hall.has_value() ) {
          res.descr = "Hall symbol \"" + hall.to_string() + "\"";
          res.settings = settingsFromHallSymbol( hall );
          if ( !res.settings.empty() )
            return res;
        }
        const auto hm = firstValue( s.hmSymbol );
        if ( hm.has_value() ) {
          res.descr = "space group symbol \"" + hm.to_string() + "\"";
          res.settings = settingsFromHMSymbol( hm );
          if ( !res.settings.empty() )
            return res;
        }
        const auto n = firstValue( s.itNumber );
        if ( n.has_value() ) {
          res.descr = "space group number \"" + n.to_string() + "\"";
          const auto nn = n.trimmed();
          Optional<std::int64_t> i;
          if ( !nn.empty() && nn.contains_only( "0123456789" ) )
            i = nn.toInt();
          if ( i.has_value() && i.value() >= 1 && i.value() <= 230 ) {
            res.descr = "space group number " + std::to_string( i.value() );
            res.settings = SpaceGroup::settingsOfNumber(
              static_cast<unsigned>( i.value() ) );
          }
        }
        return res;
      }

      bool containsUncased( StrView hay, StrView needle )
      {
        if ( needle.size() > hay.size() )
          return false;
        for ( auto i : ncrange( hay.size() - needle.size() + 1 ) )
          if ( hay.substr( i ).startswithUncased( needle ) )
            return true;
        return false;
      }

      [[noreturn]] void throwSymopError( const std::string& why,
                                         const Document& doc,
                                         bool permissiveHint )
      {
        std::string msg = ( "The symmetry operations listed in the CIF data "
                            + why + ". Please use a CIF file that either lists"
                            " all atoms of the unit cell with space group P 1,"
                            " or lists the symmetry operations of a standard"
                            " space group setting consistent with its space"
                            " group symbol." );
        if ( permissiveHint )
          msg += ( " Alternatively, try permissive mode (at your own risk)." );
        const auto& raw = doc.textData().rawData();
        if ( containsUncased( StrView( raw.begin(), raw.end() - raw.begin() ),
                              "pymatgen" ) )
          msg += ( " This file was written by pymatgen, as are Materials"
                   " Project downloads. For Materials Project entries, the"
                   " \"Conventional Standard\" CIF download from the website"
                   " works, or pass the ID as \"mpid::<number>\" instead of a"
                   " file name to have NCrystal download that itself." );
        NCRYSTAL_THROW( BadInput, msg );
      }
    }
  }
}

NC::SpaceGroup NCC::chooseSetting( const Document& doc,
                                   const StructureItems& s,
                                   CIFPermissive permissive,
                                   StrView overrideSpacegroup,
                                   const CIFWarnFct& warnfct )
{
#ifndef NDEBUG
  nc_assert_always( bool(warnfct) );
  nc_assert_always( s.block( doc ).index() < doc.nBlocks() );
#endif
  if ( overrideSpacegroup.has_value() ) {
    const auto o = overrideSpacegroup.to_string();
    warnfct( "Overriding spacegroup to \"" + o + "\" due to explicit request" );
    auto settings = settingsFromHMSymbol( overrideSpacegroup );
    if ( settings.empty() ) {
      auto p = SpaceGroup::tryParse( overrideSpacegroup,
                                     SpaceGroupParseMode::Strict );
      if ( !p.setting.has_value() ) {
        std::string why = p.errorMsg;
        if ( !why.empty() && why.front() >= 'A' && why.front() <= 'Z' )
          why.front() = static_cast<char>( why.front() - 'A' + 'a' );
        NCRYSTAL_THROW2( BadInput, "Unknown spacegroup \"" << o << "\": not a"
                         " known Hermann-Mauguin symbol, and " << why );
      }
      settings.push_back( p.setting.value() );
    }
    return resolve( settings, s, "space group \"" + o + "\"" );
  }

  //The listed operations (the first tag with values):
  std::vector<StrView> ops;
  for ( auto& o : s.operations ) {
    for ( auto i : ncrange( o.values.size() ) )
      if ( !isNull( o.values[i] ) )
        ops.push_back( o.values[i].text );
    if ( !ops.empty() )
      break;
  }
  const auto sym = symbolSettings( s );
  if ( ops.empty() )
    return resolve( sym.settings, s, sym.descr );

  std::vector<SpaceGroup> fromOps;
  std::string why;
  std::vector<SymOp> symops;
  symops.reserve( ops.size() );
  for ( auto& o : ops ) {
    auto p = SymOp::tryParse( o );
    if ( !p.op.has_value() ) {
      why = "could not be parsed (" + p.errorMsg + ")";
      break;
    }
    symops.push_back( p.op.value() );
  }
  if ( why.empty() ) {
    fromOps = settingsFromOperations( symops );
    if ( fromOps.empty() )
      why = "are not those of a tabulated space group setting";
  }
  const bool isPermissive = ( permissive == CIFPermissive::YES );
  if ( !fromOps.empty() ) {
    //(a few settings share their operations, then any is fine):
    const auto& setting = fromOps.front();
    if ( !sym.settings.empty() ) {
      bool same = false;
      for ( auto& e : sym.settings )
        same = same || ( e.number() == setting.number() );
      if ( !same ) {
        why = ( "are those of space group "
                + std::to_string( setting.number() ) + ", but the "
                + sym.descr + " is that of space group "
                + std::to_string( sym.settings.front().number() ) );
        if ( !isPermissive )
          throwSymopError( why, doc, true );
        warnfct( "The symmetry operations listed in the CIF data " + why
                 + ". The listed operations are used (accepted due to"
                 " permissive mode)" );
      }
    }
    return setting;
  }
  if ( !isPermissive )
    throwSymopError( why, doc, !sym.settings.empty() );
  auto setting = resolve( sym.settings, s, sym.descr );
  warnfct( "The symmetry operations listed in the CIF data " + why + ", so the "
           + sym.descr + " is used (accepted due to permissive mode)" );
  return setting;
}
