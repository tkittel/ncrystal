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

#include "NCrystal/internal/sgsym/NCSGIdentify.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCSGSymTable.hh"
#include <map>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {

    using SGList = std::vector<SpaceGroup>;

    //Sorted operations without duplicates, as lookup key:
    std::vector<SymOp> canonicalOps( std::vector<SymOp> ops )
    {
      std::sort( ops.begin(), ops.end() );
      ops.erase( std::unique( ops.begin(), ops.end() ), ops.end() );
      return ops;
    }

    template<class TFct>
    void forAllSettings( TFct fct )
    {
      for ( auto i : ncrange( SpaceGroupSettingIdx::max_value ) )
        fct( SpaceGroup( SpaceGroupSettingIdx{
              static_cast<std::uint16_t>( i + 1 ) } ) );
    }

    const std::map<std::vector<SymOp>,SGList>& opsMap()
    {
      static const auto s_map = []()
      {
        std::map<std::vector<SymOp>,SGList> m;
        forAllSettings( [&m]( const SpaceGroup& sg )
        {
          m[ canonicalOps( SGSymmetry::get( sg ).operations() ) ]
            .push_back( sg );
        } );
        return m;
      }();
      return s_map;
    }

    //Hall symbol with single spaces and lower case:
    std::string normaliseHall( StrView s )
    {
      std::string res;
      for ( auto& part : s.split() ) {
        if ( !res.empty() )
          res += ' ';
        res += part.to_string();
      }
      return lowerCase( res );
    }

    //Hermann-Mauguin symbol without spaces and underscores, in lower case,
    //with a trailing axes letter after a space (e.g. "R -3 m H") written as a
    //choice code (":h"):
    std::string normaliseHM( StrView s )
    {
      s = s.trimmed();
      std::string res;
      const bool trailing_axes
        = ( s.size() > 2 && s.find( ':' ) == StrView::npos
            && ( s.back() == 'H' || s.back() == 'h'
                 || s.back() == 'R' || s.back() == 'r' )
            && ( s[ s.size() - 2 ] == ' ' ) );
      for ( auto c : s ) {
        if ( c == ' ' || c == '\t' || c == '_' )
          continue;
        res += c;
      }
      if ( trailing_axes )
        res.insert( res.size() - 1, 1, ':' );
      return lowerCase( res );
    }

    //Short monoclinic symbol of a unique axis b setting (e.g. "P 21/c" for "P
    //1 21/c 1"), or an empty string:
    std::string shortMonoclinicHM( const SpaceGroup& sg )
    {
      if ( sg.number() < 3 || sg.number() > 15 || sg.choice()[0] != 'b' )
        return {};
      auto parts = StrView( sg.hmSymbol() ).split();
      nc_assert_always( parts.size() == 4 );
      return parts.at(0).to_string() + parts.at(2).to_string();
    }

    struct SymbolMaps {
      std::map<std::string,SGList> hall, hm_full, hm_base;
      std::vector<std::pair<std::string,SpaceGroup>> hm_eglide;//(full and base)
    };

    const SymbolMaps& symbolMaps()
    {
      static const auto s_maps = []()
      {
        SymbolMaps m;
        forAllSettings( [&m]( const SpaceGroup& sg )
        {
          m.hall[ normaliseHall( sg.hallSymbol() ) ].push_back( sg );
          const auto full = normaliseHM( sg.hmSymbol() );
          m.hm_full[ full ].push_back( sg );
          const auto base = full.substr( 0, full.find( ':' ) );
          m.hm_base[ base ].push_back( sg );
          const auto shortsym = shortMonoclinicHM( sg );
          if ( !shortsym.empty() )
            m.hm_base[ lowerCase( shortsym ) ].push_back( sg );
          const auto no = sg.number();
          if ( no == 39 || no == 41 || no == 64 || no == 67 || no == 68 ) {
            m.hm_eglide.emplace_back( full, sg );
            if ( base != full )
              m.hm_eglide.emplace_back( base, sg );
          }
        } );
        return m;
      }();
      return s_maps;
    }

    //Whether the symbol a, with "e" glides, matches b (in which the "e"
    //glides are written as a, b or c). Both normalised, the first character
    //(the lattice letter) must agree exactly:
    bool matchesWithEGlide( const std::string& a, const std::string& b )
    {
      if ( a.size() != b.size() || a.empty() || a.front() != b.front() )
        return false;
      for ( auto i : ncrange( std::size_t(1), a.size() ) ) {
        const char ca = a[i], cb = b[i];
        if ( ca == cb )
          continue;
        if ( ca == 'e' && ( cb == 'a' || cb == 'b' || cb == 'c' ) )
          continue;
        return false;
      }
      return true;
    }

    SGList sortedUnique( SGList v )
    {
      std::sort( v.begin(), v.end() );
      v.erase( std::unique( v.begin(), v.end() ), v.end() );
      return v;
    }
  }
}

NC::SGList NC::settingsFromOperations( const std::vector<SymOp>& ops )
{
  const auto& m = opsMap();
  auto it = m.find( canonicalOps( ops ) );
  return it == m.end() ? SGList{} : it->second;
}

NC::SGList NC::settingsFromHallSymbol( StrView s )
{
  const auto& m = symbolMaps().hall;
  auto it = m.find( normaliseHall( s ) );
  return it == m.end() ? SGList{} : it->second;
}

NC::SGList NC::settingsFromHMSymbol( StrView s )
{
  const auto& maps = symbolMaps();
  const auto key = normaliseHM( s );
  if ( key.empty() )
    return {};
  const bool has_choice = key.find( ':' ) != std::string::npos;
  const auto& m = has_choice ? maps.hm_full : maps.hm_base;
  auto it = m.find( key );
  if ( it != m.end() )
    return it->second;
  //"e" glides (never the lattice letter, which is the first character):
  if ( key.find( 'e', 1 ) == std::string::npos )
    return {};
  SGList res;
  for ( auto& e : maps.hm_eglide ) {
    if ( ( e.first.find( ':' ) != std::string::npos ) == has_choice
         && matchesWithEGlide( key, e.first ) )
      res.push_back( e.second );
  }
  return sortedUnique( std::move( res ) );
}
