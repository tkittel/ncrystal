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

#include "NCrystal/internal/query/NCQuery.hh"
#include "NCQueryImpl.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    void queryimpl_structure( std::ostream& os, const Info& info )
    {
      //["structure"]: null, or the StructureInfo as a dict (the setting as
      //its string form, e.g. "227:2", or null if unknown):
      if ( info.isMultiPhase() || !info.hasStructureInfo() ) {
        os << "null";
        return;
      }
      const StructureInfo& si = info.getStructureInfo();
      streamJSONDictEntry( os, "spacegroup", si.spacegroup,
                           JSONDictPos::FIRST );
      if ( si.spacegroup_setting.get() )
        streamJSONDictEntry( os, "spacegroup_setting",
                             SpaceGroup( si.spacegroup_setting ).toString() );
      else
        streamJSONDictEntry( os, "spacegroup_setting", json_null_t{} );
      streamJSONDictEntry( os, "a", si.lattice_a );
      streamJSONDictEntry( os, "b", si.lattice_b );
      streamJSONDictEntry( os, "c", si.lattice_c );
      streamJSONDictEntry( os, "alpha", si.alpha );
      streamJSONDictEntry( os, "beta", si.beta );
      streamJSONDictEntry( os, "gamma", si.gamma );
      streamJSONDictEntry( os, "volume", si.volume );
      streamJSONDictEntry( os, "n_atoms", si.n_atoms, JSONDictPos::LAST );
    }

    void queryimpl_sites( std::ostream& os, const Info& info )
    {
      //["sites"]: per AtomInfo (in order), the crystallographic sites as a
      //list of {position,multiplicity} dicts, or null if unknown:
      if ( info.isMultiPhase() || !info.hasAtomInfo() ) {
        os << "null";
        return;
      }
      os << '[';
      bool first_ai = true;
      for ( auto& ai : info.getAtomInfos() ) {
        if ( !first_ai )
          os << ',';
        first_ai = false;
        if ( ai.sites().empty() ) {
          os << "null";
          continue;
        }
        os << '[';
        bool first = true;
        for ( auto& site : ai.sites() ) {
          if ( !first )
            os << ',';
          first = false;
          os << "{\"position\":[";
          streamJSON( os, site.position[0] );
          os << ',';
          streamJSON( os, site.position[1] );
          os << ',';
          streamJSON( os, site.position[2] );
          os << "],\"multiplicity\":" << site.multiplicity << '}';
        }
        os << ']';
      }
      os << ']';
    }
  }
}

void NC::JSONQueryInfo( std::ostream& os, const Info& info,
                        const Query& query )
{
  QueryImpl::validateQuery( query );
  const auto& key = query.front();
  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_structure = StrView::make("structure");
  constexpr auto sv_sites = StrView::make("sites");
  if ( query.size() != 1 )
    NCRYSTAL_THROW2(BadInput, "Invalid Info JSON query (no arguments"
                    " should follow key \""<<key<<"\")");
  if ( key == sv_list ) {
    os << "[\"structure\",\"sites\"]";
  } else if ( key == sv_structure ) {
    queryimpl_structure( os, info );
  } else if ( key == sv_sites ) {
    queryimpl_sites( os, info );
  } else {
    NCRYSTAL_THROW2(BadInput, "Invalid Info JSON query key: \""<<key<<'"');
  }
}
