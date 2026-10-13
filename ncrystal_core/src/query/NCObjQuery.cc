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
      if ( !info.hasStructureInfo() ) {
        os << "null";
        return;
      }
      const StructureInfo& si = info.getStructureInfo();
      streamJSONDictEntry( os, "spacegroup", si.spacegroup,
                           JSONDictPos::FIRST );
      if ( si.spacegroup_setting.has_value() ) {
        const SpaceGroup setting( si.spacegroup_setting.value() );
        streamJSONDictEntry( os, "spacegroup_setting", setting.toString() );
      } else {
        streamJSONDictEntry( os, "spacegroup_setting", json_null_t{} );
      }
      streamJSONDictEntry( os, "a", si.lattice_a );
      streamJSONDictEntry( os, "b", si.lattice_b );
      streamJSONDictEntry( os, "c", si.lattice_c );
      streamJSONDictEntry( os, "alpha", si.alpha );
      streamJSONDictEntry( os, "beta", si.beta );
      streamJSONDictEntry( os, "gamma", si.gamma );
      streamJSONDictEntry( os, "volume", si.volume );
      streamJSONDictEntry( os, "n_atoms", si.n_atoms, JSONDictPos::LAST );
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
  if ( query.size() != 1 )
    NCRYSTAL_THROW2(BadInput, "Invalid Info JSON query (no arguments"
                    " should follow key \""<<key<<"\")");
  if ( key == sv_list ) {
    os << "[\"structure\"]";
  } else if ( key == sv_structure ) {
    queryimpl_structure( os, info );
  } else {
    NCRYSTAL_THROW2(BadInput, "Invalid Info JSON query key: \""<<key<<'"');
  }
}
