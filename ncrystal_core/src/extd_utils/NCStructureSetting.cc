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

#include "NCrystal/internal/extd_utils/NCStructureSetting.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"

namespace NC = NCrystal;

NC::Optional<NC::SpaceGroup>
NC::declaredStructureSetting( const StructureInfo& si )
{
  if ( !si.spacegroup_setting.get() )
    return NullOpt;
  si.spacegroup_setting.validate();
  const SpaceGroup setting( si.spacegroup_setting );
  if ( setting.number() != si.spacegroup )
    NCRYSTAL_THROW2(BadInput,"Declared space group setting "<<setting
                    <<" is not a setting of space group number "
                    <<si.spacegroup<<" (the settings are: "
                    <<SpaceGroup::listSettings( si.spacegroup )<<")");
  const bool rhombo = usesRhombohedralAxes( static_cast<int>( si.spacegroup ),
                                            si.alpha );
  if ( ( StrView( setting.choice() ) == "r" ) != rhombo )
    NCRYSTAL_THROW2(BadInput,"Declared space group setting "<<setting
                    <<" does not match the unit cell axes: the cell is given"
                    " in "<<( rhombo ? "rhombohedral axes (a=b=c,"
                                       " alpha=beta=gamma), for which the"
                                       " setting is "
                                     : "hexagonal axes (gamma=120), for"
                                       " which the setting is " )
                    <<si.spacegroup<<':'<<( rhombo ? 'r' : 'h' ));
  return setting;
}

void NC::completeAndCheckCell( StructureInfo& si )
{
  nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
  const auto declared = declaredStructureSetting( si );
  const SpaceGroup setting = structureSetting( si );
  const auto& cc = SGSymmetry::get( setting ).cellConstraints();
  CellParameters cp{ si.lattice_a, si.lattice_b, si.lattice_c,
                     si.alpha, si.beta, si.gamma };
  try {
    cc.complete( cp );
  } catch ( Error::BadInput& e ) {
    NCRYSTAL_THROW2(BadInput,"Unit cell not consistent with space group "
                    <<setting<<" ("<<cc
                    <<( !declared.has_value() && setting.hasMultipleSettings()
                        ? "; the setting is assumed, as none was declared"
                        : "" )
                    <<"): "<<e.what());
  }
  si.lattice_a = cp.a;
  si.lattice_b = cp.b;
  si.lattice_c = cp.c;
  si.alpha = cp.alpha;
  si.beta = cp.beta;
  si.gamma = cp.gamma;
}
