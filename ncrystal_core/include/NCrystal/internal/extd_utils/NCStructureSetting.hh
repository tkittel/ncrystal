#ifndef NCrystal_StructureSetting_hh
#define NCrystal_StructureSetting_hh

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

#include "NCrystal/interfaces/NCInfoTypes.hh"
#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  //The space group setting whose symmetry a crystal structure has: the one
  //stored on the StructureInfo (declared by the data source or detected from
  //the atom positions), else the default setting of the number (rhombohedral
  //axes from the cell angles). This is what hkl-related code (EqRefl) must
  //use. Requires spacegroup != 0:
  inline SpaceGroup structureSetting( const StructureInfo& si )
  {
    nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
    if ( si.spacegroup_setting.get() )
      return SpaceGroup( si.spacegroup_setting );
    std::string spec = std::to_string( si.spacegroup );
    if ( usesRhombohedralAxes( static_cast<int>( si.spacegroup ), si.alpha ) )
      spec += ":r";
    return SpaceGroup( spec, SpaceGroupParseMode::Relaxed );
  }

  //The setting declared by the data source (spacegroup_setting set before
  //building the Info), validated: must be a setting of the number and match
  //the hexagonal/rhombohedral axes of the cell. NullOpt if none:
  Optional<SpaceGroup> declaredStructureSetting( const StructureInfo& );

  //Fill in lattice parameters given as 0 which are implied by the symmetry
  //of the setting (e.g. b=c=a for cubic, gamma=120 for hexagonal axes) and
  //check all of them against the cell constraints of the setting (declared
  //if any, else the default of the number). Exact equality is required for
  //constrained values. Throws BadInput. Requires spacegroup != 0:
  void completeAndCheckCell( StructureInfo& );

}

#endif
