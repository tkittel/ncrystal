#ifndef NCrystal_InfoBldStructure_hh
#define NCrystal_InfoBldStructure_hh

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
#include "NCrystal/internal/infobld/NCInfoBuilder.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilderDetail {

    //Total number of atoms in the unit cell (sum of numberPerUnitCell()),
    //validating that the list is non-empty and that no entry has zero atoms:
    std::size_t totalNumberOfAtomsInUnitCell( const AtomInfoList& );

    //Validate and complete the unit cell: cell parameters (completed from
    //the symmetry of the setting), atom positions (expanded from sites if
    //given as such, remapped to [0,1), sorted and checked for duplicates),
    //the volume, the space group setting (checked against the positions),
    //and the sites (derived from the positions):
    void validateAndCompleteUnitCell( InfoBuilder::UnitCell& );

    //Verify that the atom positions have the symmetry of the declared space
    //group, in the setting declared by the data source if any
    //(spacegroup_setting set), else in some of its settings: slight
    //deviations give a warning, larger ones (or found settings with other
    //symmetry-equivalent reflections than NCrystal assumes) a BadInput
    //exception. Returns the setting (see the .cc file):
    SpaceGroupSettingIdx checkStructureSymmetry( const StructureInfo&,
                                                 const AtomInfoList& );

    //The setting declared by the data source if any, else the one implied by
    //the space group number and cell alone (R axes from the cell angles),
    //when the number has a single setting, else 0:
    SpaceGroupSettingIdx settingWithoutPositions( const StructureInfo& );

  }
}

#endif
