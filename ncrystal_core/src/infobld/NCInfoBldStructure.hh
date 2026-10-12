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
    //the space group), atom positions (remapped to [0,1), sorted, checked
    //for duplicates and for the symmetry of the space group) and the
    //volume:
    void validateAndCompleteUnitCell( InfoBuilder::UnitCell& );

    //Verify that the atom positions have the symmetry of the declared space
    //group (in some of its settings): slight deviations give a warning,
    //larger ones (or settings with other symmetry-equivalent reflections
    //than NCrystal assumes) a BadInput exception. See the .cc file:
    void checkStructureSymmetry( const StructureInfo&, const AtomInfoList& );

  }
}

#endif
