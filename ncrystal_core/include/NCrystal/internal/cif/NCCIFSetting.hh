#ifndef NCrystal_CIFSetting_hh
#define NCrystal_CIFSetting_hh

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
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace CIF {

    //The space group setting of the structure items (extracted from the
    //Document): from the listed symmetry operations (the first tag with
    //values), cross-checked against the space group symbol, or else from the
    //symbol (the Hall symbol, else the Hermann-Mauguin symbol, else the
    //number; the first value of each which is not "?" or "."), where several
    //candidate settings are resolved by the cell (lenient tolerances, with
    //the rounding errors of the printed values). An overrideSpacegroup (if
    //has_value: a setting like "70:1", a Hermann-Mauguin or a Hall symbol,
    //optionally prefixed with "HM:" or "Hall:"; unprefixed symbols whose two
    //readings differ, like "P 2", are BadInput) replaces the symmetry
    //information of the data (with a warning). Inconsistent or
    //unusable information is BadInput, unless permissive, where the listed
    //operations (or else the symbol) are used, with a warning:
    SpaceGroup chooseSetting( const Document&, const StructureItems&,
                              CIFPermissive, StrView overrideSpacegroup,
                              const CIFWarnFct& );

    //The cell of the structure items snapped to the constraints of the
    //setting (lenient tolerances, with the rounding errors of the printed
    //values, cf. SGCellConstraints::trySnap), or why it does not fit. Throws
    //BadInput if a cell value is not a number:
    SGCellConstraints::SnapResult snapCell( const Document&,
                                            const StructureItems&,
                                            const SpaceGroup& );
  }
}

#endif
