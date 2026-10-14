#ifndef NCrystal_LegacyLatticeChecks_hh
#define NCrystal_LegacyLatticeChecks_hh

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

#include "NCrystal/core/NCTypes.hh"

//The unit cell checks of NCrystal 4.4 and earlier (from NCLatticeUtils.hh,
//hand-written per crystal system), kept here as a reference for
//completeAndCheckStructureInfoCell (NCInfoUtils.hh), which uses the cell
//constraints derived from the symmetry operations.

namespace NCRYSTAL_NAMESPACE {
  namespace LegacyLattice {

    //Validate that lattice lengths are compatible with given spacegroup. For
    //space groups where a==b or a==c, it is allowed to provide b=0 or c=0, and
    //the function will then update the values of b and/or c accordingly. Any
    //inconsistencies or errors will result in BadInput exceptions thrown. It is
    //allowed to provide spacegroup==0, in which case the function won't do much:
    void checkAndCompleteLattice( unsigned spacegroup, double a, double& b, double& c );

    //Whether the space group is one of the rhombohedral (R) ones, and whether
    //a unit cell of such a space group with lattice angle alpha (degrees) is
    //described in rhombohedral axes (a=b=c, alpha=beta=gamma) rather than in
    //hexagonal axes (alpha=beta=90, gamma=120):
    bool isRhombohedralSpaceGroup( int spacegroup );
    bool usesRhombohedralAxes( int spacegroup, double alpha );

    //Same for angles, except that all three angles might be left zero in case the
    //space group defines them:
    void checkAndCompleteLatticeAngles( unsigned spacegroup, double& alpha, double& beta, double& gamma );

    //Spacegroup utils:
    enum CrystalSystem {
      Triclinic,
      Monoclinic,
      Orthorhombic,
      Tetragonal,
      Trigonal,
      Hexagonal,
      Cubic
    };
    CrystalSystem crystalSystem( int spacegroup );

  }
}

#endif
