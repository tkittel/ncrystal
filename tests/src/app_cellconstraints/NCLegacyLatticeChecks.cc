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

#include "NCLegacyLatticeChecks.hh"

namespace NC = NCrystal;

bool NC::LegacyLattice::isRhombohedralSpaceGroup( int sg )
{
  return ( sg==146 || sg==148 || sg==155 || sg==160
           || sg==161 || sg==166 || sg==167 );
}

bool NC::LegacyLattice::usesRhombohedralAxes( int sg, double alpha )
{
  return isRhombohedralSpaceGroup( sg ) && alpha != 90.0;
}

void NC::LegacyLattice::checkAndCompleteLatticeAngles( unsigned sg, double& alpha, double& beta, double& gamma )
{
  if (sg>230)
    NCRYSTAL_THROW2(BadInput,"invalid spacegroup number ("<<sg<<")");
  if ( sg < 1 )
    return;
  auto cs = crystalSystem( sg );
  switch (cs) {
  case Orthorhombic:
  case Tetragonal:
  case Cubic:
    //all 90
    if ( ( alpha>0 && alpha!=90 ) || ( beta>0 && beta!=90 ) || ( gamma>0 && gamma!=90 ) )
      NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") requires alpha=beta=gamma=90");
    alpha = beta = gamma = 90;
    return;
  case Trigonal:
  case Hexagonal:
    if ( alpha > 0 && usesRhombohedralAxes( static_cast<int>(sg), alpha ) ) {
      //Rhombohedral axes (a=b=c must be checked by caller):
      if ( ( beta>0 && beta!=alpha ) || ( gamma>0 && gamma!=alpha )
           || !( alpha < 120 ) )
        NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") in rhombohedral axes"
                        " requires alpha=beta=gamma<120");
      beta = gamma = alpha;
      return;
    }
    if ( ( alpha>0 && alpha!=90 ) || ( beta>0 && beta!=90 ) || ( gamma>0 && gamma!=120 ) )
      NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") requires alpha=beta=90 and gamma=120");
    alpha = beta = 90;
    gamma = 120;
    return;
  case Monoclinic:
    //Only the setting with unique axis b is supported (in particular by the
    //EqRefl class, which assumes it):
    if ( ( alpha>0 && alpha!=90 ) || ( gamma>0 && gamma!=90 ) )
      NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") requires alpha=gamma=90"
                      " (only monoclinic settings with unique axis b are"
                      " supported)");
    if ( !( beta>0 && beta<180 ) )
      NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") requires beta to be"
                      " set (and to a value < 180).");
    alpha = gamma = 90;
    return;
  case Triclinic:
    //Simply require all three angles to be set with values <180
    if ( !( alpha>0 && alpha<180 ) || !( beta>0 && beta<180 ) || !( gamma>0 && gamma<180 ) )
      NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<sg<<") requires all three angles to be set (and to values < 180).");
    return;
  };
}


void NC::LegacyLattice::checkAndCompleteLattice( unsigned sg, double a, double& b, double & c )
{
  if (sg>230)
    NCRYSTAL_THROW2(BadInput,"invalid spacegroup number ("<<sg<<")");
  if (sg>0) {
    //lattice_b and lattice_c does not always need to be provided if identical
    //to lattice_a by spacegroup symmetry:
    bool b_is_a = false;
    bool c_is_a = false;
    if (sg >= 75 && sg <= 194 )
      b_is_a = true;
    else if ( sg >= 195 && sg <= 230 )
      c_is_a = b_is_a = true;

    if (b_is_a) {
      if (!b) {
        b = a;
      } else if (!(b==a)) {
        NCRYSTAL_THROW2(BadInput,"lattice lengths a and b must be equal for spacegroup "<<sg);
      }
    }
    if (c_is_a) {
      if (!c) {
        c = a;
      } else if (!(c==a)) {
        NCRYSTAL_THROW2(BadInput,"lattice lengths a and c must be equal for spacegroup "<<sg);
      }
    }
  }
  if ( !(a>0.0) || !(b>0.0) || !(c>0.0) )
    NCRYSTAL_THROW(BadInput,"lattice parameters must be positive numbers");
}

NC::LegacyLattice::CrystalSystem NC::LegacyLattice::crystalSystem( int sg )
{
  if ( sg < 1 || sg > 230 )
    NCRYSTAL_THROW(BadInput,"Space group number is not in the range 1 to 230");
  if ( sg >= 143 ) {
    if ( sg >= 195 )
      return CrystalSystem::Cubic;
    return sg >= 168 ? CrystalSystem::Hexagonal : CrystalSystem::Trigonal;
  }
  if ( sg >= 75 )
    return CrystalSystem::Tetragonal;
  if ( sg >= 16 )
    return CrystalSystem::Orthorhombic;
  return sg >= 3 ? CrystalSystem::Monoclinic : CrystalSystem::Triclinic;
}

