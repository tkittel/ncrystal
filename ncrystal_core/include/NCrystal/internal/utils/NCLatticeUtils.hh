#ifndef NCrystal_LatticeUtils_hh
#define NCrystal_LatticeUtils_hh

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

#include "NCrystal/internal/utils/NCRotMatrix.hh"

namespace NCRYSTAL_NAMESPACE {

  //Construct cell rotation matrix from lattice parameters (in angstrom and
  //radians):
  RotMatrix getLatticeRot( double lattice_a, double lattice_b, double lattice_c,
                           double alpha, double beta, double gamma );

  // Validate orientation parameters (check against null-vectors, parallel
  // vectors, consistency of angles in the two frames given dirtol,
  // etc.). Throws BadInput exception in case of issues.
  void precheckLatticeOrientTol( double dirtol );
  void precheckLatticeOrientDir( const OrientDir& );
  void verifyLatticeOrientDef( const LabAxis& l1, const CrystalAxis& c1,
                               const LabAxis& l2, const CrystalAxis& c2,
                               double dirtol );

  //Same, but on OrientDir's which might include HKLPoint's, in which case the
  //validation won't be complete (due to the absence of a lattice rotation
  //matrix):
  void precheckLatticeOrientDef( const OrientDir& dir1, const OrientDir& dir2, double dirtol );

  //Fully validate orientation parameters and construct crystal rotation matrix
  //(transforms lab frame -> crystal frame). Requires a reciprocal lattice
  //rotation matrix.
  RotMatrix verifyLatticeOrientDefAndConstructCrystalRotation( const OrientDir& dir1,
                                                               const OrientDir& dir2,
                                                               double dirtol,
                                                               const RotMatrix& reci_lattice );


  //Reciprocal lattice rotation (inverted lattice transform, multiplied by 2pi):
  RotMatrix getReciprocalLatticeRot( double lattice_a, double lattice_b, double lattice_c,
                                     double alpha, double beta, double gamma );

  //Given lattice lengths, translate dcutoff into required maximum values of
  //|h|, |k| and |l| (i.e. all planes with d>=dcutoff have indices within
  //these, for any lattice angles). The result is at least 1 for each:
  struct MaxHKL{ int h, k, l; };
  MaxHKL estimateHKLRange( double dcutoff, double lattice_a,
                           double lattice_b, double lattice_c );

  //Calculator of |k|^2 = (2pi/d)^2 for Miller indices, from the textbook
  //expression
  //
  //   |k|^2 = F1*h^2 + F2*k^2 + F3*l^2 + F4*hk + F5*hl + F6*kl
  //
  //(the six factors F are the reciprocal metric of the cell). Two planes whose
  //|k|^2 are equal by the symmetry of the cell, e.g. 333 and 511 in a cubic
  //cell (both 27*F1), or 700 and 530 in hexagonal axes (both 49*F1), should
  //also give the same double, so that they compare and sort identically on
  //all platforms. This is achieved by forming the integer parts first (exact),
  //and by letting all terms whose factors are the same double share a single
  //integer sum, so that such planes are evaluated by identical operations:
  //the cubic 333 is fma(F1,27,0), and so is 511.
  struct KSqCalc {
    //Merged terms: group g has factor[g] multiplying the sum of the integer
    //parts (h^2,k^2,l^2,hk,hl,kl) of the terms t with group_of_term[t]==g:
    unsigned ngroups = 0;
    double factor[6] = {};
    unsigned group_of_term[6] = {};//(6 = dropped, zero factor)
    double ksq( int h, int k, int l ) const;//(NCFillHKL_FMA.hh has an
                                            //identical kernel)
  };
  KSqCalc getKSqCalc( double lattice_a, double lattice_b, double lattice_c,
                      double alpha, double beta, double gamma );//(radians)

  //Calculate d-spacing from Miller index and the reciprocal metric:
  double dspacingFromHKL( int h, int k, int l, const KSqCalc& );

  //Whether the space group is one of the rhombohedral (R) ones, and whether
  //a unit cell of such a space group with lattice angle alpha (degrees) is
  //described in rhombohedral axes (a=b=c, alpha=beta=gamma) rather than in
  //hexagonal axes (alpha=beta=90, gamma=120):
  bool isRhombohedralSpaceGroup( int spacegroup );
  bool usesRhombohedralAxes( int spacegroup, double alpha );

  //Try to estimate the HKL point based on normal + dspacing, essentially by hkl
  //= lattice_rot * normal / dspacing, but then returning the smallest of
  //(hkl,-hkl).  Will throw CalcError in case the normal does not actually
  //correspond to one coming from an integral hkl point.
  std::tuple<int,int,int> normalAndDSpacingToHKLIndex( const RotMatrix& lattice_rot,
                                                       double dspacing,
                                                       const Vector& normal );


}

#endif
