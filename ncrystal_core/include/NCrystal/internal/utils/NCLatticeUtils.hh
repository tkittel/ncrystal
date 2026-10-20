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

  //|k|^2 = (2pi/d)^2 and d-spacings for Miller indices in a given unit cell.
  //Planes whose |k|^2 are equal by the symmetry of the cell (e.g. 333 and 511
  //in a cubic cell) get identical values, on all platforms:
  class KSqCalc final {
  public:
    KSqCalc( double lattice_a, double lattice_b, double lattice_c,
             double alpha, double beta, double gamma );//(angles in radians)
    double ksq( int h, int k, int l ) const;
    double dspacing( int h, int k, int l ) const;//(BadInput for hkl=000)
  private:
    //For the FMA-dispatched kernel in NCFillHKL.cc:
    friend inline double detail_ksqFillHKLKernel( const KSqCalc&,
                                                  int, int, int );
    unsigned m_ngroups = 0;
    double m_factor[6] = {};
    unsigned m_group_of_term[6] = {};
  };

  //Whether the space group is one of the rhombohedral (R) ones, and whether
  //a unit cell of such a space group with lattice angle alpha (degrees) is
  //described in rhombohedral axes (a=b=c, alpha=beta=gamma) rather than in
  //hexagonal axes (alpha=beta=90, gamma=120). An unspecified alpha (0) means
  //hexagonal axes:
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
