#ifndef NCrystal_UijUtils_hh
#define NCrystal_UijUtils_hh

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

#include "NCrystal/internal/sgsym/NCSymOp.hh"

namespace NCRYSTAL_NAMESPACE {

  //Utilities for anisotropic displacement tensors, stored as the components
  //(U11,U22,U33,U12,U13,U23) of a symmetric 3x3 matrix. The CIF form U^ij is
  //relative to the normalised reciprocal axes, the fractional form is
  //U_frac,ij = U^ij a*_i a*_j (the covariance of the displacements in
  //fractional coordinates), so a symmetry operation with rotation part R maps
  //U_frac to R U_frac R^T, and the Debye-Waller factor of a reflection h is
  //exp(-2pi^2 h^T U_frac h).

  namespace UijUtils {

    using Uij = std::array<double,6>;

    //Real space metric (G_ij = a_i.a_j, row-major) and reciprocal axis
    //lengths a*_i (without 2pi factor), from cell parameters (angstrom and
    //degrees):
    struct CellMetric {
      std::array<double,9> g;
      std::array<double,3> recilen;
    };
    CellMetric cellMetric( double a, double b, double c,
                           double alpha, double beta, double gamma );

    Uij cifToFrac( const Uij&, const CellMetric& );
    Uij fracToCif( const Uij&, const CellMetric& );

    //R U_frac R^T, with R the rotation part of the operation:
    Uij transformFrac( const SymOp&, const Uij& );

    //Equivalent isotropic displacement, trace(U_cart)/3, of U_frac:
    double ueqFromFrac( const Uij&, const CellMetric& );

    //Whether the matrix is positive definite (forms related by congruence,
    //like the CIF, fractional and Cartesian forms, agree on this):
    bool isPositiveDefinite( const Uij& );

    //Largest absolute value of the components:
    double maxAbs( const Uij& );

    //Eigenvalues (ascending) of the Cartesian form of U_frac:
    std::array<double,3> cartEigenvalues( const Uij&, const CellMetric& );

    //Mean-squared displacements n^T U_cart n along a fixed set of 36
    //directions n, with weights (summing to 1), so weighted sums approximate
    //averages over all directions of smooth functions of them (6x6
    //Gauss-Legendre points in an octant of the eigenframe; relative errors
    //~1e-6 for typical tensors and ~2e-4 for eigenvalue ratios of 25 in
    //the incoherent elastic cross section):
    std::vector<PairDD> directionalMSDs( const Uij&, const CellMetric& );

  }
}

#endif
