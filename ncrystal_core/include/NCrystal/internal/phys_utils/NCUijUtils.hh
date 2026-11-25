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
  //(U11,U22,U33,U12,U13,U23) of a symmetric 3x3 matrix, in the conventions of
  //the IUCr nomenclature (Trueblood et al., Acta Cryst. A52 (1996) 770):
  //
  //  Ustd:  U^ij, the standard form relative to the normalised reciprocal axes
  //         (Aa^2), as written in NCMAT and CIF data.
  //  Ustar: U*_ij = U^ij a*_i a*_j, the covariance of the displacements in
  //         fractional coordinates. A symmetry operation with rotation part R
  //         maps it to R Ustar R^T, and the Debye-Waller factor of a
  //         reflection h is exp(-2pi^2 h^T Ustar h).
  //  Ucart: the Cartesian form (Aa^2).

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

    Uij convertUstdToUstar( const Uij&, const CellMetric& );
    Uij convertUstarToUstd( const Uij&, const CellMetric& );

    //R Ustar R^T, with R the rotation part of the operation:
    Uij applySymOpToUstar( const SymOp&, const Uij& );

    //Equivalent isotropic displacement, trace(Ucart)/3:
    double calcUeqFromUstar( const Uij&, const CellMetric& );

    //Whether the matrix is positive definite (forms related by congruence,
    //like Ustd, Ustar and Ucart, agree on this):
    bool isPositiveDefinite( const Uij& );

    //Largest absolute value of the components:
    double maxAbs( const Uij& );

  }
}

#endif
