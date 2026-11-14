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

#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"

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
    //degrees). Throws BadInput if they do not form a valid cell:
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

    //Ueq of a Ustd tensor in a cell, and an estimate of the largest error in
    //it caused by the rounding errors of the tensor components and the cell
    //parameters (exact for the components, as Ueq is linear in them, and to
    //first order for the cell, by moving each parameter by its rounding
    //error):
    PairDD calcUeqWithRoundingError( const Uij& ustd, const CellParameters&,
                                     const Uij& ustdRoundingErrors,
                                     const CellParameters& cellRoundingErrors );
    //Whether the matrix is positive definite (forms related by congruence,
    //like Ustd, Ustar and Ucart, agree on this):
    bool isPositiveDefinite( const Uij& );

    //Largest absolute value of the components:
    double maxAbs( const Uij& );

    //Eigenvalues (ascending) of Ucart:
    std::array<double,3> calcUcartEigenvalues( const Uij&, const CellMetric& );

    //Ustd of a site at the position 'from', transformed to the
    //symmetry-equivalent position 'to' (positions as fractional coordinates,
    //compared periodically within 1e-4), by the first operation of the
    //setting mapping one onto the other. Throws BadInput if there is none:
    Uij mapUstdToPosition( const SGSymmetry&, const CellMetric&,
                           const Vector& from, const Vector& to, const Uij& );

    //Ustd of a site snapped to its site symmetry, i.e. averaged over the
    //operations of the setting leaving the position invariant. Throws BadInput
    //if the input clearly breaks the site symmetry (a component changing by
    //more than relTol times the largest component, plus twice the largest
    //rounding error of the input components, if given), or if the result is
    //not positive definite. The site description (e.g. "Al (0, 0, 0)") is used
    //in error messages:
    Uij snapUstdToSiteSymmetry( const SGSymmetry&, const CellMetric&,
                                const Vector& site, const Uij&,
                                const std::string& siteDescr,
                                double relTol = 0.05,
                                const Uij& roundingErrors = Uij{} );

    namespace detail {
      //Whether the tensors invariant under the rotation parts of ops (a site
      //symmetry) have the structure which snapUstdToSiteSymmetry imposes
      //exactly: each component zero, free, or a power of two (with sign)
      //times another free component, with the free components independent.
      //For tests:
      bool uijStructureIsExact( Span<const SymOp* const> ops );
    }

    //Ustd of an orbit given at each of its positions: the one at the first
    //position snapped to its site symmetry (cf. snapUstdToSiteSymmetry, with
    //the largest rounding errors of all given ones), which must map onto
    //those given at the other positions within the same tolerance (BadInput
    //otherwise):
    Uij ustdOfOrbit( const SGSymmetry&, const CellMetric&,
                     const std::vector<Vector>& positions,
                     const std::vector<Uij>& ustds,
                     const std::string& orbitDescr, double relTol = 0.05,
                     const Uij& roundingErrors = Uij{} );

    //Mean-squared displacements n^T Ucart n (from Ustar) along a fixed set of
    //36 directions n, with weights (summing to 1), so weighted sums approximate
    //averages over all directions of smooth functions of them (6x6
    //Gauss-Legendre points in an octant of the eigenframe; relative errors
    //~1e-6 for typical tensors and ~2e-4 for eigenvalue ratios of 25 in
    //the incoherent elastic cross section):
    std::vector<PairDD> directionalMSDs( const Uij&, const CellMetric& );

  }
}

#endif
