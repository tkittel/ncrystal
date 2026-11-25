#ifndef NCrystal_CIFSites_hh
#define NCrystal_CIFSites_hh

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
#include "NCrystal/internal/phys_utils/NCUijUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace CIF {

    //An atom of the crystal structure of CIF data: one or more sites (with
    //positive occupancy) at the same positions, merged into a composition:
    struct CIFAtom {
      std::vector<std::string> labels;//of the merged sites, in data order
      //Fraction of each element (or isotope, e.g. "Al" or "H2"), sorted:
      std::vector<std::pair<double,std::string>> composition;
      double occupancy;//total, at most 1
      //Isotropic displacement (Aa^2), and the anisotropic one (Ustd, at the
      //first position, rescaled to have uiso as Ueq), absent if not given:
      Optional<double> uiso;
      Optional<UijUtils::Uij> uij;
      std::vector<Vector> positions;//orbit, symmetrised site first
      std::vector<UijUtils::Uij> uijs;//at each position (empty without uij)
      //A split site (cf. expandSiteToOrbit): number of symmetry-equivalent
      //positions near each other (0 if not split) and their distance (Aa):
      unsigned splitPositions = 0;
      double splitDistance = 0.0;
    };

    //The atoms of the structure items (extracted from the Document), in the
    //given space group setting. Messages name the data and the line of the
    //site concerned:
    //  * Sites with occupancy 0 (or below) are skipped, and the element of a
    //    site is taken from its type symbol, or else from its label checked
    //    against the chemical formula (cf. elementFromCIFTypeSymbol and
    //    elementFromCIFLabel, with a warning).
    //  * Sites are expanded into their orbits with lenient tolerances (cf.
    //    expandSiteToOrbit), with their tensors snapped to the site symmetry.
    //    Without Uiso (or Biso), the uiso is the Ueq of the tensor; with both,
    //    they must agree within the rounding errors of the printed values plus
    //    tolUij relative to the larger (else BadInput, or the Uiso is used
    //    with a warning if permissive), and the tensor is rescaled to it.
    //    Tensors with all components zero are treated as absent.
    //  * Sites at the same positions are merged (their displacements must be
    //    identical, and their total occupancy at most 1), and the total
    //    number of positions must agree with _atom_type_number_in_cell (1%)
    //    if given.
    //  * Split sites (cf. expandSiteToOrbit) are BadInput unless permissive
    //    (with a warning).
    //Anisotropic displacements given as B are ignored (with a warning):
    std::vector<CIFAtom> getAtoms( const Document&, const StructureItems&,
                                   const SpaceGroup&, CIFPermissive,
                                   double tolUij, const CIFWarnFct& );
  }
}

#endif
