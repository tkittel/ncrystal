#ifndef NCrystal_CIFDescribe_hh
#define NCrystal_CIFDescribe_hh

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

namespace NCRYSTAL_NAMESPACE {

  namespace CIF {

    //Description of CIF data from its items, for NCMAT comments: lines with
    //the title, authors ("A", "A and B", or "A, et al.") with journal and
    //year, DOI, CIF creation method, date and author, and the COD or MP
    //entry with its URL (values trimmed, multi-line values split into
    //lines). The COD and MP ids are the requested ones (those of the data
    //source, 0 if none), else those embedded in the data (0 if none), with a
    //warning if both are given and differ. The formula is
    //_chemical_formula_sum, else _cod_original_formula_sum (empty if none):
    struct CIFDescription {
      std::vector<std::string> lines;
      std::uint64_t codID = 0;
      std::uint64_t mpID = 0;
      std::string formula;
    };
    CIFDescription getDescription( const StructureItems&,
                                   const Optional<ProvenanceItems>&,
                                   std::uint64_t requestedCodID,
                                   std::uint64_t requestedMPID,
                                   const CIFWarnFct& );
  }
}

#endif
