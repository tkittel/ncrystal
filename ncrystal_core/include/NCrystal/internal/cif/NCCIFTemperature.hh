#ifndef NCrystal_CIFTemperature_hh
#define NCrystal_CIFTemperature_hh

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

    //The temperature (K) with which the Uiso values of CIF data are
    //associated, from the printed values of _diffrn_ambient_temperature and
    //_cell_measurement_temperature (in that order of preference; without
    //value if absent), or the chosen temperature (if has_value), with the
    //comment for the NCMAT data and whether it is a warning:
    //  * A chosen temperature must be in [0.1,3000] (else BadInput), and is
    //    used regardless of the fields (which the comment mentions).
    //  * Otherwise fields which are not a number in [0.1,3000] are BadInput,
    //    unless permissive (ignored with a warning), the preferred valid one
    //    is used, and another one must agree within 5% (BadInput, unless
    //    permissive, with a warning). Without valid fields, 293 K is used
    //    (with a warning):
    struct CIFUisoTemperature {
      double temperature;
      std::string comment;
      bool isWarning;
    };
    CIFUisoTemperature
    resolveUisoTemperature( StrView ambientTemperature,
                            StrView measurementTemperature,
                            Optional<double> chosen, CIFPermissive,
                            const CIFWarnFct& );
  }
}

#endif
