#ifndef NCrystal_SABLegacyFactory_hh
#define NCrystal_SABLegacyFactory_hh

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

#include "NCrystal/internal/sab/NCSABScatterHelper.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SAB {

    ///////////////////////
    // Legacy algorithms //
    ///////////////////////

    //The legacy algorithm is selected with knllux=-1. For reference and
    //validation work (the "legacy_oversample" type of the ["sab","refsample",...]
    //query), sampling artifacts can be reduced by sampling at higher energies
    //(increasing SABSampler::EGridMargin):

    enum class LegacySABAlgOpts { DEFAULT, OVERSAMPLE10 };

    //Direct factory function with no caching:
    std::unique_ptr<const SABScatterHelper>
    createScatterHelper( shared_obj<const SABData>,
                         std::shared_ptr<const VectD> energyGrid = nullptr,
                         LegacySABAlgOpts = LegacySABAlgOpts::DEFAULT );

    //Same with caching:
    shared_obj<const SABScatterHelper>
    createScatterHelperWithCache( shared_obj<const SABData>,
                                  std::shared_ptr<const VectD> egrid = nullptr,
                                  LegacySABAlgOpts = LegacySABAlgOpts::DEFAULT );

  }

}

#endif
