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

#ifndef NCrystal_TDProd_hh
#define NCrystal_TDProd_hh

#include "NCrystal/factories/NCFactImpl.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace FactImpl {

    //TextData object for path and tds, which is the previously created object
    //if an identical one is still cached (then wasCached is true):
    struct ProducedTextData {
      TextDataSP td;
      bool wasCached;
    };
    ProducedTextData produceTextDataSP_PreferPreviousObject( const TextDataPath&,
                                                             TextDataSource&& );

  }

}

#endif
