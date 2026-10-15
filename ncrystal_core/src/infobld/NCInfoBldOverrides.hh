#ifndef NCrystal_InfoBldOverrides_hh
#define NCrystal_InfoBldOverrides_hh

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

#include "NCrystal/interfaces/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilderDetail {

    //Returns orig if its overrideable fields (density, number density,
    //cfg-data, phase list) already equal the given ones. Otherwise returns an
    //Info object sharing orig's underlying data with the given fields, reusing
    //an object from the cache of that data if one with identical fields exists
    //(the cache is bounded, keeping its oldest and newest entries):
    InfoPtr overrideInfoFieldsWithCache( InfoPtr orig,
                                         const Info::OverrideableDataFields& );

    //Whether two phase lists have identical fractions and Info objects:
    bool phaseListsIdentical( const Info::PhaseList&, const Info::PhaseList& );

  }
}

#endif
