#ifndef NCrystal_InfoBldInternalState_hh
#define NCrystal_InfoBldInternalState_hh

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

  //NB: Info::InternalState is only declared fully in this private header
  //(included by NCInfoBuilder.cc and NCInfoBldOverrides.cc), so only those
  //files can directly construct Info objects!
  struct Info::InternalState {
    shared_obj<const Info::Data> data;
    std::shared_ptr<const Info::OverrideableData> oData;
  };

  inline Info::Info( InternalState&& state )
    : m_data(std::move(state.data)),
      m_oData(std::move(state.oData))
  {
    nc_assert_always(!m_data->composition.empty());
  }

}

#endif
