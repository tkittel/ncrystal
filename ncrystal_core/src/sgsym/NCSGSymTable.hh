#ifndef NCrystal_SGSymTable_hh
#define NCrystal_SGSymTable_hh

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

#include "NCrystal/internal/sgsym/NCSpaceGroupSettingIdx.hh"

namespace NCRYSTAL_NAMESPACE {

  //Access the hardwired table of the 530 space group settings,
  //i.e. the rows of Table A1.4.2.7 in:
  //
  //    Shmueli, U., Hall, S. R. & Grosse-Kunstleve, R. W. (2010).
  //    "Space-group symbols for numeric and symbolic computations."
  //    International Tables for Crystallography, Volume B: Reciprocal
  //    Space, 3rd ed., Appendix A1.4.2, pp. 122-134 (Table A1.4.2.7).
  //
  //Each entry holds the first three columns of a row.

  struct SGSymTableEntry {
    std::uint8_t number;//space group number (1..230)
    const char * choice;//setting choice code ("" if the number has a single
                        //setting, otherwise never empty)
    const char * hm;//Hermann-Mauguin symbol (ITVB's computer-entry format)
    const char * hall;//Hall symbol (ITVB's "Hall entry" format, lower-case
                      //lattice letter)
  };

  const SGSymTableEntry& sgsymTableEntry( SpaceGroupSettingIdx );

}

#endif
