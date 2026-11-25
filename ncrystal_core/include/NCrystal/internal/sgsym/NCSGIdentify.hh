#ifndef NCrystal_SGIdentify_hh
#define NCrystal_SGIdentify_hh

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

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"

//Identification of tabulated space group settings (cf. NCSpaceGroup.hh) from
//their symmetry operations or from Hall or Hermann-Mauguin symbols.

namespace NCRYSTAL_NAMESPACE {

  //The settings with exactly the given operations (any order, duplicates
  //allowed), including all centring copies, or none. Usually at most one: only
  //three pairs of origin choice 1 settings of number 68 (e.g. "C c c a:1" and
  //"C c c b:1") share their operations (and Hall symbols). The settings of all
  //530 table entries are created on the first call:
  std::vector<SpaceGroup> settingsFromOperations( const std::vector<SymOp>& );

  //The settings matching a symbol, in table order, or none. Hall symbols (e.g.
  //"-P 2ybc") are compared after normalising case and spaces. Hermann-Mauguin
  //symbols (e.g. "P 1 21/c 1", "P21/c", "R -3 m:H", "R -3 m H" or "P n n n:2")
  //are compared ignoring spaces, underscores and case. A Hermann-Mauguin symbol
  //without choice code (e.g. "R -3 m" or "P n n n") matches all settings with
  //that symbol, short monoclinic symbols (e.g. "P 21/c") imply unique axis b,
  //and "e" glides (e.g. "C m c e" for "C m c a") are recognised for the five
  //numbers using them (39, 41, 64, 67, 68). The two kinds are separate since
  //some strings are valid as both, with different meanings (e.g. "P 2"):
  std::vector<SpaceGroup> settingsFromHallSymbol( StrView );
  std::vector<SpaceGroup> settingsFromHMSymbol( StrView );

}

#endif
