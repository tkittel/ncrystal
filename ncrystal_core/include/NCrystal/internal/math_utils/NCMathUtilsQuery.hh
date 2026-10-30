#ifndef NCrystal_MathUtilsQuery_hh
#define NCrystal_MathUtilsQuery_hh

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

#include "NCrystal/core/NCSmallVector.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace MathUtils {

    // JSON queries with key "mathutils", giving access to the utilities for
    // validation and benchmarking (cf. NCQuery.hh). Arrays are passed as input
    // array placeholders (cf. registerJSONInputDblVect in NCString.hh):
    //
    // ["mathutils","fourierint",FR,FI,W,OPTIONS...]: fourierIntegral(..) with
    // FI an array or "none", and options "h=VAL" (required), "scheme=NAME"
    // (trap, trapec, lin or cubic; default lin), "eval=NAME" (direct, fftlin
    // or fftcubic; default direct) and "oversample=N" (default 8). Returns
    // {"F":[...],"cpu_seconds":VAL}.
    using Query = SmallVector<StrView,8>;
    void JSONQuery( std::ostream&, const Query& );

  }
}

#endif
