#ifndef NCrystal_SGExactOrbit_hh
#define NCrystal_SGExactOrbit_hh

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

#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"

//Internals of expandSiteToOrbit (private to the sgsym component, the detail
//namespace keeps them accessible to unit tests).

namespace NCRYSTAL_NAMESPACE {
  namespace detail {

    //Tolerances of expandSiteToOrbit (see NCSGSymmetry.hh):
    constexpr double site_merge_tol = 5e-4;
    constexpr double site_distinct_tol = 1e-2;
    [[noreturn]] void throwNearSpecial( const Vector& site, double dist );

    //Same as expandSiteToOrbit, but in exact integer arithmetic for sites
    //whose coordinates are recognised as rationals (see rationalFromDouble),
    //so the positions are the correctly rounded values of the exact images.
    //Returns false (leaving res untouched) if not applicable:
    bool expandSiteToOrbitExact( const SGSymmetry&, const Vector& site,
                                 SGSiteOrbit& res );

    //The fraction p/q with the smallest q (q<=maxDenominator) whose correctly
    //rounded double is exactly x, e.g. 1/10 for 0.1 and 49/160 for 0.30625,
    //but NullOpt for 1.0-0.9. NullOpt also if x is not finite or |x|>2^53:
    Optional<std::pair<std::int64_t,std::int64_t>>
    rationalFromDouble( double x, std::int64_t maxDenominator );

  }
}

#endif
