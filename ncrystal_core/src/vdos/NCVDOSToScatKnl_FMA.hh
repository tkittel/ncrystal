#ifndef NCrystal_VDOSToScatKnl_FMA_hh
#define NCrystal_VDOSToScatKnl_FMA_hh

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

#include "NCrystal/internal/utils/NCFMADispatch.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NCRYSTAL_NAMESPACE {
#ifndef NCRYSTAL_WIN_FMA
  namespace {
#endif

    NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
      void, vdosScatKnlAccumulate,
      ( double* ncrestrict sab, const double* ncrestrict aFact,
        double c, std::size_t n )
    )

    //sab[i] += aFact[i]*c, for i in [0,n) with explicit ncfma. Used in the
    //hottest loop of fillSABFromVDOS in NCVDOSToScatKnl.cc:
    NCRYSTAL_FMADISPATCH_DECLARATOR(void,vdosScatKnlAccumulate)
    ( double* ncrestrict sab, const double* ncrestrict aFact,
      double c, std::size_t n )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(vdosScatKnlAccumulate,(sab,aFact,c,n));
      nc_assert( buffersDisjoint( sab, n, aFact, n ) );
#ifdef __wasm__
      //ncfma is slow in wasm, so skip elements where it would return sab[i]
      //unchanged: |aFact[i]*c| below half an ulp of sab[i], with a margin for
      //rounding and underflow. Constants are 2^55 and 2^-1021:
      constexpr double k1 = 36028797018963968.0;
      constexpr double k2 = 2.0 * std::numeric_limits<double>::min();
#endif
      for ( std::size_t i = 0; i < n; ++i ) {
#ifdef __wasm__
        if ( ncabs( aFact[i] * c ) * k1 + k2 < ncabs( sab[i] ) )
          continue;
#endif
        sab[i] = ncfma( aFact[i], c, sab[i] );
      }
    }

#ifndef NCRYSTAL_WIN_FMA
  }
#endif
}

#endif
