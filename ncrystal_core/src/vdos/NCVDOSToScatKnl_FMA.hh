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

    //sab[i] += aFact[i]*c, for i in [0,n) with explicit std::fma. Used in the
    //hottest loop of fillSABFromVDOS in NCVDOSToScatKnl.cc:
    NCRYSTAL_FMADISPATCH_DECLARATOR(void,vdosScatKnlAccumulate)
    ( double* ncrestrict sab, const double* ncrestrict aFact,
      double c, std::size_t n )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(vdosScatKnlAccumulate,(sab,aFact,c,n));
      nc_assert( buffersDisjoint( sab, n, aFact, n ) );
      for ( std::size_t i = 0; i < n; ++i )
        sab[i] = std::fma( aFact[i], c, sab[i] );
    }

#ifndef NCRYSTAL_WIN_FMA
  }
#endif
}

#endif
