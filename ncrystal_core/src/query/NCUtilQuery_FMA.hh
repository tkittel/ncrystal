#ifndef NCrystal_UtilQuery_FMA_hh
#define NCrystal_UtilQuery_FMA_hh

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
      void, queryimpl_fmadiagnose_cmulDispatch,
      ( double* ncrestrict re, double* ncrestrict im,
        const double* ncrestrict cre, const double* ncrestrict cim,
        std::size_t n )
    )

    NCRYSTAL_FMADISPATCH_DECLARATOR(void,queryimpl_fmadiagnose_cmulDispatch)
    ( double* ncrestrict re, double* ncrestrict im,
      const double* ncrestrict cre, const double* ncrestrict cim,
      std::size_t n )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
        queryimpl_fmadiagnose_cmulDispatch, (re,im,cre,cim,n) );
      nc_assert( buffersDisjoint( re, n, im, n ) );
      nc_assert( buffersDisjoint( re, n, cre, n ) );
      nc_assert( buffersDisjoint( re, n, cim, n ) );
      nc_assert( buffersDisjoint( im, n, cre, n ) );
      nc_assert( buffersDisjoint( im, n, cim, n ) );
      nc_assert( buffersDisjoint( cre, n, cim, n ) );
      for ( std::size_t i = 0; i < n; ++i ) {
        double a = re[i], b = im[i], c = cre[i], d = cim[i];
        re[i] = std::fma( a, c, -(b*d) );
        im[i] = std::fma( a, d, b*c );
      }
    }

#ifndef NCRYSTAL_WIN_FMA
  }
#endif
}

#endif
