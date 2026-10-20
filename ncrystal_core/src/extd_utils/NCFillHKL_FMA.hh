#ifndef NCrystal_FillHKL_FMA_hh
#define NCrystal_FillHKL_FMA_hh

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
#include "NCrystal/internal/utils/NCVector.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace {
    NCRYSTAL_FMADISPATCH_INLINE
    double fillhkl_dot_fma( const Vector& a, const Vector& b )
    {
      return std::fma( a.x(), b.x(),
                       std::fma( a.y(), b.y(), a.z() * b.z() ) );
    }
  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, fillhkl_ksq,
    ( const double* factor, const unsigned* group_of_term, unsigned ngroups,
      int h, int k, int l )
  )
  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, fillhkl_ksq_approx,
    ( const double* m9, double h, double k, double l )
  )
  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, fillhkl_fsq,
    ( const Vector* hkl, const double* factors,
      const Vector* const* posarrs, const std::size_t* nposarr,
      std::size_t nspecies )
  )

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_fillhkl_ksq,fillhkl_ksq)
  ( const double* factor, const unsigned* group_of_term, unsigned ngroups,
    int h, int k, int l )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(fillhkl_ksq,(factor,group_of_term,ngroups,h,k,l));
    //|k|^2 from a KSqCalc (NCLatticeUtils.hh): must stay identical to
    //KSqCalc::ksq. Only the fma chain rounds (and identically everywhere):
    const std::int64_t H = h, K = k, L = l;//(64 bit: exact for |h| < 1e9)
    const std::int64_t intpart[6] = { H*H, K*K, L*L, H*K, H*L, K*L };
    std::int64_t intsum[7] = { 0, 0, 0, 0, 0, 0, 0 };
    for ( unsigned t = 0; t < 6; ++t )
      intsum[group_of_term[t]] += intpart[t];
    double acc = 0.0;
    for ( unsigned g = 0; g < ngroups; ++g )
      acc = std::fma( factor[g], static_cast<double>( intsum[g] ), acc );
    return acc;
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_fillhkl_ksq_approx,fillhkl_ksq_approx)
  ( const double* m9, double h, double k, double l )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(fillhkl_ksq_approx,(m9,h,k,l));
    //Cheap |k|^2 = (rec_lat*hkl).mag2() with rec_lat given as its 9
    //row-major elements, for the loose preselection of the box loop only
    //(within rounding of the exact kernel, which has a 1e-5 slack there):
    const double r0 = std::fma( m9[0], h, std::fma( m9[1], k, m9[2]*l ) );
    const double r1 = std::fma( m9[3], h, std::fma( m9[4], k, m9[5]*l ) );
    const double r2 = std::fma( m9[6], h, std::fma( m9[7], k, m9[8]*l ) );
    return std::fma( r0, r0, std::fma( r1, r1, r2*r2 ) );
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_fillhkl_fsq,fillhkl_fsq)
  ( const Vector* hkl, const double* factors,
    const Vector* const* posarrs, const std::size_t* nposarr,
    std::size_t nspecies )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(fillhkl_fsq,
                                        (hkl,factors,posarrs,
                                         nposarr,nspecies));

    //Time to calculate phases (hkl.dot(pos)*2pi; we speed up the expensive
    //sin+cos by a factor of 3 by shifting the phase to [0,2pi], easily done by
    //simply NOT multiplying with 2pi, and using our own fast sincos_02pi
    //through sincos_2pix.  Since typically ~99% of the hkl initialisation time
    //is spent calculating sin+cos here, that actually translates into an
    //overall speedup of a factor of 3 (measured in NCrystal v2.7.0)!
    //
    //For NCrystal 5.0.0 we additionally moved this crucial part to our newly
    //introduced FMA machinery.
    //
    //|F|^2 = real^2+imag^2 with real+i*imag = sum_i factor_i * sum_j
    //exp(2pi*i*hkl.dot(pos_ij)). Inner (atomic position) sums use StableSum
    //(additions only, so contraction-proof) of the canonical sincos_2pix
    //values; the small outer (species) sum and the final assembly are explicit
    //fma:
    double re(0.0), im(0.0);
    for ( std::size_t i = 0; i < nspecies; ++i ) {
      const double factor = factors[i];
      if ( !factor )
        continue;
      StableSum cpsum, spsum;
      const Vector* pos = posarrs[i];
      const std::size_t npos = nposarr[i];
      for ( std::size_t j = 0; j < npos; ++j ) {
        const auto spcp = sincos_2pix( fillhkl_dot_fma( *hkl, pos[j] ) );
        cpsum.add( spcp.cos );
        spsum.add( spcp.sin );
      }
      re = std::fma( factor, cpsum.sum(), re );
      im = std::fma( factor, spsum.sum(), im );
    }
    return std::fma( re, re, im*im );
  }

}

#endif
