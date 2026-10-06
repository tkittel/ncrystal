#ifndef NCrystal_SABCellInteg_FMA_hh
#define NCrystal_SABCellInteg_FMA_hh

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
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {

    //Shared by contribImpl and fillContribAtAlpha below; not itself
    //fma-dispatched (no std::fma calls), so needs no #ifdef split:
    inline double calc_bu_minus_bl_times_smiddle( bool is_bounded_by_both,
                                                  double dbpm,
                                                  double bu,
                                                  double bl,
                                                  double smiddle )
    {
      nc_assert(bu-bl > -0.01);//could be slightly negative due to
                               //numerical instabilities
      const double bumbl( is_bounded_by_both
                          ? 2.0*dbpm
                          : ncmax(0.0,bu-bl) );
      nc_assert(bumbl >= 0.0 );
      nc_assert(smiddle >= -1e-9 );
      const double res = ncmax(0.0,bumbl*smiddle);
      nc_assert(res >= 0.0);
      nc_assert(std::isfinite(res));
      return res;
    }

#ifndef NCRYSTAL_WIN_FMA
    namespace {
#endif

      NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
        void, fmaRampFill,
        ( double* ncrestrict out, unsigned i0, unsigned i1,
          double slope, double offset )
      )
      NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
        double, contribImpl,
        ( double fourE, double b1, double b2, double invdb,
          bool is_bounded_by_both,
          double alpha, double s_at_b1, double s_at_b2 )
      )
      NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
        void, fillContribAtAlpha,
        ( double* ncrestrict contrib_out, std::size_t npts,
          const double* ncrestrict Sb1_arr,
          const double* ncrestrict Sb2_arr,
          const double* ncrestrict alpha_arr,
          double foure, double E_div_kT,
          double cs_b1, double cs_b2, double invdb,
          bool is_bounded_by_betaminus,
          bool is_bounded_by_betaplus,
          bool is_bounded_on_both_sides )
      )

      NCRYSTAL_FMADISPATCH_DECLARATOR(void,fmaRampFill)
      ( double* ncrestrict out, unsigned i0, unsigned i1,
        double slope, double offset )
      {
        NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(fmaRampFill,(out,i0,i1,slope,offset));
        for ( unsigned i = i0; i < i1; ++i )
          out[i] = std::fma( slope, static_cast<double>(i), offset );
      }

      //Body of IntegrandOfA::contrib (NCSABCellInteg.cc), extracted as a free
      //function so it can be fma-dispatched:
      NCRYSTAL_FMADISPATCH_DECLARATOR(double,contribImpl)
      ( double fourE, double b1, double b2, double invdb,
        bool is_bounded_by_both,
        double alpha, double s_at_b1, double s_at_b2 )
      {
        NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
          contribImpl,
          (fourE,b1,b2,invdb,is_bounded_by_both,alpha,s_at_b1,s_at_b2) );
        nc_assert(b2>b1);
        //bl/bu via the robust getBetaMinus/getBetaPlus
        //(the naive a-+dbpm formula cancels catastrophically for a~=4E):
        const double bl = ncclamp( getBetaMinus(fourE*0.25,alpha), b1, b2 );
        const double bu = ncclamp( getBetaPlus(fourE*0.25,alpha), b1, b2 );
        const double dbpm = std::sqrt( fourE * alpha );
        //FIXME: ^^^ This calculates sqrt(4E*a) thrice instead of just once!

        //To find the contribution we integrate S(a,b) over [bl,bu]. This is
        //easy, since we always interpolate linearly in b:
        const double bmiddle( is_bounded_by_both ? alpha : (bu+bl)*0.5 );
        nc_assert(valueInInterval(-0.01,1.01,(bmiddle-b1)*invdb));
        const double rb = ncclamp((bmiddle-b1)*invdb,0.0,1.0);
        const double smiddle = nclerp( s_at_b1, s_at_b2, rb );
        return calc_bu_minus_bl_times_smiddle( is_bounded_by_both,
                                               dbpm, bu, bl, smiddle );
      }

      NCRYSTAL_FMADISPATCH_DECLARATOR(void,fillContribAtAlpha)
      ( double* ncrestrict contrib_out, std::size_t npts,
        const double* ncrestrict Sb1_arr,
        const double* ncrestrict Sb2_arr,
        const double* ncrestrict alpha_arr,
        double foure, double E_div_kT,
        double cs_b1, double cs_b2, double invdb,
        bool is_bounded_by_betaminus,
        bool is_bounded_by_betaplus,
        bool is_bounded_on_both_sides )
      {
        NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
          fillContribAtAlpha,
          (contrib_out,npts,Sb1_arr,Sb2_arr,alpha_arr,foure,E_div_kT,
           cs_b1,cs_b2,invdb,is_bounded_by_betaminus,is_bounded_by_betaplus,
           is_bounded_on_both_sides) );
        //contrib_out/Sb1_arr/Sb2_arr/alpha_arr are ncrestrict: at the call
        //site (impl_numIntRegion in NCSABCellInteg.cc) they are always
        //contrib_at_a (a local stack array) and two distinct SOfAlphaGrid
        //instances' member arrays (never overlapping, since they are separate
        //objects/members):
        //buffersDisjoint only compares addresses (contrib_out is still unset):
        // cppcheck-suppress uninitvar symbolName=contrib_out
        nc_assert( buffersDisjoint( contrib_out, npts, Sb1_arr, npts ) );
        nc_assert( buffersDisjoint( contrib_out, npts, Sb2_arr, npts ) );
        nc_assert( buffersDisjoint( contrib_out, npts, alpha_arr, npts ) );
        double bl(cs_b1), bu(cs_b2);
        for ( std::size_t i = 0; i < npts; ++i ) {
          double Sb1 = Sb1_arr[i];
          double Sb2 = Sb2_arr[i];
          double a = alpha_arr[i];
          //bl/bu via the robust getBetaMinus/getBetaPlus
          //(fixme: calculating expensive sqrt(4E*a) thrice instead of once!)
          double dbpm = std::sqrt( foure * a );//nb: expensive
          if ( is_bounded_by_betaminus )
            bl = getBetaMinus(E_div_kT,a);
          if ( is_bounded_by_betaplus )
            bu = ncmax(bl,getBetaPlus(E_div_kT,a));//ncmax as a safeguard
                                                   //against FP issues
          //To find the contribution we integrate S(a,b) over [bl,bu]. This is
          //easy, since we always interpolate linearly in b:
          const double bmiddle( is_bounded_on_both_sides ? a : (bu+bl)*0.5 );
          double rb = (bmiddle-cs_b1)*invdb;
          double smiddle = nclerp(Sb1,Sb2,rb);
          contrib_out[i] = calc_bu_minus_bl_times_smiddle( is_bounded_on_both_sides,
                                                           dbpm, bu, bl, smiddle );
        }
      }

#ifndef NCRYSTAL_WIN_FMA
    }
#endif
  }
}

#endif
