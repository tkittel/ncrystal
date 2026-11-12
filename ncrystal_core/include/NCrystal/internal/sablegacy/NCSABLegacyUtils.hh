#ifndef NCrystal_SABLegacyUtils_hh
#define NCrystal_SABLegacyUtils_hh

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

#include "NCrystal/internal/sab/NCSABUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABLegacy {

    //Kernels with fewer than 500 beta points (configurable with the
    //SAB_BETATHICKENING_MINNBETA env var) get extra beta rows, linearly
    //interpolated in S, between each pair of existing rows. This reduces
    //numerical artifacts of the legacy processing for direct kernels. Returns
    //NullOpt if the kernel is not thickened:
    Optional<SABData> thickenSparseBetaGrid( const SABData& );

    //Sample log-linear distribution on [a,b] with values fa, fb at a, b (and
    //their logs logfa, logfb):
    double sampleLogLinDist_fast(double a, double fa, double b, double fb, double rand, double logfa, double logfb);

    //Expand kinematic alpha-limits for ekin_div_kt = ekin/kT at each beta point
    //outwards to nearest points in the alpha-grid. Ranges are limited to the
    //extent of the alpha grid, and start at beta index ibeta_low. Empty ranges
    //are marked by indices both equal to alphaGrid.size() in the
    //out_alpharanges vector.
    void activeGridRanges( const SABData& data,
                           EDivKT,
                           std::vector<std::pair<uint16_t,uint16_t>>& out_alpharanges,
                           std::size_t& ibeta_low );

    //Break down how XS contributions falls over whole grid cells, as well as
    //front and back tails covering partial cells.
    struct TailedBreakdown {
      double xs_front=0.0, xs_middle=0.0, xs_back=0.0;
      unsigned imiddle_low=0, imiddle_upp=0;
      struct TailPoint {
        double alpha=0.0, sval=0.0, logsval=0.0;
      } front, back;
      bool narrow = false;//if range inside single grid bin
    };
    TailedBreakdown createTailedBreakdown( const Span<const double>& alphaGrid,
                                           const Span<const double>& sab,
                                           const Span<const double>& logsab,
                                           const Span<const double>& alphaIntegrals_cumul,
                                           double alpha_low, double alpha_upp,
                                           const unsigned aidx_low, const unsigned aidx_upp );

  }
}

////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::SABLegacy::sampleLogLinDist_fast(double a, double fa, double b, double fb, double rand, double logfa, double logfb)
{
  //fixme: new LogLinDistSampler class is better.
  nc_assert(b>a);
  nc_assert(fa>=0.);
  nc_assert(fb>=0.);
  double df = fb-fa;
  if ( fa > 0.0 && fb > 0.0 && df != 0.0 ) {
    //usual, non-degenerate case:
    double a_sub_b = a-b;
    double logfa_fb = logfb-logfa;
    if ( a_sub_b * logfa_fb != 0.0 )
      return a_sub_b*std::log(fa*std::exp(a*logfa_fb/a_sub_b)/(fa+rand*df))/logfa_fb;
    df = 0.0;
  }
  if (!df)
    return nclerp(a,b,rand);//fa=fb, select uniformly in [a,b]
  //exactly one of fa and fb is 0:
  nc_assert(fa||fb);
  double x = (b-a)*std::sqrt(rand);
  return (fa ? b-x : a + x);
}

#endif
