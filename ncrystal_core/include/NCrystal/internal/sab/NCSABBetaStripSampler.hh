#ifndef NCrystal_SABBetaStripSampler_hh
#define NCrystal_SABBetaStripSampler_hh

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

#include "NCrystal/interfaces/NCSABData.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {


    ////////////////////////////////////////////////////////////////////////////
    //                                                                        //
    // Exact sampling of (alpha,beta) from S(alpha,beta) within the           //
    // kinematically accessible region for a neutron of energy E, suitable    //
    // for production usage also in the limit E->0.                           //
    //                                                                        //
    ////////////////////////////////////////////////////////////////////////////

    // Details are given in the following, using e=E/kT, and KB(e) to denote the
    // kinematically given phasespace at that energy, namely the (alpha,beta)
    // points with (alpha-beta)^2<=4*e*alpha. A given instance of the
    // SABBetaStripSampler only works for e values up to some maximal energy,
    // emax=Emax/kT.
    //
    //When to use: At low E, where full-cell sampling is inefficient. The
    //             provided samplings are exact for all e<=emax, but efficiency
    //             degrades as emax grows, and also for e far below emax unless
    //             emax itself is small. So a solution is to use several
    //             objects, with emax values a small factor apart, using each
    //             for e just below its emax and the lowest one for all e below
    //             that. Callers can use expectedTries() (which is a calculated
    //             1/acceptance-rate) and the relative cost per try to decide at
    //             which energies full-cell sampling is better.
    //
    //Method: Rejection sampling from a carefully weighted striped overlay. More
    //        specifically, beta is split into strips (nStripsPerBetaCell equal
    //        subdivisions of each beta cell of the kernel), and each beta strip
    //        is then chopped along alpha so the resulting rectangular region
    //        fits KB(emax) as tightly as possible. The maximum S value within
    //        this rectangular area, Smax, is then cached, as is a strip weight
    //        which takes both Smax and the shape of KB(emax) across the strip
    //        into account. During sampling, a strip is first picked according
    //        to the weight, then beta is sampled uniformly in the strip and
    //        alpha uniformly in KB(e), and the chosen point is then accepted
    //        with probability S(alpha,beta)/Smax times another factor again
    //        correcting for the alpha width of KB(e) varying across the strip.

    class SABBetaStripSampler final : private MoveOnly {
    public:
      SABBetaStripSampler( shared_obj<const SABData>, double emax_div_kT,
                           unsigned nStripsPerBetaCell = 4 );

      //Sample at 0<e<=emax. Requires S>0 somewhere within KB(e), and is only
      //efficient if expectedTries(e,..) is small:
      struct Result { double alpha, beta; unsigned ntries; };
      Result sample( RNG&, double e ) const;

      double emaxDivKT() const noexcept { return m_emax; }

      //Expected tries per sample call at e, given the integral of S within
      //KB(e):
      double expectedTries( double e, double integralWithinKB ) const
      {
        nc_assert(integralWithinKB>0.0);
        return m_totalWeight * 4.0 * std::sqrt( e ) / integralWithinKB;
      }

      //Rectangle [b0,b1]x[a0,a1] and S bound of each strip used (strips
      //without S>0 are dropped), from the same calcStripBound calls as in the
      //constructor. For tests:
      struct StripBound { double b0, b1, a0, a1, smax; };
      std::vector<StripBound> stripBounds() const;

    private:
      //Strip within beta cell ib, with alpha cells ia0..ia1 covering KB(e):
      struct Strip {
        double b0, b1, smax, b1plusemax;
        std::uint32_t ib, ia0, ia1;
      };
      double evalSCell( double alpha, double beta,
                        std::size_t ia, std::size_t ib ) const;
      //Strip bound, plus the alpha cells containing a0 and a1:
      struct StripCalc { StripBound bound; std::uint32_t ia0, ia1; };
      StripCalc calcStripBound( double b0, double b1, std::uint32_t ib ) const;
      shared_obj<const SABData> m_sab;
      VectD m_logS;
      double m_emax;
      std::vector<Strip> m_strips;
      //Alias table for picking strips:
      VectD m_aliasProb;
      std::vector<std::uint32_t> m_alias;
      double m_totalWeight = 0.0;
    };
  }
}

#endif
