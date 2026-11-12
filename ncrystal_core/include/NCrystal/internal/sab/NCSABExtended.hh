#ifndef NCrystal_SABExtended_hh
#define NCrystal_SABExtended_hh

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

#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/sab/NCSABExtender.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SAB {

    class SABExtended final : private MoveOnly {

    public:

      //Class combining SABProcessor and SABExtender instances, in order to
      //extend cross-sections and sampling capabilities of an S(alpha,beta)
      //table to large energies.

      //Constructor initialises based on an extender and processor. IMPORTANT
      //NOTE: It is very important that the SABExtender was created with
      //SigmaBound=1 (fixme: redesign so this is always guaranteed?).
      SABExtended( shared_obj<const SABProcessor>,
                   shared_obj<const SAB::SABExtender> );

      //Convenience functions for initialising using a free-gas extension:
      static shared_obj<const SABExtended>
      createWithFGExtender( shared_obj<const SABProcessor> );

      static shared_obj<const SABExtended>
      createWithFGExtender( const SABCfg::Cfg&,
                            shared_obj<const SABData>,
                            std::shared_ptr<const VectD>
                            energyGrid = nullptr );

      //Same, for short-collision-time extension:
      static shared_obj<const SABExtended>
      createWithSCTExtender( shared_obj<const SABProcessor>,
                             Temperature teff,
                             const SABCfg::Cfg& );

      //Standard entry point selecting the extender from
      //cfg.extender_model: SCT when an effective temperature is
      //provided (and cfg does not select FreeGas), otherwise the
      //free-gas extender:
      static shared_obj<const SABExtended>
      create( const SABCfg::Cfg&,
              shared_obj<const SABData>,
              Optional<Temperature> teff,
              std::shared_ptr<const VectD> energyGrid = nullptr );

      //Access cross sections, always assuming SigmaBound=1barn.
      CrossSect crossSectionUnitSigmaBound( NeutronEnergy ekin ) const;

      //The evaluator of the table defining crossSectionUnitSigmaBound (see
      //NCIsoXSUtils.hh), and the kinks of that cross section (the ends of the
      //processor's energy grid):
      const IsoXSEvaluator& xsEvaluator() const noexcept { return m_xs; }
      const VectD& xsBreakpoints() const noexcept { return m_xsBreakpoints; }

      //Sample scatterings:
      ScatterOutcomeIsotropic sampleScatter( RNG&, NeutronEnergy ) const;
      PairDD sampleScatterAlphaBeta( RNG& rng, NeutronEnergy ekin ) const;

      const SABProcessor& processor() const noexcept { return m_p; }
      const SAB::SABExtender& extender() const noexcept { return m_e; }

    private:
      shared_obj<const SABProcessor> m_p;
      shared_obj<const SAB::SABExtender> m_e;
      NeutronEnergy m_emax;
      double m_kT;
      double m_invkT;
      double m_extSAtEmax;
      double m_tableSAtEmax;
      IsoXSEvaluator m_xs;
      VectD m_xsBreakpoints;
      PairDD scatABHighE( RNG&, NeutronEnergy ) const;
    };
  }
}

////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::CrossSect NCrystal::SAB::
SABExtended::crossSectionUnitSigmaBound( NeutronEnergy ekin ) const
{
  return CrossSect{ m_xs.eval( ekin.dbl() ) };
}

inline NCrystal::PairDD NCrystal::SAB::
SABExtended::sampleScatterAlphaBeta( RNG& rng, NeutronEnergy ekin ) const
{
  if ( ekin <= m_emax ) {
    auto ab = m_p->sampleScatterAlphaBeta( rng, ekin );
    return { ab.alpha, ab.beta };
  }
  return scatABHighE( rng, ekin );
}

#endif
