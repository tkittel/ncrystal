#ifndef NCrystal_SABProcessor_hh
#define NCrystal_SABProcessor_hh

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
#include "NCrystal/internal/sab/NCSABCfg.hh"
#include "NCrystal/internal/utils/NCExtraTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SAB {

    class SABProcessor final : public UniqueID {

      // Processes an S(alpha,beta) table to provide cross sections and
      // sampling. Cross sections are the integral of S (as interpolated within
      // each table cell) within the kinematic boundary, to the accuracy set in
      // SABCfg::Cfg. Sampled (alpha,beta) values follow that same S exactly.
      // Both hold up to the energy of getEMaxInfo(), above which SABExtended
      // provides the extension.
      //
      // All expensive work happens in the constructor, at a cost depending on
      // the table size and SABCfg::Cfg. Afterwards, cross sections and
      // sampling are fast at all energies. Results are reproducible: identical
      // across platforms for the same input and RNG stream.

    public:

      ~SABProcessor();

      enum class SampleSupport { YES, NO };
      enum class StoreExtraDiagnostics { YES, NO };

      //The egrid argument follows the convention of DI_ScatKnl::energyGrid()
      //(NCInfoTypes.hh), and can be null:
      SABProcessor( const SABCfg::Cfg& cfg,
                    shared_obj<const SABData>,
                    std::shared_ptr<const VectD> egrid,
                    SampleSupport = SampleSupport::YES,
                    StoreExtraDiagnostics = StoreExtraDiagnostics::NO );

      shared_obj<const SABData> sabDataPtr() const;

      double kT() const;

      //integral of S(alpha,beta) within phasespace of a given neutron energy:
      double phaseSpaceIntegral( NeutronEnergy ) const;

      //Assuming SigmaBound=1barn, get the corresponding cross section:
      CrossSect crossSectionUnitSigmaBound( NeutronEnergy ) const;

      //The same cross section as an IsoXSTable (see
      //NCIsoXSUtils.hh), agreeing to within rounding:
      IsoXSTable createXSTabulation() const;

      struct AlphaBetaOutcome {
        double alpha, beta;
        unsigned ntries;
      };
      AlphaBetaOutcome sampleScatterAlphaBeta( RNG&, NeutronEnergy ) const;

      //Below this E/kT, beta is sampled at that value and mu is isotropic
      //(both exact to ~1e-8, while mu from alpha degrades as eps/sqrt(E/kT)):
      static constexpr double lowEFloorEDivKT = 1e-16;

      bool hasSampleSupport() const;

      void toJSONFullDiagnostics( std::ostream& ) const;

      //Stream information to JSON (appropriate for usage in a Process). If an
      //optional sigma_scale is provided (which would be the elemental cross
      //section in a standard monoatomic material), it will be encoded in the
      //JSON data as well. Likewise, the short description of the method for
      //extending the table to higher energies can be provided as well,
      void toJSONProcessInfo( std::ostream&,
                              Optional<SigmaBound> sigma_scale = NullOpt,
                              Optional<std::string> extension_method
                              = NullOpt  ) const;


      //There is a maximal energy after which the SAB does not reliably cover
      //neutron physics (this is the last point in the internal energy
      //grid). Access this energy and associated information:
      struct EPtInfo final {
        NeutronEnergy ekin;
        EDivKT E_div_kT;
        double phaseSpaceIntegral;
        CrossSect crossSectionUnitSigmaBound;
      };
      EPtInfo getEMaxInfo() const;

      //For reference, access the gridded results:
      const VectD& getEDivKTGrid() const;//[units of kT]
      const VectD& getPhaseSpaceIntegralAtGrid() const;

      //Move-only:
      SABProcessor( const SABProcessor& ) = delete;
      SABProcessor& operator=( const SABProcessor& ) = delete;
      SABProcessor( SABProcessor&& ) noexcept;
      SABProcessor& operator=( SABProcessor&& ) noexcept;

    private:
      void * m_impl = nullptr;
    };

  }
}

#endif
