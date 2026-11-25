#ifndef NCrystal_VDOSDebyeGroup_hh
#define NCrystal_VDOSDebyeGroup_hh

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

#include "NCrystal/internal/sab/NCSABExtended.hh"
#include "NCrystal/internal/vdos/NCVDOSLux.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace DynInfoUtils {

    //Grouped treatment of VDOSDebye @DYNINFO entries: at vdoslux 0-5 (cf.
    //VDOSLux::allowGrouping), all entries with the same mass and temperature
    //form a group, with exact cross sections (sum of VDOSDebyeXS tables) and
    //samplings from a few Debye models. The entries are split into clusters
    //with similar msd values, and each cluster is sampled by a Debye model
    //(phonon expansion) with their mean msd. This keeps initialisation time
    //and memory usage in check for materials with many entries.
    //
    //FIXME: The cross sections of the grouped treatment (VDOSDebyeXS) and of
    //the phonon expansions (vdoslux=6 and 100-106) differ by up to ~2e-3, a
    //limitation of the expansions.
    //Changing vdoslux between 5 and 6 thus causes such jumps. Fixing the
    //expansions would solve this.

    //A VDOSDebye entry with its scale (fraction times bound cross section,
    //possibly relative to that of other entries):
    struct DebyeTempAndScale {
      double scale;
      DebyeTemperature debyeTemperature;
    };

    class VDOSDebyeGroupData final : private MoveOnly {
    public:

      //Data for a group of VDOSDebye entries with the same temperature and
      //mass, with scales normalised to a sum of 1. Cross sections are thus
      //per unit summed scale (in barn for a bound cross section of 1 barn).

      //Tolerance on the spread of log(msd) values within the clusters (see
      //clusterEntries), for a group with neff effective entries:
      //T*max(0,1-3/neff)^2, with T=0.6, 0.4, 0.25, 0.16, 0.1, 0.06
      //at vdoslux 0-5. Groups with few entries thus keep them apart (cheap),
      //while those with many (costly) are clustered more coarsely,
      //continuously in the inputs (a tiny entry changes neff very little).
      //Cross sections are exact, so this only concerns sampled
      //distributions:
      static double tolerance( VDOS::VDOSLux, double neff );

      //Whether all clusters get vdoslux and knllux two levels lower (never
      //below 1, and 0 stays 0): for neff>3, i.e. when grouping is in effect
      //(clustering is about as imprecise as the lower levels). This jumps at
      //neff=3, but by little (distributions at levels L and L-2 differ by
      //~1e-3). Since neff<=n, groups with n<=3 entries thus get neither
      //merging (except of equal msd values) nor lowering, i.e. are sampled
      //like without grouping:
      static bool lowersLevels( double neff );

      //Effective number of entries, (sum of scales)^2/(sum of scales^2):
      using VecEntries = std::vector<DebyeTempAndScale>;
      static double effectiveNumberOfEntries( const VecEntries& );

      //Clusters of entries (as indices into the input) with contiguous
      //log(msd) values (no cuts between equal values): the fewest clusters
      //with E <= tolerance^2 * sum(scale), where E = sum of
      //scale*(log(msd)-cluster mean)^2, and for that number the lowest E
      //(found by dynamic programming in O(k*n^2) for k clusters and n
      //entries). Clusters are sorted by decreasing summed scale (ties: lowest
      //msd first), their entries by increasing msd (ties: input order):
      static std::vector<std::vector<std::size_t>>
      clusterEntries( const VecEntries&, Temperature, AtomMass,
                      double tolerance );

      //A cluster of entries with its Debye model:
      struct Cluster {
        std::size_t n;//number of entries
        double scale;//summed scale of its entries
        DebyeTemperature debyeTemperature;//from their mean msd
        //Summed cross sections of its entries (thinned table, within ~5e-5
        //of the exact sum), absent if there is only one cluster (then equal
        //to xs()):
        Optional<IsoXSEvaluator> xs;
        shared_obj<const SAB::SABExtended> sabext;//(SupportXS::NO)
      };

      VDOSDebyeGroupData( Temperature, AtomMass, std::size_t nentries,
                          PairDD debyeTempRange, IsoXSTable,
                          std::vector<Cluster> );

      Temperature temperature() const noexcept { return m_temp; }
      AtomMass mass() const noexcept { return m_mass; }
      std::size_t nEntries() const noexcept { return m_n; }
      PairDD debyeTempRange() const noexcept { return m_tdrange; }

      //Cross sections of all entries (the sum of those of the clusters):
      const IsoXSEvaluator& xs() const noexcept { return m_xs; }

      //Sorted by decreasing scale:
      const std::vector<Cluster>& clusters() const noexcept { return m_clusters; }

      //Sample with a cluster picked in proportion to its cross section (no
      //random number is used for the picking with a single cluster):
      ScatterOutcomeIsotropic sampleScatter( RNG&, NeutronEnergy ) const;

    private:
      Temperature m_temp;
      AtomMass m_mass;
      std::size_t m_n;
      PairDD m_tdrange;
      IsoXSEvaluator m_xs;
      std::vector<Cluster> m_clusters;
    };

    //Creates (or fetches from a cache, cleared by clearCaches) the data for
    //entries with the given temperature and mass. Scales must be positive,
    //and only their ratios matter. Requires vdoslux.allowGrouping() and
    //knllux!=-1. Throws CalcError where VDOSDebyeXS does not support the mass
    //or temperature. Independent work is done as FactoryJobs:
    shared_obj<const VDOSDebyeGroupData>
    createVDOSDebyeGroupData( std::vector<DebyeTempAndScale>, Temperature,
                              AtomMass, VDOS::VDOSLux, int knllux );
  }
}

#endif
