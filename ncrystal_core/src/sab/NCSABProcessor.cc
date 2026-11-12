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
#include "NCrystal/internal/sab/NCSABSurveyor.hh"
#include "NCrystal/internal/sab/NCSABCellInteg.hh"
#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/sab/NCSABBetaStripSampler.hh"
#include "NCrystal/internal/sab/NCSABIdx.hh"
#include "NCrystal/internal/utils/NCMixedDataVector.hh"
#include "NCrystal/internal/utils/NCTinyVector.hh"
#include "NCrystal/internal/utils/NCFastSearch.hh"
#include "NCrystal/internal/utils/NCSpline.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SAB;

namespace NCRYSTAL_NAMESPACE {

  namespace SAB {

    namespace {

#ifndef NDEBUG
      static bool dummy_testsabidx = [](){
        namespace SI = SABIdx;
        nc_assert(SI::PackedIndex::fromAlphaIdxBetaIdx(2,6).unpackAlphaIdx()==2);
        nc_assert(SI::PackedIndex::fromAlphaIdxBetaIdx(2,6).unpackBetaIdx()==6);
        nc_assert(SI::PackedIndex::fromAlphaIdxBetaIdx(0,0).unpackAlphaIdx()==0);
        nc_assert(SI::PackedIndex::fromAlphaIdxBetaIdx(0,0).unpackBetaIdx()==0);
        nc_assert(SI::NAlpha{4}.value() == 4);
        nc_assert(SI::NAlphaCells(SI::NAlpha{4}).value() == 3);
        nc_assert(SI::NAlphaCells{4}.value() == 4);
        nc_assert(SI::NAlpha(SI::NAlphaCells{4}).value() == 5);
        return true;
      }();
#endif
      class CellMgr final : private NoCopyMove {

        void initCellIntegrals()
        {
          auto& sab = m_sab->sab();
          auto& alpha = m_sab->alphaGrid();
          auto& beta = m_sab->betaGrid();
          //Now calculate cell integrals:
          // const auto nalpha = alpha.size();
          const auto nbeta = beta.size();
          SABIdx::NAlpha nalpha(m_sab->alphaGrid());
          SABIdx::NAlphaCells nalphacells(nalpha);
          nc_assert_always(nbeta*nalpha.value()==sab.size());
          nc_assert_always(nbeta>=2&&nalpha.value()>=2);
          nc_assert( nalphacells.value()+1==alpha.size() );
          const auto nbetacells = nbeta-1;
          const auto ncells = nalphacells.value()*nbetacells;
          m_cellIntegral.reserve(ncells);
          VectD alphaint_cache;
          alphaint_cache.resize(nalphacells.value(),0.0);
          auto itS = sab.begin();
          const double * itLogS = m_logS->data();
          //NB: S(alpha[ia],beta[ib]) = sab[ ib*nalpha+ia ]
          {
            //Init alphaint_cache with cell edge integrals at ib=0:
            auto itA = alpha.begin();
            for ( auto& e : alphaint_cache ) {
              const double a_low = *itA++;
              const double s_low = *itS++;
              const double logS_low = *itLogS++;
              e = integrateLogLinOrLinWCL(a_low,s_low,*itA,*itS,
                                          logS_low, *itLogS);
            }
          }

          //fixme: check all _always in file
          for ( auto ib : ncrange(nbetacells) ) {
            //skip to beginning of beta slice:
            ++itS;
            ++itLogS;
            nc_assert( (std::size_t)std::distance(sab.begin(),itS) == (ib+1)*nalpha.value() );
            nc_assert( (std::size_t)std::distance(static_cast<const double*>(m_logS->data()),itLogS) == (ib+1)*nalpha.value() );
            nc_assert( itS == std::next(sab.begin(),(ib+1)*nalpha.value()) );
            nc_assert( itLogS == std::next(m_logS->data(),(ib+1)*nalpha.value()) );
            const double b1 = vectAt(beta,ib);//fixme: use iterators instead
            const double b2 = vectAt(beta,ib+1);
            const double half_db = 0.5 * (b2-b1);
            auto itAIntCache = alphaint_cache.begin();
            auto itAIntCacheE = alphaint_cache.end();
            auto itA = alpha.begin();
            nc_assert(alphaint_cache.size()+1==alpha.size());
            for ( ; itAIntCache!=itAIntCacheE; ++itAIntCache ) {
              nc_assert(itA<alpha.end());
#ifndef NDEBUG
              //Verify consistency of indexing:
              auto dbg_ialpha = static_cast<std::size_t>( itA-alpha.begin() );
#endif
              const double a_low = *itA++;
              const double s_low = *itS++;
              const double logS_low = *itLogS++;
              // nc_assert_always( logS_low == *itLogS++);
              double aint = integrateLogLinOrLinWCL(a_low,s_low,*itA,*itS,
                                                    logS_low, *itLogS);
              //combine linearly with previous edge integration to get full cell
              //integral:
              nc_assert( m_cellIntegral.size()
                         == SABIdx::SABCellIdx( nalphacells,
                                                dbg_ialpha, ib ).value() );
              m_cellIntegral.push_back( half_db*(aint+*itAIntCache) );
              //store this edge integration for use during next beta line:
              *itAIntCache = aint;
            }
          }
          nc_assert(m_cellIntegral.size()==ncells);
          nc_assert(std::next(itS)==sab.end());
          nc_assert(std::next(itLogS)==m_logS->data()+sab.size());
        }

      public:
        using cellidx_t = SABIdx::PackedIndex;
        const SABData& sabData() const { return *m_sab; }
        shared_obj<const SABData> sabDataPtr() const { return m_sab; }

        CellMgr( shared_obj<const SABData> sab )
          : m_sab(std::move(sab)),
            m_logS( NCS::createLogSTable( *m_sab ) )
        {
#ifndef NDEBUG
          (void)dummy_testsabidx;
#endif
          m_nalpham1 = m_sab->alphaGrid().size()-1;
          initCellIntegrals();
        }

        //Table of log(S) values, also shared with SABBetaStripSampler:
        const shared_obj<const VectD>& logSTable() const { return m_logS; }

        CellData lookupCellInfo( cellidx_t packed_idx ) const
        {
          auto ia = packed_idx.unpackAlphaIdx();
          auto ib = packed_idx.unpackBetaIdx();
          auto& sab = m_sab->sab();
          auto& alpha = m_sab->alphaGrid();
          auto& beta = m_sab->betaGrid();
          CellData res;
          res.a1 = vectAt(alpha,ia);
          res.a2 = vectAt(alpha,ia+1);
          res.b1 = vectAt(beta,ib);
          res.b2 = vectAt(beta,ib+1);
          SABIdx::NAlpha na(m_nalpham1+1);
          SABIdx::SABIdx idx_sab( na, ia, ib );
          //const auto na = alpha.size();
          //auto idx_sab = ib*na + ia;
          nc_assert( idx_sab.value() == ib*na.value() + ia );
          auto itS = std::next(sab.begin(),idx_sab.value());
          const double * itLogS = m_logS->data() + idx_sab.value();
          nc_assert(std::next(itS)<sab.end());
          nc_assert(std::next(itLogS)<m_logS->data()+sab.size());
          res.S[0] = *itS++;
          res.S[1] = *itS;
          res.logS[0] = *itLogS++;
          res.logS[1] = *itLogS;
          //Advance one row, but we already advanced one item, so only add
          //nalpha-1:
          itS += m_nalpham1;
          itLogS += m_nalpham1;
          nc_assert(std::next(itS)<sab.end());
          nc_assert(std::next(itLogS)<m_logS->data()+sab.size());
          res.S[2] = *itS++;
          res.S[3] = *itS;
          res.logS[2] = *itLogS++;
          res.logS[3] = *itLogS;
          return res;
        }

        double getCellIntegral( SABIdx::PackedIndex p_idx ) const
        {
          SABIdx::SABCellIdx sab_idx( SABIdx::NAlphaCells{m_nalpham1}, p_idx );
          return vectAt(m_cellIntegral, sab_idx.value() );
        }

      private:
        shared_obj<const SABData> m_sab;
        shared_obj<const VectD> m_logS;
        VectD m_cellIntegral;
        std::size_t m_nalpham1;
      };

      //Copy of the data with S scaled by 2^k (exact), using the smallest k>=0
      //which brings the smallest positive S up to at least 2^-478, but keeping
      //the largest S at most 2^478 (both ~1e10 inside sqrt(DBL_MIN) and
      //sqrt(DBL_MAX), leaving room for products with cell areas, sums, etc.).
      //S values still subnormal after this are set to 0, since log/exp of
      //subnormal values are unreliable (e.g. with MSVC or FTZ/DAZ modes, exp
      //can return 0 instead). Results are linear in S, so the scaling is exactly
      //undone on them with invScale=2^-k. As k depends on the data, no absolute
      //thresholds must be applied to S values or integrals of the copy:
      shared_obj<const SABData> prepareSForProcessor( const shared_obj<const SABData>& sd,
                                                      double& invScale )
      {
        constexpr int kkk = 478;
        constexpr double dblmin = std::numeric_limits<double>::min();
        double smin = kInfinity;
        double smax = 0.0;
        for ( auto v : sd->sab() ) {
          if ( v > 0.0 ) {
            smin = ncmin( smin, v );
            smax = ncmax( smax, v );
          }
        }
        int k = 0;
        if ( smax > 0.0 ) {
          int emin, emax;
          std::frexp( smin, &emin );//smin >= 2^(emin-1)
          std::frexp( smax, &emax );//smax < 2^emax
          k = ncmax( 0, ncmin( 1 - kkk - emin, kkk - emax ) );
        }
        invScale = std::ldexp( 1.0, -k );
        if ( k == 0 && !( smin < dblmin ) )
          return sd;
        VectD sab( sd->sab() );
        for ( auto& v : sab ) {
          v = std::ldexp( v, k );
          if ( v < dblmin )
            v = 0.0;
        }
        return makeSO<const SABData>( VectD( sd->alphaGrid() ),
                                      VectD( sd->betaGrid() ),
                                      std::move( sab ),
                                      sd->temperature(), sd->boundXS(),
                                      sd->elementMassAMU(),
                                      sd->suggestedEmax() );
      }

      struct SABProcImpl final : private NoCopyMove {
        using SampleSupport = SABProcessor::SampleSupport;
        using StoreExtraDiagnostics = SABProcessor::StoreExtraDiagnostics;
        SABProcImpl( const SABCfg::Cfg&,
                     shared_obj<const SABData>,
                     std::shared_ptr<const VectD> egrid,
                     SampleSupport,
                     StoreExtraDiagnostics );

        double phaseSpaceIntegral( NeutronEnergy ) const;
        using AlphaBetaOutcome = SABProcessor::AlphaBetaOutcome;
        AlphaBetaOutcome sampleAlphaBeta( RNG&,
                                          NeutronEnergy ) const;

        void toJSONFullDiagnostics( std::ostream& ) const;

        double m_kT;
        double m_invkT;
        //constant needed for low E extrapolation:
        double m_lowEExtrapolationConstant = 0.0;

        //Sorted energy grid (E/kT):
        VectD m_eGrid;
        LogGuideTable m_eGridGuide;//for fast searches in m_eGrid
        //Equivalent phase-spaced bounded integral of S for each of the energy
        //points:
        VectD m_sIntegral;

        //For smooth sIntegrals, we need:
        Optional<HymanInterp> m_sIntegralSmooth;

        //Extra diagnostics data
        VectD m_diagnostics_fcAR; //Full cell AR at at each m_eGrid pt
        //Expected tries of a strip sampler with emax at each m_eGrid pt (-1
        //if none):
        VectD m_diagnostics_stripTries;

        //Number of touched cells at each energy grid point, for full-cell
        //sampling (empty without sampling support):
        std::vector<std::uint32_t> m_nTouchedCells;
        bool canSample() const { return !m_nTouchedCells.empty(); }
        ////////////////////////////////////////////
        //Data needed for the full-cell sampling:

        //The cumulative full cell integrals, ordered in the order in which they
        //are touched by neutron phasespace as energy increases.
        RandIdxByWeightPicker m_cumulFCInt;
        //The corresponding packed cell indices:
        std::vector<SABIdx::PackedIndex> m_cumulFCInt_cellidx;

        //The original data, and the factor undoing the scaling of the S values
        //in m_cellmgr (see prepareSForProcessor) on results linear in S:
        shared_obj<const SABData> m_origSABData;
        double m_invSScale = 1.0;
        CellMgr m_cellmgr;
        //Phase-space integrals at the grid points, with the scaling undone:
        VectD m_sIntegralUnscaled;

        //Beta-strip samplers with increasing emax, and the one to use at each
        //energy grid point (-1 for none):
        std::vector<SABBetaStripSampler> m_strips;
        std::vector<std::int32_t> m_stripIdx;
        void initStripLadder( const VectD& fcAR );

        static constexpr double lowEFloorEDivKT
          = SABProcessor::lowEFloorEDivKT;
        //Strip sampler with emax at the lowest grid point, used below the
        //energy grid:
        Optional<SABBetaStripSampler> m_stripBelowGrid;
        void initStripBelowGrid();
        using ABOutcome = SABProcessor::AlphaBetaOutcome;
        ABOutcome sampleBelowFloor( RNG&, EDivKT ) const;
      private:

        //Sample a scattering event:
        struct ABRes { double a, b; std::size_t ntries; };
        ABRes sampleAlphaBetaFC( RNG&, EDivKT,
                                 std::size_t nCellsTouched ) const;

      };

      struct sIntegralAtE_Result {
        double integralWithinKB;
        double integralTouchedCells;

        //Number of cells touched (they will all be located at the beginning of
        //SABSurveyor::data()):
        std::size_t nTouchedCells;

        //Internal cache, speeding up repeated calls if they are done from low
        //to high energy values:
        struct PrevPtInfo {
          double E_div_kT = -1.0;//-1.0 means "n/a"
          StableSumKahan sumFullyCoveredCells;
          std::size_t nLastTouched = 0;
          //Cells touched but not fully covered at E_div_kT (in SABSurveyor
          //order), so later calls need not loop over all touched cells:
          std::vector<SABSurveyor::CellInfo> active;
        } prevCache;

      };

      void sIntegralAtE( EDivKT e,
                         const CellMgr& mgr,
                         const SABSurveyor& surv,
                         StdLogLinCellIntegrator::IntegrationScheme scheme,
                         sIntegralAtE_Result& result )
      {
        auto& prev = result.prevCache;
        const double prev_E = ( ( prev.E_div_kT >= 0.0
                                  && prev.E_div_kT <= e.dbl() )
                                ? prev.E_div_kT : -1.0 );


        nc_assert(e.dbl()>0&&std::isfinite(e.dbl()));
        //Find all touched cells:
        auto survCells = surv.data();
        //std::upper bound returns the first element where the comparision
        //returns true. The list must be sorted so the comparison will return
        //[false,false,...,false,true,true,...,true].

        //we find the first element coming AFTER e

        //FIXME: Why do we need the binary search? Why not just proceed over all
        //cells until e < ci.e_touch?

        auto itSearchBegin = survCells.begin();
        if ( prev_E >= 0.0 && prev.nLastTouched > 2 )
          itSearchBegin += prev.nLastTouched-1;//fixme: is this actually helping in a measurable way??

        auto itLastTouchedE
          = std::upper_bound(itSearchBegin,
                             survCells.end(),
                             e.dbl(),
                             [](double val,
                                const SABSurveyor::CellInfo& ci)
                             { return val < ci.e_touch; } );
        //upper_bound returns one past the last touched:
        if( itLastTouchedE == survCells.begin() )
          NCRYSTAL_THROW2(BadInput,"Phase-space of neutron energy E/kT="
                          <<fmt(e.dbl())
                          <<" does not touch any S(alpha,beta) cells!");
        StableSumKahan sum;
        StableSumKahan sumFullCells;
        StableSumKahan sumCoveredCells;
        result.nTouchedCells = std::distance(survCells.begin(),itLastTouchedE);
        //cells with relative max contrib below this can never contribute:
        // (fixme: sablux should modify 1e10 in range [1e6,1e17])
        const double threshold = ncclamp( 1.0/(1e10*result.nTouchedCells),
                                          1e-30,1e-6 );

        //Constants for the tapered shortcut below:
        constexpr double taperBand = 10.0;
        constexpr double invTaperBand = 1.0/taperBand;
        constexpr double logTaperBand = 2.3025850929940456840179914546843642;
        nc_assert( floateq( std::log(taperBand), logTaperBand ) );
        constexpr double twoLogTaperBand = 2.0*logTaperBand;

        //Cells fully covered at prev_E are already included in .prevCache, so
        //only visit the active ones and those newly touched (keeping the
        //SABSurveyor order, which fixes the summation order):
        auto& active = prev.active;
        if ( prev_E < 0.0 )
          active.clear();
        active.insert( active.end(),
                       survCells.begin()
                       + ( prev_E >= 0.0 ? prev.nLastTouched : 0 ),
                       itLastTouchedE );
        //Cells fully covered at e (incl. those already covered at prev_E):
        for ( auto& c : active )
          if ( e.dbl() >= c.e_cover )
            sumCoveredCells.add( mgr.getCellIntegral( c.cellidx ) );
        if ( prev_E >= 0.0 )
          sumCoveredCells.add( prev.sumFullyCoveredCells );
        //Cutoff for negligible crossed cells, from quantities at e only (so
        //the result does not depend on the previously evaluated energy):
        const double cutoffref = sumCoveredCells.sum() * threshold;
        for ( auto it = active.begin(); it != active.end(); ++it ) {
          nc_assert( e.dbl() >= it->e_touch );
          nc_assert( prev_E < it->e_cover );
          if ( e.dbl() >= it->e_cover )
            continue;//fully covered (summed above)
          const double fullCellIntegral = mgr.getCellIntegral( it->cellidx );
          double contrib;
          sumFullCells.add(fullCellIntegral);
          //We want to apply a "don't waste time on irrelevant cell" short
          //cut, but it becomes too sensitive to FFT-convolution noise in the
          //input SAB if we use a hard threshold. So we instead use a smooth
          //tapering.
          if ( fullCellIntegral < cutoffref*invTaperBand ) {
            contrib = 0.0;//safely below threshold, don't waste time
          } else {
            CellData cellData = mgr.lookupCellInfo( it->cellidx );
            StableSumKahan crossedRes;
            StdLogLinCellIntegrator::integrateWithinKB( cellData, e,
                                                        scheme, crossedRes );
            contrib = crossedRes.sum();
            if ( fullCellIntegral < cutoffref*taperBand ) {
              //fixme replace /twoLogTaperBand with *invTwoLogTaperBand
              const double t = ncclamp( ( std::log(fullCellIntegral/cutoffref)
                                          + logTaperBand )
                                        / twoLogTaperBand, 0.0, 1.0 );
              //Quintic smootherstep (Ken Perlin): 0 and 1 derivatives
              //vanish at both ends, so no kink at the band edges either.
              //FIXME: Code duplicated with NCVDOSUtils.cc, need common NCMath
              //util fct.
              //Explicit ncfma for the inner Horner steps:
              const double s = t*t*t*ncfma( t, ncfma(t,6.0,-15.0), 10.0 );
              contrib *= s;
            }
            sum.add( contrib );
          }
        }

        //Prepare cache for next:
        active.erase( std::remove_if( active.begin(), active.end(),
                                      [e]( const SABSurveyor::CellInfo& c )
                                      { return c.e_cover <= e.dbl(); } ),
                      active.end() );
        prev.sumFullyCoveredCells.set(sumCoveredCells);
        prev.E_div_kT = e.dbl();
        prev.nLastTouched = result.nTouchedCells;
        sum.add(sumCoveredCells);
        sumFullCells.add(sumCoveredCells);


        result.integralWithinKB = sum.sum();
        result.integralTouchedCells = sumFullCells.sum();
        nc_assert(std::isfinite(result.integralWithinKB));
        nc_assert(std::isfinite(result.integralTouchedCells));
      }


      //Automatically determine EMax by looking at the full cell integrals of
      //cells touched at a given energy:
      NeutronEnergy determineEMax( Temperature temp,
                                   //AtomMass,//fixme: add if needed (for FG model perhaps?)
                                   const VectD& cumulFullCellIntegral,
                                   const VectD& cellTouch_E_div_kT )
      {
        const double kT = temp.kT();
        const VectD& e = cellTouch_E_div_kT;
        const VectD& s = cumulFullCellIntegral;
        const double extra_safety = 0.98;

        //First, discard trailing entries where the integral(S) remains
        //~constant, since that is a sign that increased phase-space does not
        //add new contributions to the integral. Thus, we discard until we have
        //removed eps_emax (in logIntegral space) of the integral(S) difference:

        nc_assert_always(s.size()==e.size());
        nc_assert_always(s.size()>=16);
        const double s1(s.front());
        const double s2(s.back());
        nc_assert_always( s2>s1 && s1>=0.0 && std::isfinite(s2) );
        double tgt_s;
        constexpr double eps_emax = 1e-2;
        if ( s1==0.0 ) {
          tgt_s = s2*(1.0-eps_emax);
          } else {
          tgt_s = s1 * std::pow( s2/s1, 1.0-eps_emax );
        }
        std::size_t i = s.size()-1;
        while ( i && vectAt(s,i)>tgt_s )
          --i;
        nc_assert_always( i > 4 );//fixme nicer error

        auto max_E_div_kT = vectAt(e,i);
        nc_assert( max_E_div_kT > 0.0 );
        return NeutronEnergy{ ncmin( max_E_div_kT*kT*extra_safety,
                                     1000.0,
                                     10000*kT ) };
      }


    }//anon namespace

    namespace detail {
      //(NCRYSTAL_API and not in the anonymous namespace, for unit tests)

      //The E->0 limit of Sint(E)/sqrt(E), f0=4*int_0^inf sqrt(b)S(b,b)db, and
      //the mean of sqrt(alpha) along alpha=beta, weighted by that integrand.
      //logS is the createLogSTable table of the SABData.
      NCRYSTAL_API PairDD diagonalLowEnergyLimit( const SABData&,
                                                  const VectD& logS );

      //Find Emin/kT, the lowest point of the energy grid. Below Emin, the cross
      //section is extrapolated as 1/sqrt(E), so Emin is the lowest E where
      //f(E)=Sint(E)/sqrt(E) deviates from its E->0 value (by 1e-3 of its total
      //change up to Emax). Scans upwards in steps of 2 from the lowest E where
      //Sint is numerically reliable, so dips and peaks are not skipped.
      //diaglimit comes from diagonalLowEnergyLimit, f_of_e takes E/kT.
      NCRYSTAL_API double determineEMinDivKT( const SABCfg::Cfg&, double emax,
                                              const std::function<double(double)>& f_of_e,
                                              const PairDD& diaglimit );
    }
  }
}

NC::PairDD NCS::detail::diagonalLowEnergyLimit( const SABData& sab,
                                                const VectD& logS )
{
  //Integrates c(x)=4*sqrt(x)*S(x,x) and c(x)*sqrt(x) along the diagonal, per
  //segment within one grid cell. There S(x,x) blends two log-linear rows, so
  //it is smooth apart from its log-slope: segments are split into m equal
  //pieces, each with a log change of S of at most 0.5, where Gauss-Legendre-8
  //is accurate to ~1e-14 (3e-14 for a piece starting at x=0). Each piece is
  //integrated in t=sqrt(x), as dx=2t*dt removes the sqrt(x) singularity.
  const auto& A = sab.alphaGrid();
  const auto& B = sab.betaGrid();
  const auto& S = sab.sab();
  const std::size_t na = A.size();
  //Log change of S along alpha over [x0,x1] in cell (ia,ib), in rows ib and
  //ib+1 (none if a row has a zero corner, being linear there):
  auto logChange = [&]( std::size_t ia, std::size_t ib, double x0, double x1 )
  {
    double q = 0.0;
    for ( auto jb : { ib, ib + 1 } ) {
      const std::size_t i = SAB::calcSABIdx( na, ia, jb );
      if ( vectAt( S, i ) > 0.0 && vectAt( S, i + 1 ) > 0.0 )
        q = ncmax( q, ncabs( vectAt( logS, i + 1 ) - vectAt( logS, i ) ) );
    }
    return q * ( x1 - x0 ) / ( vectAt( A, ia + 1 ) - vectAt( A, ia ) );
  };
  StableSum tot, totsq;
  for ( auto ib : ncrange( B.size() - 1 ) ) {
    const double lo = ncmax( vectAt(B,ib), 0.0, A.front() );
    const double hi = ncmin( vectAt(B,ib+1), A.back() );
    if ( !( hi > lo ) )
      continue;
    //Split at the alpha grid points within [lo,hi]:
    std::size_t ia = static_cast<std::size_t>
      ( std::upper_bound( A.begin(), A.end(), lo ) - A.begin() ) - 1;
    double x0 = lo;
    while ( x0 < hi && ia + 1 < na ) {
      const double x1 = ncmin( hi, vectAt(A,ia+1) );
      constexpr double kkk = 1.0 / 0.5;
      const auto m = static_cast<unsigned>
        ( std::ceil( logChange( ia, ib, x0, x1 ) * kkk ) );
      const unsigned npieces = ( m > 1 ? m : 1 );
      double y0 = x0;
      for ( auto k : ncrange( 1u, npieces + 1 ) ) {
        const double y1 = ( k == npieces ? x1
                            : x0 + ( x1 - x0 ) * k / npieces );
        //S(x,x) at x=t^2 (clamped, as t^2 can round outside the piece):
        auto sAt = [&]( double t )
        {
          const double x = ncclamp( t * t, y0, y1 );
          return SAB::evalSABCellWCL( sab, logS, ia, ib, x, x );
        };
        //c(x)dx = 8*t^2*S dt, and c(x)*sqrt(x)dx = 8*t^3*S dt:
        const double t0 = std::sqrt( y0 );
        const double t1 = std::sqrt( y1 );
        tot.add( integrateGaussLegendre8( [&sAt]( double t )
                                          { return 8.0 * t * t * sAt( t ); },
                                          t0, t1 ) );
        totsq.add( integrateGaussLegendre8( [&sAt]( double t )
                                            { return 8.0*t*t*t * sAt( t ); },
                                            t0, t1 ) );
        y0 = y1;
      }
      x0 = x1;
      ++ia;
    }
  }
  const double f0 = tot.sum();
  return { f0, f0 > 0.0 ? totsq.sum() / f0 : 0.0 };
}

double NCS::detail::determineEMinDivKT( const SABCfg::Cfg& cfg, double emax,
                                        const std::function<double(double)>& f_of_e,
                                        const PairDD& diaglimit )
{
#ifndef NDEBUG
  nc_assert_always( emax > 0.0 && std::isfinite(emax) );
#endif
  const double f0 = diaglimit.first;
  const double logf2 = std::log( f_of_e( emax ) );
  if ( !( f0 > 0.0 && std::isfinite(f0) && std::isfinite(logf2) ) )
    return emax;
  //Rounding noise relative to the strip width at alpha is
  //eps*sqrt(alpha)/(4*sqrt(e)); keep it below 1e-3 of the target (f0
  //only provides the weights and a rough target here):
  constexpr double eps = std::numeric_limits<double>::epsilon();
  const double target0 = 1e-3 * ncabs( logf2 - std::log(f0) );
  if ( !( target0 > 0.0 ) )
    return emax;
  const double e1 = ncsquare( eps * diaglimit.second
                              / ( 4.0 * 1e-3 * target0 ) );
  if ( !( e1 < emax ) )
    return emax;
  const double loge1 = std::log( e1 );
  const double loge2 = std::log( emax );
  const double logref = std::log( f_of_e( e1 ) );
  const double target = 1e-3 * ncabs( logf2 - logref );
  if ( !( target > 0.0 && std::isfinite(target) ) )
    return emax;
  auto froot = [&f_of_e,logref,target,loge2,logf2]( double loge )
  {
    const double v = ( loge < loge2
                       ? std::log( f_of_e( std::exp(loge) ) )
                       : logf2 );
    return ncabs( v - logref ) - target;
  };
  constexpr double logStep = 0.69314718055994530942;//log(2)
  nc_assert( floateq( logStep, std::log(2.0) ) );
  double bracket_lo = loge1;
  double bracket_hi = loge2;
  for ( double cand = loge1 + logStep; cand < loge2; cand += logStep ) {
    if ( froot( cand ) >= 0.0 ) {
      if ( cand == loge1 + logStep )
        NCRYSTAL_THROW2(CalcError,"No low-energy regime found in the"
                        " S(alpha,beta) integral (it already changes at"
                        " E/kT="<<fmtg(std::exp(cand))<<")");
      bracket_hi = cand;
      break;
    }
    bracket_lo = cand;
  }
  const double logemin = findRoot2( froot, bracket_lo, bracket_hi,
                                    cfg.egrid_emin_accuracy );
  return ncclamp( std::exp( logemin ), e1, emax );
}

namespace NCRYSTAL_NAMESPACE {
  namespace SAB {
    namespace {
      struct SIntByETouch final : private MoveOnly {
        VectD sint, e;
      };

      SIntByETouch determineSIntegralByTouchedCells( const CellMgr& cellmgr,
                                                     const SABSurveyor& surv,
                                                     double kT )
      {
        auto it = surv.data().begin();
        auto itE = surv.data().end();
        nc_assert_always( it != itE );

        SIntByETouch res;
        res.e.reserve(4096);
        res.sint.reserve(4096);
        double eprev(0.0);
        double curr_thr(0.0);
        const double thr_eup = 1e3/kT;//1keV/kT is essentially ~infinity for our
                                      //purposes
        StableSumKahan sum;

        nc_assert_always(it!=itE);
        for ( ;it!=itE;++it ) {
          double ci = cellmgr.getCellIntegral(it->cellidx);
          if ( it->e_touch > thr_eup )
            break;
          if ( !(ci>curr_thr) )
            continue;
          sum.add(ci);
          //fixme: tune?
          curr_thr = sum.sum()*1e-26;//less than this won't ever affect result (~16
          //digits + O(1e9) cells => ~O(1e-25) is enough). (fixme: tighten?)
          if ( it->e_touch > eprev*1.001 ) {
            eprev = it->e_touch;
            res.e.push_back(it->e_touch);
            res.sint.push_back(sum.sum());
          }
        }

        return res;
      }

      //Energy grid (E/kT, while requested_egrid is in eV): either a fixed grid,
      //or (grid empty) the range [tlo,thi] of an adaptive grid:
      struct EGridSpec {
        VectD grid;
        double tlo = 0.0, thi = 0.0;
      };

      EGridSpec determineEGrid( const SABCfg::Cfg& cfg,
                            const SABData& sab,
                            const CellMgr& cellmgr,
                            const SABSurveyor& surv,
                            std::shared_ptr<const VectD> requested_egrid )
      {
        const auto scheme = cfg.integSchemeDetermineEGrid;
        //Number of points of a fixed grid requested in the data (otherwise
        //the grid is adaptive):
        std::size_t npts = 0;
        bool npts_requested = false;

        //fixme: major cleanup needed
        Optional<NeutronEnergy> suggestedEMax_by_egrid;
        Optional<NeutronEnergy> suggestedEMin_by_egrid;

        if ( requested_egrid && !requested_egrid->empty() ) {
          const auto& re = *requested_egrid;
          //fixme: more sanity checks
          const double invkT = 1.0/sab.temperature().kT();
          if ( re.size() == 1 ) {
            suggestedEMax_by_egrid = NeutronEnergy{ re.front() };
          } else if ( re.size() == 3 ) {
            //fixme: the NCMAT spec says that this should give points
            //distributed "evenly in logarithmic space". Perhaps we can solve
            //this merely by updating the language in the spec, to be a bit more
            //loose?

            // { emin or 0, emax or 0, npts or 0 }
            if ( re.at(0) != 0.0 )
              suggestedEMin_by_egrid = NeutronEnergy{ re.at(0) };
            if ( re.at(1) != 0.0 )
              suggestedEMax_by_egrid = NeutronEnergy{ re.at(1) };
            if ( re.at(2) != 0.0 ) {
              //Range check before the cast:
              const double v = re.at(2);
              if ( !( v >= 10.0 && std::floor(v) == v
                      && v < std::ldexp( 1.0, std::numeric_limits<std::size_t>
                                         ::digits ) ) )
                NCRYSTAL_THROW2(BadInput,"Invalid egrid npts (must be"
                                " integral value, at least 10): "<<re.at(2));
              npts = static_cast<std::size_t>(v);
              npts_requested = true;
            }
          } else {
            //complete egrid given directly by request, so in this case we
            //simply honour it and return it directly (we can only assume that
            //the users know what they are doing):
            if ( !( re.size() >= 10 ) )
              NCRYSTAL_THROW2(BadInput,"Invalid requested egrid (must have"
                              " length 1, 3, or >=10)");
            EGridSpec res;
            res.grid.reserve( re.size() );
            for ( auto e : re )
              res.grid.push_back( e * invkT );
            return res;
          }
        }

        {
          nc_assert_always( !npts_requested || npts >= 10 );
          if ( suggestedEMax_by_egrid.has_value() )
            suggestedEMax_by_egrid.value().validate();
          if ( suggestedEMin_by_egrid.has_value() )
            suggestedEMin_by_egrid.value().validate();
          if ( suggestedEMax_by_egrid.has_value()
               && suggestedEMin_by_egrid.has_value()
               && !( suggestedEMax_by_egrid.value()
                     > suggestedEMin_by_egrid.value() ) ) {
            NCRYSTAL_THROW2(BadInput,"egrid: emax ("
                            <<suggestedEMax_by_egrid.value()
                            <<") not greater than emin ("
                            <<suggestedEMin_by_egrid.value()<<")");
          }
        }

        auto sIntByEtouch
          = determineSIntegralByTouchedCells( cellmgr, surv,
                                              sab.temperature().kT() );
        const double kT = sab.temperature().kT();

        Optional<NeutronEnergy> suggestedEMin = suggestedEMin_by_egrid;
        Optional<NeutronEnergy> suggestedEMax;
        if ( sab.suggestedEmax() > 0.0 )
          suggestedEMax = NeutronEnergy{ sab.suggestedEmax() };
        if ( suggestedEMax_by_egrid.has_value() )
          suggestedEMax = suggestedEMax_by_egrid;

        std::pair<NeutronEnergy,NeutronEnergy> E_min_max;
        E_min_max.second = ( suggestedEMax.has_value()
                             ? suggestedEMax.value()
                             : determineEMax( sab.temperature(), sIntByEtouch.sint, sIntByEtouch.e ) );

        if ( suggestedEMin.has_value() ) {
          E_min_max.first = suggestedEMin.value();
        } else {
          sIntegralAtE_Result tmp_result;
          auto f = [&cellmgr, &surv, &tmp_result,scheme]( double E_div_kT ) {
            nc_assert(E_div_kT>0.0);
            const double se = std::sqrt(E_div_kT);
            nc_assert(se>0.0);
            sIntegralAtE( EDivKT{E_div_kT}, cellmgr, surv, scheme,
                          tmp_result );
            return tmp_result.integralWithinKB / se;
          };
          const auto diaglimit
            = detail::diagonalLowEnergyLimit( sab, *cellmgr.logSTable() );
          E_min_max.first
            = NeutronEnergy{ kT * detail::determineEMinDivKT( cfg,
                                                      E_min_max.second.get()/kT,
                                                      f, diaglimit ) };

        }
        if ( E_min_max.first >= E_min_max.second ) {
          NCRYSTAL_WARN("Lowering Emin to 10% of Emax");
          E_min_max.first = NeutronEnergy( 0.1 * E_min_max.second.get() );
        }
        if ( E_min_max.first.get() < 1e-50*E_min_max.second.get() ) {
          NCRYSTAL_WARN("Raising Emin to 1e-50*Emax");
          E_min_max.first = NeutronEnergy( 1e-50 * E_min_max.second.get() );
        }

        constexpr double lowE_factor = 1e-4;
        const double invkT = 1.0/kT;
        if ( !npts_requested ) {
          //Adaptive grid over the same range as the fixed grid below:
          EGridSpec res;
          res.tlo = E_min_max.first.get() * lowE_factor * invkT;
          res.thi = E_min_max.second.get() * invkT;
          return res;
        }

        VectD final_egrid;
        {
          //at lowE int(S)(E) = const*sqrt(E).
          //
          //So to have an acceptance rate of 0.9 between E pts, we need E points to be
          //spaced with a factor of k, given by:
          //
          //  const*sqrt(E) = 0.9*const*sqrt(E*k)
          //  <=> sqrt(k) = 1/0.9 <=> k = 1/0.9**2
          //
          // Meaning that N pts will cover a factor k^N values. So to cover En/E0=1e4,
          // one will need:
          //
          // k^N = 1e4 <=> N = ln(1e4)/lnk = ln(1e4)/ln(1/0.9**2) ~= 43.7 ~= 44
          //
          // If on the other hand we wanted to cover En/E0=1e6 with AR>=0.95, we would
          // need ~134 pts. This would also be OK, but note that the quality of
          // sampled results is the same.
          constexpr std::size_t n_lowe = 44;//fixme: we should consider a
                                            //lin-in-sqrt(E) interpolation in
                                            //the region of these cells, and a
                                            //lin-in-E above. Since we have
                                            //already determined a "lowE" cutoff
                                            //this seems sensible. For
                                            //user-provided egrid, we can for
                                            //simplicity just use lin-in-sqrt(E)
                                            //for the lower 5% of epts (or we
                                            //can detect)
          if ( npts < 80 )
            npts = 80;
          //+1 since the point at Emin is shared by both ranges:
          const std::size_t n = npts + 1 - n_lowe;
          final_egrid.reserve(n+n_lowe);
          for ( auto ee : geomspace(E_min_max.first.get()*lowE_factor,
                                    E_min_max.first.get(), n_lowe) )
            final_egrid.push_back(ee);
          final_egrid.pop_back();
          for ( auto ee : geomspace(E_min_max.first.get(),
                                    E_min_max.second.get(), n) )
            final_egrid.push_back(ee);
          nc_assert_always( final_egrid.size() == npts );
        }

        for ( auto& ee : final_egrid )
          ee *= invkT;

        EGridSpec res;
        res.grid = std::move( final_egrid );
        return res;
      }

      //Phase-space integrals at (increasing) energy grid points:
      struct EGridNode {
        double t;//E/kT
        double sint;//integral within kinematic boundary
        double touched;//integral of touched cells
        std::size_t ntouched;//number of touched cells
      };

      std::vector<EGridNode> integrateAtEGridPts( const VectD& tvals,
                                                  const CellMgr& cellmgr,
                                                  const SABSurveyor& surv,
                                                  SABCfg::IntegrationScheme scheme )
      {
        std::vector<EGridNode> res;
        res.reserve( tvals.size() );
        sIntegralAtE_Result integAtE;
        for ( auto t : tvals ) {
          nc_assert( res.empty() || t > res.back().t );
          sIntegralAtE( EDivKT{t}, cellmgr, surv, scheme, integAtE );
#ifndef NDEBUG
          //The incremental touched-cells integral must equal a direct sum
          //(checked at a few points only, to keep Debug builds fast):
          if ( res.size() % 37 == 0 || res.size() + 1 == tvals.size() ) {
            StableSumKahan touched;
            auto cells = surv.data();
            for ( auto i : ncrange( integAtE.nTouchedCells ) )
              touched.add( cellmgr.getCellIntegral( cells[i].cellidx ) );
            nc_assert_always( floateq( touched.sum(),
                                       integAtE.integralTouchedCells ) );
          }
#endif
          res.push_back( { t, integAtE.integralWithinKB,
                           integAtE.integralTouchedCells,
                           integAtE.nTouchedCells } );
        }
        return res;
      }

      //Intervals [t_i,t_{i+1}) where the interpolation of g(t)=Sint(t)/sqrt(t)
      //disagrees by more than tol (relative) at the geometric midpoint with the
      //cubic in log(t) through the 4 surrounding nodes. This is a cheap error
      //estimate, needing no new integrals:
      void addInterpDisagreementFlags( const VectD& x, const VectD& g,
                                       const HymanInterp& interp, double tol,
                                       std::vector<std::size_t>& flagged )
      {
        const std::size_t n = x.size();
        if ( n < 4 )
          return;
        VectD lx;
        lx.reserve( n );
        for ( auto t : x )
          lx.push_back( std::log( t ) );
        for ( auto i : ncrange( n - 1 ) ) {
          const double tm = std::sqrt( x[i] * x[i+1] );
          if ( !( tm > x[i] && tm < x[i+1] ) )
            continue;
          const double sm = 0.5 * ( lx[i] + lx[i+1] );
          const std::size_t j = ( i == 0 ? 0 : ncmin( i - 1, n - 4 ) );
          double lag = 0.0;
          for ( auto a : ncrange( j, j + 4 ) ) {
            double w = 1.0;
            for ( auto b : ncrange( j, j + 4 ) )
              if ( b != a )
                w *= ( sm - lx[b] ) / ( lx[a] - lx[b] );
            lag = ncfma( g[a], w, lag );
          }
          const double h = interp.eval( tm );
          if ( ncabs( h - lag ) > tol * ncmax( ncabs( lag ), ncabs( h ) ) )
            flagged.push_back( i );
        }
      }

      //Adaptive energy grid over [tlo,thi], for interpolating
      //g(t)=Sint(t)/sqrt(t) with HymanInterp (g is nearly constant at low t,
      //where Sint~sqrt(t)). Starting from cfg.egrid_n0 geometric points, it is
      //refined in passes, each computing its new integrals in one increasing
      //sweep. Intervals are refined when:
      //
      //  * A check point added in the previous pass (the geometric midpoint, or
      //    the points at 1/3 and 2/3 in log(t) with cfg.egrid_ncheck2)
      //    deviates from the interpolation through the previous nodes. Both
      //    intervals adjacent to the check point are then refined.
      //  * The interpolation disagrees with a local cubic in log(t) (see
      //    addInterpDisagreementFlags).
      //  * Sint(t_i) < 0.8*Sint(t_{i+1}), since sampling between grid points
      //    uses the sampler at t_{i+1}, so the acceptance rate is about
      //    Sint(t_i)/Sint(t_{i+1}).
      //
      //The first two use a third of cfg.egrid_tolerance, as a margin for
      //errors away from the points checked.
      std::vector<EGridNode> adaptiveEGrid( const SABCfg::Cfg& cfg,
                                            double tlo, double thi,
                                            const CellMgr& cellmgr,
                                            const SABSurveyor& surv )
      {
        nc_assert_always( tlo > 0.0 && thi > tlo && std::isfinite( thi ) );
        nc_assert_always( cfg.egrid_n0 >= 4 );
        constexpr double min_acceptance = 0.8;
        const double tol = cfg.egrid_tolerance / 3.0;
        const VectD check_fractions = ( cfg.egrid_ncheck2
                                        ? VectD{ 1.0/3.0, 2.0/3.0 }
                                        : VectD{ 0.5 } );

        VectD t0 = geomspace( tlo, thi, cfg.egrid_n0 );
        t0.front() = tlo;
        t0.back() = thi;
        auto nodes = integrateAtEGridPts( t0, cellmgr, surv,
                                          cfg.integScheme );

        //Intervals to refine, by the index of their left node (initially all):
        std::vector<std::size_t> flagged;
        for ( auto i : ncrange( nodes.size() - 1 ) )
          flagged.push_back( i );

        while ( true ) {
          VectD x, g;
          x.reserve( nodes.size() );
          g.reserve( nodes.size() );
          for ( auto& nd : nodes ) {
            x.push_back( nd.t );
            g.push_back( nd.sint / std::sqrt( nd.t ) );
          }
          HymanInterp interp( x, g );
          addInterpDisagreementFlags( x, g, interp, tol, flagged );
          for ( auto i : ncrange( nodes.size() - 1 ) )
            if ( nodes[i].sint < min_acceptance * nodes[i+1].sint )
              flagged.push_back( i );
          std::sort( flagged.begin(), flagged.end() );
          flagged.erase( std::unique( flagged.begin(), flagged.end() ),
                         flagged.end() );

          VectD tnew;
          for ( auto i : flagged ) {
            const double la = std::log( x[i] ), lb = std::log( x[i+1] );
            for ( auto w : check_fractions ) {
              const double tc = std::exp( ncfma( w, lb - la, la ) );
              if ( tc > x[i] && tc < x[i+1] )
                tnew.push_back( tc );
            }
          }
          flagged.clear();
          if ( tnew.empty() )
            break;
          if ( nodes.size() + tnew.size() > cfg.egrid_npts_max )
            NCRYSTAL_THROW2(CalcError,"SAB energy grid needs more than "
                            <<cfg.egrid_npts_max<<" points to meet its"
                            " precision target (please report it to "
                            "the NCrystal developers for investigation).");
          auto newnodes = integrateAtEGridPts( tnew, cellmgr, surv,
                                               cfg.integScheme );
          VectD failed;
          for ( auto& nd : newnodes ) {
            const double diff = ncfma( std::sqrt( nd.t ), interp.eval( nd.t ),
                                       -nd.sint );
            if ( ncabs( diff ) > tol * ncabs( nd.sint ) )
              failed.push_back( nd.t );
          }
          std::vector<EGridNode> merged;
          merged.reserve( nodes.size() + newnodes.size() );
          std::merge( nodes.begin(), nodes.end(),
                      newnodes.begin(), newnodes.end(),
                      std::back_inserter( merged ),
                      []( const EGridNode& a, const EGridNode& b )
                      { return a.t < b.t; } );
          nodes = std::move( merged );
          for ( auto t : failed ) {
            auto it = std::lower_bound( nodes.begin(), nodes.end(), t,
                                        []( const EGridNode& a, double tt )
                                        { return a.t < tt; } );
            nc_assert( it != nodes.begin() && it + 1 != nodes.end()
                       && it->t == t );
            const auto i = static_cast<std::size_t>( it - nodes.begin() );
            flagged.push_back( i - 1 );
            flagged.push_back( i );
          }
        }
        return nodes;
      }

      SABProcImpl::SABProcImpl( const SABCfg::Cfg& cfg,
                                shared_obj<const SABData> sd,
                                std::shared_ptr<const VectD> requested_egrid,
                                SampleSupport sampleSupport,
                                StoreExtraDiagnostics extraDiag )
        : m_kT( sd->temperature().kT() ),
          m_invkT( 1.0 / m_kT ),
          m_origSABData( sd ),
          m_cellmgr( prepareSForProcessor( sd, m_invSScale ) )
      {
        const bool do_sample( sampleSupport == SampleSupport::YES );
        const bool do_diag( extraDiag == StoreExtraDiagnostics::YES );
        {
          //All neutron energies down to E=0 reach the line alpha=beta (with
          //beta>=0), so the grid must cross it:
          const auto& ag = m_cellmgr.sabData().alphaGrid();
          const auto& bg = m_cellmgr.sabData().betaGrid();
          if ( !( ncmax( ag.front(), ncmax( bg.front(), 0.0 ) )
                  < ncmin( ag.back(), bg.back() ) ) )
            NCRYSTAL_THROW2( BadInput, "The S(alpha,beta) grid (alpha in ["
                             << fmt( ag.front() ) << "," << fmt( ag.back() )
                             << "], beta in [" << fmt( bg.front() ) << ","
                             << fmt( bg.back() ) << "]) does not cover the"
                             " line alpha=beta (beta>=0), which neutrons of"
                             " all energies down to E=0 can reach." );
        }
        SABSurveyor surv(m_cellmgr.sabData());

        std::vector<EGridNode> egridNodes;
        {
          auto egspec = determineEGrid( cfg, m_cellmgr.sabData(), m_cellmgr,
                                        surv, requested_egrid );
          egridNodes = ( egspec.grid.empty()
                         ? adaptiveEGrid( cfg, egspec.tlo, egspec.thi,
                                          m_cellmgr, surv )
                         : integrateAtEGridPts( egspec.grid, m_cellmgr, surv,
                                                cfg.integScheme ) );
        }
        m_eGrid.reserve( egridNodes.size() );
        for ( auto& nd : egridNodes )
          m_eGrid.push_back( nd.t );
        nc_assert_always( !m_eGrid.empty() && m_eGrid.front() > 0.0 );
        m_eGridGuide = LogGuideTable( m_eGrid.data(), m_eGrid.size() );

        std::size_t nTouchedCellsMax(0);
        VectD fcAR;//full-cell acceptance rate at each energy grid point
        {
          static_assert
            ( -static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::min())
              >= static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max()),
              "" );
          constexpr auto n_cellmax
            = static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());

          if ( surv.data().size() > n_cellmax )
            NCRYSTAL_THROW2(BadInput,"Too many SAB cells (does not fit in int32)");
          m_sIntegral.reserve( m_eGrid.size() );
          if ( do_sample ) {
            m_nTouchedCells.reserve( m_eGrid.size() );
            fcAR.reserve( m_eGrid.size() );
          }
          if ( do_diag )
            m_diagnostics_fcAR.reserve( m_eGrid.size() );
          for ( auto& nd : egridNodes ) {
            m_sIntegral.push_back( nd.sint );
            if ( do_diag )
              m_diagnostics_fcAR.push_back( nd.touched
                                            ? ( nd.sint / nd.touched )
                                            : 1.0 );
            if (!do_sample)
              continue;

            fcAR.push_back( nd.touched > 0.0 ? ( nd.sint / nd.touched ) : 0.0 );
            //Full-cell sampling, unless initStripLadder picks a strip sampler:
            nc_assert( nd.ntouched < n_cellmax );
            nc_assert( nd.ntouched > 0 );
            m_nTouchedCells.push_back( static_cast<std::uint32_t>
                                       ( nd.ntouched ) );
            nTouchedCellsMax = std::max<std::size_t>( nTouchedCellsMax,
                                                      nd.ntouched );
          }
        }

        if ( do_sample) {
          auto itCellInfo = surv.data().begin();
          VectD cumulFCInt;
          cumulFCInt.reserve(nTouchedCellsMax);
          m_cumulFCInt_cellidx.reserve(nTouchedCellsMax);
          StableSumKahan sum_cumulFCInt;
          for ( std::size_t i = 0; i < nTouchedCellsMax; ++i, ++itCellInfo ) {
            nc_assert(itCellInfo != surv.data().end());
            sum_cumulFCInt.add( m_cellmgr.getCellIntegral( itCellInfo->cellidx ) );
            cumulFCInt.push_back(sum_cumulFCInt.sum());
            m_cumulFCInt_cellidx.push_back(itCellInfo->cellidx);
          }
          m_cumulFCInt = RandIdxByWeightPicker( std::move( cumulFCInt ) );
          nc_assert( m_cumulFCInt.size() == nTouchedCellsMax );
          nc_assert( m_cumulFCInt_cellidx.size() == nTouchedCellsMax );
          m_stripIdx.assign( m_eGrid.size(), -1 );
          initStripLadder( fcAR );
          if ( do_diag ) {
            m_diagnostics_stripTries.reserve( m_eGrid.size() );
            for ( auto k : ncrange( m_eGrid.size() ) ) {
              const EDivKT e{ vectAt( m_eGrid, k ) };
              const double integ = vectAt( m_sIntegral, k );
              SABBetaStripSampler st( m_cellmgr.sabDataPtr(),
                                      m_cellmgr.logSTable(), e );
              const bool ok = st.hasStrips() && integ > 0.0;
              m_diagnostics_stripTries.push_back( ok
                                                  ? st.expectedTries( e, integ )
                                                  : -1.0 );
            }
          }
        }


        m_lowEExtrapolationConstant = ( m_sIntegral.front()
                                        * std::sqrt(1.0/m_eGrid.front()) );
        if ( do_sample )
          initStripBelowGrid();

        {
          //Sint can only increase with E, as KB(E) grows. Closely spaced grid
          //points can however show tiny decreases from the integration noise:
          //clamp those (anything below the precision target), and treat larger
          //decreases as errors:
          for ( auto i : ncrange( std::size_t(1), m_sIntegral.size() ) ) {
            const double prev = vectAt( m_sIntegral, i - 1 );
            double& cur = m_sIntegral[i];
            if ( !( cur < prev ) )
              continue;
            if ( prev - cur > cfg.egrid_tolerance * prev )
              NCRYSTAL_THROW2(CalcError,"Integral(S(alpha,beta)) over"
                              " available phasespace decreases with E at"
                              " E/kT=" << fmtg( vectAt( m_eGrid, i ) ) );
            cur = prev;
          }
        }
        {
          //Interpolation of g(t)=Sint(t)/sqrt(t) with HymanInterp (cubic
          //Hermite, fourth-order slopes, no overshoot where the data are
          //monotone). g is nearly constant at low t where Sint~sqrt(t), which is
          //checked at the first point:
          VectD gvals;
          gvals.reserve( m_eGrid.size() );
          for ( auto i : ncrange( m_eGrid.size() ) )
            gvals.push_back( m_sIntegral[i] / std::sqrt( m_eGrid[i] ) );
          m_sIntegralSmooth.emplace( m_eGrid, gvals );
          //The slope of Sint=sqrt(t)*g at the first point must be within 10% of
          //that of sqrt(t)*g(t0), i.e. |2*t0*g'(t0)/g(t0)| <= 0.1 (not very
          //tight, just a sanity check):
          const double t0 = m_eGrid.front();
          const double gslope = m_sIntegralSmooth.value().nodeSlopes().front();
          if ( !( ncabs( 2.0 * t0 * gslope ) <= 0.1 * ncabs( gvals.front() ) ) ) {
            const double expected = 0.5 * m_sIntegral.front() / t0;
            const double actual = ncfma( std::sqrt( t0 ), gslope, expected );
            NCRYSTAL_THROW2(CalcError,"Integral(S(alpha,beta)) over available"
                            " phasespace does not appear to decrease"
                            " as sqrt(E) at low E (actual slope was "
                            <<fmtg(actual)
                            <<" while expected value was around "
                            <<fmtg(expected)<<".");
          }
          //No constraint on the slope at Emax (the grid ends where enough of
          //the kernel is covered, and SABExtended takes over above).
        }
        m_sIntegralUnscaled.reserve( m_sIntegral.size() );
        for ( auto v : m_sIntegral )
          m_sIntegralUnscaled.push_back( v * m_invSScale );
      }

      double SABProcImpl::phaseSpaceIntegral( NeutronEnergy ekin ) const
      {
        //fixme: just a temporary implementation, we can certainly improve and
        //for instance cache some of the 1.0/(x1-x0) factors.
        const double t = ekin.dbl() * m_invkT;
        auto& x = m_eGrid;
        auto& y = m_sIntegral;
        nc_assert( x.size() == y.size() );
        nc_assert( x.size() >= 2 );
        nc_assert( x.front() >= 0.0 );
        nc_assert( t >= 0.0 );
        nc_assert( std::isfinite(t) );

        if (t >= x.back())
          return y.back();
        if (t <= x.front()) {
          return m_lowEExtrapolationConstant * std::sqrt(t);
        }
        nc_assert( m_sIntegralSmooth.has_value() );
        //Sint=sqrt(t)*g(t), clipped at 0.0 just in case (HymanInterp keeps the
        //sign of the node values):
        return ncmax( 0.0, std::sqrt( t ) * m_sIntegralSmooth.value().eval( t ) );
      }



      SABProcessor::AlphaBetaOutcome
      SABProcImpl::sampleAlphaBeta( RNG& rng, NeutronEnergy ekin ) const
      {
        nc_assert( canSample() );
        //Find the index of the first value in the energy grid whose value exceeds
        //ekin (or the last value if ekin>=m_eGrid.back()):
        const EDivKT e{ ekin.dbl()*m_invkT };
        if ( m_stripBelowGrid.has_value() ) {
          if ( e.dbl() < lowEFloorEDivKT )
            return sampleBelowFloor( rng, e );
          if ( e.dbl() < m_eGrid.front() ) {
            auto r = m_stripBelowGrid.value().sample( rng, e );
            return { r.alpha, r.beta, r.ntries };
          }
        }
        nc_assert(m_eGrid.size()>=2);
        std::size_t idx_E_overlay
          = m_eGridGuide.lowerBoundIdx(m_eGrid.data(), m_eGrid.size(), e.dbl());
        if ( idx_E_overlay == m_eGrid.size() )
          --idx_E_overlay;
        nc_assert(m_nTouchedCells.size()==m_eGrid.size());
        if ( e.dbl() > 0.0 && e.dbl() <= m_eGrid.back() ) {
          const std::int32_t istrip = vectAt( m_stripIdx, idx_E_overlay );
          if ( istrip >= 0 ) {
            auto r = vectAt( m_strips, static_cast<std::size_t>( istrip ) )
              .sample( rng, e );
            return { r.alpha, r.beta, r.ntries };
          }
        }
        const std::size_t nCellsTouched = vectAt(m_nTouchedCells,idx_E_overlay);
        ABRes alphabeta = sampleAlphaBetaFC( rng, e, nCellsTouched );
        nc_assert( ncsquare( alphabeta.a - alphabeta.b )
                   <= e.dbl() * 4.0 * alphabeta.a );
        SABProcessor::AlphaBetaOutcome res;
        res.alpha = alphabeta.a;
        res.beta = alphabeta.b;
        res.ntries = static_cast<unsigned>( alphabeta.ntries );
        return res;
      }

      void SABProcImpl::initStripLadder( const VectD& fcAR )
      {
        //Walking up the energy grid, a strip sampler is built with emax at the
        //highest grid point within ladderFactor of the lowest one not yet
        //covered, and used at the grid points where stripTryCost times its
        //expected tries is below the expected full-cell tries. Stops at the
        //first sampler used nowhere, if the full-cell acceptance rate is at
        //least minFullCellAR at all its grid points:
        constexpr double ladderFactor = 5.0;
        constexpr double stripTryCost = 1.0;
        constexpr double minFullCellAR = 0.4;
        const std::size_t n = m_eGrid.size();
        nc_assert_always( fcAR.size() == n && m_stripIdx.size() == n );
        std::size_t h = 0;
        while ( h < n ) {
          std::size_t j = h;
          while ( j + 1 < n
                  && ( vectAt( m_eGrid, j + 1 )
                       <= ladderFactor * vectAt( m_eGrid, h ) ) )
            ++j;
          SABBetaStripSampler st( m_cellmgr.sabDataPtr(),
                                  m_cellmgr.logSTable(),
                                  EDivKT{ vectAt( m_eGrid, j ) } );
          if ( !st.hasStrips() ) {
            h = j + 1;
            continue;//S=0 within KB(emax), nothing to sample
          }
          bool used = false;
          double minAR = 1.0;
          for ( std::size_t k = h; k <= j; ++k ) {
            const double ar = vectAt( fcAR, k );
            const double integ = vectAt( m_sIntegral, k );
            minAR = ncmin( minAR, ar );
            if ( !( integ > 0.0 && ar > 0.0 ) )
              continue;
            const EDivKT e{ vectAt( m_eGrid, k ) };
            if ( stripTryCost * st.expectedTries( e, integ ) < 1.0 / ar ) {
              vectAt( m_stripIdx, k )
                = static_cast<std::int32_t>( m_strips.size() );
              used = true;
            }
          }
          if ( used )
            m_strips.push_back( std::move( st ) );
          else if ( minAR >= minFullCellAR )
            break;
          h = j + 1;
        }
      }

      void SABProcImpl::initStripBelowGrid()
      {
        //Efficient at all e below its emax (tries independent of e), so a
        //single sampler covers all energies below the grid:
        if ( !( m_lowEExtrapolationConstant > 0.0 ) )
          return;//S=0 within KB(e) at low e: never sampled there
        SABBetaStripSampler st( m_cellmgr.sabDataPtr(), m_cellmgr.logSTable(),
                                EDivKT{ m_eGrid.front() } );
        if ( !st.hasStrips() )
          return;
        m_stripBelowGrid = std::move( st );
      }

      SABProcessor::AlphaBetaOutcome
      SABProcImpl::sampleBelowFloor( RNG& rng, EDivKT e ) const
      {
        //Beta sampled at e=lowEFloorEDivKT (resampling the rare
        //beta<-E/kT), alpha uniformly within KB(E):
        nc_assert( e.dbl() >= 0.0 && e.dbl() < lowEFloorEDivKT );
        const auto& st = m_stripBelowGrid.value();
        nc_assert( st.emaxDivKT().dbl() >= lowEFloorEDivKT );
        SABProcessor::AlphaBetaOutcome res;
        res.ntries = 0;
        while ( true ) {
          auto r = st.sample( rng, EDivKT{ lowEFloorEDivKT } );
          res.ntries += r.ntries;
          if ( r.beta >= -e.dbl() ) {
            res.beta = r.beta;
            break;
          }
        }
        const auto alim = getAlphaLimits( e, res.beta );
        res.alpha = randInterval( rng, alim.first, alim.second );
        return res;
      }

      SABProcImpl::ABRes
      SABProcImpl::sampleAlphaBetaFC( RNG& rng,
                                      EDivKT e,
                                      std::size_t nCellsTouched ) const
      {
        ABRes res;
        res.ntries = 1;
        const double foure = 4.0*e.dbl();
        while ( true ) {
          std::size_t randidx = m_cumulFCInt.pick( rng, nCellsTouched );
          auto cellidx = vectAt(m_cumulFCInt_cellidx,randidx);
          auto cellData = m_cellmgr.lookupCellInfo( cellidx );
          auto ab = FullCellSampler::sampleOneAlphaBeta(cellData,rng);
          if ( ncsquare(ab.first-ab.second)<=foure*ab.first ) {
            res.a = ab.first;
            res.b = ab.second;
            return res;
          }
          ++res.ntries;
          //Guard against endless loops in pathological cases: warn at 1e4
          //tries, abort at 1e5:
          if ( res.ntries == 10000 ) {
            static Msg::WarnFirstTime wf;
            wf("S(alpha,beta) full-cell sampling needed 10000 tries.");
          }
          nc_assert_always( res.ntries < 100000 );
        }
      }

      void SABProcImpl::toJSONFullDiagnostics( std::ostream& os ) const
      {
        //Lots of diagnostics: Things related to egrid search (like how many
        //calls needed for root finding, how emin/emax were determined, how many
        //points added below emin ).
        //
        //Things related to: sampling method at each egrid

        nc_assert_always( m_eGrid.size()
                          < std::numeric_limits<std::uint32_t>::max() );
        os <<"{\"egrid\":";
        streamJSONHugeDblVect(os, VectD( m_eGrid ));
        os <<",\"sintegral\":";
        streamJSONHugeDblVect(os, VectD( m_sIntegralUnscaled ));

        os <<",\"kT\":";
        streamJSON(os,m_kT);

        const bool can_sample = canSample();

        if (!m_diagnostics_fcAR.empty()) {
          os << ",\"diagnostics\":{\"full_cell_AR\":";
          streamJSON(os,m_diagnostics_fcAR);
          if ( !m_diagnostics_stripTries.empty() ) {
            os << ",\"dedicated_strip_expected_tries\":";
            streamJSON(os,m_diagnostics_stripTries);
          }
          os << '}';
        }

        if ( can_sample ) {
          os << ",\"sample_info\":{";
          nc_assert_always( m_nTouchedCells.size() == m_eGrid.size() );
          const std::uint32_t n
            = static_cast<std::uint32_t>(m_nTouchedCells.size());
          //Beta-strip sampling: grid points and expected tries:
          std::vector<std::uint32_t> epts_with_strip;
          VectD strip_tries;
          for ( std::uint32_t i = 0; i < n; ++i ) {
            const auto is = vectAt(m_stripIdx,i);
            if ( is < 0 )
              continue;
            epts_with_strip.push_back(i);
            strip_tries.push_back( vectAt( m_strips,
                                           static_cast<std::size_t>(is) )
                                   .expectedTries( EDivKT{vectAt(m_eGrid,i)},
                                                   vectAt(m_sIntegral,i) ) );
          }
          os << "\"epts_with_beta_strip_sampling\":";
          streamJSON(os,epts_with_strip);
          os << ",\"beta_strip_expected_tries\":";
          streamJSON(os,strip_tries);
          os << '}';
        }
        os << '}';
      }

      inline const SABProcImpl* sp_cimpl( const void* implptr ) noexcept
      {
        assert(implptr!=nullptr);//no exceptions
        return static_cast<const SABProcImpl*>(implptr);
      }

    }//anon namespace
  }
}

NCS::SABProcessor::~SABProcessor()
{
  if ( m_impl ) {
    SABProcImpl* sp = static_cast<SABProcImpl*>(m_impl);
    delete sp;
    m_impl = nullptr;
  }
}

double NCS::SABProcessor::phaseSpaceIntegral( NeutronEnergy ekin ) const
{
  auto sp = sp_cimpl(m_impl);
  return sp->phaseSpaceIntegral(ekin) * sp->m_invSScale;
}

NC::CrossSect
NCS::SABProcessor::crossSectionUnitSigmaBound( NeutronEnergy ekin ) const
{
  //Fixme: optimise (most likely by caching in new vector in which we can
  //directly interpolate results as per the old SABXSProvider?)
  //fixme: keep synchronized with getEMaxInfo()
  auto sp = sp_cimpl(m_impl);
  if ( !( ekin.dbl() * sp->m_invkT > sp->m_eGrid.front() ) ) {
    //Below the energy grid Sint=k*sqrt(E/kT), so xs=k*sqrt(kT)/(4*sqrt(E)),
    //computed in a way which avoid zero division and infinities in the result.
    const double e = ncmax( ekin.dbl(),
                            std::numeric_limits<double>::denorm_min() );
    return CrossSect{ ncmin( sp->m_lowEExtrapolationConstant * sp->m_invSScale
                             * 0.25 * std::sqrt( sp->m_kT ) / std::sqrt( e ),
                             std::numeric_limits<double>::max() ) };
  }
  const double Sint = sp->phaseSpaceIntegral(ekin) * sp->m_invSScale;
  const double kTdiv4 = sp->m_kT * 0.25;//fixme: cache?
  return CrossSect{ Sint * kTdiv4 / ekin.dbl() };
}

NC::IsoXSTable NCS::SABProcessor::createXSTabulation() const
{
  //On the grid, e*xs = (kT/4)*Sint(t) = sqrt(e)*(sqrt(kT)/4)*g(t), with g the
  //HymanInterp interpolation (a cubic Hermite in u, the same u in t and e):
  auto sp = sp_cimpl(m_impl);
  const auto& x = sp->m_eGrid;
  const auto& y = sp->m_sIntegral;
  nc_assert_always( x.size() >= 2 && x.size() == y.size() );
  nc_assert_always( sp->m_sIntegralSmooth.has_value() );
  const auto m = sp->m_sIntegralSmooth.value().nodeSlopes();
  nc_assert_always( m.size() == x.size() );
  const double kT = sp->m_kT;
  //Also undoing the scaling of the S values (see prepareSForProcessor):
  const double c = 0.25 * std::sqrt( kT ) * sp->m_invSScale;
  IsoXSTable res;
  res.E.reserve( x.size() );
  for ( auto t : x )
    res.E.push_back( t * kT );
  res.q.reserve( 4 * ( x.size() - 1 ) );
  for ( auto i : ncrange( x.size() - 1 ) ) {
    const double g0 = y[i] / std::sqrt( x[i] );
    const double g1 = y[i+1] / std::sqrt( x[i+1] );
    const double h = x[i+1] - x[i];
    const double d0 = h * m[i], d1 = h * m[i+1];
    res.q.insert( res.q.end(), { c * g0, c * d0,
                                 c * ncfma( 3.0, g1 - g0,
                                            -ncfma( 2.0, d0, d1 ) ),
                                 c * ncfma( 2.0, g0 - g1, d0 + d1 ) } );
  }
  //Below the grid xs=K*sqrt(kT)/(4*sqrt(e)), above it xs=(kT/4)*Sint(Emax)/e
  //(cf. crossSectionUnitSigmaBound):
  res.tail[0] = sp->m_lowEExtrapolationConstant * c;
  res.tail[2] = 0.25 * kT * y.back() * sp->m_invSScale;
  return res;
}

NCS::SABProcessor::EPtInfo NCS::SABProcessor::getEMaxInfo() const
{
  auto sp = sp_cimpl(m_impl);
  EPtInfo res;
  res.E_div_kT = EDivKT{ sp->m_eGrid.back() };
  res.phaseSpaceIntegral = sp->m_sIntegralUnscaled.back();
  res.ekin = NeutronEnergy{ res.E_div_kT.dbl() * sp->m_kT };
  res.crossSectionUnitSigmaBound
    = CrossSect{ res.phaseSpaceIntegral * 0.25 / res.E_div_kT.dbl() };
  return res;
}

NCS::SABProcessor::AlphaBetaOutcome
NCS::SABProcessor::sampleScatterAlphaBeta( RNG& rng, NeutronEnergy ekin ) const
{
  return sp_cimpl(m_impl)->sampleAlphaBeta(rng,ekin);
}

NCS::SABProcessor::SABProcessor( const SABCfg::Cfg& cfg,
                                 shared_obj<const SABData> sd,
                                 std::shared_ptr<const VectD> req_egrid,
                                 SampleSupport sampleSupport,
                                 StoreExtraDiagnostics extraDiag )
  : m_impl( static_cast<void*>
            ( new SABProcImpl( cfg, std::move(sd),std::move(req_egrid),
                               sampleSupport,extraDiag ) ) )
{
}

NCS::SABProcessor::SABProcessor( SABProcessor&& o ) noexcept
{
  std::swap(m_impl,o.m_impl);
}

NCS::SABProcessor& NCS::SABProcessor::operator=( SABProcessor&& o ) noexcept
{
  std::swap(m_impl,o.m_impl);
  return *this;
}

void NCS::SABProcessor::toJSONFullDiagnostics( std::ostream& os ) const
{
  return sp_cimpl(m_impl)->toJSONFullDiagnostics(os);
}

bool NCS::SABProcessor::hasSampleSupport() const
{
  return sp_cimpl(m_impl)->canSample();
}

double NCS::SABProcessor::kT() const
{
  return sp_cimpl(m_impl)->m_kT;
}

const NC::VectD& NCS::SABProcessor::getEDivKTGrid() const
{
  return sp_cimpl(m_impl)->m_eGrid;
}

const NC::VectD& NCS::SABProcessor::getPhaseSpaceIntegralAtGrid() const
{
  return sp_cimpl(m_impl)->m_sIntegralUnscaled;
}

NC::shared_obj<const NC::SABData> NCS::SABProcessor::sabDataPtr() const
{
  return sp_cimpl(m_impl)->m_origSABData;
}

void NCS::SABProcessor::toJSONProcessInfo( std::ostream& os,
                                           Optional<SigmaBound> sigma_scale,
                                           Optional<std::string>
                                           extension_method ) const
{
  const auto sp = sp_cimpl(m_impl);
  const auto& sab = sp->m_cellmgr.sabData();
  const auto emax = getEMaxInfo().ekin;
  const auto negrid = sp->m_eGrid.size();
  const auto emin = NeutronEnergy{ sp->m_eGrid.front()*sp->m_kT };
  {
    std::ostringstream tmp;
    tmp << "nalpha="<<sab.alphaGrid().size()<<";nbeta="<<sab.betaGrid().size();
    tmp << ";Emax="<<emax;
    tmp << ";T="<<sab.temperature();
    tmp << ";M="<<sab.elementMassAMU();
    if ( extension_method.has_value() && extension_method.value() != "freegas")
      tmp << ";extend="<<extension_method.value();
    if ( sigma_scale.has_value() )
      tmp << ";sigma_free="
          <<sigma_scale.value().free(sab.elementMassAMU());


    streamJSONDictEntry( os, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( os, "Emax", emax.dbl()  );
  streamJSONDictEntry( os, "Emin", emin.dbl()  );

  streamJSONDictEntry( os, "negrid", negrid  );
  streamJSONDictEntry( os, "T", sab.temperature().dbl()  );
  streamJSONDictEntry( os, "M", sab.elementMassAMU().dbl()  );
  if ( sigma_scale.has_value() ) {
    streamJSONDictEntry( os, "sigma_bound=",
                         sigma_scale.value().dbl() );
    streamJSONDictEntry( os, "sigma_free=",
                         sigma_scale.value().free(sab.elementMassAMU()).dbl() );
  }
  streamJSONDictEntry( os, "nalpha", sab.alphaGrid().size()  );
  streamJSONDictEntry( os, "nbeta", sab.betaGrid().size()  );
  if ( extension_method.has_value() )
    streamJSONDictEntry( os, "extension_method", extension_method.value() );
  streamJSONDictEntry( os, "sabprocessor_uid",
                       getUniqueID().value, JSONDictPos::LAST );
}
