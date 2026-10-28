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
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCSpline.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <fstream>

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

//#define NCRYSTAL_SABPROC_SMOOTH_INTERP_LOGLOG

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

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
              e = integrateAlphaInterval_fast(a_low,s_low,*itA,*itS,
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
              double aint = integrateAlphaInterval_fast(a_low,s_low,*itA,*itS,
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
        ScatterOutcomeIsotropic
        ab2scatout( RNG&, NeutronEnergy, const AlphaBetaOutcome& ) const;

        void toJSON( std::ostream& ) const;

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
        Optional<PCHIPInterp> m_sIntegralSmooth;

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

        CellMgr m_cellmgr;

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
        ABOutcome sampleBelowFloor( RNG&, double E_div_kT ) const;
      private:

        //Sample a scattering event:
        struct ABRes { double a, b; std::size_t ntries; };
        ABRes sampleAlphaBetaFC( RNG&, double E_div_kT,
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

      void sIntegralAtE( double E_div_kT,
                         const CellMgr& mgr,
                         const SABSurveyor& surv,
                         StdLogLinCellIntegrator::IntegrationScheme scheme,
                         sIntegralAtE_Result& result )
      {
        auto& prev = result.prevCache;
        const double prev_E = ( ( prev.E_div_kT >= 0.0
                                  && prev.E_div_kT <= E_div_kT )
                                ? prev.E_div_kT : -1.0 );
        const double prev_contrib( prev_E >= 0.0
                                   ? prev.sumFullyCoveredCells.sum() : 0.0 );


        nc_assert(E_div_kT>0&&std::isfinite(E_div_kT));
        //Find all touched cells:
        auto survCells = surv.data();
        //std::upper bound returns the first element where the comparision
        //returns true. The list must be sorted so the comparison will return
        //[false,false,...,false,true,true,...,true].

        //we find the first element coming AFTER E_div_kT

        //FIXME: Why do we need the binary search? Why not just proceed over all
        //cells until e < ci.e_touch?

        auto itSearchBegin = survCells.begin();
        if ( prev_E >= 0.0 && prev.nLastTouched > 2 )
          itSearchBegin += prev.nLastTouched-1;//fixme: is this actually helping in a measurable way??

        auto itLastTouchedE
          = std::upper_bound(itSearchBegin,
                             survCells.end(),
                             E_div_kT,
                             [](double e,
                                const SABSurveyor::CellInfo& ci)
                             { return e < ci.e_touch; } );
        //upper_bound returns one past the last touched:
        if( itLastTouchedE == survCells.begin() )
          NCRYSTAL_THROW2(BadInput,"Phase-space of neutron energy E/kT="
                          <<fmt(E_div_kT)
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
        for ( auto it = active.begin(); it != active.end(); ++it ) {
          nc_assert( E_div_kT >= it->e_touch );
          nc_assert( prev_E < it->e_cover );

          const double fullCellIntegral = mgr.getCellIntegral( it->cellidx );
          if ( E_div_kT >= it->e_cover ) {
            //fully covers cells:
            sumCoveredCells.add(fullCellIntegral);
          } else {
            double contrib;
            const double cutoffref = prev_contrib * threshold;
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
              StdLogLinCellIntegrator::integrateWithinKB( cellData, E_div_kT,
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
                //Explicit std::fma for the inner Horner steps:
                const double s = t*t*t*std::fma( t, std::fma(t,6.0,-15.0), 10.0 );
                contrib *= s;
              }
              sum.add( contrib );
            }
          }
        }
        //Handle interactions with prev_E cache:
        if ( prev_E >= 0.0 )
          sumCoveredCells.add(prev.sumFullyCoveredCells);

        //Prepare cache for next:
        active.erase( std::remove_if( active.begin(), active.end(),
                                      [E_div_kT]( const SABSurveyor::CellInfo& c )
                                      { return c.e_cover <= E_div_kT; } ),
                      active.end() );
        prev.sumFullyCoveredCells.set(sumCoveredCells);
        prev.E_div_kT = E_div_kT;
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


      //Given f(e) = Sintegral(e)/sqrt(e), and a search range [e1,e2] find in as
      //few f(e) evaluations as possible, the highest possible emin value where
      //f(emin)~=f(e1). It is assumed that [e1,e2] is already narrow enough that
      //f(e) is a monotonic function over the interval.
      double determineEMinDivKT( const SABCfg::Cfg& cfg,
                                 const PairDD& e_search_range,
                                 const std::function<double(double)>& f_of_e )
      {
        //fixme: this function could use some luxury parameters
        const double root_acc = cfg.egrid_emin_accuracy;

        const double e1 = e_search_range.first;
        const double e2 = e_search_range.second;
        const double logf1 = std::log(f_of_e( e1 ));
        const double logf2 = std::log(f_of_e( e2 ));
        //f usually increases, but can decrease (e.g. for kernels where S falls
        //off quickly along the diagonal), so the target is a small change in
        //either direction (flat or non-finite: no low-energy regime found):
        const double logS_raise_target = 1e-3*(logf2-logf1);//fixme: lux?
        if ( !( logS_raise_target != 0.0 && std::isfinite(logS_raise_target) ) )
          return e2;
        const double sgn = ( logS_raise_target > 0.0 ? 1.0 : -1.0 );
        const double logf1pluseps = logf1 + logS_raise_target;
        const double loge1 = std::log(e1);
        const double loge2 = std::log(e2);
        unsigned long ncalls = 2;
        auto froot = [&f_of_e,logf1pluseps,&ncalls,sgn,
                      logf1,logf2,loge1,loge2]( double loge )
        {
          ++ncalls;
          nc_assert(std::isfinite(loge));
          double val;
          if ( !(loge>loge1) ) {
            val = logf1;
          } else if ( !(loge<loge2) ) {
            val = logf2;
          } else {
            val = std::log(f_of_e(std::exp(loge)));
          }
          return sgn * ( val - logf1pluseps );
        };
        if ( false ) {//fixme
          std::string fn = "ncrystal_logemin_info.txt";
          NCRYSTAL_WARN("Writing to "<<fn
                        <<" if it does not already exist");
          if (!file_exists(fn)) {
            std::ofstream ofs(fn.c_str(), std::ofstream::out);
            ofs << "#ncrystal_xycurve\n";
            ofs << "#colnames = loge;deltalogf\n";
            for ( auto loge : linspace(loge1,loge2,100000) )
              ofs<<fmt(loge)<<" "<<fmt(froot(loge))<<"\n";
            ofs.close();
          }
        }

        //Rough narrowing with exponential steps before the fine root refinement
        //below (the general root finder can not know like we do that the root
        //often sits very near one edge):
        double bracket_lo = loge1, bracket_hi = loge2;
        {
          //fixme: stepFactor of 10.0 might be better?
          constexpr double stepFactor = 100.0;//~2 decades per step
          const double logStep = std::log(stepFactor);
          double cand = bracket_lo;
          while ( cand + logStep < bracket_hi ) {
            cand += logStep;
            if ( froot(cand) >= 0.0 ) {
              bracket_hi = cand;
              break;
            }
            bracket_lo = cand;
          }
        }

        double logemin;
        try {
          logemin = findRoot2(froot,bracket_lo,bracket_hi, root_acc);
        } catch ( NC::Error::CalcError& e ) {
          logemin = -1.0;
          //Write fct to file for debugging:
          std::string fn = "ncrystal_logemin_fail_info.txt";
          NCRYSTAL_WARN("Failed to find emin via root-finding. Will attempt"
                        " to write debugging info to "<<fn
                        <<" if it does not already exist");
          if (!file_exists(fn)) {
            std::ofstream ofs(fn.c_str(), std::ofstream::out);
            ofs << "#ncrystal_xycurve\n";
            ofs << "#colnames = loge;deltalogf\n";
            for ( auto loge : linspace(loge1,loge2,100000) )
              ofs<<fmt(loge)<<" "<<fmt(froot(loge))<<"\n";
            ofs.close();
          }
          throw;
        }

        return ncclamp(std::exp(logemin),e1,e2);
      }

      //Based on full-cell integrals, determine range in which a more refined
      //search for EMin/kT should take place.
      PairDD determineEMinDivKTminRange( const VectD& cumulFullCellIntegral,
                                         const VectD& cellTouch_E_div_kT,
                                         const double EMax_div_kT )
      {
        const VectD& e = cellTouch_E_div_kT;
        const VectD& s = cumulFullCellIntegral;

        //We check for the point where the S-integral has increased some amount
        //from the initial value.

        nc_assert_always(s.size()==e.size());
        nc_assert_always(s.size()>=16);
        const double s1(s.front());
        const double s2(s.back());
        nc_assert_always( s2>s1 && s1>=0.0 && std::isfinite(s2) );
        const std::size_t ns = s.size();

        std::size_t i1 = 0;
        {
          constexpr double eps_emin1 = 1e-6;
          const double tgt_s1 = (s1?s1*std::pow(s2/s1,eps_emin1):s2*eps_emin1);
          while ( i1<ns && vectAt(s,i1)<tgt_s1 )
            ++i1;
          nc_assert_always( i1+8 < ns );
        }

        //Upper reference point: use EMax_div_kT directly rather than a
        //touched-cell-array crossing index (near-degenerate cell "touch"
        //energies can swap order under FMA-contraction noise, perturbing the
        //Emin found downstream). EMax_div_kT is most often a single
        //well-defined energy, and the real S-integral is still evaluated there
        //via f_of_e:
        constexpr double safety1 = 0.001;
        PairDD res( vectAt(e,i1)*safety1, EMax_div_kT );
        nc_assert_always( std::isfinite(res.first) );
        nc_assert_always( std::isfinite(res.second) );
        nc_assert_always( std::isfinite(res.first) );
        nc_assert_always( res.second >= res.first );
        return res;
      }

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

      //fixme: clarify that this returns in E/kT, but requested_egrid is in eV
      VectD determineEGrid( const SABCfg::Cfg& cfg,
                            const SABData& sab,
                            const CellMgr& cellmgr,
                            const SABSurveyor& surv,
                            std::shared_ptr<const VectD> requested_egrid )
      {
        const auto scheme = cfg.integSchemeDetermineEGrid;
        static_assert( std::is_same<decltype(cfg.egrid_npts),unsigned>
                       ::value, "" );
        auto npts = static_cast<std::size_t>( cfg.egrid_npts );

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
            }
          } else {
            //complete egrid given directly by request, so in this case we
            //simply honour it and return it directly (we can only assume that
            //the users know what they are doing):
            if ( !( re.size() >= 10 ) )
              NCRYSTAL_THROW2(BadInput,"Invalid requested egrid (must have"
                              " length 1, 3, or >=10)");
            VectD res;
            res.reserve( re.size() );
            for ( auto e : re )
              res.push_back( e * invkT );
            return res;
          }
        }

        {
          nc_assert_always( npts >= 10 );
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
          auto r = determineEMinDivKTminRange( sIntByEtouch.sint, sIntByEtouch.e,
                                               E_min_max.second.get()/kT );

          sIntegralAtE_Result tmp_result;
          auto f = [&cellmgr, &surv, &tmp_result,scheme]( double E_div_kT ) {
            nc_assert(E_div_kT>0.0);
            const double se = std::sqrt(E_div_kT);
            nc_assert(se>0.0);
            sIntegralAtE( E_div_kT, cellmgr, surv, scheme, tmp_result );
            return tmp_result.integralWithinKB / se;
          };
          E_min_max.first = NeutronEnergy{ kT*determineEMinDivKT( cfg, r, f ) };

        }
        if ( E_min_max.first >= E_min_max.second ) {
          NCRYSTAL_WARN("Lowering Emin to 10% of Emax");
          E_min_max.first = NeutronEnergy( 0.1 * E_min_max.second.get() );
        }
        if ( E_min_max.first.get() < 1e-50*E_min_max.second.get() ) {
          NCRYSTAL_WARN("Raising Emin to 1e-50*Emax");
          E_min_max.first = NeutronEnergy( 1e-50 * E_min_max.second.get() );
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
          constexpr double lowE_factor = 1e-4;
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

        const double invkT = 1.0/kT;
        for ( auto& ee : final_egrid )
          ee *= invkT;

        return final_egrid;
      }

      SABProcImpl::SABProcImpl( const SABCfg::Cfg& cfg,
                                shared_obj<const SABData> sd,
                                std::shared_ptr<const VectD> requested_egrid,
                                SampleSupport sampleSupport,
                                StoreExtraDiagnostics extraDiag )
        : m_kT( sd->temperature().kT() ),
          m_invkT( 1.0 / m_kT ),
          m_cellmgr(sd)
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

        m_eGrid = determineEGrid( cfg, m_cellmgr.sabData(), m_cellmgr,
                                  surv, requested_egrid );
        nc_assert_always( !m_eGrid.empty() && m_eGrid.front() > 0.0 );
        m_eGridGuide = LogGuideTable( m_eGrid.data(), m_eGrid.size() );

        const auto scheme = cfg.integScheme;

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
          sIntegralAtE_Result integAtE;
          if ( do_sample ) {
            m_nTouchedCells.reserve( m_eGrid.size() );
            fcAR.reserve( m_eGrid.size() );
          }
          if ( do_diag )
            m_diagnostics_fcAR.reserve( m_eGrid.size() );
          for ( auto& E_div_kT : m_eGrid ) {
            sIntegralAtE( E_div_kT, m_cellmgr, surv, scheme, integAtE );
            m_sIntegral.push_back( integAtE.integralWithinKB );
#ifndef NDEBUG
            //The incremental touched-cells integral must equal a direct sum
            //(checked at a few energy grid points only, to keep Debug builds
            //fast):
            const std::size_t ie = m_sIntegral.size() - 1;
            if ( ie % 37 == 0 || ie + 1 == m_eGrid.size() ) {
              StableSumKahan touched;
              auto cells = surv.data();
              for ( auto i : ncrange( integAtE.nTouchedCells ) )
                touched.add( m_cellmgr.getCellIntegral( cells[i].cellidx ) );
              nc_assert_always( floateq( touched.sum(),
                                         integAtE.integralTouchedCells ) );
            }
#endif
            if ( do_diag )
              m_diagnostics_fcAR.push_back( integAtE.integralTouchedCells
                                            ? ( integAtE.integralWithinKB
                                                /integAtE.integralTouchedCells )
                                            : 1.0 );
            if (!do_sample)
              continue;

            fcAR.push_back( integAtE.integralTouchedCells > 0.0
                            ? ( integAtE.integralWithinKB
                                / integAtE.integralTouchedCells )
                            : 0.0 );
            //Full-cell sampling, unless initStripLadder picks a strip sampler:
            nc_assert( integAtE.nTouchedCells < n_cellmax );
            nc_assert( integAtE.nTouchedCells > 0 );
            m_nTouchedCells.push_back( static_cast<std::uint32_t>
                                       ( integAtE.nTouchedCells ) );
            nTouchedCellsMax = std::max<std::size_t>( nTouchedCellsMax,
                                                      integAtE.nTouchedCells );
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
              const double e = vectAt( m_eGrid, k );
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

        if ( cfg.sIntegralInterp == SABCfg::SIntegralInterpolation::Smooth ) {
          //We set up PCHIP interpolation. For more "proper" PCHIP behaviour, we
          //let the PCHIP class determine endpoint slopes automatically, but do
          //some slope sanity checks afterwards.
          VectD xvals(m_eGrid);
#ifdef NCRYSTAL_SABPROC_SMOOTH_INTERP_LOGLOG
          for ( auto& e : xvals ) {
            nc_assert( e > 0.0 && std::isfinite(e) );
            e = std::log(e);
          }
#endif
          VectD yvals(m_sIntegral);
#ifdef NCRYSTAL_SABPROC_SMOOTH_INTERP_LOGLOG
          for ( auto& e : yvals ) {
            nc_assert( e > 0.0 && std::isfinite(e) );
            e = std::log(e);
          }
#endif
#ifdef NCRYSTAL_SABPROC_SMOOTH_INTERP_LOGLOG
          //f(e)=k*sqrt(e) => log(f(e)) = log(k) + 0.5*log(e)
          //               => d(logf)/d(loge) = 0.5
          const double leftSlope = 0.5;
#else
          //f(e)=k*sqrt(e) => f'(e)=k/(2*sqrt(e))
          //And with k= s(emin)/sqrt(emin) we get f'(emin)=s(emin)/(2emin)
          const auto leftSlope = 0.5 * m_sIntegral.front() / m_eGrid.front();
#endif
          //Always autodetect slope at Emax, since we can not be sure to what
          //degree the limiting behaviour has been reached.
          m_sIntegralSmooth.emplace( xvals, yvals,
                                     NullOpt, NullOpt,
                                     PCHIPInterp::OptimiseForLogSpacedX::YES );
          auto slopes = m_sIntegralSmooth.value().endPointSlopes();
          constexpr double leftSlopeTol = 0.1;//not very tight, just a sanity
                                              //check
          if ( !floateq(leftSlope,slopes.first,leftSlopeTol,1e-10) ) {
            NCRYSTAL_THROW2(CalcError,"Integral(S(alpha,beta)) over available"
                            " phasespace does not appear to decrease"
                            " as sqrt(E) at low E (actual slope was "
                            <<fmtg(slopes.first)
                            <<" while expected value was around "
                            <<fmtg(leftSlope)<<".");
          }
          //No constraint on the slope at Emax (the grid ends where enough of
          //the kernel is covered, and SABExtended takes over above), but the
          //S integral can only increase with E, as KB(E) grows. With such
          //values, the PCHIP interpolation never decreases either:
          for ( auto i : ncrange( std::size_t(1), m_sIntegral.size() ) )
            if ( vectAt( m_sIntegral, i ) < vectAt( m_sIntegral, i - 1 ) )
              NCRYSTAL_THROW2(CalcError,"Integral(S(alpha,beta)) over"
                              " available phasespace decreases with E at"
                              " E/kT=" << fmtg( vectAt( m_eGrid, i ) ) );
        }
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
        if ( m_sIntegralSmooth.has_value() ) {
          if ( t == 0.0 )
            return kInfinity;
#ifdef NCRYSTAL_SABPROC_SMOOTH_INTERP_LOGLOG
          return std::exp(m_sIntegralSmooth.value().eval(std::log(t)));
#else
          //Clip at 0.0 just in case the interpolation in some rare cases could
          //give values below 0.0 (should hopefully almost never happen in
          //practice).
          return ncmax(0.0,m_sIntegralSmooth.value().eval(t));
#endif
        }
        //linear interpolation:
        std::size_t i = m_eGridGuide.lowerBoundIdx(x.data(), x.size(), t); // 1..n-1
        nc_assert( i>0 );
        nc_assert( i<x.size() );
        double x0 = x[i-1], x1 = x[i];
        double y0 = y[i-1], y1 = y[i];
        double r = (t - x0) / (x1 - x0);//fixme: cache more! (if not in smooth mode)
        return nclerp( y0, y1, r );
      }

      ScatterOutcomeIsotropic
      SABProcImpl::ab2scatout( RNG& rng,
                               NeutronEnergy ekin,
                               const AlphaBetaOutcome& alphabeta ) const
      {
        auto r = convertABToDeltaEMuWithIsotropicFallback( alphabeta.alpha,
                                                           alphabeta.beta,
                                                           ekin, m_kT, rng,
                                                           lowEFloorEDivKT );
        return { NeutronEnergy{ ncmax( 0.0, ekin.dbl() + r.deltaE ) },
                 CosineScatAngle{ r.mu } };
      }



      SABProcessor::AlphaBetaOutcome
      SABProcImpl::sampleAlphaBeta( RNG& rng, NeutronEnergy ekin ) const
      {
        nc_assert( canSample() );
        //Find the index of the first value in the energy grid whose value exceeds
        //ekin (or the last value if ekin>=m_eGrid.back()):
        const double E_div_kT = ekin.dbl()*m_invkT;
        if ( m_stripBelowGrid.has_value() ) {
          if ( E_div_kT < lowEFloorEDivKT )
            return sampleBelowFloor( rng, E_div_kT );
          if ( E_div_kT < m_eGrid.front() ) {
            auto r = m_stripBelowGrid.value().sample( rng, E_div_kT );
            return { r.alpha, r.beta, r.ntries };
          }
        }
        nc_assert(m_eGrid.size()>=2);
        std::size_t idx_E_overlay
          = m_eGridGuide.lowerBoundIdx(m_eGrid.data(), m_eGrid.size(), E_div_kT);
        if ( idx_E_overlay == m_eGrid.size() )
          --idx_E_overlay;
        nc_assert(m_nTouchedCells.size()==m_eGrid.size());
        if ( E_div_kT > 0.0 && E_div_kT <= m_eGrid.back() ) {
          const std::int32_t istrip = vectAt( m_stripIdx, idx_E_overlay );
          if ( istrip >= 0 ) {
            auto r = vectAt( m_strips, static_cast<std::size_t>( istrip ) )
              .sample( rng, E_div_kT );
            return { r.alpha, r.beta, r.ntries };
          }
        }
        const std::size_t nCellsTouched = vectAt(m_nTouchedCells,idx_E_overlay);
        ABRes alphabeta = sampleAlphaBetaFC( rng, E_div_kT, nCellsTouched );
        nc_assert( ncsquare( alphabeta.a - alphabeta.b )
                   <= E_div_kT * 4.0 * alphabeta.a );
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
                                  vectAt( m_eGrid, j ) );
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
            if ( stripTryCost * st.expectedTries( vectAt( m_eGrid, k ), integ )
                 < 1.0 / ar ) {
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
                                m_eGrid.front() );
        if ( !st.hasStrips() )
          return;
        m_stripBelowGrid = std::move( st );
      }

      SABProcessor::AlphaBetaOutcome
      SABProcImpl::sampleBelowFloor( RNG& rng, double E_div_kT ) const
      {
        //Beta sampled at e=lowEFloorEDivKT (resampling the rare
        //beta<-E/kT), alpha uniformly within KB(E):
        nc_assert( E_div_kT >= 0.0 && E_div_kT < lowEFloorEDivKT );
        const auto& st = m_stripBelowGrid.value();
        nc_assert( st.emaxDivKT() >= lowEFloorEDivKT );
        SABProcessor::AlphaBetaOutcome res;
        res.ntries = 0;
        while ( true ) {
          auto r = st.sample( rng, lowEFloorEDivKT );
          res.ntries += r.ntries;
          if ( r.beta >= -E_div_kT ) {
            res.beta = r.beta;
            break;
          }
        }
        const auto alim = getAlphaLimits( E_div_kT, res.beta );
        res.alpha = randInterval( rng, alim.first, alim.second );
        return res;
      }

      SABProcImpl::ABRes
      SABProcImpl::sampleAlphaBetaFC( RNG& rng,
                                      double E_div_kT,
                                      std::size_t nCellsTouched ) const
      {
        ABRes res;
        res.ntries = 1;
        const double foure = 4.0*E_div_kT;
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

      void SABProcImpl::toJSON( std::ostream& os ) const
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
        streamJSONHugeDblVect(os, VectD( m_sIntegral ));
        os <<",\"egrid_interpolation\":";
        streamJSON(os,(m_sIntegralSmooth.has_value()?"smooth":"linear"));
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
                                   .expectedTries( vectAt(m_eGrid,i),
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
  return sp_cimpl(m_impl)->phaseSpaceIntegral(ekin);
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
    return CrossSect{ ncmin( sp->m_lowEExtrapolationConstant * 0.25
                             * std::sqrt( sp->m_kT ) / std::sqrt( e ),
                             std::numeric_limits<double>::max() ) };
  }
  const double Sint = sp->phaseSpaceIntegral(ekin);
  const double kTdiv4 = sp->m_kT * 0.25;//fixme: cache?
  return CrossSect{ Sint * kTdiv4 / ekin.dbl() };
}

NCS::SABProcessor::EPtInfo NCS::SABProcessor::getEMaxInfo() const
{
  auto sp = sp_cimpl(m_impl);
  EPtInfo res;
  res.E_div_kT = sp->m_eGrid.back();
  res.phaseSpaceIntegral = sp_cimpl(m_impl)->m_sIntegral.back();
  res.ekin = NeutronEnergy{ res.E_div_kT * sp->m_kT };
  res.crossSectionUnitSigmaBound
    = CrossSect{ res.phaseSpaceIntegral * 0.25 / res.E_div_kT };
  return res;
}

NC::ScatterOutcomeIsotropic
NCS::SABProcessor::sampleScatter( RNG& rng, NeutronEnergy ekin ) const
{
  auto alphabeta = sp_cimpl(m_impl)->sampleAlphaBeta(rng,ekin);
  return sp_cimpl(m_impl)->ab2scatout( rng, ekin, alphabeta );
}

NCS::SABProcessor::ScatOutcomeDiag
NCS::SABProcessor::sampleScatterDiag( RNG& rng, NeutronEnergy ekin ) const
{
  ScatOutcomeDiag res;
  res.alphabeta = sp_cimpl(m_impl)->sampleAlphaBeta(rng,ekin);
  res.ekinmu = sp_cimpl(m_impl)->ab2scatout( rng, ekin, res.alphabeta );
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

void NCS::SABProcessor::toJSON( std::ostream& os ) const
{
  return sp_cimpl(m_impl)->toJSON(os);
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
  return sp_cimpl(m_impl)->m_sIntegral;
}

NC::shared_obj<const NC::SABData> NCS::SABProcessor::sabDataPtr() const
{
  return sp_cimpl(m_impl)->m_cellmgr.sabDataPtr();
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
    if ( !sp->m_sIntegralSmooth.has_value() )
      tmp << ";Einterp=linear";

    streamJSONDictEntry( os, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( os, "Emax", emax.dbl()  );
  streamJSONDictEntry( os, "Emin", emin.dbl()  );
  streamJSONDictEntry( os, "Einterp", (sp->m_sIntegralSmooth.has_value()
                                       ?"smooth":"linear")  );
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
