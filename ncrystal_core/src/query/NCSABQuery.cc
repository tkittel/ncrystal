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

#include "NCrystal/internal/sab/NCSABRefEval.hh"
#include "NCrystal/internal/sab/NCSABKBCellSmpl.hh"
#include "NCrystal/internal/sab/NCSABSurveyor.hh"
#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/sab/NCSABCellInteg.hh"
#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/sab/NCSABRefSampler.hh"
#include "NCrystal/internal/sab/NCSABExtended.hh"
#include "NCrystal/internal/sab/NCSABExtender.hh"
#include "NCrystal/internal/extd_utils/NCSABAnalyser.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/sab/NCSABIntegrator.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/core/NCSmallVector.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCSABQuery.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {

    namespace {

      struct ExtractedSabData {
        shared_obj<const SABData> sab;
        //As the factories would use it (incl. auto-detection for
        //direct kernels):
        Optional<Temperature> teff;
        //Strictly VDOS-derived reference values:
        Optional<Temperature> teff_vdos;
        Optional<double> msd_vdos;//[Aa^2]
        //Energy grid and fraction of the atom in the material:
        DI_ScatKnl::EGridShPtr egrid;
        double fraction;
      };

      ExtractedSabData
      query_impl_extractSabData( const std::string& matcfgstr,
                                 Optional<std::string> atomdsplbl )
      {
        MatCfg cfg(matcfgstr);
        auto info = FactImpl::createInfo(cfg);
        const DynamicInfo* di = nullptr;
        if ( !atomdsplbl.has_value() ) {
          if ( info->getDynamicInfoList().size() > 1 )
            NCRYSTAL_THROW(BadInput, "ATOMDISPLAYLABEL required for"
                           " polyatomic materials");
          di = info->getDynamicInfoList().front().get();
        } else {
          di = InfoUtils::findDynInfo( info, atomdsplbl );
        }
        auto di_knl = dynamic_cast<const DI_ScatKnl*>(di);
        if ( !di_knl )
          NCRYSTAL_THROW(BadInput,"Indicated DynInfo object does not provide"
                         " S(alpha,beta) kernels.");
        Optional<Temperature> teff_vdos;
        Optional<double> msd_vdos;
        if ( auto di_vdos = dynamic_cast<const DI_VDOS*>( di_knl ) ) {
          VDOSEval ve( di_vdos->vdosData() );
          teff_vdos = Temperature{ DoValidate,
                                   ve.calcEffectiveTemperature() };
          msd_vdos = ve.getMSD( ve.calcGamma0() );
        }
        const VDOS::VDOSLux vdoslux( cfg.get_vdoslux() );
        return { DynInfoUtils::extractSABDataFromDynInfo( di_knl, vdoslux ),
                 DynInfoUtils::extractTeffFromDynInfo( di_knl ),
                 teff_vdos,
                 msd_vdos,
                 di_knl->energyGrid(),
                 di_knl->fraction() };
      }

      struct SampleResult {
        VectD a;
        VectD b;
        std::uint64_t ntries = 0;
        double prob1 = -1.0;//probability of edge @ beta=b1
      };
      SampleResult testFullCellOverlaySample( EDivKT e,
                                              RNG& rng,
                                              const CellData& c,
                                              std::uint64_t nsample,
                                              std::uint64_t ntries_max )
      {
        FullCellSampler fcsampler( &c );
        const double foure = 4*e.dbl();
        SampleResult res;
        res.prob1 = fcsampler.probabilityEdge1();
        res.a.reserve( countToSizeT( nsample ) );
        res.b.reserve( countToSizeT( nsample ) );
        while( res.a.size() < nsample ) {
          if ( res.ntries > ntries_max ) {
            NCRYSTAL_WARN("Early abort of FullCellOverlaySample due to bad AR.");
            break;
          }
          ++res.ntries;
          auto ab = fcsampler.sampleAlphaBeta( rng );
          if ( ncsquare(ab.first-ab.second) <= ab.first*foure ) {
            res.a.push_back( ab.first );
            res.b.push_back( ab.second );
          }
        }
        return res;
      }

      SampleResult produceRefCellSamples( EDivKT e,
                                          RNG& rng,
                                          const CellData& c,
                                          std::uint64_t nsample,
                                          std::uint64_t ntries_max )
      {
        RefCellSampler refsampler( c, e );
        const double foure = 4*e.dbl();
        SampleResult res;
        res.prob1 = -1.0;//not available
        res.a.reserve( countToSizeT( nsample ) );
        res.b.reserve( countToSizeT( nsample ) );
        while( res.a.size() < nsample ) {
          if ( res.ntries > ntries_max )
            NCRYSTAL_THROW(CalcError,"RefCellSampler too inefficient.");
          auto ab = refsampler.sampleAlphaBeta( rng );
          res.ntries += ab.ntries;
          if ( !( ncsquare(ab.alpha-ab.beta) <= ab.alpha*foure ) )
            NCRYSTAL_THROW(CalcError,
                           "RefCellSampler produced inaccessible point.");
          res.a.push_back( ab.alpha );
          res.b.push_back( ab.beta );
        }
        return res;
      }

      void query_impl_analyse( std::ostream& os,
                               const ExtractedSabData& xd,
                               const SABAnalyser::Options& opt )
      {
        auto r = SABAnalyser::analyse( *xd.sab, opt );
        auto stream_opt_dbl = [&os]( const Optional<double>& v )
        {
          if ( v.has_value() )
            streamJSON( os, v.value() );
          else
            streamJSON( os, json_null_t{} );
        };
        os << "{\"sabanalyse\":{\"teff\":";
        stream_opt_dbl( r.teff.has_value()
                        ? Optional<double>{ r.teff.value().dbl() }
                        : Optional<double>{} );
        os << ",\"teff_relspread\":";
        streamJSON( os, r.teff_relspread );
        os << ",\"teff_nrows\":";
        streamJSON( os, static_cast<std::uint64_t>(r.teff_nrows) );
        os << ",\"msd\":";
        stream_opt_dbl( r.msd );
        os << ",\"msd_relspread\":";
        streamJSON( os, r.msd_relspread );
        os << ",\"msd_nrows\":";
        streamJSON( os, static_cast<std::uint64_t>(r.msd_nrows) );
        os << ",\"max_row_m0\":";
        streamJSON( os, r.max_row_m0 );
        os << ",\"recoil_center_ratio\":";
        stream_opt_dbl( r.recoil_center_ratio );
        os << ",\"teff_extender\":";
        //What the kernel-extension factories would actually use
        //(explicit NCMAT v8 value, or estimate):
        stream_opt_dbl( xd.teff.has_value()
                        ? Optional<double>{ xd.teff.value().dbl() }
                        : Optional<double>{} );
        os << ",\"teff_vdos\":";
        stream_opt_dbl( xd.teff_vdos.has_value()
                        ? Optional<double>{ xd.teff_vdos.value().dbl() }
                        : Optional<double>{} );
        os << ",\"msd_vdos\":";
        stream_opt_dbl( xd.msd_vdos );
        //What the turn-key consumer (acceptance policy applied,
        //default options, cached) provides for this kernel:
        auto auto_tm = SABAnalyser::estimateTeffMSD( *xd.sab );
        os << ",\"auto\":{\"teff\":";
        stream_opt_dbl( auto_tm.effectiveTemperature.has_value()
                        ? Optional<double>{ auto_tm.effectiveTemperature
                                            .value().dbl() }
                        : Optional<double>{} );
        os << ",\"msd\":";
        stream_opt_dbl( auto_tm.msd );
        os << "},\"temperature\":";
        streamJSON( os, xd.sab->temperature().dbl() );
        os << ",\"elementmassamu\":";
        streamJSON( os, xd.sab->elementMassAMU().dbl() );
        os << ",\"diagnostics\":";
        if ( r.diagnostics == nullptr ) {
          streamJSON( os, json_null_t{} );
        } else {
          auto& d = *r.diagnostics;
          auto stream_statuses
            = [&os]( const std::vector<SABAnalyser::RowStatus>& v )
          {
            os << '[';
            for ( std::size_t i = 0; i < v.size(); ++i ) {
              if ( i )
                os << ',';
              os << static_cast<unsigned>( vectAt( v, i ) );
            }
            os << ']';
          };
          os << "{\"floor_value\":";
          streamJSON( os, d.floor_value );
          os << ",\"alpha\":";
          streamJSONHugeDblVect( os, VectD(d.alpha) );
          os << ",\"m0\":";
          streamJSONHugeDblVect( os, VectD(d.m0) );
          os << ",\"mean\":";
          streamJSONHugeDblVect( os, VectD(d.mean) );
          os << ",\"variance\":";
          streamJSONHugeDblVect( os, VectD(d.variance) );
          os << ",\"teff_row\":";
          streamJSONHugeDblVect( os, VectD(d.teff_row) );
          os << ",\"msd_row\":";
          streamJSONHugeDblVect( os, VectD(d.msd_row) );
          os << ",\"teff_row_status\":";
          stream_statuses( d.teff_row_status );
          os << ",\"msd_row_status\":";
          stream_statuses( d.msd_row_status );
          os << '}';
        }
        os << "}}";
      }

      void query_impl_xs( std::ostream& os,
                          ExtractedSabData xd,
                          SABCfg::ExtenderModel extender,
                          const VectD& energies,
                          int knllux )
      {
        //The cross section contributed to the material, like the factories
        //would provide it, but with the given extender:
        if ( knllux < 0 )
          knllux = SABCfg::sablux_default_luxury;
        auto cfg = SABCfg::createConfig(knllux);
        cfg.extender_model = extender;
        if ( extender == SABCfg::ExtenderModel::SCT && !xd.teff.has_value() )
          NCRYSTAL_THROW(BadInput,"The SCT extender requires an effective"
                         " temperature, which is not available.");
        const double sb = xd.sab->boundXS().dbl();
        auto sabext = SABExtended::create( cfg, xd.sab, xd.teff,
                                           std::move( xd.egrid ) );
        VectD xs;
        xs.reserve( energies.size() );
        for ( auto e : energies )
          xs.push_back( xd.fraction
                        * sabext->xsEvaluator().evalScaled( e, 1.0, sb ) );
        os << "{\"extension_method\":";
        streamJSON( os, sabext->extender().shortName().value_or("custom") );
        os << ",\"xs\":";
        streamJSON( os, xs );
        os << '}';
      }

      void query_impl_proc( std::ostream& os,
                            shared_obj<const SABData> sab,
                            Optional<Temperature> teff,
                            VectD egrid,
                            std::uint64_t nsample,
                            std::uint64_t seed,
                            const std::vector<NeutronEnergy>& sample_ekins,
                            int knllux )
      {
        if ( knllux < 0 )
          knllux = SABCfg::sablux_default_luxury;
        const auto cfg = SABCfg::createConfig(knllux);

        //Create with diagnostics:
        auto sp = makeSO<SABProcessor>
          ( cfg,
            sab,
            std::make_shared<VectD>( std::move(egrid) ),
            SABProcessor::SampleSupport::YES,
            SABProcessor::StoreExtraDiagnostics::YES );
        //Extender per cfg, like the actual factories (SCT needs teff):
        const bool use_sct
          = ( cfg.extender_model == SABCfg::ExtenderModel::SCT
              && teff.has_value() );
        auto spe = ( use_sct
                     ? SABExtended::createWithSCTExtender( sp, teff.value(),
                                                           cfg )
                     : SABExtended::createWithFGExtender( sp ) );

        os << "{\"extender\":";
        streamJSON( os, use_sct ? "SCT" : "FreeGas" );
        os << ",\"sabproc\":";
        spe->processor().toJSONFullDiagnostics(os);
        os << ",\"sample\":";
        if ( !nsample ) {
          streamJSON(os,json_null_t{});
          os << '}';
          return;
        }
        //One sample per energy (a list if more than one), each with a new RNG
        //from the seed:
        const bool many = sample_ekins.size() > 1;
        if ( many )
          os << '[';
        for ( auto i : ncrange( sample_ekins.size() ) ) {
          auto rng = createBuiltinRNG( seed );
          VectD a, b;
          a.reserve( countToSizeT( nsample ) );
          b.reserve( countToSizeT( nsample ) );
          for ( std::uint64_t n = 0; n < nsample; ++n ) {
            auto evt = spe->sampleScatterAlphaBeta( rng, sample_ekins.at(i) );
            a.push_back( evt.first );
            b.push_back( evt.second );
          }
          os << ( i ? "," : "" ) << "{\"alpha\":";
          streamJSONHugeDblVect(os,std::move(a));
          os << ",\"beta\":";
          streamJSONHugeDblVect(os,std::move(b));
          os << '}';
        }
        if ( many )
          os << ']';
        os << '}';
      }

      //Processor JSON (with diagnostics) and the mean number of tries of
      //nperpt samplings of (alpha,beta) at each energy grid point:
      void query_impl_proctries( std::ostream& os,
                                 shared_obj<const SABData> sab,
                                 std::uint64_t nperpt,
                                 std::uint64_t seed,
                                 int knllux )
      {
        if ( knllux < 0 )
          knllux = SABCfg::sablux_default_luxury;
        SABProcessor sp( SABCfg::createConfig(knllux), std::move(sab),
                         std::make_shared<VectD>(),
                         SABProcessor::SampleSupport::YES,
                         SABProcessor::StoreExtraDiagnostics::YES );
        os << "{\"sabproc\":";
        sp.toJSONFullDiagnostics(os);
        auto rng = createBuiltinRNG( seed );
        VectD meantries;
        for ( auto e : sp.getEDivKTGrid() ) {
          //Just below the grid point (e*kT/kT might round up to the next one):
          const NeutronEnergy ekin{ e * sp.kT() * ( 1.0 - 1e-12 ) };
          StableSum sum;
          for ( std::uint64_t i = 0; i < nperpt; ++i )
            sum.add( sp.sampleScatterAlphaBeta( rng, ekin ).ntries );
          meantries.push_back( nperpt ? sum.sum() / nperpt : 0.0 );
        }
        os << ",\"mean_tries\":";
        streamJSON(os,meantries);
        os << '}';
      }

      void query_impl_refsample( std::ostream& os,
                                 shared_obj<const SABData> sab,
                                 Optional<Temperature> teff,
                                 std::uint64_t nsample,
                                 std::uint64_t seed,
                                 const std::vector<NeutronEnergy>& sample_ekins,
                                 int knllux )
      {
        //Extend beyond Emax per cfg, like the actual factories (the
        //default RefSampleExtension extender is free-gas):
        if ( knllux < 0 )
          knllux = SABCfg::sablux_default_luxury;
        const auto cfg = SABCfg::createConfig(knllux);
        SABRef::RefSampleExtension ext;
        const bool use_sct
          = ( cfg.extender_model == SABCfg::ExtenderModel::SCT
              && teff.has_value() );
        if ( use_sct )
          ext.extender = makeSO<SAB::SABSCTExtender>( sab->temperature(),
                                                      teff.value(),
                                                      sab->elementMassAMU(),
                                                      SigmaBound{1.0},
                                                      cfg.sct_table_npts );
        //Expensive, so once for all energies:
        ext.emax = SABRef::autodetectRefSampleEmax( sab );
        //One sample per energy (a list if more than one), each with a new RNG
        //from the seed:
        const bool many = sample_ekins.size() > 1;
        os << "{\"refsample\":";
        if ( many )
          os << '[';
        for ( auto i : ncrange( sample_ekins.size() ) ) {
          const NeutronEnergy sample_ekin = sample_ekins.at(i);
          auto rng = createBuiltinRNG( seed );
          const double E_div_kT = sample_ekin.dbl() / sab->temperature().kT();
          auto ab = SABRef::refSampleAlphaBeta( rng, sab, EDivKT{E_div_kT},
                                                nsample, ext );
          os << ( i ? "," : "" ) << "{\"extender\":";
          streamJSON( os, use_sct ? "SCT" : "FreeGas" );
          os << ",\"E\":";
          streamJSON(os,sample_ekin);
          os << ",\"E_div_kT\":";
          streamJSON(os,E_div_kT);
          os << ",\"kT\":";
          streamJSON(os,sab->temperature().kT());
          os << ",\"seed\":";
          streamJSON(os,seed);
          os << ",\"nsample\":";
          streamJSON(os,nsample);
          os << ",\"alpha\":";
          streamJSONHugeDblVect(os,std::move(ab.first));
          os << ",\"beta\":";
          streamJSONHugeDblVect(os,std::move(ab.second));
          os << '}';
        }
        if ( many )
          os << ']';
        os << '}';
      }

      void query_impl_legacysample( std::ostream& os,
                                    shared_obj<const SABData> sab,
                                    std::uint64_t nsample,
                                    std::uint64_t seed,
                                    NeutronEnergy sample_ekin,
                                    bool legacy_oversample )
      {
        auto rng = createBuiltinRNG( seed );
        auto scathelper = SAB::
          createScatterHelperWithCache( sab, nullptr,
                                        ( legacy_oversample
                                          ? SAB::LegacySABAlgOpts::OVERSAMPLE10
                                          : SAB::LegacySABAlgOpts::DEFAULT ) );
        auto& sampler = scathelper->sampler;
        const double E_div_kT = sample_ekin.dbl() / sab->temperature().kT();
        VectD alpha, beta;
        alpha.reserve( countToSizeT( nsample ) );
        beta.reserve( countToSizeT( nsample ) );
        for ( std::uint64_t i = 0; i < nsample; ++i ) {
          auto ab = sampler.sampleAlphaBeta( sample_ekin, rng );
          alpha.push_back( ab.first );
          beta.push_back( ab.second );
        }
        os << "{\"legacysample\":{\"E\":";
        streamJSON(os,sample_ekin);
        os << ",\"E_div_kT\":";
        streamJSON(os,E_div_kT);
        os << ",\"kT\":";
        streamJSON(os,sab->temperature().kT());
        os << ",\"seed\":";
        streamJSON(os,seed);
        os << ",\"nsample\":";
        streamJSON(os,nsample);
        os << ",\"alpha\":";
        streamJSONHugeDblVect(os,std::move(alpha));
        os << ",\"beta\":";
        streamJSONHugeDblVect(os,std::move(beta));
        os << "}}";
      }


      void query_impl_sglcell( std::ostream& os, EDivKT e,
                               double a1, double a2, double b1, double b2,
                               double s11, double s12, double s21, double s22,
                               std::uint64_t nsample, std::uint64_t seed )
      {
        const std::uint64_t ntries_max = nsample*100;
        PairDD alpha(a1,a2), beta(b1,b2);
        VectD alpha_v = {a1,a2};
        VectD beta_v = {b1,b2};
        SABSurveyor surv( alpha_v, beta_v );
        SABCellSurvey cellsurv( a1, a2, b1, b2, e );
        double cellinteg_full;
        SmallVector<std::pair<StrView,double>,32> cellinteg_pb_list;
        double cellinteg_pb_chosen_for_fcsample = -1.0;
        CellData cell;//fixme: use this in all relevant interfaces?
        {
          cell.a1 = a1;
          cell.a2 = a2;
          cell.b1 = b1;
          cell.b2 = b2;
          cell.S[0] = s11;
          cell.S[1] = s12;
          cell.S[2] = s21;
          cell.S[3] = s22;
          for ( auto i : ncrange(4) )
            cell.logS[i] = NC::SABUtils::safeLogOrElse( cell.S[i], 0.0 );
        }
        {
          StableSumKahan sum_full;
          namespace SCI = StdLogLinCellIntegrator;
          SCI::integrateFullCell(cell,sum_full);
          cellinteg_full = sum_full.sum();

          StrView sv_allschemes(SABCfg::allIntegSchemesAsStr());
          for ( auto sv_scheme : sv_allschemes.split(';') ) {
            auto scheme = SABCfg::str2IntegScheme( sv_scheme );
            StableSumKahan sum_pb;
            SCI::integrateWithinKB( cell, e, scheme, sum_pb);
            cellinteg_pb_list.emplace_back(sv_scheme,sum_pb.sum());
            if ( scheme == SABCfg::IntegrationScheme::Flex5 )
              cellinteg_pb_chosen_for_fcsample = sum_pb.sum();
          }
        }

        const double predicted_fc_ar
          = ( cellinteg_full > 0.0
              ? cellinteg_pb_chosen_for_fcsample/cellinteg_full
              : -1.0 );

        nc_assert_always(surv.data().size()==1);
        const double E_div_kT_touch = surv.data().front().e_touch;
        SampleResult samples_fc, samples_ref;
        Optional<double> prob1_fc;
        auto rng = createBuiltinRNG( seed );
        if ( nsample>0 && e.dbl()>E_div_kT_touch ) {
          samples_ref = produceRefCellSamples( e, rng, cell,
                                               nsample, ntries_max );
          nc_assert_always( samples_ref.a.size() == nsample );
          nc_assert_always( samples_ref.b.size() == nsample );
          //Full-cell sampling only when not too inefficient:
          if ( predicted_fc_ar > 0.001 )
            samples_fc = testFullCellOverlaySample( e, rng, cell,
                                                    nsample, ntries_max );
          prob1_fc = samples_fc.prob1;
        }
        os<<"{\"alpha\":";
        streamJSON(os,alpha);
        os<<",\"beta\":";
        streamJSON(os,beta);
        os<<",\"S\":";
        SmallVector<double,4> s_v = {s11, s12, s21, s22};
        streamJSON(os,s_v);
        os<<",\"surveyor\":{\"E_div_kT_touch\":";
        streamJSON(os,surv.data().front().e_touch);
        os<<",\"E_div_kT_cover\":";
        streamJSON(os,surv.data().front().e_cover);

        {
          //Fixme: SABCellEval is not trustworthy -> replace eventually with new
          //code!
          double s[4] = {s11, s12, s21, s22};//fixme check order
          using CellEval = SABCellEval<InterpolationScheme::LOGLIN,
                                       SABInterpolationOrder::ALPHA_FIRST>;
          CellEval ce( alpha, beta, s );
          os<<"},\"celleval_OBSOLETE\":{\"full_integral\":";//fixme: remove this obsolete class
          streamJSON(os,ce.integral());
          os<<",\"phasespace_integral\":";
          streamJSON(os,ce.integralWithinKinematicBounds( e ));
          os<<",\"phasespace_E_div_kT\":";
          streamJSON(os,e.dbl());
        }
        os<<"},\"cellintegral\":{\"full_integral\":";
        streamJSON(os,cellinteg_full);
        os<<",\"phasespace_E_div_kT\":";
        streamJSON(os,e.dbl());
        os <<",\"integration_regions\":";
        cellsurv.toJSON(os);
        os<<",\"phasespace_integral_format\":[\"scheme\",\"value\"]";
        os<<",\"phasespace_integral\":[";
        bool first=true;
        for ( auto& sv : cellinteg_pb_list ) {
          if (!first)
            os<<',';
          first = false;
          os<<'[';
          streamJSON(os,sv.first);
          os<<',';
          streamJSON(os,sv.second);
          os<<']';
        }
        os <<"]},\"sampling\":{";
        os <<"\"fc_predicted_AR\":";
        streamJSON(os,predicted_fc_ar);
        os <<",\"fc_sampled_alpha\":";
        streamJSONHugeDblVect(os,std::move(samples_fc.a));
        os <<",\"fc_sampled_beta\":";
        streamJSONHugeDblVect(os,std::move(samples_fc.b));
        os <<",\"fc_sampled_ntries\":";
        streamJSON(os,samples_fc.ntries);
        os <<",\"fc_prob_edge_1\":";
        streamJSON(os,prob1_fc);
        os <<",\"ref_sampled_alpha\":";
        streamJSONHugeDblVect(os,std::move(samples_ref.a));
        os <<",\"ref_sampled_beta\":";
        streamJSONHugeDblVect(os,std::move(samples_ref.b));
        os <<",\"ref_sampled_ntries\":";
        streamJSON(os,samples_ref.ntries);
        os <<"}}";
      }

      void query_impl_surveyor( std::ostream& os,
                                const VectD& alpha,
                                const VectD& beta )
      {
        SABSurveyor surv( alpha, beta );
        os<<"{\"alpha\":";
        streamJSON(os,alpha);
        os<<",\"beta\":";
        streamJSON(os,beta);
        //Each entry in a CellList is (E/kT,cell index):
        using CellList = std::vector<std::pair<double,SABSurveyor::cellidx_t>>;
        auto sortCellList = [](CellList& cl)
        {
          using E = CellList::value_type;
          std::sort( cl.begin(), cl.end(),
                     []( const E& a, const E& b )
                     {
                       if ( a.first != b.first )
                         return a.first < b.first;
                       return a.second.val < b.second.val;
                     });
        };

        auto streamCellList = [&os](const CellList& cl)
        {
          nc_assert_always(cl.size()>=1);
          os << '[';
          bool first = true;
          for ( auto& e : cl ) {
            if ( first )
              first = false;
            else
              os << ',';
            os<<"[";
            streamJSON(os,e.first);
            os<<',';
            streamJSON(os,e.second.unpackAlphaIdx());
            os<<',';
            streamJSON(os,e.second.unpackBetaIdx());
            os<<']';
          }
          os << ']';
        };
        os<<",\"list_format\":[\"E_div_kT\",\"cell_ialpha\",\"cell_ibeta\"]";
        CellList touchList, coverList;
        {
          touchList.reserve( surv.data().size() );
          for ( auto& cc : surv.data() )
            touchList.emplace_back( cc.e_touch, cc.cellidx );
          sortCellList( touchList );

          coverList.reserve( surv.data().size() );
          for ( auto& cc : surv.data() )
            coverList.emplace_back( cc.e_cover, cc.cellidx );
          sortCellList( coverList );
        }
        os<<",\"touch_list\":";
        streamCellList(touchList);
        os<<",\"cover_list\":";
        streamCellList(coverList);
        os<<'}';
      }

      void query_impl_SABRefEval( std::ostream& os,
                                  const std::string& matcfg,
                                  NeutronEnergy eval,
                                  Optional<std::string> atomDisplayLabel,
                                  std::uint64_t nsample )
      {
        //Fixme: is this obsolete? Should we remove ths option + the
        //NCSABRefEval header again?
        auto sabdata = query_impl_extractSabData( matcfg,
                                                  atomDisplayLabel ).sab;
        SABRefEval<> refeval( sabdata, eval );
        RNG* rngptr = nullptr;
        std::shared_ptr<RNGStream> rngholder;
        if ( nsample > 0 ) {
          rngholder = getRNG();
          rngptr = rngholder.get();
        }
        refeval.toJSON( os, nsample, rngptr );
      }

      class RNGSpy final : public NC::RNGStream {
        shared_obj<RNGStream> m_realrng;
        std::uint64_t m_count = 0;
      public:
        RNGSpy(shared_obj<RNGStream> rs) : m_realrng( std::move(rs) ) {}
        std::uint64_t count() const { return m_count; }
      protected:
        double actualGenerate() override {
          ++m_count;
          return m_realrng->generate(); }
      };

      void query_impl_samplexyparabolicband( std::ostream& os,
                                             double x0, double y0,
                                             double x1, double y1,
                                             std::uint64_t nsample )
      {
        RNGSpy rng( getRNG() );

        ParabolicBandBoxSampler sampler( x0, y0, x1, y1 );
        if (!sampler.canSample())
          nsample = 0;

        nc_assert_always(x0>=0);
        nc_assert_always(y0>=0);
        nc_assert_always(x1>x0);
        nc_assert_always(y1>y0);

        os << "{\"x0\":";
        streamJSON(os,x0);
        os << ",\"y0\":";
        streamJSON(os,y0);
        os << ",\"x1\":";
        streamJSON(os,x1);
        os << ",\"y1\":";
        streamJSON(os,y1);
        os << ",\"samples\":[";
        for ( std::uint64_t i = 0; i < nsample; ++i ) {
          auto res = sampler.sample(rng);
          if ( i )
            os << ',';
          streamJSON(os,res);
        }
        os << "],\"rng_per_sample\":";
        if ( nsample )
          streamJSON(os, double(rng.count())/nsample );
        else
          streamJSON(os, json_null_t{} );
        os << ",\"sampler_details\":";
        sampler.toJSON(os);
        os<<'}';
      }
    }
  }
}

void NC::SABUtils::JSONQuery( std::ostream& os, const Query& query )
{
  auto invalid = [&query](const char * reason){
    std::ostringstream ss;
    ss << "Invalid SAB JSON query: ";
    streamJSON( ss, query );//trick: using streamJSON for easy format.
    if ( reason )
      ss<< " ("<<reason<<')';
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  constexpr auto sv_sab = StrView::make("sab");
  if ( query.size() < 2 || query.front() != sv_sab ) {
    invalid(nullptr);
    return;
  }
  //shift off the "sab" and key entries:
  const auto& key = query.at(1).trimmed();
  auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
  auto argstr = [&arg]( std::size_t i ) { return arg(i).to_string(); };
  const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_refeval = StrView::make("refeval");
  constexpr auto sv_samplepb = StrView::make("samplepb");
  constexpr auto sv_surveyor = StrView::make("surveyor");
  constexpr auto sv_sglcell = StrView::make("sglcell");
  constexpr auto sv_integschemes = StrView::make("integschemes");
  constexpr auto sv_proc = StrView::make("proc");
  constexpr auto sv_proctries = StrView::make("proctries");
  constexpr auto sv_refsample = StrView::make("refsample");
  constexpr auto sv_analyse = StrView::make("analyse");
  constexpr auto sv_xs = StrView::make("xs");

  if ( key == sv_refeval ) {
    //query like: ncrystal_query sab refeval 0.025 'bla.ncmat' 1000 ['Al']
    if ( nargs != 4 && nargs != 3 )
      invalid("correct usage: [\"sab\",\"refeval\","
              "ENERGY_EV,MATCFGSTR,NSAMPLE,ATOMDISPLAYLABEL]"
              " (the ATOMDISPLAYLABEL can be left out for monoatomic materials)");
    auto opt_eval = arg(0).toDbl();
    if ( !opt_eval.has_value() || !(opt_eval.value()>0.0) )
      invalid("invalid ENERGY_EV value");
    auto opt_nsample = arg(2).toUInt64();
    if ( !opt_nsample.has_value() )
      invalid("invalid NSAMPLE value");
    Optional<std::string> atomDisplayLabel;
    if ( nargs == 4 )
      atomDisplayLabel = argstr(3);
    query_impl_SABRefEval( os,
                           argstr(1),
                           NeutronEnergy{ DoValidate, opt_eval.value() },
                           atomDisplayLabel,
                           opt_nsample.value() );
  } else if ( key == sv_samplepb ) {
    //sab samplepb x0 y0 x1 y1 nsample
    if ( nargs != 5 )
      invalid("correct usage: [\"sab\",\"samplepb\",X0,Y0,X1,Y1,NSAMPLE");
    auto opt_nsample = arg(4).toUInt64();
    if ( !opt_nsample.has_value() )
      invalid("invalid NSAMPLE value");
    auto s_x0 = arg(0).toDbl();
    auto s_y0 = arg(1).toDbl();
    auto s_x1 = arg(2).toDbl();
    auto s_y1 = arg(3).toDbl();
    if ( !s_x0.has_value() || !s_y0.has_value()
         || !s_x1.has_value() || !s_y1.has_value() )
      invalid("bad coordinate value");
    query_impl_samplexyparabolicband( os, s_x0.value(), s_y0.value(),
                                      s_x1.value(), s_y1.value(),
                                      opt_nsample.value() );
  } else if ( key == sv_surveyor ) {
    const char * usage = ( "correct usage: [\"sab\",\"surveyor\",ALPHAGRID,"
                           "BETAGRID], with grids encoded as '@val1@..@valn'" );
    if ( nargs != 2 )
      invalid(usage);
    auto decodeGrid = [](StrView v) -> VectD
    {
      if (!v.startswith('@'))
        return {};
      auto vs = v.substr(1).splitTrimmed('@');
      VectD res;
      res.reserve(vs.size());
      for ( auto& e : vs ) {
        auto val = e.toDbl();
        if (!val.has_value())
          return {};
        res.push_back(val.value());
      }
      return res;
    };
    const VectD alpha = decodeGrid( arg(0) );
    const VectD beta = decodeGrid( arg(1) );
    if ( beta.size() < 2 || alpha.size() < 2)
      invalid(usage);
    query_impl_surveyor( os, alpha, beta );
  } else if ( key == sv_sglcell ) {
    const char * usage = ( "correct usage: [\"sab\",\"sglcell\",EDIVKT,ALPHA1,"
                           "ALPHA2,BETA1,BETA2,SA1B1,SA2B1,SA1B2,SA2B2,"
                           "NSAMPLE,SEED], with negative beta values prefixed"
                           " with '@' and NSAMPLE+SEED being optional." );
    std::uint64_t nsample = 0;
    std::uint64_t seed = 123456;
    if ( nargs>=10 ) {
      auto optns = arg(9).toUInt64();
      if ( !optns.has_value() )
        invalid(usage);
      nsample = optns.value();
    }
    if ( nargs>=11 ) {
      auto optsd = arg(10).toUInt64();
      if ( !optsd.has_value() )
        invalid(usage);
      seed = optsd.value();
    }
    if (!arg(3).startswith('@')||!arg(4).startswith('@'))
      invalid(usage);
    SmallVector<Optional<double>,9> v;
    for ( auto i : ncrange(9)) {
      v.push_back( arg(i).startswith('@')
                   ? arg(i).substr(1).toDbl()
                   : arg(i).toDbl() );
    }
    for ( auto e : v ) {
      if (!e.has_value()||ncisnan(e.value())||!std::isfinite(e.value()))
        invalid(usage);
    }
    const double E_div_kT = v.at(0).value();
    const double a1 = v.at(1).value();
    const double a2 = v.at(2).value();
    const double b1 = v.at(3).value();
    const double b2 = v.at(4).value();
    const double s11 = v.at(5).value();
    const double s12 = v.at(6).value();
    const double s21 = v.at(7).value();
    const double s22 = v.at(8).value();
    if ( !(E_div_kT>0.0) || !(a1<a2) || !(a1>=0.0) || !(b2>b1)
         || !(s11>=0.0) || !(s12>=0.0)
         || !(s21>=0.0) || !(s22>=0.0) )
      invalid(usage);
    query_impl_sglcell( os, EDivKT{E_div_kT}, a1, a2, b1, b2,
                        s11, s12, s21, s22, nsample, seed );
  } else if ( key == sv_analyse ) {
    //query like: ["sab","analyse",MATCFGSTR,ATOMDISPLAYLABEL,"diag"]
    const char * usage
      = ( "correct usage: [\"sab\",\"analyse\",MATCFGSTR,"
          "ATOMDISPLAYLABEL,\"diag\"] (the ATOMDISPLAYLABEL can be left"
          " as an empty string for monoatomic materials, and the optional"
          " \"diag\" enables extra diagnostics)" );
    if ( nargs < 1 || nargs > 3 )
      invalid(usage);
    Optional<std::string> atomdsplbl;
    if ( nargs >= 2 && !arg(1).empty() )
      atomdsplbl = argstr(1);
    SABAnalyser::Options opt;
    if ( nargs == 3 ) {
      if ( arg(2) != "diag" )
        invalid(usage);
      opt.collect_diagnostics = true;
    }
    auto xd = query_impl_extractSabData( argstr(0), atomdsplbl );
    query_impl_analyse( os, xd, opt );
  } else if ( key == sv_proc || key == sv_refsample ) {
    //Almost same usage+parsing of proc/refsample:
    const bool is_proc = key==sv_proc;
    const char * usage
      = ( is_proc
          ? ( "correct usage: [\"sab\",\"proc\",MATCFGSTR,ATOMDISPLAYLABEL,"
              "NSAMPLE,SEED,SAMPLE_EVAL,EGRID0,..,EGRIDN]"
              " (the EGRID pts are optional, SAMPLE_EVAL can be a comma"
              " separated list, and the ATOMDISPLAYLABEL can"
              " be left as an empty string for monoatomic materials)" )
          : ( "correct usage: [\"sab\",\"refsample\",MATCFGSTR,"
              "ATOMDISPLAYLABEL,NSAMPLE,SEED,SAMPLE_EVAL,TYPE] where TYPE is one"
              " of \"legacy\", \"legacy_oversample\", and \"ref\" (SAMPLE_EVAL"
              " can be a comma separated list for \"ref\", and the"
              " ATOMDISPLAYLABEL can be left as an empty string for monoatomic"
              " materials)" ) );
    if ( nargs < 5 )
      invalid(usage);

    Optional<std::string> atomdsplbl;
    if ( !arg(1).empty() )
      atomdsplbl = argstr(1);
    auto extracted = query_impl_extractSabData( argstr(0), atomdsplbl );
    auto sabdata = std::move( extracted.sab );
    const auto teff = extracted.teff;
    const int knllux = MatCfg(argstr(0)).get_knllux();

    auto optns = arg(2).toUInt64();
    if ( !optns.has_value() )
      invalid(usage);
    std::uint64_t nsample = optns.value();
    auto optsd = arg(3).toUInt64();
    if ( !optsd.has_value() )
      invalid(usage);
    std::uint64_t seed = optsd.value();
    //SAMPLE_EVAL can be a comma separated list of energies:
    std::vector<NeutronEnergy> sample_ekins;
    for ( auto& part : arg(4).split(',') ) {
      auto v = part.toDbl();
      if ( !v.has_value() || !( v.value() >= 0.0 ) )
        invalid(usage);
      sample_ekins.emplace_back( DoValidate_t{}, v.value() );
    }
    if ( sample_ekins.empty() )
      invalid(usage);
    VectD egrid;
    bool refsample_type_is_legacy(false);
    bool legacy_oversample(false);
    if ( !is_proc ) {
      if ( nargs != 6 )
        invalid(usage);
      //parse TYPE
      if ( arg(5) == "ref" ) {
        refsample_type_is_legacy = false;
      } else if ( arg(5) == "legacy" ) {
        refsample_type_is_legacy = true;
      } else if ( arg(5) == "legacy_oversample" ) {
        refsample_type_is_legacy = true;
        legacy_oversample = true;
      } else {
        invalid(usage);
      }
    } else if ( nargs > 5 ) {
      if (!is_proc)
        invalid(usage);
      egrid.reserve(static_cast<std::size_t>(nargs-5));
      for ( std::size_t i = 5; i < nargs; ++i ) {
        auto val = arg(i).toDbl();
        if (!val.has_value()||!(val.value()>0.0))
          invalid(usage);
        egrid.push_back(val.value());
      }
    }
    if ( key == sv_proc ) {
      query_impl_proc( os, std::move(sabdata), teff, std::move(egrid),
                       nsample, seed, sample_ekins, knllux );
    } else {
      if ( refsample_type_is_legacy ) {
        if ( sample_ekins.size() != 1 )
          invalid(usage);
        query_impl_legacysample( os, std::move(sabdata), nsample, seed,
                                 sample_ekins.front(), legacy_oversample );
      } else {
        query_impl_refsample( os, std::move(sabdata), teff, nsample, seed,
                              sample_ekins, knllux );
      }
    }
  } else if ( key == sv_proctries ) {
    const char * usage
      = ( "correct usage: [\"sab\",\"proctries\",MATCFGSTR,"
          "ATOMDISPLAYLABEL,NPERPT,SEED] (the ATOMDISPLAYLABEL can be left as"
          " an empty string for monoatomic materials)" );
    if ( nargs != 4 )
      invalid(usage);
    Optional<std::string> atomdsplbl;
    if ( !arg(1).empty() )
      atomdsplbl = argstr(1);
    auto extracted = query_impl_extractSabData( argstr(0), atomdsplbl );
    auto optn = arg(2).toUInt64();
    auto optsd = arg(3).toUInt64();
    if ( !optn.has_value() || !optsd.has_value() )
      invalid(usage);
    query_impl_proctries( os, std::move( extracted.sab ), optn.value(),
                          optsd.value(), MatCfg(argstr(0)).get_knllux() );
  } else if ( key == sv_xs ) {
    const char * usage
      = ( "correct usage: [\"sab\",\"xs\",MATCFGSTR,ATOMDISPLAYLABEL,"
          "EXTENDER,E0,..,EN] where EXTENDER is \"sct\" or \"freegas\""
          " and the E values are neutron energies in eV (the ATOMDISPLAYLABEL"
          " can be left as an empty string for monoatomic materials)" );
    if ( nargs < 4 )
      invalid(usage);
    SABCfg::ExtenderModel extender;
    if ( arg(2) == "sct" )
      extender = SABCfg::ExtenderModel::SCT;
    else if ( arg(2) == "freegas" )
      extender = SABCfg::ExtenderModel::FreeGas;
    else
      invalid(usage);
    VectD energies;
    energies.reserve( static_cast<std::size_t>( nargs - 3 ) );
    for ( std::size_t i = 3; i < nargs; ++i ) {
      auto val = arg(i).toDbl();
      if ( !val.has_value() || !( val.value() >= 0.0 ) )
        invalid(usage);
      energies.push_back( val.value() );
    }
    Optional<std::string> atomdsplbl;
    if ( !arg(1).empty() )
      atomdsplbl = argstr(1);
    query_impl_xs( os, query_impl_extractSabData( argstr(0), atomdsplbl ),
                   extender, energies, MatCfg(argstr(0)).get_knllux() );
  } else if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"sab\",\"list\"]");
    streamJSON( os, std::array<StrView,10>{ sv_analyse,
                                            sv_integschemes,
                                            sv_proc,
                                            sv_proctries,
                                            sv_refeval,
                                            sv_refsample,
                                            sv_samplepb,
                                            sv_sglcell,
                                            sv_surveyor,
                                            sv_xs } );
  } else {
    invalid(nullptr);
  }
}
