
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

#include "NCrystal/internal/elincscatter/NCElIncScatter.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/phys_utils/NCElIncXS.hh"
#include "NCrystal/internal/phys_utils/NCDebyeMSD.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCStableDbl.hh"
#include <sstream>

namespace NC = NCrystal;

NC::ElIncScatter::~ElIncScatter() = default;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    struct ElIncInfo {
      VectD msd, bixs, scale;
    };

    enum class ForCapabilityCheckOnly{ Yes, No };
    Optional<ElIncInfo> extractInfo( const Info& info,
                                     const ElIncScatterCfg& cfg,
                                     ForCapabilityCheckOnly fcco = ForCapabilityCheckOnly::No )
    {
      const bool capOnly = ( fcco == ForCapabilityCheckOnly::Yes );
      if ( capOnly ) {
        if ( cfg.scale_factor == 0.0 )
          return NullOpt;
        if ( ! ( cfg.use_sigma_incoherent || cfg.use_sigma_coherent ) )
          return NullOpt;
        if ( !info.hasTemperature() )
          return NullOpt;
      } else {
        nc_assert_always( !(cfg.scale_factor<=0.0) );
        nc_assert_always( cfg.use_sigma_incoherent || cfg.use_sigma_coherent );
        if ( !info.hasTemperature() )
          NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter lacks temperature.");
      }

      ElIncInfo res;
      SmallVector<double,6> msd;
      SmallVector<double,6> bixs;
      SmallVector<double,6> scale;

      auto getSigma = [&cfg](const AtomData& ad)
      {
        return (  cfg.use_sigma_incoherent ? ad.incoherentXS().get() : 0.0 )
          + ( cfg.use_sigma_coherent ? ad.coherentXS().get() : 0.0 );
      };


      //Prefer initialising via atom infos (in this case we require all to contribute):
      if ( info.hasAtomInfo() ) {

        //Average number of atoms per unit cell (positions*occupancy):
        const double ntot = info.getStructureInfo().n_atoms_average;

        for ( auto& ai : info.getAtomInfos() ) {
          //scale factor + cross section (per actual atom). Randomly placed
          //vacancies make part of the coherent scattering incoherent: the
          //disorder term 4pi*o(1-o)*b^2 per position, i.e. (1-o)*sigma_coh
          //per atom. Not added if the full sigma_coh is already included
          //(use_sigma_coherent, i.e. coherent elastic scattering treated in
          //the incoherent approximation, which already contains it):
          const double o = ai.occupancy();
          scale.push_back( ( StableDbl( double(ai.numberOfPositions()) ) * o
                             * cfg.scale_factor ).value() / ntot );
          const double disorder
            = ( cfg.use_sigma_incoherent && !cfg.use_sigma_coherent
                ? ( ( StableDbl( 1.0 ) - StableDbl( o ) )
                    * ai.atomData().coherentXS().get() ).value()
                : 0.0 );
          bixs.push_back( getSigma(ai.atomData()) + disorder );
          //msd:
          if ( ai.msd().has_value() ) {
            msd.push_back( ai.msd().value() );
          } else {
            //Fall-back to calculating MSDs from the isotropic Debye model. Eventually
            //we would like to avoid this here, and make sure this is done on the Info
            //object itself.
            if ( !ai.debyeTemp().has_value() ) {
              if ( capOnly )
                return NullOpt;
              else
                NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter has AtomInfo object without "
                               "mean-square-displacements (MSD), and there is not enough information to"
                               " estimate one (a Debye temperature + material temperature is required).");
            }
            auto debyeTemp = ai.debyeTemp().value();
            auto temperature = info.getTemperature();
            auto atomMass = ai.atomData().averageMassAMU();
            nc_assert(debyeTemp.get()>0.0&&temperature.get()>0.0&&atomMass.get()>0.0);
            msd.push_back( debyeIsotropicMSD( debyeTemp, temperature, atomMass ) );
          }
        }
      } else {
        //Try to initialise via dyninfo sections (ok if some, but not all, have missing info):
        std::string missing_msd_elems;
        for ( auto& di : info.getDynamicInfoList() ) {
          auto di_knl = dynamic_cast<const DI_ScatKnl*>(di.get());
          if ( !di_knl )
            continue;//freegas/sterile: never an msd carrier
          Optional<double> msd_value = extractMSDFromDynInfo( di_knl );
          if ( msd_value.has_value() ) {
            msd.push_back( msd_value.value() );
            scale.push_back( di->fraction() * cfg.scale_factor );
            bixs.push_back( getSigma( di->atomData() ) );
          } else {
            if ( !missing_msd_elems.empty() )
              missing_msd_elems += ", ";
            missing_msd_elems += di->atomData().description(false);
          }
        }
        if ( !capOnly && !msd.empty() && !missing_msd_elems.empty() ) {
          //Partial contribution for a solid: arguably better than refusing flat
          //out, but warn the user.
          NCRYSTAL_WARN("Incoherent-elastic scattering of solid material"
                        " will be missing the contribution of component(s)"
                        " without available mean-squared displacements: "
                        <<missing_msd_elems);
        }
        if ( msd.empty() ) {
          if ( capOnly )
            return NullOpt;
          else
            NCRYSTAL_THROW(MissingInfo,
                           "Info object passed to ElIncScatter lacks"
                           " information to create Debye-Waller factors.");
        }
      }

      if ( capOnly ) {
        //Return non-NullOpt if at least one component has actual contribution:
        for ( auto i : ncrange(bixs.size()) )
          if ( bixs.at(i) * scale.at(i) > 0 )
            return res;
        return NullOpt;
      } else {
        //Copy and return actual results.
        res.msd.reserve( msd.size() );
        res.bixs.reserve( bixs.size() );
        res.scale.reserve( scale.size() );
        res.msd.insert(res.msd.end(), msd.begin(), msd.end() );
        res.bixs.insert(res.bixs.end(), bixs.begin(), bixs.end() );
        res.scale.insert(res.scale.end(), scale.begin(), scale.end() );
        return res;
      }
    }

    std::unique_ptr<ElIncXS> createElIncXS( const Info& info,
                                            const ElIncScatterCfg& cfg )
    {
      auto res = extractInfo( info, cfg );
      if (!res.has_value())
        NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter lacks information to create Debye-Waller factors.");
      return ncmake_unique<ElIncXS>( std::move(res.value().msd),
                                     std::move(res.value().bixs),
                                     std::move(res.value().scale) );
    }

    //Sum of the tabulations of the elements:
    IsoXSTable tabulateElInc( const ElIncXS& ei )
    {
      const auto elms = ei.elementData();
      if ( elms.empty() )
        return ElIncXS::tabulateMonoAtomic( 0.0, SigmaBound{ 0.0 } );
      if ( elms.size() == 1 )
        return ElIncXS::tabulateMonoAtomic( elms[0].first,
                                            SigmaBound{ elms[0].second } );
      std::vector<IsoXSTable> tables;
      tables.reserve( elms.size() );
      std::vector<std::pair<double,const IsoXSTable*>> inputs;
      inputs.reserve( elms.size() );
      for ( auto& e : elms )
        tables.push_back( ElIncXS::tabulateMonoAtomic( e.first,
                                                       SigmaBound{ e.second } ) );
      for ( auto& t : tables )
        inputs.emplace_back( 1.0, &t );
      return mergeIsoXSTables( inputs );
    }

    std::vector<IsoXSEvaluator> elementEvaluators( const ElIncXS& ei )
    {
      std::vector<IsoXSEvaluator> res;
      const auto elms = ei.elementData();
      if ( elms.size() < 2 )
        return res;
      res.reserve( elms.size() );
      for ( auto& e : elms )
        res.emplace_back( ElIncXS::tabulateMonoAtomic( e.first,
                                                       SigmaBound{ e.second } ) );
      return res;
    }

  }
}

bool NC::ElIncScatter::hasSufficientInfo( const Info& info, const ElIncScatterCfg& cfg )
{
  auto res = extractInfo( info, cfg, ForCapabilityCheckOnly::Yes );
  return res.has_value();
}


NC::ElIncScatter::ElIncScatter( const Info& info, const ElIncScatterCfg& cfg )
  : ElIncScatter( createElIncXS( info, cfg ) )
{
}

NC::ElIncScatter::ElIncScatter( const VectD& elements_meanSqDisp,
                                const VectD& elements_boundincohxs,
                                const VectD& elements_scale )
  : ElIncScatter( ncmake_unique<ElIncXS>( elements_meanSqDisp,
                                          elements_boundincohxs,
                                          elements_scale ) )
{
}

NC::CrossSect NC::ElIncScatter::crossSectionIsotropic( CachePtr&, NeutronEnergy ekin ) const
{
  return CrossSect{ m_xs.eval( ekin.dbl() ) };
}

NC::CosineScatAngle NC::ElIncScatter::sampleMu( RNG& rng, NeutronEnergy ekin ) const
{
  //Pick element by the values of their tabulations:
  const auto elms = m_elincxs->elementData();
  nc_assert( m_elmxs.size() == ( elms.size() > 1 ? elms.size() : 0 ) );
  std::size_t idx = 0;
  if ( m_elmxs.size() > 1 ) {
    SmallVector<double,32> commul;
    StableSumKahan sum;
    for ( auto& ev : m_elmxs ) {
      sum.add( ev.eval( ekin.dbl() ) );
      commul.push_back( sum.sum() );
    }
    idx = RandIdxByWeightPicker::oneShotPick( rng, commul );
  }
  nc_assert( idx < elms.size() );
  return ElIncXS::sampleMuMonoAtomic( rng, ekin, elms[idx].first );
}

NC::Optional<NC::IsoXSTable>
NC::ElIncScatter::createIsoXSTable() const
{
  return cloneIsoXSTable( m_xs.table() );
}

NC::ScatterOutcomeIsotropic NC::ElIncScatter::sampleScatterIsotropic( CachePtr&, RNG& rng, NeutronEnergy ekin ) const
{
  return { ekin, sampleMu( rng, ekin ) };
}

NC::ElIncScatter::ElIncScatter( std::unique_ptr<ElIncXS> p )
  : m_elincxs( std::move(p) ),
    m_xs( tabulateElInc( *m_elincxs ) ),
    m_elmxs( elementEvaluators( *m_elincxs ) )
{
}

std::shared_ptr<NC::ProcImpl::Process> NC::ElIncScatter::createMerged( const Process& oraw,
                                                                       double scale_self,
                                                                       double scale_other ) const
{
  auto optr = dynamic_cast<const ElIncScatter*>(&oraw);
  if (!optr)
    return nullptr;
  auto& o = *optr;
  nc_assert( m_elincxs != nullptr );
  nc_assert( o.m_elincxs != nullptr );
  return std::make_shared<ElIncScatter>(ncmake_unique<ElIncXS>( *m_elincxs, scale_self,
                                                                *o.m_elincxs, scale_other ));
}

NC::Optional<std::string> NC::ElIncScatter::specificJSONDescription() const
{
  std::ostringstream ss;
  auto xs_lowE = m_elincxs->evaluate( NeutronEnergy{ 0.0 } );
  auto nelem = m_elincxs == nullptr ? 0 : m_elincxs->nElements();
  {
    std::ostringstream tmp;
    tmp << "nelements="<<nelem;
    tmp << ";max_contrib="<<xs_lowE;
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "sigma_lowE_limit", xs_lowE.dbl() );
  streamJSONDictEntry( ss, "nelements", nelem, JSONDictPos::LAST );
  return ss.str();
}

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
void NC::ElIncScatter::evalManyXSIsotropic( CachePtr&, const double* ekin, std::size_t N, double* out_xs ) const
{
  for ( auto i : ncrange( N ) )
    out_xs[i] = m_xs.eval( ekin[i] );
}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NC::ElIncScatter::evalXSAndSampleScatterIsotropic( CachePtr&, RNG& rng, NeutronEnergy ekin ) const
{
  return { CrossSect{ m_xs.eval( ekin.dbl() ) }, { ekin, sampleMu( rng, ekin ) } };
}

#endif
