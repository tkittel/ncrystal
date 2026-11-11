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

#include "NCrystal/internal/sabscatter/NCSABScatter.hh"
#include "NCrystal/internal/sab/NCSABExtended.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
namespace NC = NCrystal;

NC::SABScatter::~SABScatter() = default;

NC::SABScatter::SABScatter( shared_obj<const SABExtended> sh, double scale )
  : m_sh(std::move(sh)), m_scale(scale)
{
  nc_assert( m_sh!=nullptr);
  if ( !( m_scale > 0.0 && std::isfinite( m_scale ) ) )
    NCRYSTAL_THROW2(BadInput,"SABScatter: invalid scale (must be positive"
                    " and finite): "<<m_scale);
}

NC::SABScatter::SABScatter( shared_obj<const SABExtended> sh,
                            SigmaBound sb, double scale )
  : SABScatter( std::move(sh), sb.dbl() * scale )
{
}

NC::CrossSect
NC::SABScatter::crossSectionIsotropic( CachePtr&, NeutronEnergy ekin ) const
{
  return CrossSect{ m_sh->xsEvaluator().evalScaled( ekin.dbl(), 1.0, m_scale ) };
}

NC::Optional<NC::IsoXSTable>
NC::SABScatter::createIsoXSTable() const
{
  return scaleIsoXSTable( m_sh->xsEvaluator().table(), 1.0, m_scale );
}

NC::VectD NC::SABScatter::createIsoXSBreakpoints() const
{
  return m_sh->xsBreakpoints();
}

NC::ScatterOutcomeIsotropic
NC::SABScatter::sampleScatterIsotropic( CachePtr&,
                                        RNG& rng, NeutronEnergy ekin ) const
{
  return m_sh->sampleScatter( rng, ekin );
}

std::shared_ptr<NC::ProcImpl::Process>
NC::SABScatter::createMerged( const Process& oraw,
                              double scale_self,
                              double scale_other ) const
{
  nc_assert(scale_other>0.0);
  nc_assert(scale_self>0.0);
  auto optr = dynamic_cast<const SABScatter*>(&oraw);
  if (!optr)
    return nullptr;
  auto& o = *optr;

  if ( m_sh != o.m_sh )//fixme: it is actually enough that SABProcessor+SABExtender are compatible
    return nullptr;

  const double newscale = scale_self*m_scale + scale_other*o.m_scale;
  return std::make_shared<SABScatter>( m_sh, newscale );
}

NC::Optional<std::string> NC::SABScatter::specificJSONDescription() const
{
  std::string extension_method
    = m_sh->extender().shortName().value_or("custom");
  std::ostringstream ss;
  m_sh->processor().toJSONProcessInfo( ss,
                                       SigmaBound{ m_scale },
                                       extension_method );
  return ss.str();
}
