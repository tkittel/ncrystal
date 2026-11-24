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

#include "NCrystal/internal/extd_utils/NCRotatedScatter.hh"

namespace NC = NCrystal;

NC::RotatedScatter::RotatedScatter( ProcImpl::ProcPtr proc,
                                    const RotMatrix& rot )
  : m_proc( std::move(proc) ),
    m_rot( Matrix( MatrixAllowCopy, rot ) ),
    m_invrot( ~rot )
{
  nc_assert_always( m_proc->processType() == ProcessType::Scatter );
  nc_assert_always( ncabs( rot.determinant() - 1.0 ) < 1e-6 );
}

NC::EnergyDomain NC::RotatedScatter::domain() const noexcept
{
  return m_proc->domain();
}

NC::CrossSect
NC::RotatedScatter::crossSection( CachePtr& cp, NeutronEnergy ekin,
                                  const NeutronDirection& dir ) const
{
  return m_proc->crossSection( cp, ekin, ( m_invrot * dir.as<Vector>()
                                           ).as<NeutronDirection>() );
}

NC::ScatterOutcome
NC::RotatedScatter::sampleScatter( CachePtr& cp, RNG& rng, NeutronEnergy ekin,
                                   const NeutronDirection& dir ) const
{
  const auto d = ( m_invrot * dir.as<Vector>() ).as<NeutronDirection>();
  auto res = m_proc->sampleScatter( cp, rng, ekin, d );
  if ( res.direction == d )
    res.direction = dir;//unscattered, preserve exactly
  else
    res.direction = ( m_rot * res.direction.as<Vector>()
                      ).as<NeutronDirection>();
  return res;
}

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
bool NC::RotatedScatter::isPureElasticScatter() const
{
  return m_proc->isPureElasticScatter();
}
#endif

NC::Optional<std::string> NC::RotatedScatter::specificJSONDescription() const
{
  //Same layout as for ProcComposition, so the wrapped process shows up as a
  //sub-component in process summaries.
  std::ostringstream ss;
  ss << "{\"components\":[[1.0," << m_proc->jsonDescription() << "]]}";
  return ss.str();
}
