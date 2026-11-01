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

#include "NCrystal/internal/phys_utils/NCFreeGasSampler.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCFreeGasSampler_FMA.hh"

namespace NC = NCrystal;

NC::FreeGasSampler::FreeGasSampler( Temperature t, AtomMass m )
  : m_kT( t.kT() ),
    m_invkT( 1.0 / m_kT ),
    m_A( m.relativeToNeutronMass() ),
    m_k( 1.0 / ( 1.0 + m_A ) )
{
  t.validate();
  m.validate();
  nc_assert_always( m_kT > 0.0 && m_A > 0.0 );
}

NC::PairDD NC::FreeGasSampler::sampleAlphaBeta( RNG& rng,
                                                NeutronEnergy ekin ) const
{
  //Sampling in NCFreeGasSampler_FMA.hh (for FMA dispatch):
  const EDivKT e{ ekin.dbl() * m_invkT };
  double res[2];
  freeGasSampleColl( rng, m_A, m_k, e.dbl(), false, res );
  const double alpha = res[0], beta = res[1];
  //Guard against rounding outside the kinematic limits [a-,a+], with
  //(a-2c-beta)^2=4*c*cp (c=E/kT, cp=c+beta) at a=a+-. Skip the costly (sqrt)
  //clamp if alpha is clearly inside (margin far above rounding errors,
  //cp>=1e-6*c keeps them small), as the clamp is then a no-op:
  const double cp = e.dbl() + beta;
  if ( cp >= 1e-6 * e.dbl() ) {
    const double d = alpha - ( 2.0 * e.dbl() + beta );//no fma (2*c exact)
    if ( d * d < ( 4.0 * e.dbl() ) * cp * ( 1.0 - 1e-9 ) ) {
#ifndef NDEBUG
      auto alim = getAlphaLimits( e, beta );
      nc_assert_always( alpha > alim.first && alpha < alim.second );
#endif
      return { alpha, beta };
    }
  }
  auto alim = getAlphaLimits( e, beta );
  return { ncclamp( alpha, alim.first, alim.second ), beta };
}

NC::PairDD NC::FreeGasSampler::sampleDeltaEMu( RNG& rng,
                                               NeutronEnergy ekin ) const
{
  double res[3];
  freeGasSampleColl( rng, m_A, m_k, ekin.dbl() * m_invkT, true, res );
  return { res[1] * m_kT, res[2] };
}
