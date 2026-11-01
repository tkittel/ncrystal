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

NC::FreeGasSampler::Coll
NC::FreeGasSampler::sampleColl( RNG& rng, double c ) const
{
  //Work with velocities in units of sqrt(2kT/M): neutron n=(0,0,y), target
  //T=x*(sT,0,muT) (azimuth irrelevant). Sample (x,muT) with density
  //~|n-T|*x^2*exp(-x^2), then rotate the relative velocity r=n-T isotropically
  //(CM frame). The neutron velocity change is then A*k*g with g=|r|*Om-r,
  //giving alpha=A*k^2*g^2, beta=2*y*k*g_z+alpha (no large cancellations).
  //
  //Contraction safety: every a*b+c is an explicit ncfma (or exact, like
  //1-2*u), so results do not depend on compiler FP-contraction settings.
  nc_assert( c >= 0.0 );
  constexpr double halfsqrtpi = 0.5 * kSqrtPi;
  const double y = std::sqrt( m_A * c );
  const double p_gamma2 = 1.0 / ncfma( halfsqrtpi, y, 1.0 );
  const double fourY = 4.0 * y;
  double x, u, r2;
  while ( true ) {
    double x2;
    if ( rng.generate() < p_gamma2 ) {
      x2 = -std::log( rng.generate() * rng.generate() );//Gamma(2)
    } else {
      const double g = randNorm( rng );
      x2 = ncfma( 0.5 * g, g, randExp( rng ) );//Gamma(3/2)
    }
    x = std::sqrt( x2 );
    u = rng.generate();//muT = 1-2u
    const double ymx = y - x;
    r2 = ncfma( ymx, ymx, ( fourY * x ) * u );//|n-T|^2 (stable form)
    if ( ncsquare( rng.generate() * ( y + x ) ) <= r2 )
      break;
  }
  const double muT = 1.0 - 2.0 * u;//exact
  const double sT = 2.0 * std::sqrt( u * ( 1.0 - u ) );
  const double r = std::sqrt( r2 );
  const Vector om = randIsotropicDirection( rng );
  const double gx = ncfma( r, om.x(), x * sT );
  const double gy = r * om.y();
  const double gz = ncfma( r, om.z(), -ncfma( -x, muT, y ) );
  const double gxy2 = ncfma( gx, gx, gy * gy );
  const double Ak = m_A * m_k;
  Coll res;
  res.alpha = ( Ak * m_k ) * ncfma( gz, gz, gxy2 );
  res.beta = ncmax( -c, ncfma( ( 2.0 * y ) * m_k, gz, res.alpha ) );
  const double nz = ncfma( Ak, gz, y );//final n velocity: (Ak*gx,Ak*gy,nz)
  const double nmag2 = ncfma( nz, nz, ncsquare( Ak ) * gxy2 );
  res.mu = ( nmag2 > 0.0
             ? ncclamp( nz / std::sqrt( nmag2 ), -1.0, 1.0 )
             : 2.0 * rng.generate() - 1.0 );//exact
  return res;
}

NC::PairDD NC::FreeGasSampler::sampleAlphaBeta( RNG& rng,
                                                NeutronEnergy ekin ) const
{
  const EDivKT e{ ekin.dbl() * m_invkT };
  auto s = sampleColl( rng, e.dbl() );
  //Guard against rounding outside kinematic limits:
  auto alim = getAlphaLimits( e, s.beta );
  return { ncclamp( s.alpha, alim.first, alim.second ), s.beta };
}

NC::PairDD NC::FreeGasSampler::sampleDeltaEMu( RNG& rng,
                                               NeutronEnergy ekin ) const
{
  auto s = sampleColl( rng, ekin.dbl() * m_invkT );
  return { s.beta * m_kT, s.mu };
}
