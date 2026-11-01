#ifndef NCrystal_FreeGasSampler_FMA_hh
#define NCrystal_FreeGasSampler_FMA_hh

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

#include "NCrystal/internal/utils/NCFMADispatch.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NCRYSTAL_NAMESPACE {
#ifndef NCRYSTAL_WIN_FMA
  namespace {
#endif

    NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
      void, freeGasSampleColl,
      ( RNG& rng, double A, double k, double c, bool wantmu, double* res )
    )

    //Sample a free-gas collision (cf. FreeGasSampler), for target mass A (in
    //neutron masses), k=1/(1+A) and c=E/kT. Fills res[0]=alpha, res[1]=beta
    //and (only if wantmu) res[2]=mu=cos(theta_scat), so res must have room
    //for 2 or 3 values.
    NCRYSTAL_FMADISPATCH_DECLARATOR(void,freeGasSampleColl)
    ( RNG& rng, double A, double k, double c, bool wantmu, double* res )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
        freeGasSampleColl, (rng,A,k,c,wantmu,res) );
      //Work with velocities in units of sqrt(2kT/M): neutron n=(0,0,y),
      //target T=x*(sT,0,muT) (azimuth irrelevant). Sample (x,muT) with
      //density ~|n-T|*x^2*exp(-x^2), then rotate the relative velocity r=n-T
      //isotropically (CM frame). The neutron velocity change is then A*k*g
      //with g=|r|*Om-r, giving alpha=A*k^2*g^2, beta=2*y*k*g_z+alpha (no
      //large cancellations).
      //
      //Contraction safety: every a*b+c is an explicit ncfma (or exact, like
      //1-2*u), so results do not depend on FP-contraction settings (and this
      //function is safe for NCRYSTAL_FMADISPATCH_ATTR).
      nc_assert( c >= 0.0 );
      constexpr double halfsqrtpi = 0.5 * kSqrtPi;
      const double y = std::sqrt( A * c );
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
      const double Ak = A * k;
      const double alpha = ( Ak * k ) * ncfma( gz, gz, gxy2 );
      res[0] = alpha;
      res[1] = ncmax( -c, ncfma( ( 2.0 * y ) * k, gz, alpha ) );
      if ( !wantmu )
        return;
      const double nz = ncfma( Ak, gz, y );//n final: (Ak*gx,Ak*gy,nz)
      const double nmag2 = ncfma( nz, nz, ncsquare( Ak ) * gxy2 );
      res[2] = ( nmag2 > 0.0
                 ? ncclamp( nz / std::sqrt( nmag2 ), -1.0, 1.0 )
                 : 2.0 * rng.generate() - 1.0 );//exact
    }

#ifndef NCRYSTAL_WIN_FMA
  }
#endif
}

#endif
