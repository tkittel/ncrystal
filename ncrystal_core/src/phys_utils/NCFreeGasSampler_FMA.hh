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
      bool, freeGasCollStep,
      ( double y, bool gamma2, double a, double b, double u, double ua,
        double* xr2 )
    )
    NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
      bool, freeGasCollFinish,
      ( double A, double k, double c, double y, double x, double u,
        double r2, double omx, double omy, double omz, bool wantmu,
        double* res )
    )

    //Free-gas collision sampling (cf. FreeGasSampler) for target mass A (in
    //neutron masses), k=1/(1+A) and c=E/kT, split in two fma-dispatched
    //steps without any RNG access (it might throw), with velocities in units
    //of sqrt(2kT/M): neutron n=(0,0,y) with y=sqrt(A*c), target
    //T=x*(sT,0,muT) (azimuth irrelevant). (x,muT) is sampled with density
    //~|n-T|*x^2*exp(-x^2), then the relative velocity r=n-T is rotated
    //isotropically (CM frame). The neutron velocity change is then A*k*g with
    //g=|r|*Om-r, giving alpha=A*k^2*g^2, beta=2*y*k*g_z+alpha (no large
    //cancellations).
    //
    //Contraction safety: every a*b+c is an explicit ncfma (or exact, like
    //1-2*u), so results do not depend on FP-contraction settings (and these
    //functions are safe for NCRYSTAL_FMADISPATCH_ATTR).

    //One rejection-loop trial: x^2 from Gamma(2) (gamma2, a and b uniform)
    //or Gamma(3/2) (a normal, b exponential), u (muT=1-2u) and ua uniform.
    //Sets xr2[0]=x, xr2[1]=|n-T|^2 and returns whether the trial is accepted.
    NCRYSTAL_FMADISPATCH_DECLARATOR(bool,freeGasCollStep)
    ( double y, bool gamma2, double a, double b, double u, double ua,
      double* xr2 )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
        freeGasCollStep, (y,gamma2,a,b,u,ua,xr2) );
      const double x2 = ( gamma2
                          ? -std::log( a * b )
                          : ncfma( 0.5 * a, a, b ) );
      const double x = std::sqrt( x2 );
      const double ymx = y - x;
      const double r2 = ncfma( ymx, ymx, ( ( 4.0 * y ) * x ) * u );
      xr2[0] = x;
      xr2[1] = r2;
      return ncsquare( ua * ( y + x ) ) <= r2;
    }

    //Kinematics of an accepted trial, with (omx,omy,omz) an isotropic unit
    //vector. Fills res[0]=alpha, res[1]=beta and (if wantmu)
    //res[2]=mu=cos(theta_scat). Returns false if mu is undefined (zero final
    //neutron velocity), in which case res[2] must be sampled by the caller.
    NCRYSTAL_FMADISPATCH_DECLARATOR(bool,freeGasCollFinish)
    ( double A, double k, double c, double y, double x, double u,
      double r2, double omx, double omy, double omz, bool wantmu,
      double* res )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
        freeGasCollFinish, (A,k,c,y,x,u,r2,omx,omy,omz,wantmu,res) );
      const double muT = 1.0 - 2.0 * u;//exact
      const double sT = 2.0 * std::sqrt( u * ( 1.0 - u ) );
      const double r = std::sqrt( r2 );
      const double gx = ncfma( r, omx, x * sT );
      const double gy = r * omy;
      const double gz = ncfma( r, omz, -ncfma( -x, muT, y ) );
      const double gxy2 = ncfma( gx, gx, gy * gy );
      const double Ak = A * k;
      const double alpha = ( Ak * k ) * ncfma( gz, gz, gxy2 );
      res[0] = alpha;
      res[1] = ncmax( -c, ncfma( ( 2.0 * y ) * k, gz, alpha ) );
      if ( !wantmu )
        return true;
      const double nz = ncfma( Ak, gz, y );//n final: (Ak*gx,Ak*gy,nz)
      const double nmag2 = ncfma( nz, nz, ncsquare( Ak ) * gxy2 );
      if ( !( nmag2 > 0.0 ) )
        return false;
      res[2] = ncclamp( nz / std::sqrt( nmag2 ), -1.0, 1.0 );
      return true;
    }

#ifndef NCRYSTAL_WIN_FMA
  }
#endif
}

#endif
