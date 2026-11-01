#ifndef NCrystal_FreeGasSampler_hh
#define NCrystal_FreeGasSampler_hh

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

#include "NCrystal/core/NCTypes.hh"

////////////////////////////////////////////////////////////////////////////////
// Free-gas scattering sampled as an explicit elastic two-body collision with //
// a target whose velocity is drawn from a Maxwellian, weighted by the        //
// relative speed (constant sigma_free means collision rate ~ |v_n-V|).       //
// For constant sigma_free this is exactly equivalent to the closed-form      //
// free-gas S(alpha,beta) (whose cross sections FreeGasXSProvider provides),  //
// and needs no E-dependent setup. Target speed x (in units of                //
// sqrt(2kT/M)) and cosine muT wrt. the neutron are sampled with the classic  //
// rejection scheme: x^2 from Gamma(2) or Gamma(3/2) with relative weights    //
// 1 : y*sqrt(pi)/2 (y=sqrt(A*E/kT) the reduced neutron speed), muT flat, and //
// acceptance |v_rel|/(v_n+V). The outgoing direction is then isotropic in    //
// the CM frame. Refs: R.R. Coveyou, R.R. Bate, R.K. Osborn, J. Nucl. Energy  //
// 2 (1956) 153-167, doi:10.1016/0891-3919(55)90030-9, eqs. 5-7 (method there //
// credited to H. Kahn, Rand Corp., private communication). Modern write-up:  //
// P.K. Romano, J.A. Walsh, Ann. Nucl. Energy 114 (2018) 318-324,             //
// doi:10.1016/j.anucene.2017.12.044, sec. 2.1 and alg. 1.                    //
//                                                                            //
// Results use the usual S(alpha,beta) conventions (alpha in units of         //
// hbar^2Q^2/(2*m_n*kT), beta=(E'-E)/kT).                                     //
////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  class FreeGasSampler final {
  public:
    FreeGasSampler( Temperature, AtomMass );

    PairDD sampleAlphaBeta( RNG&, NeutronEnergy ) const;
    PairDD sampleDeltaEMu( RNG&, NeutronEnergy ) const;//(E'-E,mu)

  private:
    struct Coll { double alpha, beta, mu; };
    Coll sampleColl( RNG&, double c ) const;//c=E/kT
    double m_kT, m_invkT, m_A, m_k;//k=1/(1+A)
  };

}

#endif
