#ifndef NCrystal_FreeGasUtils_hh
#define NCrystal_FreeGasUtils_hh

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
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"

////////////////////////////////////////////////////////////////////////////////
// Cross sections in the free-gas model (targets with no internal states, in  //
// thermal equilibrium, constant sigma_free), from integrating the closed-    //
// form free-gas S(alpha,beta) (eq. 19 of the NCrystal sampling paper,        //
// 10.1016/j.jcp.2018.11.043). Scatterings are sampled with                   //
// FreeGasSampler (NCFreeGasSampler.hh).                                      //
////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  class FreeGasXSProvider final {
  public:
    FreeGasXSProvider( Temperature, AtomMass, SigmaFree );
    FreeGasXSProvider( Temperature, AtomMass, SigmaBound );
    ~FreeGasXSProvider();

    //The cross section from the formula:
    CrossSect crossSectionUntabulated( NeutronEnergy ekin ) const;

    //Evaluate (1+1/(2a^2))*erf(a)+exp(-a^2)/(sqrt(pi)*a) (used internally, but
    //exposed here for testing):
    static double evalXSShapeASq(double a_squared);

    SigmaFree sigmaFree() const { return SigmaFree{ m_sigmaFree }; }

    //High energy form, sigma_free*(1+kT/(2*A*E)) above A*E/kT=36:
    IsoXSHighEForm highEForm() const;

    //The cross section from a universal table (built once and shared), with a
    //relative deviation below 1e-10 from crossSectionUntabulated:
    CrossSect crossSectionTabulated( NeutronEnergy ) const;

    //The same as an IsoXSTable (see NCIsoXSUtils.hh), agreeing
    //with crossSectionTabulated to within rounding:
    IsoXSTable createXSTabulation() const;
  private:
    double m_sigmaFree, m_ca;
  };

}

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  inline CrossSect FreeGasXSProvider::crossSectionUntabulated( NeutronEnergy ekin ) const
  {
    return CrossSect{ m_sigmaFree * evalXSShapeASq( m_ca * ekin.dbl() ) };
  }

}

#endif
