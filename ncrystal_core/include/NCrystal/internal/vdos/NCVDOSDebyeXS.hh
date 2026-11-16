#ifndef NCrystal_VDOSDebyeXS_hh
#define NCrystal_VDOSDebyeXS_hh

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

namespace NCRYSTAL_NAMESPACE {

  namespace VDOS {

    class VDOSDebyeXS final : private MoveOnly {
    public:

      //Inelastic cross section of the idealised Debye VDOS (density ~E^2 up
      //to k*TDebye), per unit bound scattering cross section, without any
      //VDOS expansion: the one-phonon part is computed exactly, the
      //multi-phonon part comes from a fit (NCVDOSDebyeXSData.hh, made by
      //devel/simplebuild/static/NCDev/scripts/fit_vdosdebye_xs) with a max
      //relative error on the total below 1e-3 for T <= 50*TDebye, and below
      //1e-2 for the unphysical T in [50,1000]*TDebye (above melting for any
      //real solid). The high energy limit is exact. Throws CalcError for T >
      //maxTOverTDebye()*TDebye or for masses outside massRange():
      VDOSDebyeXS( DebyeTemperature, Temperature, AtomMass );

      //Supported range:
      static std::pair<AtomMass,AtomMass> massRange();//[0.95,260]*m_neutron
      static double maxTOverTDebye();//1000

      //sigma/sigma_bound, and its one- and multi-phonon parts:
      double sigma( NeutronEnergy ) const;
      PairDD sigmaParts( NeutronEnergy ) const;

      //sigma/sigma_bound as an IsoXSTable (see NCIsoXSUtils.hh), agreeing
      //with sigma(..) to ~1e-4 (exact high-energy limit):
      IsoXSTable createXSTabulation() const;

    private:
      double m_x, m_A, m_kT, m_ktheta, m_lam, m_sigfree, m_umax, m_sigUmax;
      //Per fit piece: weight in x and u-series coefficients (x and A parts
      //already summed):
      struct Series { double wx; VectD c; };
      std::array<Series,6> m_series;
      double sigma1( double e ) const;
      double sigmaMP( double u ) const;
    };

  }
}

#endif
