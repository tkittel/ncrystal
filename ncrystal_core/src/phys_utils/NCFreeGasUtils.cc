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

#include "NCrystal/internal/phys_utils/NCFreeGasUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCSpline.hh"
namespace NC=NCrystal;

NC::FreeGasXSProvider::FreeGasXSProvider( Temperature temp_kelvin,
                                          AtomMass target_mass_amu,
                                          SigmaFree sigma )
{
  //Evaluate cross-sections via sampling paper (10.1016/j.jcp.2018.11.043)
  //eq. 20, with explicit prefactor of sigma_free=sigma_bound*(A/(1+A))^2" (see
  //also discussion after eq. 20 in
  //https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927245/). And discussion on
  //DGSW-399.
  //
  //Sampling paper eq. 20 does not include A=M/mn in "a", but here we scale a by
  //sqrt(A). Not sure where we got this from, but it was heavily validated using
  //numerical integration of sampling paper eq. 19.
  temp_kelvin.validate();
  target_mass_amu.validate();
  sigma.validate();

  double A = target_mass_amu.relativeToNeutronMass();//mass_nuclide / mass_neutron
  m_sigmaFree = sigma.get();
  nc_assert_always( temp_kelvin.get() > 0.0);
  nc_assert_always( target_mass_amu.get() > 0.0);
  nc_assert_always( m_sigmaFree > 0.0);

  nc_assert( temp_kelvin.kT() > 0.0 );
  m_ca = A / temp_kelvin.kT();//a^2 = A*E/kT = ca*E with ca = A/kT [units 1/eV]:
}

NC::FreeGasXSProvider::FreeGasXSProvider( Temperature temp_kelvin,
                                          AtomMass target_mass_amu,
                                          SigmaBound sb )
  : FreeGasXSProvider( temp_kelvin, target_mass_amu, sb.free(target_mass_amu) )
{
}

NC::FreeGasXSProvider::~FreeGasXSProvider() = default;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //With x=a^2=ca*E, xs=sigma_free*G(x)/x, where
    //G(x)=(x+1/2)*erf(sqrt(x))+sqrt(x)*exp(-x)/sqrt(pi)=sqrt(x)*H(x), with H
    //smooth in x. So xs=sigma_free*U(ca*E), where U is the cross section for
    //sigma_free=1 and ca=1, i.e. with E=x:
    //
    //  x*U(x) = sqrt(x) * H(x)
    //
    //The universal table of U tabulates H in x (the Q part), for x in
    //[1e-16,36]: below it U=2/sqrt(pi*x) is good to x/3<4e-17 (so it joins the
    //table to rounding), and above it the highEForm U=1+1/(2x) applies.
    double freeGasH( double x )
    {
      const double s = std::sqrt( x );
      return stable_exp( -x ) * kInvSqrtPi + ( s + 0.5 / s ) * ncerf( s );
    }

    IsoXSTable tabulateFreeGasUniversal()
    {
      constexpr double xmax = 36.0;//must match highEForm
      auto fit = fitPiecewiseCubic( freeGasH, 1e-16, xmax, 1e-10, 4.0 );
      nc_assert_always( fit.x.back() == xmax );
      IsoXSTable t;
      t.E = std::move( fit.x );
      t.q = std::move( fit.c );
      t.tail[0] = 2.0 * kInvSqrtPi;
      t.tail[2] = 0.5;
      t.tail[3] = 1.0;
      return t;
    }

    const IsoXSEvaluator& freeGasUniversal()
    {
      static const IsoXSEvaluator ev( tabulateFreeGasUniversal() );
      return ev;
    }
  }
}

NC::CrossSect NC::FreeGasXSProvider::crossSectionTabulated( NeutronEnergy ekin ) const
{
  return CrossSect{ freeGasUniversal().evalScaled( ekin.dbl(), m_ca, m_sigmaFree ) };
}

NC::IsoXSTable NC::FreeGasXSProvider::createXSTabulation() const
{
  nc_assert( m_ca > 0.0 );
  return scaleIsoXSTable( freeGasUniversal().table(), m_ca, m_sigmaFree );
}

NC::IsoXSHighEForm NC::FreeGasXSProvider::highEForm() const
{
  //evalXSShapeASq is 1+0.5/a^2 above a^2=36 (neglecting terms of order
  //exp(-36)/a), with a^2=m_ca*E:
  nc_assert( m_ca > 0.0 );
  return { m_sigmaFree, 0.5 * m_sigmaFree / m_ca, NeutronEnergy{ 36.0 / m_ca } };
}

double NC::FreeGasXSProvider::evalXSShapeASq(double a_squared)
{
  //a^2 = A*E/kT.
  if (a_squared>36.0)
    return 1.0 + 0.5 / a_squared;
  const double a = std::sqrt(a_squared);
  if (a<0.1) {
    if (a==0.0)
      return kInfinity;
    constexpr double c1 = 2.0/3.0;
    constexpr double c2 = 1.0/15.0;
    constexpr double c3 = 1.0/105.0;
    constexpr double c4 = 1.0/756.0;
    constexpr double c5 = 1.0/5940.0;
    const double a2 = a_squared;
    return kInvSqrtPi * ( 2.0 / a +  a *( c1- a2*(c2-a2*(c3-a2*(c4-a2*c5)))));
  };
  //intermediate region, fall-back to full formula (slow):
  const double inva = 1.0 / a;
  return ( 1.0 + 0.5*inva*inva ) * ncerf(a) + kInvSqrtPi * std::exp(-a_squared)*inva;
}
