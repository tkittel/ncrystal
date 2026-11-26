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

#include "NCrystal/internal/scbragg/NCSCBragg.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NC = NCrystal;

int main(int , char**)
{
  NC::SCOrientation orient;
  using Basis = NC::OrientDir::CrysDir::Basis;
  orient.setPrimaryDirection( { { Basis::Cartesian, { 1, 0, 0 } },
                                { 1, 0, 0 } } );
  orient.setSecondaryDirection( { { Basis::Cartesian, { 0, 1, 0 } },
                                  { 0, 1, 0 } }, 1e-6 );
  const NC::MosaicityFWHM mosaicity{ 1 * NCrystal::kDeg };

  //focusing on the  0   2  -2 family
  NCrystal::MatCfg cfg("Ge_sg227.ncmat;dcutoff=2;dcutoffup=2.1;inelas=0;incoh_elas=0");
  auto ci = NC::FactImpl::createInfo(cfg);
  auto sc = NC::makeSO<NC::SCBragg>(ci,orient,mosaicity);
  auto pc = NC::FactImpl::createScatter(cfg);

  auto rng  = NC::getRNG();

  //  NC::NeutronDirection incident{-1,0.1,0};

  nc_assert_always(!ci->hklList().empty());
  double thr = 2 * ci->hklList().front().dspacing;
  auto wlt_vec = NC::linspace(thr*0.9,thr*0.99,10);
  wlt_vec.insert(wlt_vec.begin(),thr*0.6);
  wlt_vec.push_back(thr*0.99999);
  wlt_vec.push_back(thr*1.00001);
  unsigned nLoop = 100000;
  //unit test of xs
  NC::CachePtr cp_sc, cp_pc;

  for(auto wl: wlt_vec)
  {
    double xss=0., xsp=0.;
    NC::NeutronEnergy ekin = NC::NeutronWavelength{wl};
    for(unsigned i=0;i<nLoop;i++)
    {
      xss += sc->crossSection( cp_sc, ekin,
                               NC::randIsotropicNeutronDirection(rng) ).get();
    }
    xss/=nLoop;
    xsp=pc->crossSectionIsotropic(cp_pc, ekin ).get();

    printf("wl %g, pc %g, sgl_cry %g, error %g%% \n", wl , xsp, xss, (xsp?(xss-xsp)/xsp*100.:(xss?9999.:0.)));
    if (wl>thr&&xss>0.0) {
      printf("Forbidden xs beyond threshold seen!\n");
      return 1;
    }
  }

  //Neutrons exactly along plane normals near back-scattering: cross sections
  //must be as for slightly tilted neutrons, and scatterings must be Bragg
  //reflections.
  for ( double wlfact : { 0.9999, 0.99995, 0.99999 } ) {
    const NC::NeutronEnergy ekin = NC::NeutronWavelength{ thr*wlfact };
    const double mu_bragg = 1.0 - 2.0*wlfact*wlfact;
    for ( auto& n : { NC::Vector(1,0,0), NC::Vector(-1,0,0),
                      NC::Vector(1,1,0), NC::Vector(-1,-1,0),
                      NC::Vector(0,1,-1), NC::Vector(1,0,1),
                      NC::Vector(-1,0,1), NC::Vector(0,-1,-1) } ) {
      const NC::Vector v = n.unit();
      const NC::Vector vt = NC::Vector( v.x()+1e-6, v.y()-2e-6,
                                        v.z()+1e-6 ).unit();
      const double xs = sc->crossSection( cp_sc, ekin,
                                          v.as<NC::NeutronDirection>() ).get();
      const double xst = sc->crossSection( cp_sc, ekin,
                                           vt.as<NC::NeutronDirection>() ).get();
      if ( n.mag2() == 1.0 ) {
        //Not a normal of the 220 family.
        nc_assert_always( xs == 0.0 && xst == 0.0 );
        continue;
      }
      nc_assert_always( xs > 0.0 );
      nc_assert_always( NC::floateq( xs, xst, 1e-3, 0.0 ) );
      for ( unsigned i = 0; i < 100; ++i ) {
        auto out = sc->sampleScatter( cp_sc, *rng, ekin,
                                      v.as<NC::NeutronDirection>() );
        nc_assert_always( out.ekin.get() == ekin.get() );
        nc_assert_always( NC::floateq( out.direction.as<NC::Vector>().dot(v),
                                       mu_bragg ) );
      }
    }
  }

  return 0;
}
