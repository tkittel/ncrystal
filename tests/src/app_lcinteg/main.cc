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

////////////////////////////////////////////////////////////////////////////////
// Checks the precision of LCStdFrame::calcXSIntegral, by comparing with a
// brute-force integration of LCStdFrame::calcXS.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/extd_utils/NCLCUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/interfaces/NCRNG.hh"

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  void test( NC::RNG& rng, double mosdeg, double prec )
  {
    NC::LCStdFrame frame( NC::MosaicityFWHM{ mosdeg * NC::kDeg }, prec );
    const auto& gm = frame.gaussMos();
    const double cta = gm.mosaicityCosTruncationAngle();
    const double sta = gm.mosaicitySinTruncationAngle();

    //Planesets with various d-spacings and angles to the lcaxis, sorted by
    //d-spacing:
    std::vector<NC::LCPlaneSet> planes;
    double dsp = 3.5;
    for ( auto i : NC::ncrange(40) ) {
      const double alpha = ( i%10==0 ? NC::kPiHalf
                             : ( 0.02 + 0.98*rng.generate() ) * NC::kPiHalf );
      planes.emplace_back( dsp, alpha, gm.mosaicityTruncationAngle(),
                           0.1 + rng.generate() );
      dsp *= ( 1.0 - 0.09*rng.generate() );
    }

    std::size_t nroi = 0;
    for ( auto i : NC::ncrange(40) ) {
      const double wl = ( i%2 ? 0.1 + 0.9*rng.generate()
                          : 0.5 + 6.0*rng.generate() );
      const double c3 = ( i%7 ? rng.generate() : 0.001*rng.generate() );
      NC::LCROIFinder finder( wl, c3, cta, sta );
      std::vector<NC::LCROI> rois;
      for ( auto& ps : planes ) {
        if ( wl > ps.twodsp )
          break;
        finder.findROIs( &ps, rois );
      }
      NC::LCStdFrame::NeutronPars neutron( wl, c3, std::sqrt(1.0-c3*c3) );
      double sum(0.0), sum_ref(0.0);
      std::vector<NC::PairDD> vals;
      for ( auto& roi : rois ) {
        if ( roi.isDegenerate() )
          continue;
        NC::LCStdFrame::NormalPars normal( roi.planeset, roi.normal_sign );
        const double a(roi.rotmin), b(roi.rotmax);
        //Brute-force reference, with the trapezoidal rule:
        constexpr unsigned n = 512;
        const double h = (b-a)/n;
        NC::StableSum ss;
        ss.add( 0.5*frame.calcXS( neutron, normal, std::cos(a) ) );
        ss.add( 0.5*frame.calcXS( neutron, normal, std::cos(b) ) );
        for ( unsigned j = 1; j < n; ++j )
          ss.add( frame.calcXS( neutron, normal, std::cos(a+j*h) ) );
        const double ref = ss.sum() * h;
        const double val = frame.calcXSIntegral( neutron, normal, a, b );
        REQUIRE( ref > 0.0 && val > 0.0 );
        vals.emplace_back( val, ref );
        sum += val;
        sum_ref += ref;
        ++nroi;
      }
      //Check the sum, i.e. the cross section. For the individual integrals,
      //the brute-force reference is not precise enough for the negligible
      //ones (those are dominated by the non-smooth truncation edges).
      REQUIRE( NC::ncabs( sum - sum_ref ) <= prec*sum_ref );
      for ( auto& e : vals )
        REQUIRE( NC::ncabs( e.first - e.second )
                 <= 2.0*prec*e.second + 1e-3*prec*sum_ref );
    }
    REQUIRE( nroi > 200 );
  }
}

int main()
{
  auto rng = NC::createBuiltinRNG( 12345 );
  for ( double mosdeg : { 0.01, 0.5, 5.0 } )
    for ( double prec : { 1e-3, 1e-5 } )
      test( rng, mosdeg, prec );
  return 0;
}
