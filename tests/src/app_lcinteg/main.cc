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
// Checks the precision of LCStdFrame::calcXSIntegral and the two methods it
// uses (series expansion and numerical integration), against independent
// high precision reference values from devel/misc/gen_mosaicrefs.py.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/extd_utils/NCLCUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCTestUtils/NCTestFindData.hh"
#include <fstream>
#include <iostream>
#include <map>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  struct RefCase {
    double mosdeg, ntrunc, prec, dsp, alpha, fsq, wl, c3, sign, ref;
  };

  std::vector<RefCase> loadRefs()
  {
    std::vector<RefCase> res;
    for ( auto& v : nctest::load_test_data_numbers( "mosaicrefs",
                                                    "lcinteg.txt", 10 ) )
      res.push_back( { v[0], v[1], v[2], v[3], v[4],
                       v[5], v[6], v[7], v[8], v[9] } );
    REQUIRE( res.size() >= 100 );
    return res;
  }

  double relerr( double val, double ref )
  {
    return NC::ncabs( val - ref ) / ref;
  }
}

int main()
{
  std::map<std::tuple<double,double,double>,
           std::unique_ptr<NC::LCStdFrame>> frames;
  unsigned nseries(0), nnoseries(0);
  double worst_series(0.0), worst_num(0.0);//relative errors over prec
  double bias_num(0.0);//sum of signed relative errors over prec
  for ( auto& rc : loadRefs() ) {
    auto& frameptr = frames[ std::make_tuple( rc.mosdeg, rc.ntrunc, rc.prec ) ];
    if ( !frameptr )
      frameptr = std::make_unique<NC::LCStdFrame>
        ( NC::MosaicityFWHM{ rc.mosdeg * NC::kDeg }, rc.prec, rc.ntrunc );
    const auto& frame = *frameptr;
    const auto& gm = frame.gaussMos();
    const NC::LCPlaneSet ps( rc.dsp, rc.alpha, gm.mosaicityTruncationAngle(),
                             rc.fsq );
    const NC::LCStdFrame::NeutronPars neutron( rc.wl, rc.c3,
                                               std::sqrt( 1.0-rc.c3*rc.c3 ) );
    const NC::LCStdFrame::NormalPars normal( &ps, rc.sign );
    //Regions of interest in phi, from the ROI finder:
    NC::LCROIFinder finder( rc.wl, rc.c3, gm.mosaicityCosTruncationAngle(),
                            gm.mosaicitySinTruncationAngle() );
    std::vector<NC::LCROI> rois;
    finder.findROIs( &ps, rois );
    double val(0.0), num(0.0);
    for ( auto& roi : rois ) {
      if ( roi.isDegenerate() || roi.normal_sign != rc.sign )
        continue;
      val += frame.calcXSIntegral( neutron, normal, roi.rotmin, roi.rotmax );
      num += frame.calcXSIntegralNumerical( neutron, normal,
                                            roi.rotmin, roi.rotmax );
    }
    const auto series = frame.calcXSIntegralSeries( neutron, normal );
    //Tolerances: The series is only used when its error bound is below
    //0.1*prec, and the numerical integration (also relying on approximations
    //in GaussOnSphere) aims at a precision of prec, but in practice does much
    //better (and without a systematic bias, checked below):
    REQUIRE( relerr( num, rc.ref ) <= 0.5 * rc.prec );
    REQUIRE( relerr( val, rc.ref ) <= 0.5 * rc.prec );
    worst_num = NC::ncmax( worst_num, relerr( num, rc.ref ) / rc.prec );
    bias_num += ( num / rc.ref - 1.0 ) / rc.prec;
    if ( series.has_value() ) {
      ++nseries;
      REQUIRE( relerr( series.value(), rc.ref ) <= 0.1 * rc.prec );
      REQUIRE( val == series.value() );//used when available
      worst_series = NC::ncmax( worst_series,
                                relerr( series.value(), rc.ref ) / rc.prec );
    } else {
      ++nnoseries;
      REQUIRE( val == num );
    }
  }
  std::cout << "Series used for " << nseries << " cases, and numerical"
            << " integration for " << nnoseries << " cases." << std::endl;
  std::cout << "Worst relative errors over prec: series "
            << NC::fmt( worst_series ) << ", numerical "
            << NC::fmt( worst_num ) << std::endl;
  bias_num /= ( nseries + nnoseries );
  std::cout << "Mean signed relative error over prec, numerical: "
            << NC::fmt( bias_num ) << std::endl;
  REQUIRE( NC::ncabs( bias_num ) <= 0.05 );
  //Both methods must be tested:
  REQUIRE( nseries >= 20 && nnoseries >= 20 );
  return 0;
}
