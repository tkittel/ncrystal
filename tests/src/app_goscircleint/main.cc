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
// Checks the precision of GaussOnSphere::circleIntegral where it is evaluated
// numerically, by comparison with a brute-force integration of the exact
// density along the circle.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/phys_utils/NCGaussOnSphere.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/interfaces/NCRNG.hh"

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  unsigned s_nconfigs = 0;

  //Brute-force: Simpson's rule over the angle t along the circle, with the
  //density from the angle beta to the centre of the Gaussian. Small values of
  //w=1-cos(beta) are found as for gamma-alpha.
  double refIntegral( const NC::GaussOnSphere& gos, double gamma, double alpha )
  {
    const double ta = gos.getTruncangle();
    const double sa = std::sin(alpha), sg = std::sin(gamma);
    const double sasg = sa * sg;
    const double sd = std::sin(0.5*(gamma-alpha));
    const double w0 = 2.0*sd*sd;
    const double st = std::sin(0.5*ta);
    const double y = ( 2.0*st*st - w0 ) / ( 2.0*sasg );
    if ( !( y > 0.0 ) )
      return 0.0;
    const double tmax = ( y >= 1.0 ? NC::kPi : 2.0*std::asin(std::sqrt(y)) );
    constexpr unsigned n = 20000;
    const double h = tmax / n;
    NC::StableSum sum;
    for ( unsigned i = 0; i <= n; ++i ) {
      const double sh = std::sin( 0.5 * i * h );
      const double w = w0 + 2.0*sasg*sh*sh;
      const double beta = 2.0*std::asin(std::sqrt(0.5*w));
      const double wt = ( i==0||i==n ? 1.0 : ( i%2 ? 4.0 : 2.0 ) );
      sum.add( wt * gos.evalXExact( NC::ncmin( beta, ta ) ) );
    }
    return 2.0 * sa * sum.sum() * h / 3.0;
  }

  void test( NC::RNG& rng, double mosdeg, double ntrunc, double prec )
  {
    const double sigma = mosdeg*NC::kDeg/2.354820045030949;
    //Only test where the closed form is never used (it is valid for
    //truncation angles x where x^4/24 < prec/2, roughly):
    if ( !( ntrunc*sigma > 1.2 * std::pow( 12.0*prec, 0.25 )
            && ntrunc*sigma < NC::kPiHalf ) )
      return;
    ++s_nconfigs;
    NC::GaussOnSphere gos( sigma, ntrunc*sigma, prec );
    const double ta = gos.getTruncangle();
    //Required precision, with some room for the reference:
    const double acc = NC::ncclamp( prec*0.1, 1e-7, 1e-4 ) + 1e-9;
    unsigned ntested = 0;
    for ( auto i : NC::ncrange(400) ) {
      //Circles of all sizes, also smaller than the truncation radius:
      const double alpha = ( i%5==0 ? 2.0*ta*rng.generate()
                             : NC::kPiHalf*rng.generate() );
      const double gamma = alpha + (2.0*rng.generate()-1.0)*NC::ncmin(ta,5*sigma);
      if ( !( alpha > 0.0 && gamma > 0.0 && gamma < NC::kPi ) )
        continue;
      const double ref = refIntegral( gos, gamma, alpha );
      const double val = gos.circleIntegral( std::cos(gamma), std::sin(gamma),
                                             std::cos(alpha), std::sin(alpha) );
      //Ignore negligible values:
      if ( !( ref > 1e-3 * gos.getNormFactor() * sigma ) )
        continue;
      ++ntested;
      REQUIRE( NC::ncabs( val - ref ) <= acc * ref );
    }
    REQUIRE( ntested > 100 );
  }
}

int main()
{
  auto rng = NC::createBuiltinRNG( 12345 );
  for ( double mosdeg : { 2.0, 5.0, 20.0, 40.0 } )
    for ( double prec : { 1e-2, 1e-4, 1e-5, 1e-7 } )
      test( *rng, mosdeg, NC::GaussOnSphere::estimateNTruncFromPrec(prec), prec );
  //Large truncations (several segments), at small mosaicities:
  for ( double ntrunc : { 9.0, 15.0, 25.0 } )
    test( *rng, 1.0, ntrunc, 1e-6 );
  REQUIRE( s_nconfigs >= 10 );
  return 0;
}
