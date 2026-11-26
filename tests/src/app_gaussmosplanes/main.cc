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
// Checks that GaussMos::calcCrossSections with a PlaneList (which preselects
// planes in single precision) gives exactly the same as a straightforward
// loop over all families and normals.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/phys_utils/NCGaussMos.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/interfaces/NCRNG.hh"

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  struct Family {
    double inv2dsp, xsfact;
    std::vector<NC::Vector> normals;
  };
  using ScatCaches = std::vector<NC::GaussMos::ScatCache>;

  //Reference: all planes, without any preselection.
  void calcRef( const NC::GaussMos& gm, const std::vector<Family>& fams,
                double wl, const NC::Vector& indir,
                ScatCaches& cache, NC::VectD& xs_commul )
  {
    const double cta = gm.mosaicityCosTruncationAngle();
    const double inv2dcutoff
      = (1.0-2*std::numeric_limits<double>::epsilon())/wl;
    NC::GaussMos::InteractionPars ip;
    for ( auto& fam : fams ) {
      if ( fam.inv2dsp >= inv2dcutoff
           || NC::GaussMos::cacheRound(fam.inv2dsp) * NC::GaussMos::cacheRound(wl) >= 1.0 )
        break;
      ip.set( wl, fam.inv2dsp, fam.xsfact );
      const double inv2dsp = NC::GaussMos::cacheRound(fam.inv2dsp);
      const double sinth = NC::GaussMos::cacheRound(wl) * inv2dsp;
      const double xsoffset = xs_commul.empty() ? 0.0 : xs_commul.back();
      double xssum(0.0);
      for ( auto& normal : fam.normals ) {
        const double dot = normal.dot(indir);
        const double sdotcptsq = (1.0 - dot * dot)*(1 - sinth*sinth);
        const double ds = dot * sinth;
        const double Am = cta - ds;
        if ( Am < 0.0 || sdotcptsq > Am*Am ) {
          const double xs = gm.calcRawCrossSectionValue(ip, dot );
          if (xs) {
            xs_commul.push_back(xsoffset + (xssum += xs));
            cache.emplace_back(-normal, inv2dsp );
          }
        }
        const double Ap = cta + ds;
        if ( Ap < 0.0 || sdotcptsq > Ap*Ap ) {
          const double xs = gm.calcRawCrossSectionValue(ip, -dot );
          if (xs) {
            xs_commul.push_back(xsoffset + (xssum += xs));
            cache.emplace_back(normal, inv2dsp );
          }
        }
      }
    }
  }

  //Compare for a given neutron, returning number of contributions.
  std::size_t testNeutron( const NC::GaussMos& gm,
                           const std::vector<Family>& fams,
                           const NC::GaussMos::PlaneList& pl,
                           double wl, const NC::Vector& indir )
  {
    //Start with existing content, which must be left alone.
    ScatCaches c_ref, c_new;
    NC::VectD xs_ref = { 0.25 };
    NC::VectD xs_new = { 0.25 };
    c_ref.emplace_back( NC::Vector{ 0., 0., 1. }, 1.0 );
    c_new.emplace_back( NC::Vector{ 0., 0., 1. }, 1.0 );
    calcRef( gm, fams, wl, indir, c_ref, xs_ref );
    gm.calcCrossSections( pl, wl, indir, c_new, xs_new );
    REQUIRE( xs_ref.size() == xs_new.size() );
    REQUIRE( c_ref.size() == c_new.size() );
    REQUIRE( c_ref.size() == xs_ref.size() );
    for ( auto i : NC::ncrange( xs_ref.size() ) ) {
      //Exact equality (also when infinite):
      REQUIRE( xs_ref.at(i) == xs_new.at(i) );
      REQUIRE( c_ref.at(i).plane_inv2d() == c_new.at(i).plane_inv2d() );
      auto& n1 = c_ref.at(i).plane_normal();
      auto& n2 = c_new.at(i).plane_normal();
      REQUIRE( n1.x()==n2.x() && n1.y()==n2.y() && n1.z()==n2.z() );
    }
    return xs_ref.size() - 1;
  }

  //Returns the number of contributions:
  std::size_t testPlanes( NC::RNG& rng, std::size_t nplanes_target,
                          NC::MosaicityFWHM mos, double prec )
  {
    //Random families in order of increasing inv2dsp. Some share inv2dsp, and
    //some normals are along or perpendicular to the axes.
    std::vector<Family> fams;
    NC::GaussMos::PlaneList pl;
    double inv2dsp = 0.1;
    std::size_t nplanes = 0;
    while ( nplanes < nplanes_target ) {
      if ( fams.empty() || rng.generate() > 0.2 )
        inv2dsp *= ( 1.0 + rng.generate() * 30.0 / (nplanes_target+10.0) );
      Family fam{ inv2dsp, 0.01 + rng.generate(), {} };
      const std::size_t n = std::min<std::size_t>
        ( nplanes_target - nplanes, 1 + std::size_t(rng.generate()*24.0) );
      for ( auto i : NC::ncrange(n) ) {
        (void)i;
        const double r = rng.generate();
        if ( r < 0.02 )
          fam.normals.emplace_back( 1., 0., 0. );
        else if ( r < 0.04 )
          fam.normals.emplace_back( 0., -1., 0. );
        else
          fam.normals.push_back( NC::randIsotropicDirection(rng) );
      }
      nplanes += n;
      pl.addFamily( fam.inv2dsp, fam.xsfact, fam.normals );
      fams.push_back( std::move(fam) );
    }
    REQUIRE( pl.nPlanes() == nplanes_target );
    REQUIRE( pl.nFamilies() == fams.size() );
    REQUIRE( fams.empty() || pl.minInv2dsp() == fams.front().inv2dsp );
    REQUIRE( fams.empty() || pl.maxInv2dsp() == fams.back().inv2dsp );

    NC::GaussMos gm( mos, prec );
    std::size_t ncontrib = 0;
    const double wlmax = fams.empty() ? 1.0 : 1.0/fams.front().inv2dsp;
    const double wlmin = fams.empty() ? 1.0 : 1.0/fams.back().inv2dsp;
    for ( auto i : NC::ncrange(300) ) {
      //Random wavelengths (some beyond or exactly at Bragg thresholds):
      double wl = rng.generate() * wlmax * 1.01;
      if ( i%10 == 0 )
        wl *= 1e-3;
      else if ( i%10 == 1 && !fams.empty() )
        wl = 1.0 / fams.at( fams.size()/3 ).inv2dsp;
      else if ( i%10 == 2 )
        wl = wlmin * ( 0.95 + 0.1*rng.generate() );
      else if ( i%10 == 3 && !fams.empty() )//just below a Bragg threshold
        wl = ( 1.0 / fams.at( i % fams.size() ).inv2dsp )
          * ( 1.0 - ( i % 4 ) * std::numeric_limits<double>::epsilon() );
      //Random directions, or along axes (i.e. along some normals):
      NC::Vector dir = NC::randIsotropicDirection(rng);
      if ( i%7 == 0 )
        dir.set( 1., 0., 0. );
      else if ( i%7 == 1 )
        dir.set( 0., 1., 0. );
      else if ( i%7 == 2 )
        dir.set( 0., 0., -1. );
      else if ( i%7 == 3 && !fams.empty() ) {
        //Exactly fulfilling the Bragg condition for a plane:
        const auto& fam = fams.at( i % fams.size() );
        const auto& n = fam.normals.front();
        const double sinth = wl * fam.inv2dsp;
        if ( sinth < 1.0 ) {
          const NC::Vector perp = n.cross( NC::randIsotropicDirection(rng) );
          dir = n * (-sinth) + perp.unit() * std::sqrt( 1.0 - sinth*sinth );
          dir.normalise();
        }
      }
      ncontrib += testNeutron( gm, fams, pl, wl, dir );
    }
    return ncontrib;
  }
  void testThreshold()
  {
    //Independently of the internal rounding, planes never contribute when
    //wl>=2d, and cross sections are always finite (neutrons along the normal,
    //i.e. back-scattering, which contributes for wl just below 2d):
    const NC::GaussMos gm( NC::MosaicityFWHM{ 0.5 * NC::kDeg } );
    const NC::Vector n( 0.3, -0.4, std::sqrt( 0.75 ) );
    for ( double inv2dsp : { 0.1, 0.3333333333333333, 0.7 } ) {
      NC::GaussMos::PlaneList pl;
      pl.addFamily( inv2dsp, 0.5, std::vector<NC::Vector>{ n } );
      for ( int k = -8; k <= 8; ++k ) {
        const double wl = ( 1.0 / inv2dsp )
          * ( 1.0 + k * std::numeric_limits<double>::epsilon() );
        ScatCaches cache;
        NC::VectD xs;
        gm.calcCrossSections( pl, wl, -n, cache, xs );
        REQUIRE( xs.size() == cache.size() );
        for ( auto x : xs )
          REQUIRE( std::isfinite( x ) && x > 0.0 );
        if ( wl * inv2dsp >= 1.0 )
          REQUIRE( xs.empty() );
      }
      ScatCaches cache;
      NC::VectD xs;
      gm.calcCrossSections( pl, 0.99999 / inv2dsp, -n, cache, xs );
      REQUIRE( !xs.empty() );
    }
  }
}

int main()
{
  testThreshold();
  auto rng = NC::createBuiltinRNG( 12345 );
  //Number of planes below, at, and above multiples of the preselection chunk
  //sizes (16 and 64), and a wide range of mosaicities (incl. extreme values):
  //Also with non-default precision parameters (affecting the truncation):
  for ( std::size_t nplanes : { 0, 1, 2, 15, 16, 17, 31, 32, 33, 47, 48, 49,
                                63, 64, 65, 79, 80, 81, 111, 112, 113, 127,
                                128, 129, 191, 192, 193, 255, 256, 257, 1000,
                                4096, 5001 } ) {
    for ( double prec : { 1e-3, 1e-2, 1e-5, 1e-7 } ) {
      if ( prec != 1e-3 && !NC::isOneOf( nplanes, 1u, 17u, 65u, 129u ) )
        continue;//save time
      for ( double mosdeg : { 0.0003, 0.02, 0.5, 5.0, 30.0 } ) {
        const auto n = testPlanes( rng, nplanes,
                                   NC::MosaicityFWHM{ mosdeg * NC::kDeg },
                                   prec );
        //Non-trivial, thanks to neutrons fulfilling the Bragg condition:
        REQUIRE( n > 0 || nplanes == 0 );
      }
    }
  }
  return 0;
}
