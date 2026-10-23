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


#include "NCrystal/internal/sab/NCSABBetaStripSampler.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>
#include <map>
#include <array>

//Tests the strips of SABBetaStripSampler against references computed
//independently of the class: S is interpolated by its own implementation,
//alpha limits of KB(e) come from cancellation-free closed forms, and the
//extrema over strips are found numerically.

namespace NC = NCrystal;
namespace NS = NC::SABUtils;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  NC::VectD geom( double a, double b, unsigned n )
  {
    NC::VectD v;
    for ( auto i : NC::ncrange( n ) )
      v.push_back( a * std::pow( b / a, double(i) / ( n - 1 ) ) );
    return v;
  }

  struct Kernel {
    std::string name;
    NC::VectD ag, bg, sab;//sab[ib*na+ia]
    NC::shared_obj<const NC::SABData> sd;
  };

  Kernel makeKernel( std::string name, NC::VectD bg, NC::VectD ag,
                     NC::VectD sab )
  {
    Kernel k{ std::move(name), ag, bg, sab,
              NC::makeSO<const NC::SABData>( std::move( ag ), std::move( bg ),
                                             std::move( sab ),
                                             NC::Temperature{ 300.0 },
                                             NC::SigmaBound{ 1.0 },
                                             NC::AtomMass{ 1.0 } ) };
    return k;
  }

  //Smooth S, zero in the beta rows listed in zerorows:
  Kernel makeSmoothKernel( std::string name, NC::VectD ag, NC::VectD bg,
                           std::vector<unsigned> zerorows = {} )
  {
    NC::VectD sab;
    for ( auto ib : NC::ncrange( bg.size() ) ) {
      const bool zero = std::find( zerorows.begin(), zerorows.end(),
                                   ib ) != zerorows.end();
      for ( auto a : ag ) {
        const double b = bg.at( ib );
        sab.push_back( zero ? 0.0 : std::exp( -a - b * b / ( 1.0 + a ) )
                       / std::sqrt( 1.0 + a ) );
      }
    }
    return makeKernel( std::move(name), std::move(bg), std::move(ag),
                       std::move(sab) );
  }

  //Reference interpolation of S (log-linear along alpha, or linear if a value
  //is 0, then linear along beta; 0 outside the grid):
  double refS( const Kernel& k, double alpha, double beta )
  {
    const auto& ag = k.ag;
    const auto& bg = k.bg;
    if ( !( alpha >= ag.front() && alpha <= ag.back()
            && beta >= bg.front() && beta <= bg.back() ) )
      return 0.0;
    std::size_t ia = 0, ib = 0;
    while ( ia + 2 < ag.size() && alpha > ag.at( ia + 1 ) )
      ++ia;
    while ( ib + 2 < bg.size() && beta > bg.at( ib + 1 ) )
      ++ib;
    const std::size_t na = ag.size();
    const double ta = ( alpha - ag.at( ia ) )
                      / ( ag.at( ia + 1 ) - ag.at( ia ) );
    auto along_alpha = [&]( std::size_t jb ) {
      const double s1 = k.sab.at( jb * na + ia );
      const double s2 = k.sab.at( jb * na + ia + 1 );
      if ( s1 > 0.0 && s2 > 0.0 )
        return std::exp( std::log( s1 ) + ( std::log( s2 ) - std::log( s1 ) )
                         * ta );
      return s1 + ( s2 - s1 ) * ta;
    };
    const double tb = ( beta - bg.at( ib ) )
                      / ( bg.at( ib + 1 ) - bg.at( ib ) );
    return along_alpha( ib ) * ( 1.0 - tb ) + along_alpha( ib + 1 ) * tb;
  }

  //Alpha limits of KB(e) at beta>=-e, from (alpha-beta)^2<=4*e*alpha, in
  //cancellation-free forms:
  double refAlphaMinus( double e, double beta )
  {
    return NC::ncsquare( beta / ( std::sqrt( e + beta ) + std::sqrt( e ) ) );
  }
  double refAlphaPlus( double e, double beta )
  {
    return NC::ncsquare( std::sqrt( e + beta ) + std::sqrt( e ) );
  }

  //Minimum of f over [x0,x1] by a scan followed by golden-section refinement
  //around the best scan point:
  template<class TFct>
  double numMin( TFct f, double x0, double x1 )
  {
    const unsigned n = 200;
    double best = NC::kInfinity;
    unsigned ibest = 0;
    auto xi = [&]( unsigned i ) {
      return i == n ? x1 : x0 + ( x1 - x0 ) * i / n;
    };
    for ( auto i : NC::ncrange( n + 1 ) ) {
      const double v = f( xi( i ) );
      if ( v < best ) {
        best = v;
        ibest = i;
      }
    }
    double lo = xi( ibest ? ibest - 1 : 0 );
    double hi = xi( ibest < n ? ibest + 1 : n );
    const double g = 0.5 * ( std::sqrt( 5.0 ) - 1.0 );
    for ( auto it : NC::ncrange( 120 ) ) {
      (void)it;
      const double m1 = hi - g * ( hi - lo );
      const double m2 = lo + g * ( hi - lo );
      if ( f( m1 ) < f( m2 ) )
        hi = m2;
      else
        lo = m1;
      best = NC::ncmin( best, NC::ncmin( f( m1 ), f( m2 ) ) );
    }
    return best;
  }

  //Smallest alpha range containing KB(e) for beta in [b0,b1] and all e<=emax
  //(KB(e) grows with e, so this is the alpha range of KB(emax) in the strip):
  NC::PairDD refAlphaRange( double emax, double b0, double b1 )
  {
    const double bb0 = NC::ncmax( b0, -emax );
    REQUIRE( bb0 < b1 );
    return { numMin( [emax]( double b ) { return refAlphaMinus( emax, b ); },
                     bb0, b1 ),
             -numMin( [emax]( double b ) { return -refAlphaPlus( emax, b ); },
                      bb0, b1 ) };
  }

  //Largest beta where KB(emax) has alpha below amax (bisection):
  double refBetaCap( double emax, double amax )
  {
    double lo = 0.0;
    double hi = 1.0;
    while ( refAlphaMinus( emax, hi ) <= amax )
      hi *= 2.0;
    for ( auto it : NC::ncrange( 200 ) ) {
      (void)it;
      const double m = 0.5 * ( lo + hi );
      ( refAlphaMinus( emax, m ) <= amax ? lo : hi ) = m;
    }
    return lo;
  }

  //Maximum of S in a rectangle, on a grid including corners, alpha grid
  //points inside, and interior points:
  double refSMax( const Kernel& k, double a0, double a1, double b0, double b1 )
  {
    NC::VectD agrid;
    const unsigned na = 10;
    for ( auto i : NC::ncrange( na + 1 ) )
      agrid.push_back( i == na ? a1 : a0 + ( a1 - a0 ) * i / na );
    for ( auto a : k.ag )
      if ( a > a0 && a < a1 )
        agrid.push_back( a );
    double smax = 0.0;
    const unsigned nb = 4;
    for ( auto j : NC::ncrange( nb + 1 ) ) {
      const double b = ( j == nb ? b1 : b0 + ( b1 - b0 ) * j / nb );
      for ( auto a : agrid )
        smax = NC::ncmax( smax, refS( k, a, b ) );
    }
    return smax;
  }

  bool approxEq( double a, double b, double scale )
  {
    return NC::ncabs( a - b ) <= 1e-12 * scale;
  }

  struct Counts { unsigned nconf = 0, nused = 0, ndropped = 0, nthrow = 0; };

  //Reference strips (independent of S): the beta cells in the beta range of
  //KB(emax) on the kernel grid, each split in nstrips, with the smallest alpha
  //range containing KB(e) for all e<=emax:
  struct RefStrip { double b0, b1, a0, a1; };
  std::vector<RefStrip> refStrips( const NC::VectD& ag, const NC::VectD& bg,
                                   double emax, unsigned nstrips )
  {
    std::vector<RefStrip> res;
    const double blo = NC::ncmax( -emax, bg.front() );
    const double bhi = NC::ncmin( bg.back(), refBetaCap( emax, ag.back() ) );
    if ( !( bhi > blo ) )
      return res;
    NC::VectD cellpts{ blo };
    for ( auto b : bg )
      if ( b > blo && b < bhi )
        cellpts.push_back( b );
    cellpts.push_back( bhi );
    for ( auto i : NC::ncrange( cellpts.size() - 1 ) ) {
      const double c0 = cellpts.at( i );
      const double c1 = cellpts.at( i + 1 );
      for ( auto j : NC::ncrange( nstrips ) ) {
        RefStrip rs;
        rs.b0 = c0 + ( c1 - c0 ) * j / nstrips;
        rs.b1 = ( j + 1 == nstrips ? c1
                  : c0 + ( c1 - c0 ) * ( j + 1 ) / nstrips );
        const auto ar = refAlphaRange( emax, rs.b0, rs.b1 );
        rs.a0 = ar.first;
        rs.a1 = ar.second;
        res.push_back( rs );
      }
    }
    return res;
  }

  //Whether S>0 somewhere in the reference rectangle, shrunk by 1e-12 (the
  //accuracy of the reference rectangles):
  bool refHasS( const Kernel& k, const RefStrip& rs )
  {
    const double da = 1e-12 * ( rs.a1 - rs.a0 );
    const double db = 1e-12 * ( rs.b1 - rs.b0 );
    return refSMax( k, rs.a0 + da, rs.a1 - da, rs.b0 + db, rs.b1 - db ) > 0.0;
  }

  //Checks the strips of a SABBetaStripSampler on kernel k against the
  //reference strips: same beta ranges and rectangles, and the reference S
  //maximum in the rectangle. Strips must only be dropped where S=0, and a
  //sampler without any strips must refuse construction:
  void checkStrips( const Kernel& k, double emax, unsigned nstrips,
                    const std::vector<RefStrip>& expected, Counts& counts )
  {
    ++counts.nconf;
    NC::Optional<NS::SABBetaStripSampler> sampler;
    try {
      sampler.emplace( k.sd, emax, nstrips );
    } catch ( NC::Error::BadInput& ) {
      for ( auto& rs : expected )
        REQUIRE( !refHasS( k, rs ) );
      ++counts.nthrow;
      return;
    }
    const auto strips = sampler.value().stripBounds();
    std::size_t istrip = 0;
    for ( auto& rs : expected ) {
      const bool used = ( istrip < strips.size()
                          && approxEq( strips.at( istrip ).b0, rs.b0,
                                       NC::ncmax( 1.0, NC::ncabs( rs.b0 ) ) ) );
      if ( !used ) {
        REQUIRE( !refHasS( k, rs ) );
        ++counts.ndropped;
        continue;
      }
      const auto& s = strips.at( istrip++ );
      ++counts.nused;
      const double ascale = NC::ncmax( rs.a1, 1e-300 );
      REQUIRE( approxEq( s.b1, rs.b1, NC::ncmax( 1.0, NC::ncabs( rs.b1 ) ) ) );
      REQUIRE( approxEq( s.a0, rs.a0, ascale ) );
      REQUIRE( approxEq( s.a1, rs.a1, ascale ) );
      REQUIRE( s.smax > 0.0 );
      REQUIRE( NC::floateq( s.smax, refSMax( k, s.a0, s.a1, s.b0, s.b1 ),
                            1e-10, 0.0 ) );
    }
    REQUIRE( istrip == strips.size() );
  }

  void checkStrips( const Kernel& k, double emax, unsigned nstrips,
                    Counts& counts )
  {
    checkStrips( k, emax, nstrips, refStrips( k.ag, k.bg, emax, nstrips ),
                 counts );
  }

  //Position of v relative to the sorted special values s1<s2: 0 below s1,
  //1 on s1, 2 between, 3 on s2, 4 above s2:
  unsigned position( double v, double s1, double s2 )
  {
    return v < s1 ? 0 : v == s1 ? 1 : v < s2 ? 2 : v == s2 ? 3 : 4;
  }
}

int main()
{
  //1) Grid cells in all topologies relative to the special points of KB(e)
  //(beta=-e, beta=0, alpha=0, alpha=4e), with every cell edge below, on or
  //above each of them, and S at the cell corners in all combinations of 0,
  //1e-300, 1 and 1e300. Grids are built from these candidate points, where
  //emax=0.5 puts candidates exactly on -e=-0.5, 0 and 4e=2 (grids must start
  //below beta=0, so cells above get an extra cell in front with S=1):
  const NC::VectD bcand = { -2.0, -1.0, -0.5, -0.25, -0.1, 0.0, 0.5, 3.0 };
  const NC::VectD acand = { 0.0, 0.25, 1.0, 2.0, 3.0, 6.0 };
  struct Topo { NC::VectD bg, ag; };
  std::map<std::array<unsigned,4>,Topo> topos;
  for ( auto i : NC::ncrange( bcand.size() ) ) {
    for ( auto j : NC::ncrange( i + 1, bcand.size() ) ) {
      const double b0 = bcand.at( i );
      const double b1 = bcand.at( j );
      NC::VectD bgr{ b0, b1 };
      if ( !( b0 < 0.0 ) )
        bgr.insert( bgr.begin(), -2.0 );
      for ( auto l : NC::ncrange( acand.size() ) ) {
        for ( auto m : NC::ncrange( l + 1, acand.size() ) ) {
          const double a0 = acand.at( l );
          const double a1 = acand.at( m );
          topos.emplace( std::array<unsigned,4>{ position( b0, -0.5, 0.0 ),
                                                 position( b1, -0.5, 0.0 ),
                                                 position( a0, 0.0, 2.0 ),
                                                 position( a1, 0.0, 2.0 ) },
                         Topo{ bgr, { a0, a1 } } );
        }
      }
    }
  }
  //All topologies: beta edge pairs from positions 0..4 (13 pairs, as both
  //edges can not be on the same point), alpha edge pairs from positions 1..4
  //(8 pairs):
  REQUIRE( topos.size() == 13 * 8 );
  Counts ctopo;
  const NC::VectD svals = { 0.0, 1e-300, 1.0, 1e300 };
  const std::vector<std::pair<double,unsigned>> otherconfs
    = { {0.3,1}, {0.3,3}, {0.8,1}, {0.8,3}, {1e-3,1}, {1e-3,3} };
  for ( auto& tp : topos ) {
    const auto& bgr = tp.second.bg;
    const auto& agr = tp.second.ag;
    const auto expected = refStrips( agr, bgr, 0.5, 3 );
    std::vector<std::vector<RefStrip>> expected_other;
    for ( auto& oc : otherconfs )
      expected_other.push_back( refStrips( agr, bgr, oc.first, oc.second ) );
    for ( auto ip : NC::ncrange( 256 ) ) {
      //Corner values (last two beta rows), S=1 in any extra row in front:
      NC::VectD sab( bgr.size() * 2, 1.0 );
      for ( auto c : NC::ncrange( 4 ) )
        sab.at( sab.size() - 4 + c ) = svals.at( ( ip >> ( 2 * c ) ) & 3 );
      const auto k = makeKernel( "topo", bgr, agr, sab );
      checkStrips( k, 0.5, 3, expected, ctopo );
      //Also with the special points elsewhere, for corner S of 0 and 1:
      if ( ( ip & 0x55 ) == 0 )
        for ( auto i : NC::ncrange( otherconfs.size() ) )
          checkStrips( k, otherconfs.at( i ).first, otherconfs.at( i ).second,
                       expected_other.at( i ), ctopo );
    }
  }
  std::cout << "Cell topologies: " << topos.size() << " with "
            << svals.size() << "^4 corner S combinations, " << ctopo.nconf
            << " configurations, " << ctopo.nused << " strips as expected, "
            << ctopo.ndropped << " dropped with S=0, " << ctopo.nthrow
            << " samplers refused without any S>0" << std::endl;

  //2) Larger smooth kernels over a wide range of emax:
  const NC::VectD ag = geom( 1e-3, 50.0, 25 );
  NC::VectD bsym = geom( 1e-3, 20.0, 12 );
  {
    NC::VectD neg;
    for ( auto b : bsym )
      neg.push_back( -b );
    std::reverse( neg.begin(), neg.end() );
    neg.push_back( 0.0 );
    bsym.insert( bsym.begin(), neg.begin(), neg.end() );
  }
  //Asymmetric, with few points below 0:
  NC::VectD basym = geom( 1e-3, 20.0, 12 );
  basym.insert( basym.begin(), { -0.7, -2e-3, 0.0 } );
  std::vector<Kernel> kernels = {
    makeSmoothKernel( "symmetric beta grid", ag, bsym ),
    makeSmoothKernel( "asymmetric beta grid", ag, basym ),
    makeSmoothKernel( "beta grid without 0", ag,
                      { -7.3, -2.2, -0.61, -0.05, 0.13, 0.9, 3.3, 11.0 } ),
    makeSmoothKernel( "short alpha grid", geom( 1e-3, 0.5, 10 ),
                      { -1.0, 0.0, 0.3, 2.0, 8.0, 30.0 } ),
    makeSmoothKernel( "zero beta rows", ag, bsym,
                      { 0, 1, 2, 3, 15, 16, 17 } ),
  };
  for ( auto& k : kernels ) {
    Counts c;
    for ( double emax : { 1e-10, 1e-5, 1e-3, 0.05, 0.7, 3.0, 15.0 } )
      for ( unsigned nstrips : { 1, 3, 4 } )
        checkStrips( k, emax, nstrips, c );
    std::cout << k.name << ": " << c.nconf << " configurations, " << c.nused
              << " strips as expected, " << c.ndropped
              << " dropped with S=0" << std::endl;
  }
  return 0;
}
