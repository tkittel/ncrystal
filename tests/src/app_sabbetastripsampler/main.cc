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
#include "NCrystal/internal/utils/NCStableDbl.hh"
#include "NCrystal/interfaces/NCRNG.hh"
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

//With 0, the sampling of the cell topologies uses fewer samples and only
//e=emax, for a faster test. FIXME: Consider disabling this by default once
//the beta-strip sampler code has settled:
#define NCRYSTAL_MORE_THOROUGH_TESTS 1

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

  struct Counts { unsigned nconf = 0, nused = 0, ndropped = 0, nrefused = 0; };

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
  //sampler without any strips must not be created:
  void checkStrips( const Kernel& k, double emax, unsigned nstrips,
                    const std::vector<RefStrip>& expected, Counts& counts )
  {
    ++counts.nconf;
    NS::SABBetaStripSampler sampler( k.sd, emax, nstrips );
    if ( !sampler.hasStrips() ) {
      for ( auto& rs : expected )
        REQUIRE( !refHasS( k, rs ) );
      ++counts.nrefused;
      return;
    }
    const auto strips = sampler.stripBounds();
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

  //Integral of the reference S over alpha in [x0,x1] at fixed beta, in
  //closed form: along each cell edge S is s1*exp(L*u) (L=ln(s2/s1), u the
  //relative position in the alpha cell) or linear, and it is linear in beta
  //between the two edges. Exponentials are evaluated in log space, so even
  //values like 1e-300 and 1e300 on one edge can not overflow:
  double refAlphaIntegralClosed( const Kernel& k, double beta,
                                 double x0, double x1 )
  {
    const auto& ag = k.ag;
    const auto& bg = k.bg;
    if ( !( beta >= bg.front() && beta <= bg.back() ) )
      return 0.0;
    x0 = NC::ncmax( x0, ag.front() );
    x1 = NC::ncmin( x1, ag.back() );
    if ( !( x1 > x0 ) )
      return 0.0;
    std::size_t ib = 0;
    while ( ib + 2 < bg.size() && beta > bg.at( ib + 1 ) )
      ++ib;
    const double tb = ( beta - bg.at( ib ) )
                      / ( bg.at( ib + 1 ) - bg.at( ib ) );
    const std::size_t na = ag.size();
    NC::StableDbl sum;
    //Alpha cells overlapping [x0,x1]:
    const std::size_t ia_begin = static_cast<std::size_t>(
      std::upper_bound( ag.begin(), ag.end(), x0 ) - ag.begin() ) - 1;
    for ( std::size_t ia = ia_begin; ia + 1 < na && ag.at( ia ) < x1; ++ia ) {
      const double c0 = ag.at( ia );
      const double c1 = ag.at( ia + 1 );
      const double p0 = NC::ncmax( x0, c0 );
      const double p1 = NC::ncmin( x1, c1 );
      if ( !( p1 > p0 ) )
        continue;
      const double d = c1 - c0;
      const double u0 = ( p0 - c0 ) / d;
      const double u1 = ( p1 - c0 ) / d;
      for ( auto jb : { ib, ib + 1 } ) {
        const double wb = ( jb == ib ? 1.0 - tb : tb );
        const double s1 = k.sab.at( jb * na + ia );
        const double s2 = k.sab.at( jb * na + ia + 1 );
        double v;
        if ( s1 > 0.0 && s2 > 0.0 ) {
          //integral of s1*exp(L*u)*d over [u0,u1] is d*(S(u1)-S(u0))/L with
          //S(u)=exp(ln(s1)+L*u)<=max(s1,s2). For |L*du|<1 use instead the
          //form d*S(u0)*expm1(L*du)/L, which avoids cancellation:
          const double ls1 = std::log( s1 );
          const double L = std::log( s2 ) - ls1;
          const double du = u1 - u0;
          const double x = L * du;
          if ( NC::ncabs( x ) < 1.0 )
            v = d * std::exp( ls1 + L * u0 )
              * ( NC::ncabs( x ) < 1e-12 ? du : std::expm1( x ) / L );
          else
            v = d * ( std::exp( ls1 + L * u1 ) - std::exp( ls1 + L * u0 ) )
              / L;
        } else {
          v = d * ( s1 * ( u1 - u0 )
                    + ( s2 - s1 ) * 0.5 * ( u1 * u1 - u0 * u0 ) );
        }
        sum += NC::StableDbl( wb * v );
      }
    }
    return sum.value();
  }

  //Reference integral of S over the part of KB(e) with beta in [b0,b1] and
  //relative alpha position u=(alpha-alpha-)/(alpha+-alpha-) in [f0,f1].
  //
  //With t=sqrt(beta+e) we have alpha+-=(t+-sqrt(e))^2 and dbeta=2t*dt, so the
  //u-bin edges x_f(t)=(t-sqrt(e))^2+4*f*sqrt(e)*t are quadratic in t. The
  //integrand in t (the closed-form alpha integral above) is smooth and varies
  //by less than a factor ~exp(1) on each piece, when splitting at: beta grid
  //points, and where a u-bin edge crosses an alpha grid point or one of the
  //alpha values dividing each alpha cell into parts where S changes by at
  //most a factor exp(1) along the edges (all in closed form). Each piece is
  //integrated with fixed-order Romberg integration (33 points), verified
  //against 17 points.
  double refBinIntegral( const Kernel& k, double e, double b0, double b1,
                         double f0, double f1 )
  {
    const double sqrte = std::sqrt( e );
    const double t0 = std::sqrt( NC::ncmax( 0.0, b0 + e ) );
    const double t1 = std::sqrt( NC::ncmax( 0.0, b1 + e ) );
    if ( !( t1 > t0 ) )
      return 0.0;
    //alpha values where the integrand may change behaviour:
    const auto& ag = k.ag;
    const std::size_t na = ag.size();
    NC::VectD alevels( ag.begin(), ag.end() );
    for ( auto ia : NC::ncrange( na - 1 ) ) {
      double lmax = 0.0;
      for ( auto jb : NC::ncrange( k.bg.size() ) ) {
        const double s1 = k.sab.at( jb * na + ia );
        const double s2 = k.sab.at( jb * na + ia + 1 );
        if ( s1 > 0.0 && s2 > 0.0 )
          lmax = NC::ncmax( lmax, NC::ncabs( std::log( s2 )
                                             - std::log( s1 ) ) );
      }
      const unsigned m = static_cast<unsigned>( std::ceil( lmax ) );
      for ( unsigned j = 1; j < m; ++j )
        alevels.push_back( ag.at( ia ) + ( ag.at( ia + 1 ) - ag.at( ia ) )
                           * j / m );
    }
    NC::VectD tbreaks{ t0, t1 };
    auto addBreak = [&tbreaks,t0,t1]( double tt )
    {
      if ( tt > t0 && tt < t1 )
        tbreaks.push_back( tt );
    };
    for ( auto b : k.bg )
      if ( b > -e )
        addBreak( std::sqrt( b + e ) );
    for ( double f : { f0, f1 } ) {
      for ( auto a : alevels ) {
        //x_f(t)=a <=> t = (1-2f)*sqrt(e) +- sqrt(a-4f(1-f)e):
        const double d = a - 4.0 * f * ( 1.0 - f ) * e;
        if ( d >= 0.0 ) {
          addBreak( ( 1.0 - 2.0 * f ) * sqrte + std::sqrt( d ) );
          addBreak( ( 1.0 - 2.0 * f ) * sqrte - std::sqrt( d ) );
        }
      }
    }
    std::sort( tbreaks.begin(), tbreaks.end() );
    tbreaks.erase( std::unique( tbreaks.begin(), tbreaks.end() ),
                   tbreaks.end() );
    //Pieces are inside or outside the beta grid. Beta is clamped to
    //[b0,b1] (which never contains beta grid points other than at its ends),
    //so rounding of t*t-e at the piece ends can not pick a value from a
    //neighbouring beta cell or from outside the grid:
    const double bgmin = k.bg.front();
    const double bgmax = k.bg.back();
    auto integrand = [&k,e,sqrte,f0,f1,b0,b1]( double tt )
    {
      const double am = NC::ncsquare( tt - sqrte );
      const double w = 4.0 * sqrte * tt;
      return 2.0 * tt
        * refAlphaIntegralClosed( k, NC::ncclamp( tt * tt - e, b0, b1 ),
                                  am + f0 * w, am + f1 * w );
    };
    NC::StableDbl sum33, sum17;
    for ( auto i : NC::ncrange( tbreaks.size() - 1 ) ) {
      const double ta = tbreaks.at( i );
      const double tb = tbreaks.at( i + 1 );
      const double tmid = 0.5 * ( ta + tb );
      if ( !( tmid * tmid - e >= bgmin && tmid * tmid - e <= bgmax ) )
        continue;//outside the beta grid
      sum33 += NC::StableDbl( NC::integrateRomberg33( integrand, ta, tb ) );
      sum17 += NC::StableDbl( NC::integrateRomberg17( integrand, ta, tb ) );
    }
    //Relative tolerance, or an absolute one at the scale of the integral (for
    //results which are essentially 0):
    const double scale = *std::max_element( k.sab.begin(), k.sab.end() )
      * ( b1 - b0 ) * 4.0 * std::sqrt( e * ( b1 + e ) );
    REQUIRE( NC::floateq( sum33.value(), sum17.value(), 1e-6,
                          1e-10 * scale ) );
    return sum33.value();
  }

  //Verifies refAlphaIntegralClosed against numerical quadrature of refS,
  //split at alpha grid points and into parts where S changes by at most a
  //factor exp(1):
  void checkClosedAlphaIntegral( const Kernel& k, double beta,
                                 double x0, double x1 )
  {
    NC::VectD pts{ x0 };
    for ( auto a : k.ag )
      if ( a > x0 && a < x1 )
        pts.push_back( a );
    pts.push_back( x1 );
    NC::StableDbl num33, num17;
    //Parts are inside or outside the alpha grid (where S=0). For those
    //inside, alpha is clamped to the grid against rounding at the ends:
    const double agmin = k.ag.front();
    const double agmax = k.ag.back();
    auto f = [&k,beta,agmin,agmax]( double a )
    {
      return refS( k, NC::ncclamp( a, agmin, agmax ), beta );
    };
    for ( auto ip : NC::ncrange( pts.size() - 1 ) ) {
      const double a0 = pts.at( ip );
      const double a1 = pts.at( ip + 1 );
      if ( !( a0 >= agmin && a1 <= agmax ) )
        continue;//outside the alpha grid
      const double s0 = f( a0 );
      const double s1 = f( a1 );
      unsigned n = 1;
      if ( s0 > 0.0 && s1 > 0.0 )
        n += static_cast<unsigned>( NC::ncabs( std::log( s1 )
                                               - std::log( s0 ) ) );
      for ( auto j : NC::ncrange( n ) ) {
        const double pa = a0 + ( a1 - a0 ) * j / n;
        const double pb = ( j + 1 == n ? a1
                            : a0 + ( a1 - a0 ) * ( j + 1 ) / n );
        num33 += NC::StableDbl( NC::integrateRomberg33( f, pa, pb ) );
        num17 += NC::StableDbl( NC::integrateRomberg17( f, pa, pb ) );
      }
    }
    //Fixed-order Romberg (33 points), verified against 17 points:
    REQUIRE( NC::floateq( num33.value(), num17.value(), 1e-6, 0.0 ) );
    REQUIRE( NC::floateq( refAlphaIntegralClosed( k, beta, x0, x1 ),
                          num33.value(), 1e-9, 0.0 ) );
  }

  //Samples (alpha,beta) at e and compares with the distribution of the
  //reference S within KB(e), in bins of beta (fine bins grouped to have
  //similar expected counts) and of the relative alpha position u in
  //[alpha-,alpha+], with a chi-square test. Also checks that the mean number
  //of tries agrees with expectedTries. Sampling is skipped where it is not
  //possible (no S>0 within KB(e)) or very inefficient (expectedTries>200).
  //Returns whether samples were checked.
  bool checkSampling( const Kernel& kernel, double emax, double e,
                      unsigned nstrips, unsigned nsamples,
                      unsigned nfine_per_cell = 4 )
  {
    NS::SABBetaStripSampler sampler( kernel.sd, emax, nstrips );
    if ( !sampler.hasStrips() )
      return false;//no strips at all (tested elsewhere)
    const auto& ag = kernel.ag;
    const auto& bg = kernel.bg;
    constexpr unsigned nu = 4;

    //Fine beta bins, splitting beta cells (and beta=0) within the beta range
    //of KB(e) (S=0 outside):
    const double blo = NC::ncmax( -e, bg.front() );
    const double bhi = NC::ncmin( bg.back(), refBetaCap( e, ag.back() ) );
    if ( !( bhi > blo ) )
      return false;
    NC::VectD cellpts{ blo };
    for ( auto b : bg )
      if ( b > blo && b < bhi )
        cellpts.push_back( b );
    if ( blo < 0.0 && bhi > 0.0
         && std::find( cellpts.begin(), cellpts.end(), 0.0 ) == cellpts.end() )
      cellpts.push_back( 0.0 );
    cellpts.push_back( bhi );
    std::sort( cellpts.begin(), cellpts.end() );
    NC::VectD fine;
    for ( auto i : NC::ncrange( cellpts.size() - 1 ) )
      for ( auto j : NC::ncrange( nfine_per_cell ) )
        fine.push_back( cellpts.at( i ) + ( cellpts.at( i + 1 )
                                            - cellpts.at( i ) )
                        * j / nfine_per_cell );
    fine.push_back( bhi );

    //Reference: integral of S over KB(e) in each (fine beta bin, u bin):
    const std::size_t nfine = fine.size() - 1;
    std::vector<std::array<double,nu>> ref( nfine );
    NC::StableSum reftot;
    for ( auto i : NC::ncrange( nfine ) ) {
      for ( auto iu : NC::ncrange( nu ) ) {
        const double f0 = double( iu ) / nu;
        const double f1 = double( iu + 1 ) / nu;
        const double v = refBinIntegral( kernel, e, fine.at( i ),
                                         fine.at( i + 1 ), f0, f1 );
        ref.at( i ).at( iu ) = v;
        reftot.add( v );
      }
    }
    const double integralWithinKB = reftot.sum();
    if ( !( integralWithinKB > 0.0 ) )
      return false;//nothing to sample
    const double expTries = sampler.expectedTries( e, integralWithinKB );
    if ( expTries > 200.0 )
      return false;//S bounds of strips far above S within KB(e)

    //Sample:
    auto rng = NC::getRNG();
    std::vector<std::array<double,nu>> counts( nfine );
    for ( auto& c : counts )
      c.fill( 0.0 );
    unsigned long ntries = 0;
    for ( auto i : NC::ncrange( nsamples ) ) {
      (void)i;
      const auto r = sampler.sample( *rng, e );
      REQUIRE( r.ntries >= 1 );
      ntries += r.ntries;
      REQUIRE( r.beta >= blo && r.beta <= bhi );
      const NC::PairDD al{ refAlphaMinus( e, r.beta ),
                           refAlphaPlus( e, r.beta ) };
      const double atol = 1e-12 * al.second;
      REQUIRE( r.alpha >= al.first - atol && r.alpha <= al.second + atol );
      REQUIRE( r.alpha >= ag.front() && r.alpha <= ag.back() );
      REQUIRE( refS( kernel, r.alpha, r.beta ) > 0.0 );
      const auto ib = std::upper_bound( fine.begin(), fine.end(), r.beta )
                      - fine.begin();
      const std::size_t ifine = std::min<std::size_t>(
        nfine - 1, std::size_t( std::max<std::ptrdiff_t>( ib, 1 ) - 1 ) );
      const double u = ( al.second > al.first
                         ? ( r.alpha - al.first ) / ( al.second - al.first )
                         : 0.5 );
      const unsigned iu = NC::ncmin( nu - 1, unsigned( u * nu ) );
      counts.at( ifine ).at( iu ) += 1.0;
    }

    //Group fine beta bins to have >=~nsamples/15 expected counts each, and
    //compute chi-square over all (group,u) bins, pooling bins with fewer than
    //5 expected counts into one:
    const double target = integralWithinKB / 15.0;
    double chi2 = 0.0;
    unsigned nbins = 0;
    double pool_exp = 0.0, pool_cnt = 0.0;
    std::array<double,nu> gref, gcnt;
    gref.fill( 0.0 );
    gcnt.fill( 0.0 );
    double gsum = 0.0;
    auto addBin = [&]( double expected, double cnt )
    {
      chi2 += NC::ncsquare( cnt - expected ) / expected;
      ++nbins;
    };
    auto flush = [&]()
    {
      for ( auto iu : NC::ncrange( nu ) ) {
        const double expected = gref.at( iu ) / integralWithinKB * nsamples;
        if ( expected >= 5.0 ) {
          addBin( expected, gcnt.at( iu ) );
        } else {
          pool_exp += expected;
          pool_cnt += gcnt.at( iu );
        }
      }
      gref.fill( 0.0 );
      gcnt.fill( 0.0 );
      gsum = 0.0;
    };
    for ( auto i : NC::ncrange( nfine ) ) {
      for ( auto iu : NC::ncrange( nu ) ) {
        gref.at( iu ) += ref.at( i ).at( iu );
        gcnt.at( iu ) += counts.at( i ).at( iu );
        gsum += ref.at( i ).at( iu );
      }
      if ( gsum >= target || i + 1 == nfine )
        flush();
    }
    if ( pool_exp >= 5.0 )
      addBin( pool_exp, pool_cnt );
    else
      REQUIRE( pool_cnt <= pool_exp + 5.0 * std::sqrt( pool_exp ) + 5.0 );
    //Generous threshold (p-value below ~1e-9), as many configurations are
    //tested:
    const double ndf = NC::ncmax( 1.0, nbins - 1.0 );
    if ( !( chi2 < ndf + 10.0 * std::sqrt( 2.0 * ndf ) + 20.0 ) ) {
      std::cout << "chi2 test failed: chi2/ndf = " << chi2 << "/" << ndf
                << std::endl;
      REQUIRE( false );
    }

    //Mean number of tries (geometric distribution) vs. expectedTries:
    const double meanTries = double( ntries ) / nsamples;
    const double sdev = std::sqrt( expTries * ( expTries - 1.0 )
                                   / nsamples );
    if ( !( NC::ncabs( meanTries - expTries )
            < 6.0 * sdev + 2e-3 * expTries ) ) {
      std::cout << "tries test failed: mean " << meanTries << " expected "
                << expTries << std::endl;
      REQUIRE( false );
    }
    return true;
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
            << ctopo.ndropped << " dropped with S=0, " << ctopo.nrefused
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

  //The closed-form reference alpha integral, against numerical quadrature,
  //for the smooth kernels and for a single cell with extreme corner S values:
  for ( auto& k : kernels )
    for ( double beta : { -0.37, 0.0, 0.004, 1.3, 7.7 } )
      for ( auto xr : { NC::PairDD{ 0.0, 60.0 }, NC::PairDD{ 1e-3, 0.37 },
                        NC::PairDD{ 0.02, 0.0201 }, NC::PairDD{ 3.3, 17.0 } } )
        checkClosedAlphaIntegral( k, beta, xr.first, xr.second );
  for ( auto& pat : std::vector<std::array<double,4>>{
      { 1.0, 1e300, 1e-300, 1.0 }, { 1e-300, 1e300, 0.0, 1e-300 },
      { 0.0, 1e300, 1e300, 0.0 }, { 1e-300, 1e-300, 1e300, 1e300 } } ) {
    const auto k = makeKernel( "cell", { -1.0, 1.0 }, { 0.5, 2.0 },
                               { pat.at(0), pat.at(1), pat.at(2), pat.at(3) } );
    for ( double beta : { -1.0, -0.3, 0.9, 1.0 } )
      for ( auto xr : { NC::PairDD{ 0.0, 3.0 }, NC::PairDD{ 0.6, 1.9 },
                        NC::PairDD{ 1.99, 2.0 } } )
        checkClosedAlphaIntegral( k, beta, xr.first, xr.second );
  }
  std::cout << "Closed-form reference alpha integrals agree with"
            << " numerical quadrature" << std::endl;

  //Sampled distributions (nstrips=1 and 4) at e=emax and far below:
  for ( auto& k : kernels ) {
    for ( double emax : { 1e-4, 0.05 } ) {
      for ( double efrac : { 1.0, 0.01 } ) {
        for ( unsigned nstrips : { 1, 4 } ) {
          REQUIRE( checkSampling( k, emax, efrac * emax, nstrips, 30000 ) );
          std::cout << k.name << ": sampling at emax=" << emax
                    << ", e/emax=" << efrac << ", nstrips=" << nstrips
                    << " consistent with S" << std::endl;
        }
      }
    }
  }

  //Cell with an edge needing large-k treatment (S=1e-30 and 1) next to a
  //linear one (S=0 and 1), with alpha1 just below alpha-(1,0.5):
  {
    const auto k = makeKernel( "large-k cell", { -0.25, 0.5, 1.0 },
                               { 0.0505102572167269, 1.0 },
                               { 1.0, 1.0, 1e-30, 1.0, 0.0, 1.0 } );
    for ( unsigned nstrips : { 1, 4 } ) {
      REQUIRE( checkSampling( k, 1.0, 1.0, nstrips, 30000 ) );
      std::cout << k.name << ": sampling at emax=1, e/emax=1, nstrips="
                << nstrips << " consistent with S" << std::endl;
    }
  }

  //The cell topologies from above, with corner S values of 1, or with one
  //corner 0, also at e<emax, and with all corners 1e300 or 1e-300 (the
  //geometry checks above also cover single extreme corners):
  {
    std::vector<std::pair<std::array<double,4>,bool>> patterns;//(S,also e<emax)
    patterns.push_back( { { 1.0, 1.0, 1.0, 1.0 }, true } );
    for ( auto c : NC::ncrange( 4 ) ) {
      std::array<double,4> pat{ 1.0, 1.0, 1.0, 1.0 };
      pat.at( c ) = 0.0;
      patterns.push_back( { pat, true } );
    }
    for ( double v : { 1e300, 1e-300 } )
      patterns.push_back( { { v, v, v, v }, false } );
    constexpr bool thorough = ( NCRYSTAL_MORE_THOROUGH_TESTS != 0 );
    constexpr unsigned nsamples = ( thorough ? 1000 : 500 );
    unsigned nconf = 0, nsampled = 0;
    for ( auto& tp : topos ) {
      const auto& bgr = tp.second.bg;
      for ( auto& pat : patterns ) {
        NC::VectD sab( bgr.size() * 2, 1.0 );
        for ( auto c : NC::ncrange( 4 ) )
          sab.at( sab.size() - 4 + c ) = pat.first.at( c );
        const auto k = makeKernel( "topo", bgr, tp.second.ag, sab );
        nsampled += checkSampling( k, 0.5, 0.5, 3, nsamples, 1 ) ? 1 : 0;
        ++nconf;
        if ( thorough && pat.second ) {
          nsampled += checkSampling( k, 0.5, 0.1, 3, nsamples, 1 ) ? 1 : 0;
          ++nconf;
        }
      }
    }
    //Many configurations can not be sampled (no S>0 within KB(e)) or are
    //skipped as very inefficient, but most must have been checked (744 and 503
    //when this was written, with some margin for platform differences):
    REQUIRE( nsampled >= ( thorough ? 700u : 470u ) );
    std::cout << "Cell topologies: " << nconf << " sampling configurations"
              << " consistent with S" << std::endl;
  }
  return 0;
}
