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


#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include <iostream>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  static_assert( !std::is_copy_constructible<NC::IsoXSDescr>::value, "" );
  static_assert( !std::is_copy_constructible<NC::IsoXSDescr::Tabulation>::value, "" );
  static_assert( std::is_move_constructible<NC::IsoXSDescr>::value, "" );

  std::uint64_t g_nchecks = 0;

  //Straightforward reference implementation of the formulas in the
  //NCIsoXSUtils.hh comments:
  double refXS( const NC::IsoXSDescr::Tabulation& p, double e )
  {
    if ( e < p.E.front() ) {
      const double ee = ( e == 0.0 ? std::numeric_limits<double>::denorm_min() : e );
      return ( p.tail[0] ? p.tail[0] / std::sqrt( ee ) : 0.0 ) + p.tail[1];
    }
    if ( !( e < p.E.back() ) )
      return p.tail[3] + ( p.tail[2] + p.tail[4] / e ) / e;
    std::size_t i = static_cast<std::size_t>( std::upper_bound( p.E.begin(), p.E.end(), e ) - p.E.begin() ) - 1;
    const double u = ( e - p.E[i] ) * ( 1.0 / ( p.E[i+1] - p.E[i] ) );
    auto hor = [u]( const double * c ) { return NC::ncfma( u, NC::ncfma( u, NC::ncfma( u, c[3], c[2] ), c[1] ), c[0] ); };
    double exs = p.p.empty() ? 0.0 : hor( p.p.data() + 4*i );
    if ( !p.q.empty() )
      exs = NC::ncfma( std::sqrt( e ), hor( p.q.data() + 4*i ), exs );
    if ( p.r.empty() )
      return exs / e;
    const double inv = 1.0 / e;
    return NC::ncfma( p.r[i], inv, exs ) * inv;
  }

  bool same( double a, double b )
  {
    return a == b || ( std::isnan(a) && std::isnan(b) );
  }

  std::uint64_t dblBits( double x ) { std::uint64_t b; std::memcpy( &b, &x, sizeof(b) ); return b; }
  double fromBits( std::uint64_t b ) { double x; std::memcpy( &x, &b, sizeof(x) ); return x; }

  //Coefficient styles: 0: uniform in [0,1); 1: random signs and magnitudes
  //1e-30..1e30; 2: like 0, but q and r zero in a random subset of intervals.
  void checkPieces( const char * label, NC::VectD E, bool withq, bool withr,
                    int coefstyle, NC::RandXRSRImpl& rng )
  {
    auto coef = [&rng,coefstyle]() {
      if ( coefstyle != 1 )
        return rng.generate();
      const double sign = ( rng.generate() < 0.5 ? -1.0 : 1.0 );
      return sign * std::pow( 10.0, -30.0 + 60.0 * rng.generate() );
    };
    NC::IsoXSDescr d;
    d.breakpts = { E.front(), E.back() };
    NC::IsoXSDescr::Tabulation p;
    p.E = std::move(E);
    for ( std::size_t i = 0; i + 1 < p.E.size(); ++i ) {
      const bool zeroqr = ( coefstyle == 2 && rng.generate() < 0.5 );
      for ( int k = 0; k < 4; ++k ) {
        p.p.push_back( coef() );
        if ( withq )
          p.q.push_back( zeroqr ? 0.0 : coef() );
      }
      if ( withr )
        p.r.push_back( zeroqr ? 0.0 : coef() * p.E[i] );
    }
    for ( auto& c : p.tail )
      c = coef();
    d.table = std::move(p);
    const NC::IsoXSEvaluator ev( std::move(d) );
    const auto& pp = ev.descr().table.value();
    {
      //cloneIsoXSDescr gives an identical description:
      const auto cp = NC::cloneIsoXSDescr( ev.descr() );
      const auto& ct = cp.table.value();
      REQUIRE( cp.breakpts == ev.descr().breakpts );
      REQUIRE( ct.E == pp.E && ct.p == pp.p && ct.q == pp.q && ct.r == pp.r );
      for ( int k = 0; k < 5; ++k )
        REQUIRE( dblBits( ct.tail[k] ) == dblBits( pp.tail[k] ) );
      REQUIRE( !NC::cloneIsoXSDescr( NC::IsoXSDescr() ).table.has_value() );
    }
    auto check = [&]( double e ) {
      if ( !( e >= 0.0 ) )
        return;
      const double ref = refXS( pp, e );
      REQUIRE( same( NC::evalIsoXS( pp, e ), ref ) );
      REQUIRE( same( ev.eval( e ), ref ) );
      ++g_nchecks;
    };
    const NC::VectD& E2 = pp.E;
    //Breakpoints, their neighbours, and points inside each interval:
    for ( auto e : E2 ) {
      check( e );
      check( std::nextafter( e, 0.0 ) );
      check( std::nextafter( e, NC::kInfinity ) );
    }
    for ( std::size_t i = 0; i + 1 < E2.size(); ++i )
      for ( double u : { 0.25, 0.5, 0.75 } )
        check( E2[i] + u * ( E2[i+1] - E2[i] ) );
    //Special energies:
    for ( double e : { 0.0, std::numeric_limits<double>::denorm_min(),
                       std::numeric_limits<double>::min(), E2.front() * 1e-3,
                       E2.back() * 1e3, std::numeric_limits<double>::max(),
                       NC::kInfinity } )
      check( e );
    //The lookup table's bucket starts (as computed by LogGuideTable with 8
    //buckets per point) and their neighbours:
    {
      const std::uint64_t base = dblBits( E2.front() );
      const std::uint64_t range = dblBits( E2.back() ) - base;
      unsigned shift = 0;
      while ( ( range >> shift ) >= 8 * static_cast<std::uint64_t>( E2.size() ) )
        ++shift;
      const std::uint64_t nb = ( range >> shift ) + 1;
      const std::uint64_t step = std::max<std::uint64_t>( 1, nb / 5000 );
      for ( std::uint64_t b = 0; b < nb; b += step ) {
        const double v = fromBits( base + ( b << shift ) );
        check( v );
        check( std::nextafter( v, 0.0 ) );
        check( std::nextafter( v, NC::kInfinity ) );
      }
    }
    //Random energies, uniform in log(e) and uniform in the bit pattern:
    for ( int k = 0; k < 5000; ++k )
      check( E2.front() * 1e-2 * std::pow( 1e4 * E2.back() / E2.front(), rng.generate() ) );
    {
      const std::uint64_t lo = dblBits( E2.front() ), hi = dblBits( E2.back() );
      for ( int k = 0; k < 5000; ++k )
        check( fromBits( lo + static_cast<std::uint64_t>( rng.generate() * double( hi - lo + 1 ) ) ) );
    }
    std::cout << "OK: " << label << " (n=" << E2.size() << ( withq ? ", with q" : "" )
              << ( withr ? ", with r" : "" ) << ", coefficient style " << coefstyle << ")" << std::endl;
  }

  NC::VectD randomGrid( NC::RandXRSRImpl& rng )
  {
    //2..3000 points, each gap a random factor between 1 ulp and 1e6, below
    //1e100 (e*xs must stay finite, see NCIsoXSUtils.hh):
    const std::size_t n = 2 + static_cast<std::size_t>( rng.generate() * 2999 );
    NC::VectD x{ std::pow( 10.0, -12.0 + 14.0 * rng.generate() ) };
    while ( x.size() < n ) {
      const double r = rng.generate();
      double v;
      if ( r < 0.2 )
        v = std::nextafter( x.back(), NC::kInfinity );
      else
        v = x.back() * ( 1.0 + std::pow( 10.0, -15.0 + 21.0 * rng.generate() ) );
      if ( !( v > x.back() ) || !std::isfinite( v ) )
        v = std::nextafter( x.back(), NC::kInfinity );
      if ( !std::isfinite( v ) || v > 1e100 )
        break;
      x.push_back( v );
    }
    return x;
  }
  //Random tabulation on grid E (styles as in checkPieces):
  NC::IsoXSDescr makeDescr( NC::VectD E, bool withq, bool withr, int coefstyle,
                            NC::RandXRSRImpl& rng )
  {
    auto coef = [&rng,coefstyle]() {
      if ( coefstyle != 1 )
        return rng.generate();
      const double sign = ( rng.generate() < 0.5 ? -1.0 : 1.0 );
      return sign * std::pow( 10.0, -30.0 + 60.0 * rng.generate() );
    };
    NC::IsoXSDescr d;
    d.breakpts = { E.front(), E.back() };
    NC::IsoXSDescr::Tabulation p;
    p.E = std::move(E);
    for ( std::size_t i = 0; i + 1 < p.E.size(); ++i ) {
      for ( int k = 0; k < 4; ++k ) {
        p.p.push_back( coef() );
        if ( withq )
          p.q.push_back( coef() );
      }
      if ( withr )
        p.r.push_back( coef() * p.E[i] );
    }
    for ( auto& c : p.tail )
      c = coef();
    d.table = std::move(p);
    return d;
  }

  std::uint64_t g_nmergechecks = 0;
  double g_maxmergeulp = 0.0;

  //Merge the inputs and check the merged cross section against the scaled sum
  //of the inputs' own cross sections everywhere. Coefficients can have mixed
  //signs and cancel within a cubic, so deviations are measured relative to the
  //same sum evaluated with the absolute values of all coefficients (the
  //condition of the evaluation):
  void checkMerge( const char * label,
                   const std::vector<std::pair<double,const NC::IsoXSDescr*>>& in,
                   double tolulps, NC::RandXRSRImpl& rng )
  {
    NC::IsoXSDescr m = NC::mergeIsoXSDescr( in );
    REQUIRE( m.table.has_value() );
    REQUIRE( std::is_sorted( m.breakpts.begin(), m.breakpts.end() ) );
    for ( auto& e : in )
      for ( double b : e.second->breakpts )
        if ( e.first > 0.0 )
          REQUIRE( std::binary_search( m.breakpts.begin(), m.breakpts.end(), b ) );
    const NC::IsoXSEvaluator ev( std::move(m) );
    const auto& t = ev.descr().table.value();
    std::vector<NC::IsoXSDescr::Tabulation> absin;
    for ( auto& c : in ) {
      const auto& ti = c.second->table.value();
      NC::IsoXSDescr::Tabulation a;
      a.E = ti.E; a.p = ti.p; a.q = ti.q; a.r = ti.r;
      for ( auto* v : { &a.p, &a.q, &a.r } )
        for ( auto& x : *v )
          x = std::fabs( x );
      for ( int k = 0; k < 5; ++k )
        a.tail[k] = std::fabs( ti.tail[k] );
      absin.push_back( std::move(a) );
    }
    double maxulp = 0.0;
    auto check = [&]( double e ) {
      if ( !( e > 0.0 ) || !std::isfinite( e ) )
        return;
      double sum = 0.0, abssum = 0.0;
      for ( std::size_t ic = 0; ic < in.size(); ++ic ) {
        const auto& c = in[ic];
        if ( c.first == 0.0 )
          continue;
        sum += c.first * NC::evalIsoXS( c.second->table.value(), e );
        abssum += c.first * NC::evalIsoXS( absin[ic], e );
      }
      const double v = ev.eval( e );
      REQUIRE( v == NC::evalIsoXS( t, e ) );
      if ( !std::isfinite( abssum ) || abssum == 0.0 ) {
        if ( abssum == 0.0 )
          REQUIRE( v == 0.0 );
        return;
      }
      //(e*xs) itself is what the table holds, so also allow for cancellation
      //between the sqrt and 1/e terms by using the absolute evaluation above.
      const double ulps = std::fabs( v - sum ) / ( abssum * std::numeric_limits<double>::epsilon() );
      maxulp = std::max( maxulp, ulps );
      if ( !( ulps <= tolulps ) ) {
        std::cout << "MISMATCH " << label << " e=" << e << " merged=" << v
                  << " sum=" << sum << " (" << ulps << " ulp of the absolute evaluation)" << std::endl;
        REQUIRE( false );
      }
      ++g_nmergechecks;
    };
    const NC::VectD& E = t.E;
    for ( auto e : E ) {
      check( e );
      check( std::nextafter( e, 0.0 ) );
      check( std::nextafter( e, NC::kInfinity ) );
    }
    for ( std::size_t i = 0; i + 1 < E.size(); ++i )
      for ( double u : { 0.1, 0.5, 0.9 } )
        check( E[i] + u * ( E[i+1] - E[i] ) );
    for ( double e : { std::numeric_limits<double>::denorm_min(), E.front() * 1e-3, E.back() * 1e3 } )
      check( e );
    for ( int k = 0; k < 3000; ++k )
      check( E.front() * 1e-2 * std::pow( 1e4 * E.back() / E.front(), rng.generate() ) );
    g_maxmergeulp = std::max( g_maxmergeulp, maxulp );
    (void)label;
  }

  void testMerge( NC::RandXRSRImpl& rng )
  {
    //Structural cases:
    {
      const std::vector<std::pair<double,const NC::IsoXSDescr*>> none;
      NC::IsoXSDescr m = NC::mergeIsoXSDescr( none );
      REQUIRE( m.breakpts.empty() && !m.table.has_value() );
    }
    {
      //A single input with scale 1 merges to itself, bit for bit:
      for ( int mode : { 0, 1, 2, 3 } ) {
        NC::IsoXSDescr d = makeDescr( randomGrid( rng ), mode & 1, mode & 2, 1, rng );
        const std::vector<std::pair<double,const NC::IsoXSDescr*>> in{ { 1.0, &d } };
        NC::IsoXSDescr m = NC::mergeIsoXSDescr( in );
        const auto& a = d.table.value();
        const auto& b = m.table.value();
        REQUIRE( a.E == b.E && a.p == b.p && a.q == b.q && a.r == b.r );
        for ( int k = 0; k < 5; ++k )
          REQUIRE( a.tail[k] == b.tail[k] );
        REQUIRE( m.breakpts == d.breakpts );
      }
    }
    {
      //No table if any (non-zero scale) input lacks one; scale 0 inputs ignored:
      NC::IsoXSDescr a = makeDescr( { 1.0, 2.0 }, false, false, 0, rng );
      NC::IsoXSDescr b;
      b.breakpts = { 0.5, 3.0 };
      std::vector<std::pair<double,const NC::IsoXSDescr*>> in{ { 1.0, &a }, { 2.0, &b } };
      NC::IsoXSDescr m = NC::mergeIsoXSDescr( in );
      REQUIRE( !m.table.has_value() );
      REQUIRE( ( m.breakpts == NC::VectD{ 0.5, 1.0, 2.0, 3.0 } ) );
      in[1].first = 0.0;
      m = NC::mergeIsoXSDescr( in );
      REQUIRE( m.table.has_value() && ( m.breakpts == NC::VectD{ 1.0, 2.0 } ) );
    }
    {
      //Zero coefficients merge to zero (with empty p, q and r):
      NC::IsoXSDescr a = makeDescr( { 1.0, 2.0 }, false, false, 0, rng );
      NC::IsoXSDescr b = makeDescr( { 1.5, 3.0 }, false, false, 0, rng );
      for ( auto* d : { &a, &b } ) {
        auto& t = d->table.value();
        for ( auto& c : t.p ) c = 0.0;
        for ( auto& c : t.tail ) c = 0.0;
      }
      const std::vector<std::pair<double,const NC::IsoXSDescr*>> in{ { 1.0, &a }, { 3.0, &b } };
      NC::IsoXSDescr m = NC::mergeIsoXSDescr( in );
      REQUIRE( m.table.value().p.empty() && m.table.value().q.empty() && m.table.value().r.empty() );
      for ( double e : { 0.0, 0.5, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 10.0 } )
        REQUIRE( NC::evalIsoXS( m.table.value(), e ) == 0.0 );
    }
    {
      //Inputs with only a Q part (empty p, tails giving no P terms) merge to a
      //table with empty p, agreeing with the scaled sum:
      NC::RandXRSRImpl rng2( 4242 );
      NC::IsoXSDescr a = makeDescr( randomGrid( rng2 ), true, false, 0, rng2 );
      NC::IsoXSDescr b = makeDescr( randomGrid( rng2 ), true, false, 0, rng2 );
      for ( auto* d : { &a, &b } ) {
        auto& t = d->table.value();
        t.p.clear();
        for ( int k : { 1, 2, 3, 4 } )
          t.tail[k] = 0.0;
      }
      const std::vector<std::pair<double,const NC::IsoXSDescr*>> in{ { 0.7, &a }, { 1.3, &b } };
      NC::IsoXSDescr m = NC::mergeIsoXSDescr( in );
      REQUIRE( m.table.value().p.empty() );
      for ( int k = 0; k < 2000; ++k ) {
        const double e = 1e-12 * std::pow( 1e14, rng2.generate() );
        const double ref = 0.7 * NC::evalIsoXS( a.table.value(), e ) + 1.3 * NC::evalIsoXS( b.table.value(), e );
        REQUIRE( NC::floateq( NC::evalIsoXS( m.table.value(), e ), ref, 1e-12, 0.0 ) );
      }
    }
    //Numerical cases: 2..8 inputs on random grids (overlapping, disjoint,
    //sharing breakpoints, or offset by 1 ulp), random scales:
    for ( int trial = 0; trial < 60; ++trial ) {
      const int coefstyle = trial % 2 ? 1 : 0;
      const std::size_t K = 2 + static_cast<std::size_t>( rng.generate() * 7 );
      std::vector<NC::IsoXSDescr> ds;
      NC::VectD shared = randomGrid( rng );
      for ( std::size_t k = 0; k < K; ++k ) {
        NC::VectD E;
        const double r = rng.generate();
        if ( r < 0.25 ) {
          E = shared;//identical breakpoints
        } else if ( r < 0.5 ) {
          E = shared;//breakpoints offset by 1 ulp
          for ( auto& e : E )
            e = std::nextafter( e, NC::kInfinity );
          for ( std::size_t i = 1; i < E.size(); ++i )
            if ( !( E[i] > E[i-1] ) )
              E[i] = std::nextafter( E[i-1], NC::kInfinity );
        } else {
          E = randomGrid( rng );
        }
        const int mode = static_cast<int>( rng.generate() * 4 );
        ds.push_back( makeDescr( E, mode & 1, mode & 2, coefstyle, rng ) );
      }
      std::vector<std::pair<double,const NC::IsoXSDescr*>> in;
      for ( auto& d : ds ) {
        const double r = rng.generate();
        const double scale = ( r < 0.1 ? 0.0 : ( r < 0.3 ? 1.0 : std::pow( 10.0, -3.0 + 6.0 * rng.generate() ) ) );
        in.emplace_back( scale, &d );
      }
      if ( in.front().first == 0.0 )
        in.front().first = 1.0;
      checkMerge( "random inputs", in, 64.0, rng );
    }
  }

}

namespace {
  std::uint64_t g_nscalechecks = 0;

  //scaleIsoXSDescr against IsoXSEvaluator::evalScaled, on random tables
  //(positive coefficients) and scales. The random tables jump at their nodes,
  //evalIsoXSDerivative agrees with central differences of evalIsoXS, at
  //random energies inside the pieces of random tables (away from the nodes,
  //so the differences stay within a piece), and in the tails:
  std::uint64_t g_nderivchecks = 0;
  void testDerivative( NC::RandXRSRImpl& rng )
  {
    for ( int k = 0; k < 40; ++k ) {
      const int mode = k % 4;
      NC::IsoXSDescr d = makeDescr( randomGrid( rng ), mode & 1, mode & 2, 0, rng );
      auto& t = d.table.value();
      const double tails[5] = { 0.3, 0.2, 0.5, 0.7, 0.11 };
      std::copy( tails, tails + 5, t.tail );
      const NC::VectD& E = t.E;
      auto check = [&t]( double e, double h )
      {
        //(dividing by the difference of the rounded energies):
        const double e1 = e - h;
        const double e2 = e + h;
        const double num = ( ( NC::evalIsoXS( t, e2 ) - NC::evalIsoXS( t, e1 ) )
                             / ( e2 - e1 ) );
        const double ana = NC::evalIsoXSDerivative( t, e );
        //Scale of the derivative, for the comparison:
        const double scale = std::max( NC::ncabs( ana ),
                                       NC::ncabs( NC::evalIsoXS( t, e ) ) / e );
        REQUIRE( NC::ncabs( ana - num ) <= 1e-6 * scale );
        ++g_nderivchecks;
      };
      for ( std::size_t i = 0; i + 1 < E.size(); ++i ) {
        const double w = E[i+1] - E[i];
        if ( !( w > 1e-8 * E[i+1] ) )
          continue;//interval too narrow for central differences
        for ( int j = 0; j < 3; ++j )
          check( E[i] + ( 0.1 + 0.8 * rng.generate() ) * w, 1e-5 * w );
      }
      check( 0.5 * E.front(), 1e-5 * E.front() );
      check( 2.0 * E.back(), 1e-5 * E.back() );
    }
  }

  //where E/escale and escale*e might round to opposite sides, so the
  //neighbours of the nodes are not checked:
  void testScale( NC::RandXRSRImpl& rng )
  {
    for ( int trial = 0; trial < 40; ++trial ) {
      const std::size_t n = 2 + static_cast<std::size_t>( rng.generate() * 200 );
      NC::VectD E = NC::geomspace( 1e-6 * std::pow( 10.0, rng.generate() ),
                                   10.0 * std::pow( 10.0, rng.generate() ), n );
      const int mode = trial % 4;
      NC::IsoXSDescr d = makeDescr( std::move(E), mode & 1, mode & 2, 0, rng );
      if ( trial % 5 == 0 )
        d.table.value().p.clear();
      const double escale = ( trial % 7 == 0 ? 1.0 : std::pow( 10.0, -3.0 + 8.0 * rng.generate() ) );
      const double xsscale = ( trial % 9 == 0 ? 0.0 : std::pow( 10.0, -3.0 + 6.0 * rng.generate() ) );
      const NC::IsoXSDescr sd = NC::scaleIsoXSDescr( d, escale, xsscale );
      const NC::IsoXSEvaluator ev( std::move(d) );
      const auto& st = sd.table.value();
      REQUIRE( sd.breakpts.size() == ev.descr().breakpts.size() );
      auto check = [&]( double e ) {
        if ( !( e >= 0.0 ) )
          return;
        const double a = ev.evalScaled( e, escale, xsscale );
        const double b = NC::evalIsoXS( st, e );
        REQUIRE( std::isfinite( a ) && std::isfinite( b ) );
        //Allowed: 8 ulp, plus the change of the cross section when moving the
        //energy by up to 4 ulp (the energy transformation rounds):
        double elo = e, ehi = e;
        for ( int k = 0; k < 4; ++k ) {
          elo = std::nextafter( elo, 0.0 );
          ehi = std::nextafter( ehi, NC::kInfinity );
        }
        const double slope = std::fabs( ev.evalScaled( ehi, escale, xsscale )
                                        - ev.evalScaled( elo, escale, xsscale ) );
        const double tol = 8.0 * std::numeric_limits<double>::epsilon() * std::fabs( a ) + slope;
        const double dev = std::fabs( a - b );
        REQUIRE( dev <= tol );
        ++g_nscalechecks;
      };
      for ( std::size_t i = 0; i + 1 < st.E.size(); ++i )
        for ( double u : { 0.001, 0.25, 0.5, 0.75, 0.999 } )
          check( st.E[i] + u * ( st.E[i+1] - st.E[i] ) );
      for ( double e : { 0.0, std::numeric_limits<double>::denorm_min(), 1e-320, 1e-310,
                         std::numeric_limits<double>::min(), st.E.front() * 1e-3, st.E.back() * 1e3 } )
        check( e );
      for ( int k = 0; k < 2000; ++k )
        check( st.E.front() * 1e-2 * std::pow( 1e4 * st.E.back() / st.E.front(), rng.generate() ) );
    }
  }
}

int main()
{
  NC::RandXRSRImpl rng( 54321 );
  std::vector<std::pair<std::string,NC::VectD>> grids;
  grids.emplace_back( "two points", NC::VectD{ 1e-3, 1.0 } );
  grids.emplace_back( "two points, ratio 1e400", NC::VectD{ 1e-200, 1e200 } );
  grids.emplace_back( "geometric 1e-9..10", NC::geomspace( 1e-9, 10.0, 300 ) );
  grids.emplace_back( "geometric 1e-300..1e300", NC::geomspace( 1e-300, 1e300, 50 ) );
  {
    NC::VectD x;
    for ( int k = -40; k <= 40; ++k )
      x.push_back( std::ldexp( 1.0, k ) );
    grids.emplace_back( "powers of two", x );
  }
  {
    NC::VectD x;
    for ( auto v : NC::geomspace( 1e-3, 1e3, 40 ) ) {
      x.push_back( v );
      x.push_back( std::nextafter( v * 1.01, NC::kInfinity ) );
      x.push_back( v * 1.01 );
    }
    std::sort( x.begin(), x.end() );
    grids.emplace_back( "with 1-ulp intervals", x );
  }
  {
    NC::VectD x{ 1e-3 };
    double v = 1.0;
    for ( int k = 0; k < 2000; ++k ) {
      x.push_back( v );
      v = std::nextafter( v, NC::kInfinity );
    }
    x.push_back( 1e3 );
    grids.emplace_back( "2000 consecutive doubles", x );
  }
  grids.emplace_back( "subnormal start", NC::VectD{ std::numeric_limits<double>::denorm_min(),
                                                    1e-320, 1e-310, 1e-300, 1e-10, 1.0 } );
  for ( int k = 0; k < 10; ++k )
    grids.emplace_back( "random grid", randomGrid( rng ) );
  for ( auto& g : grids )
    for ( int mode : { 0, 1, 2, 3 } )
      for ( int coefstyle : { 0, 1, 2 } )
        checkPieces( g.first.c_str(), g.second, ( mode & 1 ), ( mode & 2 ), coefstyle, rng );
  {
    //All coefficients zero (e.g. like Bragg below its first edge): zero.
    NC::IsoXSDescr d;
    d.table.emplace();
    d.table.value().E = { 1.0, 2.0 };
    d.table.value().p = { 0.0, 0.0, 0.0, 0.0 };
    const NC::IsoXSEvaluator ev( std::move(d) );
    const auto& t = ev.descr().table.value();
    for ( double e : { 0.0, 0.5, 1.5, 5.0 } )
      REQUIRE( NC::evalIsoXS( t, e ) == 0.0 && ev.eval( e ) == 0.0 );
  }
  {
    //Empty p means zero (here with only a Q part, like free-gas tables):
    NC::IsoXSDescr::Tabulation a, b;
    a.E = b.E = { 0.5, 1.0, 3.0 };
    a.q = b.q = { 1.0, -0.5, 0.25, 0.1, 2.0, 0.3, -0.2, 0.05 };
    b.p.assign( 8, 0.0 );
    a.tail[0] = b.tail[0] = 1.5;
    a.tail[3] = b.tail[3] = 0.7;
    NC::IsoXSDescr da;
    da.table = std::move( a );
    const NC::IsoXSEvaluator ev( std::move( da ) );
    const auto& ta = ev.descr().table.value();
    for ( double e : { 0.0, 0.1, 0.5, 0.75, 1.0, 2.0, 2.999, 3.0, 10.0 } ) {
      const double ref = refXS( b, e );
      REQUIRE( NC::evalIsoXS( ta, e ) == ref && NC::evalIsoXS( b, e ) == ref && ev.eval( e ) == ref );
    }
  }
  {
    //At e=0, tail[0]/sqrt(e) uses the smallest positive double (finite):
    NC::IsoXSDescr d;
    d.table.emplace();
    d.table.value().E = { 1.0, 2.0 };
    d.table.value().p = { 1.0, 1.0, 0.0, 0.0 };
    d.table.value().tail[0] = 3.0;
    d.table.value().tail[1] = 0.5;
    const NC::IsoXSEvaluator ev( std::move(d) );
    const auto& t = ev.descr().table.value();
    const double ref = 3.0 / std::sqrt( std::numeric_limits<double>::denorm_min() ) + 0.5;
    REQUIRE( std::isfinite( ref ) );
    REQUIRE( NC::evalIsoXS( t, 0.0 ) == ref && ev.eval( 0.0 ) == ref );
  }
  std::cout << "evalIsoXS and IsoXSEvaluator agree with the reference in " << g_nchecks << " checks" << std::endl;
  testMerge( rng );
  std::cout << "mergeIsoXSDescr agrees with the sum of its inputs in " << g_nmergechecks
            << " checks (max deviation " << ( g_maxmergeulp < 8 ? "below 8" : "8 or more" )
            << " ulp of the absolute evaluation)" << std::endl;
  {
    NC::RandXRSRImpl rng3( 777 );
    testScale( rng3 );
  }
  std::cout << "scaleIsoXSDescr agrees with IsoXSEvaluator::evalScaled in " << g_nscalechecks
            << " checks (within rounding of the energy)" << std::endl;
  {
    NC::RandXRSRImpl rng4( 4242 );
    testDerivative( rng4 );
  }
  std::cout << "evalIsoXSDerivative agrees with differences of evalIsoXS in "
            << g_nderivchecks << " checks" << std::endl;
  return 0;
}
