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
    if ( e < p.E.front() )
      return ( p.tail[0] ? p.tail[0] / std::sqrt( e ) : 0.0 ) + p.tail[1];
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
    //2..3000 points, each gap a random factor between 1 ulp and 1e6:
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
      if ( !std::isfinite( v ) )
        break;
      x.push_back( v );
    }
    return x;
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
  std::cout << "evalIsoXS and IsoXSEvaluator agree with the reference in " << g_nchecks << " checks" << std::endl;
  return 0;
}
