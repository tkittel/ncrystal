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

  void checkPieces( const char * label, NC::VectD E, bool withq, bool withr,
                    NC::RandXRSRImpl& rng )
  {
    NC::IsoXSDescr d;
    d.breakpts = { E.front(), E.back() };
    NC::IsoXSDescr::Tabulation p;
    p.E = std::move(E);
    for ( std::size_t i = 0; i + 1 < p.E.size(); ++i ) {
      for ( int k = 0; k < 4; ++k ) {
        p.p.push_back( rng.generate() );
        if ( withq )
          p.q.push_back( rng.generate() );
      }
      if ( withr )
        p.r.push_back( rng.generate() * p.E[i] );
    }
    for ( auto& c : p.tail )
      c = rng.generate();
    d.table = std::move(p);
    const auto& pp = d.table.value();
    auto check = [&]( double e ) {
      if ( !( e >= 0.0 ) )
        return;
      const double ref = refXS( pp, e );
      const double v = NC::evalIsoXS( pp, e );
      REQUIRE( v == ref || ( std::isinf(v) && std::isinf(ref) ) );
      ++g_nchecks;
    };
    for ( auto e : pp.E ) {
      check( e );
      check( std::nextafter( e, 0.0 ) );
      check( std::nextafter( e, NC::kInfinity ) );
    }
    for ( double e : { 0.0, pp.E.front() * 1e-3, pp.E.back() * 1e3 } )
      check( e );
    for ( int k = 0; k < 20000; ++k )
      check( pp.E.front() * 1e-2 * std::pow( 1e4 * pp.E.back() / pp.E.front(), rng.generate() ) );
    std::cout << "OK: " << label << " (n=" << pp.E.size() << ( withq ? ", with q" : "" )
              << ( withr ? ", with r" : "" ) << ")" << std::endl;
  }

}

int main()
{
  NC::RandXRSRImpl rng( 54321 );
  for ( int mode : { 0, 1, 2, 3 } ) {
    const bool withq = ( mode & 1 ), withr = ( mode & 2 );
    checkPieces( "two points", { 1e-3, 1.0 }, withq, withr, rng );
    checkPieces( "geometric 1e-9..10", NC::geomspace( 1e-9, 10.0, 300 ), withq, withr, rng );
    NC::VectD x;
    for ( auto v : NC::geomspace( 1e-3, 1e3, 40 ) ) {
      x.push_back( v );
      x.push_back( std::nextafter( v * 1.01, NC::kInfinity ) );
      x.push_back( v * 1.01 );
    }
    std::sort( x.begin(), x.end() );
    checkPieces( "with 1-ulp intervals", x, withq, withr, rng );
  }
  {
    //All coefficients zero (e.g. like Bragg below its first edge): zero.
    NC::IsoXSDescr::Tabulation t;
    t.E = { 1.0, 2.0 };
    t.p = { 0.0, 0.0, 0.0, 0.0 };
    for ( double e : { 0.0, 0.5, 1.5, 5.0 } )
      REQUIRE( NC::evalIsoXS( t, e ) == 0.0 );
  }
  {
    //Empty p means zero (here with only a Q part, like free-gas tables):
    NC::IsoXSDescr::Tabulation a, b;
    a.E = b.E = { 0.5, 1.0, 3.0 };
    a.q = b.q = { 1.0, -0.5, 0.25, 0.1, 2.0, 0.3, -0.2, 0.05 };
    b.p.assign( 8, 0.0 );
    a.tail[0] = b.tail[0] = 1.5;
    a.tail[3] = b.tail[3] = 0.7;
    for ( double e : { 0.0, 0.1, 0.5, 0.75, 1.0, 2.0, 2.999, 3.0, 10.0 } )
      REQUIRE( NC::evalIsoXS( a, e ) == NC::evalIsoXS( b, e ) && NC::evalIsoXS( a, e ) == refXS( b, e ) );
  }
  std::cout << "evalIsoXS agrees with the reference in " << g_nchecks << " checks" << std::endl;
  return 0;
}
