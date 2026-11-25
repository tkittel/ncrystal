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
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCrystal/internal/vdos/NCVDOSDebyeXS.hh"
#include <iostream>
#include <functional>

//Tests of thinIsoXSTable.

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  using Fct = std::function<double(double)>;

  //Table with nodes at E (q-only cubic Hermite of g=sqrt(e)*xs, from xs and
  //its derivative at the nodes, evaluated just inside each interval so kinks
  //and jumps at nodes are represented), and the given tails:
  NC::IsoXSTable tabulate( const NC::VectD& E, Fct xs, Fct dxs,
                           double t0, double t1, double t2, double t3 )
  {
    NC::IsoXSTable t;
    t.E = E;
    for ( auto i : NC::ncrange( E.size() - 1 ) ) {
      const double a = E[i], b = E[i+1], w = b - a;
      const double ai = std::nextafter( a, b ), bi = std::nextafter( b, a );
      auto dg = [&xs,&dxs]( double e, double ei )
      {
        return 0.5 * xs( ei ) / std::sqrt( e ) + std::sqrt( e ) * dxs( ei );
      };
      const double ga = std::sqrt( a ) * xs( ai ), gb = std::sqrt( b ) * xs( bi );
      const double d0 = w * dg( a, ai ), d1 = w * dg( b, bi );
      t.q.insert( t.q.end(), { ga, d0, 3.0 * ( gb - ga ) - 2.0 * d0 - d1,
                               2.0 * ( ga - gb ) + d0 + d1 } );
    }
    t.tail[0] = t0;
    t.tail[1] = t1;
    t.tail[2] = t2;
    t.tail[3] = t3;
    return t;
  }

  bool isSubset( const NC::VectD& sub, const NC::VectD& full )
  {
    for ( auto e : sub )
      if ( !std::binary_search( full.begin(), full.end(), e ) )
        return false;
    return true;
  }

  //Max relative deviation of b from a, sampled densely within each interval
  //of a (and just around its nodes, from both sides):
  double maxRelDev( const NC::IsoXSTable& a, const NC::IsoXSTable& b )
  {
    double res = 0.0;
    auto upd = [&res,&a,&b]( double e )
    {
      const double xa = NC::evalIsoXS( a, e ), xb = NC::evalIsoXS( b, e );
      if ( xa == 0.0 )
        res = NC::ncmax( res, xb == 0.0 ? 0.0 : NC::kInfinity );
      else
        res = NC::ncmax( res, NC::ncabs( xb / xa - 1.0 ) );
    };
    for ( auto i : NC::ncrange( a.E.size() - 1 ) ) {
      const double lo = a.E[i], w = a.E[i+1] - lo;
      for ( auto k : NC::ncrange( 1, 12 ) )
        upd( lo + w * k / 12.0 );
      upd( lo );
      upd( std::nextafter( a.E[i+1], 0.0 ) );
    }
    upd( a.E.front() * 0.5 );
    upd( a.E.back() * 2.0 );
    return res;
  }

  void checkThinned( const NC::IsoXSTable& t, double tol, const char * descr )
  {
    auto th = NC::thinIsoXSTable( t, tol );
    REQUIRE( th.E.front() == t.E.front() && th.E.back() == t.E.back() );
    REQUIRE( isSubset( th.E, t.E ) );
    for ( auto k : NC::ncrange( 5 ) )
      REQUIRE( th.tail[k] == t.tail[k] );
    //Checked at nodes and midpoints, so allow a bit more in between:
    const double dev = maxRelDev( t, th );
    REQUIRE( dev <= 2.0 * tol + 1e-13 );
    std::cout << descr << " tol=" << tol << ": OK" << std::endl;
  }

  void testSmooth()
  {
    //xs = 2 + 0.3/sqrt(e) + 0.5*exp(-e), heavily oversampled:
    Fct f = []( double e )
    {
      return 2.0 + 0.3 / std::sqrt( e ) + 0.5 * std::exp( -e );
    };
    Fct df = []( double e )
    {
      return -0.15 / ( e * std::sqrt( e ) ) - 0.5 * std::exp( -e );
    };
    auto t = tabulate( NC::geomspace( 1e-5, 10.0, 4000 ), f, df,
                       0.3, 2.5, 0.1, 2.0 );
    std::size_t nprev = t.E.size() + 1;
    for ( double tol : { 0.0, 1e-10, 1e-7, 5e-5, 1e-3, 1e-1 } ) {
      checkThinned( t, tol, "smooth" );
      if ( tol == 0.0 ) {
        //Identical (no node can be dropped without any deviation):
        const auto th = NC::thinIsoXSTable( t, tol );
        REQUIRE( th.E == t.E && th.p == t.p && th.q == t.q && th.r == t.r );
      }
      const auto n = NC::thinIsoXSTable( t, tol ).E.size();
      REQUIRE( n <= nprev );
      nprev = n;
    }
    //Strong reduction for a smooth oversampled table:
    REQUIRE( NC::thinIsoXSTable( t, 5e-5 ).E.size() * 10 < t.E.size() );
  }

  void testKinkAndJump()
  {
    //Kink at e=1 (a node):
    Fct f = []( double e ) { return 1.0 + NC::ncabs( e - 1.0 ); };
    Fct df = []( double e ) { return e < 1.0 ? -1.0 : 1.0; };
    auto E = NC::linspace( 0.5, 1.5, 201 );
    E[100] = 1.0;
    auto t = tabulate( E, f, df, 0.0, 1.5, 0.0, 1.5 );
    for ( double tol : { 1e-6, 1e-3 } ) {
      checkThinned( t, tol, "kink" );
      const auto th = NC::thinIsoXSTable( t, tol );
      REQUIRE( std::binary_search( th.E.begin(), th.E.end(), 1.0 ) );
    }
    //Jump at e=2 (a node): xs=1 below, 3 above:
    auto E2 = NC::linspace( 1.0, 3.0, 101 );
    REQUIRE( E2[50] == 2.0 );
    auto tj = tabulate( E2, []( double e ) { return e < 2.0 ? 1.0 : 3.0; },
                        []( double ) { return 0.0; }, 0.0, 1.0, 0.0, 3.0 );
    for ( double tol : { 1e-6, 1e-2 } ) {
      checkThinned( tj, tol, "jump" );
      const auto th = NC::thinIsoXSTable( tj, tol );
      REQUIRE( std::binary_search( th.E.begin(), th.E.end(), 2.0 ) );
      //(the Hermite pieces are in sqrt(e)*xs, so a constant xs needs a few
      //nodes for tight tolerances):
      if ( tol >= 1e-2 )
        REQUIRE( th.E.size() <= 4 );
    }
  }

  void testMixedAndDebye()
  {
    //A sum of Debye model tables (q parts, and p parts from the tails of
    //tables ending inside the others), as in the grouped VDOSDebye
    //treatment:
    std::vector<NC::IsoXSTable> tbls;
    for ( double td : { 150.0, 300.0, 310.0, 600.0 } )
      tbls.push_back( NC::VDOS::VDOSDebyeXS( NC::DebyeTemperature{ td },
                                             NC::Temperature{ 293.15 },
                                             NC::AtomMass{ 15.999 } )
                      .createXSTabulation() );
    std::vector<std::pair<double,const NC::IsoXSTable*>> parts;
    for ( auto& t : tbls )
      parts.emplace_back( 0.25, &t );
    auto sum = NC::mergeIsoXSTables( parts );
    REQUIRE( !sum.p.empty() && !sum.q.empty() );
    for ( double tol : { 0.0, 5e-5, 1e-3 } )
      checkThinned( sum, tol, "Debye sum" );
    REQUIRE( NC::thinIsoXSTable( sum, 5e-5 ).E.size() * 2 < sum.E.size() );
  }

  void testEdgeCases()
  {
    //Zero cross section in part of the range:
    Fct f = []( double e ) { return e < 1.0 ? 0.0 : ( e - 1.0 ) * ( e - 1.0 ); };
    Fct df = []( double e ) { return e < 1.0 ? 0.0 : 2.0 * ( e - 1.0 ); };
    auto tz = tabulate( NC::linspace( 0.5, 2.0, 61 ), f, df, 0.0, 0.0, 0.0, 1.0 );
    for ( double tol : { 0.0, 1e-4 } )
      checkThinned( tz, tol, "zero range" );
    //Two nodes:
    auto t2 = tabulate( { 1.0, 2.0 }, []( double e ) { return 1.0 / e; },
                        []( double e ) { return -1.0 / ( e * e ); },
                        0.0, 1.0, 0.5, 0.0 );
    checkThinned( t2, 0.1, "two nodes" );
    REQUIRE( NC::thinIsoXSTable( t2, 0.1 ).E.size() == 2 );
    //Invalid tolerances:
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for ( double tol : { -1e-3, NC::kInfinity, nan } ) {
      bool threw = false;
      try {
        NC::thinIsoXSTable( t2, tol );
      } catch ( NC::Error::BadInput& ) {
        threw = true;
      }
      REQUIRE( threw );
    }
    std::cout << "invalid tolerances: OK" << std::endl;
  }
}

int main()
{
  testSmooth();
  testKinkAndJump();
  testMixedAndDebye();
  testEdgeCases();
  return 0;
}
