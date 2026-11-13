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

#include "NCrystal/internal/phys_utils/NCUijUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

//Tests of the UijUtils functions for anisotropic displacement tensors.

namespace NC = NCrystal;
namespace U = NC::UijUtils;
#define REQUIRE(x) nc_assert_always(x)

namespace {

  //(values below 1e-15 printed as 0, for platform-independent output):
  std::string str( const U::Uij& u )
  {
    std::ostringstream ss;
    ss << '(';
    for ( auto i : NC::ncrange( 6 ) )
      ss << ( i ? ", " : "" )
         << NC::fmt( NC::ncabs( u[i] ) < 1e-15 ? 0.0 : u[i], "%.10g" );
    ss << ')';
    return ss.str();
  }

  bool uijeq( const U::Uij& a, const U::Uij& b, double tol = 1e-12 )
  {
    for ( auto i : NC::ncrange( 6 ) )
      if ( !NC::floateq( a[i], b[i], tol, tol ) )
        return false;
    return true;
  }

  const NC::SGSymmetry& sgsym( const char * s )
  {
    return NC::SGSymmetry::get( NC::SpaceGroup( s,
                                                NC::SpaceGroupParseMode::Strict ) );
  }

  //Run f, which must throw BadInput, and print the message:
  template<class TFct>
  void expectBadInput( TFct f )
  {
    try {
      f();
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  BadInput: " << e.what() << std::endl;
      return;
    }
    NCRYSTAL_THROW(LogicError,"Expected BadInput exception not thrown");
  }

  void testConversions()
  {
    std::cout << "Conversions:" << std::endl;
    const U::Uij u{ 0.012, 0.009, 0.015, 0.002, 0.003, 0.001 };
    //In a cubic cell, Ustar = Ustd/a^2 and Ueq = trace/3:
    const auto mc = U::cellMetric( 4.0, 4.0, 4.0, 90.0, 90.0, 90.0 );
    const auto us = U::convertUstdToUstar( u, mc );
    for ( auto i : NC::ncrange( 6 ) )
      REQUIRE( NC::floateq( us[i], u[i] / 16.0, 1e-14, 0.0 ) );
    REQUIRE( NC::floateq( U::calcUeqFromUstar( us, mc ), 0.012, 1e-14, 0.0 ) );
    //Round trips in a triclinic cell:
    const auto mt = U::cellMetric( 5.0, 6.0, 7.0, 80.0, 95.0, 100.0 );
    const auto ust = U::convertUstdToUstar( u, mt );
    REQUIRE( uijeq( U::convertUstarToUstd( ust, mt ), u ) );
    REQUIRE( uijeq( U::convertUstdToUstar( U::convertUstarToUstd( ust, mt ),
                                           mt ), ust ) );
    std::cout << "  Ueq (triclinic): "
              << NC::fmt( U::calcUeqFromUstar( ust, mt ), "%.10g" ) << std::endl;
    //Ueq is the mean of the eigenvalues of Ucart:
    const auto ev = U::calcUcartEigenvalues( ust, mt );
    REQUIRE( ev[0] <= ev[1] && ev[1] <= ev[2] && ev[0] > 0.0 );
    REQUIRE( NC::floateq( ( ev[0] + ev[1] + ev[2] ) / 3.0,
                          U::calcUeqFromUstar( ust, mt ), 1e-12, 0.0 ) );
    //Known eigenvalues in a cubic cell (where Ucart=Ustd):
    const auto ev2 = U::calcUcartEigenvalues(
      U::convertUstdToUstar( { 1.0, 1.0, 1.0, 0.5, 0.0, 0.0 }, mc ), mc );
    std::cout << "  eigenvalues: " << NC::fmt( ev2[0], "%.10g" ) << ", "
              << NC::fmt( ev2[1], "%.10g" ) << ", "
              << NC::fmt( ev2[2], "%.10g" ) << std::endl;
    REQUIRE( NC::floateq( ev2[0], 0.5, 1e-12, 0.0 )
             && NC::floateq( ev2[1], 1.0, 1e-12, 0.0 )
             && NC::floateq( ev2[2], 1.5, 1e-12, 0.0 ) );
  }

  void testBasics()
  {
    std::cout << "Basics:" << std::endl;
    const U::Uij u{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    //y,x,-z swaps U11/U22 and maps (U13,U23) to (-U23,-U13):
    const auto r = U::applySymOpToUstar( NC::SymOp( std::string( "y,x,-z" ) ),
                                         u );
    std::cout << "  y,x,-z: " << str( u ) << " -> " << str( r ) << std::endl;
    REQUIRE( r == U::Uij( { 2.0, 1.0, 3.0, 4.0, -6.0, -5.0 } ) );
    REQUIRE( U::applySymOpToUstar( NC::SymOp(), u ) == u );
    REQUIRE( U::applySymOpToUstar( NC::SymOp( std::string( "-x,-y,-z" ) ),
                                   u ) == u );
    REQUIRE( U::maxAbs( { 0.1, -0.3, 0.2, 0.0, 0.0, 0.0 } ) == 0.3 );
    REQUIRE( U::maxAbs( U::Uij{} ) == 0.0 );
    REQUIRE( U::isPositiveDefinite( { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 } ) );
    REQUIRE( U::isPositiveDefinite( { 1.0, 1.0, 1.0, 0.9, 0.0, 0.0 } ) );
    REQUIRE( !U::isPositiveDefinite( { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0 } ) );
    REQUIRE( U::isPositiveDefinite( { 1.0, 1.0, 1.0, 0.7, 0.7, 0.7 } ) );
    REQUIRE( !U::isPositiveDefinite( { 1.0, 1.0, 1.0, -0.7, -0.7, -0.7 } ) );
    REQUIRE( !U::isPositiveDefinite( { 1.0, -1.0, 1.0, 0.0, 0.0, 0.0 } ) );
    REQUIRE( !U::isPositiveDefinite( U::Uij{} ) );
  }

  void testDirectionalMSDs()
  {
    std::cout << "Directional MSDs:" << std::endl;
    const auto mt = U::cellMetric( 5.0, 6.0, 7.0, 80.0, 95.0, 100.0 );
    const auto us = U::convertUstdToUstar( { 0.012, 0.009, 0.015,
                                             0.002, 0.003, 0.001 }, mt );
    const auto d = U::directionalMSDs( us, mt );
    std::cout << "  " << d.size() << " directions" << std::endl;
    NC::StableSum sw, swm;
    const auto ev = U::calcUcartEigenvalues( us, mt );
    for ( auto& e : d ) {
      REQUIRE( e.second > 0.0 );
      REQUIRE( e.first >= ev[0] * ( 1.0 - 1e-12 )
               && e.first <= ev[2] * ( 1.0 + 1e-12 ) );
      sw.add( e.second );
      swm.add( e.second * e.first );
    }
    REQUIRE( NC::floateq( sw.sum(), 1.0, 1e-12, 0.0 ) );
    //The average of n^T Ucart n over all directions is Ueq:
    REQUIRE( NC::floateq( swm.sum(), U::calcUeqFromUstar( us, mt ),
                          1e-12, 0.0 ) );
    //Isotropic: the same in all directions:
    const auto mc = U::cellMetric( 4.0, 4.0, 4.0, 90.0, 90.0, 90.0 );
    for ( auto& e : U::directionalMSDs(
            U::convertUstdToUstar( { 0.01, 0.01, 0.01, 0.0, 0.0, 0.0 }, mc ),
            mc ) )
      REQUIRE( NC::floateq( e.first, 0.01, 1e-12, 0.0 ) );
  }

  void testSiteSymmetry()
  {
    std::cout << "Site symmetry:" << std::endl;
    const auto& pm = sgsym( "6:b" );//Pm, mirror plane y=0
    const auto m = U::cellMetric( 5.0, 5.4, 6.0, 90.0, 95.0, 90.0 );
    const U::Uij u{ 0.012, 0.009, 0.015, 0.002, 0.003, 0.001 };
    const U::Uij umirror{ 0.012, 0.009, 0.015, 0.0, 0.003, 0.0 };
    const NC::Vector onplane( 0.1, 0.0, 0.2 );
    //On the mirror plane, U12 and U23 must vanish:
    REQUIRE( uijeq( U::snapUstdToSiteSymmetry( pm, m, onplane, umirror, "O" ),
                    umirror ) );
    U::Uij u2 = umirror;
    u2[3] = 0.0003;//(2.5% of the largest component)
    const auto s2 = U::snapUstdToSiteSymmetry( pm, m, onplane, u2, "O" );
    std::cout << "  snapped " << str( u2 ) << " -> " << str( s2 ) << std::endl;
    REQUIRE( uijeq( s2, umirror ) );
    u2[3] = 0.0012;//(10%)
    expectBadInput( [&pm,&m,&onplane,&u2](){
      U::snapUstdToSiteSymmetry( pm, m, onplane, u2, "O" ); } );
    //Accepted within the rounding errors of the input:
    const U::Uij err{ 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005 };
    REQUIRE( uijeq( U::snapUstdToSiteSymmetry( pm, m, onplane, u2, "O",
                                               0.05, err ), umirror ) );
    //Also the relative tolerance can be changed:
    REQUIRE( uijeq( U::snapUstdToSiteSymmetry( pm, m, onplane, u2, "O",
                                               0.2 ), umirror ) );
    expectBadInput( [&pm,&m,&onplane](){
      U::snapUstdToSiteSymmetry( pm, m, onplane,
                                 { 0.012, -0.009, 0.015, 0.0, 0.003, 0.0 },
                                 "O" ); } );
    //General positions are unconstrained, also when an image differs by more
    //than one cell in a coordinate (here y: 0.8 -> -0.8, a difference of -1.6,
    //which a min(d,1-d) check on unwrapped differences mistakes for being on
    //the mirror plane):
    REQUIRE( uijeq( U::snapUstdToSiteSymmetry( pm, m, { 0.1, 0.3, 0.2 }, u, "O" ),
                    u ) );
    REQUIRE( uijeq( U::snapUstdToSiteSymmetry( pm, m, { 0.1, 0.8, 0.2 }, u, "O" ),
                    u ) );
    //Cubic site symmetry (m-3m at the origin of Fm-3m): isotropic only:
    const auto& fm3m = sgsym( "225" );
    const auto mc = U::cellMetric( 4.0, 4.0, 4.0, 90.0, 90.0, 90.0 );
    const U::Uij ucub{ 0.0100, 0.0101, 0.0099, 0.0002, 0.0, 0.0 };
    const auto sc = U::snapUstdToSiteSymmetry( fm3m, mc, { 0.0, 0.0, 0.0 },
                                               ucub, "Al" );
    std::cout << "  snapped " << str( ucub ) << " -> " << str( sc ) << std::endl;
    REQUIRE( uijeq( sc, { 0.01, 0.01, 0.01, 0.0, 0.0, 0.0 } ) );
    expectBadInput( [&fm3m,&mc](){
      U::snapUstdToSiteSymmetry( fm3m, mc, { 0.0, 0.0, 0.0 },
                                 { 0.010, 0.012, 0.010, 0.0, 0.0, 0.0 },
                                 "Al" ); } );
  }

  void testOrbits()
  {
    std::cout << "Orbits:" << std::endl;
    const auto& pm = sgsym( "6:b" );
    const auto m = U::cellMetric( 5.0, 5.4, 6.0, 90.0, 95.0, 90.0 );
    const U::Uij u{ 0.012, 0.009, 0.015, 0.002, 0.003, 0.001 };
    const U::Uij uimg{ 0.012, 0.009, 0.015, -0.002, 0.003, -0.001 };
    //The mirror image (y -> -y) flips the signs of U12 and U23:
    const auto r = U::mapUstdToPosition( pm, m, { 0.1, 0.3, 0.2 },
                                         { 0.1, 0.7, 0.2 }, u );
    std::cout << "  mapped " << str( u ) << " -> " << str( r ) << std::endl;
    REQUIRE( uijeq( r, uimg ) );
    REQUIRE( uijeq( U::mapUstdToPosition( pm, m, { 0.1, 0.7, 0.2 },
                                          { 0.1, 0.3, 0.2 }, r ), u ) );
    //Also across cell boundaries:
    REQUIRE( uijeq( U::mapUstdToPosition( pm, m, { 0.1, 0.8, 0.2 },
                                          { 1.1, 0.2, -0.8 }, u ), uimg ) );
    REQUIRE( uijeq( U::mapUstdToPosition( pm, m, { 0.1, 0.3, 0.2 },
                                          { 0.1, 1.3, 0.2 }, u ), u ) );
    expectBadInput( [&pm,&m,&u](){
      U::mapUstdToPosition( pm, m, { 0.1, 0.3, 0.2 }, { 0.1, 0.4, 0.2 }, u );
    } );
    //An orbit given at each of its positions:
    const std::vector<NC::Vector> pos{ { 0.1, 0.3, 0.2 }, { 0.1, 0.7, 0.2 } };
    REQUIRE( uijeq( U::ustdOfOrbit( pm, m, pos, { u, uimg }, "O" ), u ) );
    expectBadInput( [&pm,&m,&pos,&u](){
      U::ustdOfOrbit( pm, m, pos, { u, u }, "O" ); } );
    //Small deviations are accepted (the first tensor is the result):
    U::Uij uimg2 = uimg;
    uimg2[3] = -0.0021;
    REQUIRE( uijeq( U::ustdOfOrbit( pm, m, pos, { u, uimg2 }, "O" ), u ) );
    //The first is also snapped to its site symmetry:
    const std::vector<NC::Vector> pos1{ { 0.1, 0.0, 0.2 } };
    U::Uij u3{ 0.012, 0.009, 0.015, 0.0003, 0.003, 0.0 };
    REQUIRE( uijeq( U::ustdOfOrbit( pm, m, pos1, { u3 }, "O" ),
                    { 0.012, 0.009, 0.015, 0.0, 0.003, 0.0 } ) );
  }
}

int main()
{
  testConversions();
  testBasics();
  testDirectionalMSDs();
  testSiteSymmetry();
  testOrbits();
  return 0;
}
