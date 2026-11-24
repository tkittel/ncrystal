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
// Checks of the conversion of crystal directions (Cartesian, hkl plane
// normals, lattice directions) and of the crystal rotations based on them.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  using NC::Vector;
  using NC::CrystalAxis;
  using NC::HKLPoint;
  using NC::UVWDir;
  using NC::LabAxis;
  using NC::OrientDir;

  bool close( const Vector& a, const Vector& b )
  {
    return ( a - b ).mag2() < 1e-24;
  }

  template<class TFct>
  void requireBadInput( TFct&& fct )
  {
    bool ok = false;
    try { fct(); } catch ( NC::Error::BadInput& ) { ok = true; }
    REQUIRE( ok );
  }

  struct Cell { double a, b, c, alpha, beta, gamma; };//Aa and degrees

  void testCell( const Cell& cell )
  {
    const auto lat = NC::getLatticeRot( cell.a, cell.b, cell.c,
                                        cell.alpha*NC::kDeg,
                                        cell.beta*NC::kDeg,
                                        cell.gamma*NC::kDeg );
    const auto reci = NC::getReciprocalLatticeRot( cell.a, cell.b, cell.c,
                                                   cell.alpha*NC::kDeg,
                                                   cell.beta*NC::kDeg,
                                                   cell.gamma*NC::kDeg );
    //Lattice vectors are the rows of lat, with a along x and b in the
    //xy-plane:
    const Vector va( lat[0][0], lat[0][1], lat[0][2] );
    const Vector vb( lat[1][0], lat[1][1], lat[1][2] );
    const Vector vc( lat[2][0], lat[2][1], lat[2][2] );
    REQUIRE( va.x() > 0.0 && va.y() == 0.0 && va.z() == 0.0 );
    REQUIRE( vb.y() > 0.0 && vb.z() == 0.0 && vc.z() > 0.0 );

    auto cart = [&reci]( const OrientDir::CrysDir& dir )
    {
      auto v = NC::crysDirToCartesian( dir, reci );
      REQUIRE( v.isUnitVector() );
      return v;
    };

    //The Cartesian crystal frame has x along a and z along c*:
    REQUIRE( close( cart( UVWDir{1.,0.,0.} ), Vector{1.,0.,0.} ) );
    REQUIRE( close( cart( HKLPoint{0.,0.,1.} ), Vector{0.,0.,1.} ) );
    REQUIRE( close( cart( CrystalAxis{0.,3.,4.} ), Vector{0.,0.6,0.8} ) );

    const std::vector<NC::ThreeVector> idx = {
      {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}, {1.,1.,0.}, {1.,-1.,0.},
      {0.,1.,1.}, {1.,2.,-1.}, {-3.,1.,2.}, {5.,1.,1.}, {0.5,-0.25,2.5}
    };
    for ( auto& e : idx ) {
      const UVWDir uvw{ e[0], e[1], e[2] };
      const HKLPoint hkl{ e[0], e[1], e[2] };
      //Lattice directions follow the lattice vectors:
      REQUIRE( close( cart( uvw ),
                      ( va*e[0] + vb*e[1] + vc*e[2] ).unit() ) );
      //Plane normals follow the reciprocal lattice vectors, defined by
      //a_i*b_j = 2pi*delta_ij:
      const Vector g = reci * Vector( e[0], e[1], e[2] );
      REQUIRE( close( cart( hkl ), g.unit() ) );
      REQUIRE( NC::floateq( g.dot(va), NC::k2Pi*e[0], 1e-12, 1e-12 ) );
      REQUIRE( NC::floateq( g.dot(vb), NC::k2Pi*e[1], 1e-12, 1e-12 ) );
      REQUIRE( NC::floateq( g.dot(vc), NC::k2Pi*e[2], 1e-12, 1e-12 ) );
      //Length is irrelevant:
      REQUIRE( close( cart( uvw ), cart( UVWDir{ 7*e[0], 7*e[1], 7*e[2] } ) ) );
      REQUIRE( close( cart( hkl ), cart( HKLPoint{ 7*e[0], 7*e[1], 7*e[2] } ) ) );
      //Zone law: [uvw] is in the plane (hkl) when hu+kv+lw=0.
      for ( auto& e2 : idx ) {
        const double zone = e[0]*e2[0] + e[1]*e2[1] + e[2]*e2[2];
        const double dot = cart( uvw ).dot( cart( HKLPoint{ e2[0], e2[1],
                                                            e2[2] } ) );
        REQUIRE( zone == 0.0 ? NC::ncabs(dot) < 1e-12 : dot*zone > 0.0 );
      }
    }

    //Crystal rotations take crystal directions to the lab frame, with the
    //primary direction exactly as requested and the secondary direction as
    //close as possible:
    for ( auto& e : idx ) {
      for ( auto& e2 : idx ) {
        const OrientDir::CrysDir c1 = UVWDir{ e[0], e[1], e[2] };
        const OrientDir::CrysDir c2 = HKLPoint{ e2[0], e2[1], e2[2] };
        if ( cart(c1).isParallel( cart(c2), 1e-3 ) )
          continue;
        const Vector l1( 0.3, -1.2, 0.7 ), l2( 1., 1., 1. );
        for ( bool swap : { false, true } ) {
          const OrientDir dir1{ swap ? c2 : c1, l1.as<LabAxis>() };
          const OrientDir dir2{ swap ? c1 : c2, l2.as<LabAxis>() };
          const auto rot
            = NC::verifyLatticeOrientDefAndConstructCrystalRotation
            ( dir1, dir2, NC::kPi, reci );
          REQUIRE( NC::floateq( rot.determinant(), 1.0, 1e-12, 1e-12 ) );
          REQUIRE( close( rot * cart( dir1.crystal ), l1.unit() ) );
          const Vector v2 = rot * cart( dir2.crystal );
          REQUIRE( NC::ncabs( v2.dot( l1.cross(l2) ) ) < 1e-12 );
          REQUIRE( v2.dot( l2 - l1 * ( l2.dot(l1) / l1.mag2() ) ) > 0.0 );
          //With consistent angles, the tolerance is not needed:
          const OrientDir dir2b{ dir2.crystal, v2.as<LabAxis>() };
          const auto rot2
            = NC::verifyLatticeOrientDefAndConstructCrystalRotation
            ( dir1, dir2b, 1e-10, reci );
          REQUIRE( close( rot2 * cart( dir1.crystal ), l1.unit() ) );
          REQUIRE( close( rot2 * cart( dir2.crystal ), v2 ) );
        }
      }
    }
  }

  void testBadInput()
  {
    const LabAxis l1{ 0., 0., 1. }, l2{ 0., 1., 0. };
    //Null vectors:
    requireBadInput( [&]()
    { NC::precheckLatticeOrientDir( { UVWDir{0.,0.,0.}, l1 } ); } );
    requireBadInput( [&]()
    { NC::precheckLatticeOrientDir( { UVWDir{1.,0.,0.}, {0.,0.,0.} } ); } );
    //Parallel directions:
    requireBadInput( [&]()
    {
      NC::precheckLatticeOrientDef( { UVWDir{1.,2.,3.}, l1 },
                                    { UVWDir{-2.,-4.,-6.}, l2 }, 1e-4 );
    } );
    requireBadInput( [&]()
    {
      NC::precheckLatticeOrientDef( { UVWDir{1.,2.,3.}, l1 },
                                    { UVWDir{1.,2.,4.}, l1 }, 1e-4 );
    } );
    //Not parallel (and mixed types can not be compared without a lattice):
    NC::precheckLatticeOrientDef( { UVWDir{1.,2.,3.}, l1 },
                                  { UVWDir{1.,2.,4.}, l2 }, 1e-4 );
    NC::precheckLatticeOrientDef( { UVWDir{1.,2.,3.}, l1 },
                                  { HKLPoint{1.,2.,3.}, l2 }, 1e-4 );
    //Types and comparisons:
    const OrientDir::CrysDir uvw = UVWDir{1.,2.,3.};
    const OrientDir::CrysDir hkl = HKLPoint{1.,2.,3.};
    const OrientDir::CrysDir crys = CrystalAxis{1.,2.,3.};
    REQUIRE( uvw.has_value<UVWDir>() && !uvw.has_value<HKLPoint>()
             && !uvw.has_value<CrystalAxis>() );
    REQUIRE( hkl.has_value<HKLPoint>() && crys.has_value<CrystalAxis>() );
    REQUIRE( uvw.get<UVWDir>() == ( UVWDir{1.,2.,3.} ) );
    REQUIRE( !( uvw == hkl ) && !( uvw == crys ) && !( hkl == crys ) );
    REQUIRE( uvw == OrientDir::CrysDir( UVWDir{1.,2.,3.} ) );
    REQUIRE( !( uvw == OrientDir::CrysDir( UVWDir{1.,2.,4.} ) ) );
  }
}

int main()
{
  testCell( { 5.6575, 5.6575, 5.6575, 90., 90., 90. } );//cubic
  testCell( { 4.2, 4.2, 6.1, 90., 90., 90. } );//tetragonal
  testCell( { 3.1, 4.2, 5.3, 90., 90., 90. } );//orthorhombic
  testCell( { 2.46, 2.46, 6.71, 90., 90., 120. } );//hexagonal
  testCell( { 5.1, 5.1, 5.1, 57.3, 57.3, 57.3 } );//rhombohedral
  testCell( { 5.0, 6.0, 7.0, 90., 105., 90. } );//monoclinic
  testCell( { 4.1, 5.2, 6.3, 75., 82., 101. } );//triclinic
  testBadInput();
  return 0;
}
