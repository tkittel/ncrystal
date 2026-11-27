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

// Tests of site orbits, the grouping of positions into orbits, and structure
// and cell checks of NCSGSymmetry.hh.

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include <iostream>
#include <sstream>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace {
  //A cell with lengths of about len (Aa) consistent with the constraints of
  //the setting (20 Aa by default, so general positions are not too close to
  //each other even with 192 of them):
  NC::CellParameters testCell( const NC::SGSymmetry& sym, double len = 20.0 )
  {
    using AC = NC::SGCellConstraints::AngleConstraint;
    using LC = NC::SGCellConstraints::LengthConstraint;
    const auto& cc = sym.cellConstraints();
    NC::CellParameters cp{ len, cc.b == LC::Free ? 1.1 * len : 0.0,
                           cc.c == LC::Free ? 1.2 * len : 0.0,
                           cc.alpha == AC::Free
                           ? ( cc.beta == AC::EqualToAlpha ? 80.0 : 90.0 ) : 0.0,
                           cc.beta == AC::Free ? 90.0 : 0.0,
                           cc.gamma == AC::Free ? 90.0 : 0.0 };
    cc.complete( cp );
    return cp;
  }

  //Orbit in the test cell, with strict tolerances by default:
  NC::SGSiteOrbit orbitOf( const NC::SGSymmetry& sym, const NC::Vector& site,
                           double occupancy = 1.0,
                           const NC::SGSiteTolerances& tol
                           = NC::SGSiteTolerances::strict(),
                           double len = 20.0 )
  {
    return NC::expandSiteToOrbit( sym, testCell( sym, len ), site, occupancy,
                                  tol );
  }

  std::vector<std::vector<unsigned>>
  groupOrbits( const NC::SGSymmetry& sym, const std::vector<NC::Vector>& pos )
  {
    return NC::groupPositionsIntoOrbits( sym, testCell( sym ), pos, 1.0,
                                         NC::SGSiteTolerances::strict() );
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    //Declared in the private header src/sgsym/NCSGExactOrbit.hh:
    Optional<std::pair<std::int64_t,std::int64_t>>
    rationalFromDouble( double x, std::int64_t maxDenominator );
  }
}

////////////////////////////////////////////////////////////////////////////////
// Expansion of sites into orbits
////////////////////////////////////////////////////////////////////////////////

//Tests expandSiteToOrbit: general and special positions in all 530 space
//group settings, symmetrisation of rounded special positions, rejection of
//sites near (but not on) special positions, and wrapping of coordinates.

namespace test_orbit {

namespace {
  double pdist( const NC::Vector& a, const NC::Vector& b )
  {
    double m = 0.0;
    const double d[3] = { a.x() - b.x(), a.y() - b.y(), a.z() - b.z() };
    for ( double x : d )
      m = NC::ncmax( m, NC::ncabs( x - std::round( x ) ) );
    return m;
  }

  //Value as recovered by rationalFromDouble, e.g. "1/10" or "0":
  std::string fracStr( double x )
  {
    const auto r = NC::detail::rationalFromDouble( x, 100 );
    nc_assert_always( r.has_value() );
    std::ostringstream ss;
    ss << r.value().first;
    if ( r.value().second != 1 )
      ss << '/' << r.value().second;
    return ss.str();
  }

  bool inUnitCell( const NC::Vector& v )
  {
    return ( v.x() >= 0.0 && v.x() < 1.0 && v.y() >= 0.0 && v.y() < 1.0
             && v.z() >= 0.0 && v.z() < 1.0 );
  }

  //(in a cell of about 5 Aa):
  NC::SGSiteOrbit expand( const char * sg, const NC::Vector& site,
                          double occupancy = 1.0,
                          const NC::SGSiteTolerances& tol
                          = NC::SGSiteTolerances::strict() )
  {
    return orbitOf( NC::SGSymmetry::get( NC::SpaceGroup( sg, NC::SpaceGroupParseMode::Relaxed ) ),
                    site, occupancy, tol, 5.0 );
  }

  void show( const char * sg, const NC::Vector& site, unsigned expected_mult )
  {
    const auto orbit = expand( sg, site );
    std::cout << "  " << sg << " site (" << NC::fmt( site.x(), "%.6g" )
              << ", " << NC::fmt( site.y(), "%.6g" ) << ", "
              << NC::fmt( site.z(), "%.6g" ) << "): "
              << orbit.positions.size() << " positions, site symmetry order "
              << orbit.siteSymmetryOrder << std::endl;
    REQUIRE( orbit.positions.size() == expected_mult );
  }

  void expectBad( const char * sg, const NC::Vector& site,
                  double occupancy = 1.0,
                  const NC::SGSiteTolerances& tol
                  = NC::SGSiteTolerances::strict() )
  {
    try {
      expand( sg, site, occupancy, tol );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  " << sg << " -> BadInput: " << e.what() << std::endl;
      return;
    }
    REQUIRE( false );
  }
}

void run()
{
  //General position in all settings:
  //(at least 1.4 Aa from its images in all settings, in the test cells):
  const NC::Vector generic( 0.34, 0.4643, 0.0873 );
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
    const auto& sym = NC::SGSymmetry::get( sg );
    const auto orbit = orbitOf( sym, generic );
    REQUIRE( orbit.positions.size() == sym.order() );
    REQUIRE( orbit.siteSymmetryOrder == 1 );
    REQUIRE( pdist( orbit.positions.front(), generic ) == 0.0 );
    for ( std::size_t i = 0; i < orbit.positions.size(); ++i ) {
      REQUIRE( inUnitCell( orbit.positions[i] ) );
      for ( std::size_t j = 0; j < i; ++j )
        REQUIRE( pdist( orbit.positions[i], orbit.positions[j] ) > 1e-2 );
    }
  }
  std::cout << "General positions OK in all 530 settings" << std::endl;

  std::cout << "Special positions:" << std::endl;
  show( "225", NC::Vector( 0.0, 0.0, 0.0 ), 4 );
  show( "225", NC::Vector( 0.25, 0.25, 0.25 ), 8 );
  show( "227:1", NC::Vector( 0.0, 0.0, 0.0 ), 8 );
  show( "227:2", NC::Vector( 0.125, 0.125, 0.125 ), 8 );
  show( "194", NC::Vector( 1.0/3, 2.0/3, 0.25 ), 2 );
  show( "166:h", NC::Vector( 0.0, 0.0, 0.0 ), 3 );
  show( "166:r", NC::Vector( 0.0, 0.0, 0.0 ), 1 );
  show( "14:b1", NC::Vector( 0.0, 0.0, 0.0 ), 2 );
  show( "62", NC::Vector( 0.1, 0.25, 0.2 ), 4 );
  show( "221", NC::Vector( 0.5, 0.5, 0.5 ), 1 );

  std::cout << "Rounded special positions:" << std::endl;
  //4 digits are symmetrised to the exact special position:
  for ( auto s : { NC::Vector( 0.3333, 0.6667, 0.25 ),
                   NC::Vector( 0.3333, 0.6666, 0.25 ),
                   NC::Vector( 0.33333, 0.66667, 0.25 ) } ) {
    const auto orbit = expand( "194", s );
    REQUIRE( orbit.positions.size() == 2 );
    REQUIRE( orbit.siteSymmetryOrder == 12 );
    REQUIRE( orbit.positions.front() == NC::Vector( 1.0/3, 2.0/3, 0.25 ) );
  }
  std::cout << "  194 site (0.3333, 0.6667, 0.25) etc.: symmetrised to"
            << " (1/3, 2/3, 1/4)" << std::endl;
  //Also with more digits than the exact expansion handles (z off by 1e-10,
  //within the tolerances), still giving exactly the doubles of 1/3, 2/3 and
  //1/4 (as all coordinates fixed by the site symmetry):
  for ( auto s : { NC::Vector( 0.333333333, 0.666666667, 0.2500000001 ),
                   NC::Vector( 0.3333333333333, 0.6666666666667, 0.25 ) } ) {
    const auto orbit = expand( "194", s );
    REQUIRE( orbit.positions.size() == 2 );
    REQUIRE( orbit.positions.front() == NC::Vector( 1.0/3, 2.0/3, 0.25 ) );
    REQUIRE( orbit.positions.back() == NC::Vector( 2.0/3, 1.0/3, 0.75 ) );
  }
  std::cout << "  194 site (0.333333333, 0.666666667, 0.2500000001) etc.:"
            << " symmetrised to (1/3, 2/3, 1/4)" << std::endl;
  //In all settings, sites with coordinates k/24 (many of them special
  //positions) and the same sites moved by about 1e-11 (expanded in doubles,
  //not exactly) give the same values for all coordinates fixed by the site
  //symmetry, while the others move by about 1e-11:
  {
    std::size_t nfixed = 0, nfree = 0;
    std::uint64_t lcg = 1;
    auto rnd24 = [&lcg]()
    {
      lcg = lcg * 6364136223846793005ULL + 1442695040888963407ULL;
      return static_cast<double>( ( lcg >> 33 ) % 24 ) / 24.0;
    };
    for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
      const auto& sym = NC::SGSymmetry::get(
                                NC::SpaceGroup{ NC::SpaceGroupSettingIdx{ hn } } );
      for ( int isite = 0; isite < 10; ++isite ) {
        const NC::Vector s( rnd24(), rnd24(), rnd24() );
        const NC::Vector s2( s.x() + 1.1e-11, s.y() + 3.7e-11,
                             s.z() + 8.3e-11 );
        const auto o1 = orbitOf( sym, s );
        const auto o2 = orbitOf( sym, s2 );
        REQUIRE( o1.positions.size() == o2.positions.size() );
        for ( auto& p : o1.positions ) {
          const NC::Vector* q = nullptr;
          for ( auto& e : o2.positions )
            if ( !q || pdist( p, e ) < pdist( p, *q ) )
              q = &e;
          const double a[3] = { p.x(), p.y(), p.z() };
          const double b[3] = { q->x(), q->y(), q->z() };
          for ( auto i : NC::ncrange( 3 ) ) {
            const double d = NC::ncabs( b[i] - a[i]
                                        - std::round( b[i] - a[i] ) );
            REQUIRE( d == 0.0 || ( d > 1e-13 && d < 1e-9 ) );
            ++( d == 0.0 ? nfixed : nfree );
          }
        }
      }
    }
    REQUIRE( nfixed > 1000 && nfree > 10000 );
    std::cout << "  sites k/24 in all settings: fixed coordinates identical"
              << " when expanded in doubles" << std::endl;
  }
  //3 digits are rejected (copies about 0.005 Aa apart) unless the occupancy
  //allows the copies to coexist:
  expectBad( "194", NC::Vector( 0.333, 0.667, 0.25 ) );
  expectBad( "227:2", NC::Vector( 0.125, 0.125, 0.12 ) );
  expectBad( "227:2", NC::Vector( 0.125, 0.125, 0.12 ), 0.3 );

  std::cout << "Sites near a mirror plane (Pm with y->-y, b=5.5 Aa):"
            << std::endl;
  {
    using Tol = NC::SGSiteTolerances;
    //Averaged (distance 2*0.00005*5.5=5.5e-4 Aa, below strict 1e-3 Aa):
    const auto o1 = expand( "6:b", NC::Vector( 0.1, 0.00005, 0.2 ) );
    REQUIRE( o1.positions.size() == 1 && o1.siteSymmetryOrder == 2 );
    REQUIRE( o1.positions.front().y() == 0.0 && o1.nearCopies == 0 );
    //Distance 5.5e-3 Aa: averaged with lenient tolerances, but with strict
    //ones a copy too close for occupancy 1:
    const auto o2 = expand( "6:b", NC::Vector( 0.1, 0.0005, 0.2 ), 1.0,
                            Tol::lenient() );
    REQUIRE( o2.positions.size() == 1 && o2.positions.front().y() == 0.0 );
    expectBad( "6:b", NC::Vector( 0.1, 0.0005, 0.2 ) );
    //Distance 0.066 Aa: a split site, allowed for occupancy 1/2:
    expectBad( "6:b", NC::Vector( 0.1, 0.006, 0.2 ), 0.51, Tol::lenient() );
    const auto o3 = expand( "6:b", NC::Vector( 0.1, 0.006, 0.2 ), 0.5 );
    REQUIRE( o3.positions.size() == 2 && o3.siteSymmetryOrder == 1 );
    REQUIRE( o3.nearCopies == 1 );
    REQUIRE( NC::ncabs( o3.nearestCopyDistance - 0.066 ) < 1e-12 );
    //Distance 0.66 Aa: separate positions, for any occupancy:
    const auto o4 = expand( "6:b", NC::Vector( 0.1, 0.06, 0.2 ) );
    REQUIRE( o4.positions.size() == 2 && o4.nearCopies == 0 );
    std::cout << "  y=0.00005: averaged (y=0); y=0.0005: averaged only with"
              << " lenient tolerances; y=0.006: split site (1 copy at "
              << NC::fmt( o3.nearestCopyDistance, "%.3g" ) << " Aa, needs occupancy"
              << " <= 1/2); y=0.06: two positions" << std::endl;
    //Across the periodic boundary (y=0.99995 is close to y=0.00005):
    const auto o5 = expand( "6:b", NC::Vector( 0.1, 0.99995, 0.2 ) );
    REQUIRE( o5.positions.size() == 1 && o5.positions.front().y() == 0.0 );
  }

  std::cout << "Wrapping and invalid input:" << std::endl;
  {
    const auto o = expand( "1", NC::Vector( 1.25, -0.5, 2.0 ) );
    REQUIRE( o.positions.size() == 1 );
    REQUIRE( o.positions.front().x() == 0.25 && o.positions.front().y() == 0.5
             && o.positions.front().z() == 0.0 );
    const auto o2 = expand( "1", NC::Vector( -0.0, -1e-30, 0.0 ) );
    REQUIRE( inUnitCell( o2.positions.front() ) );
    std::cout << "  (1.25, -0.5, 2.0) -> (0.25, 0.5, 0)" << std::endl;
  }
  expectBad( "1", NC::Vector( std::numeric_limits<double>::quiet_NaN(),
                              0.0, 0.0 ) );
  expectBad( "1", NC::Vector( 0.0, NC::kInfinity, 0.0 ) );

  //Rational coordinates are expanded exactly, so all positions are the
  //correctly rounded values of their fractions (e.g. 1-1/10 gives exactly the
  //double of 9/10, unlike 1.0-0.1 in double arithmetic):
  std::cout << "Exact expansion of rational coordinates:" << std::endl;
  {
    //Images of (1/10,7/20,7/10) have denominators dividing lcm(20,24)=120:
    std::size_t ncoords = 0;
    for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
      const NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
      const auto orbit = orbitOf( NC::SGSymmetry::get( sg ),
                                                NC::Vector( 0.1, 0.35, 0.7 ) );
      for ( auto& p : orbit.positions ) {
        for ( double c : { p.x(), p.y(), p.z() } ) {
          const auto r = NC::detail::rationalFromDouble( c, 1000000 );
          REQUIRE( r.has_value() && 120 % r.value().second == 0 );
          ++ncoords;
        }
      }
    }
    std::cout << "  (0.1, 0.35, 0.7) in all 530 settings: all " << ncoords
              << " coordinates exact" << std::endl;
    //Printing the positions:
    const auto o = expand( "62", NC::Vector( 0.1, 0.25, 0.2 ) );
    std::cout << "  62 site (0.1, 0.25, 0.2):";
    for ( auto& p : o.positions ) {
      std::cout << " (" << fracStr( p.x() ) << ", " << fracStr( p.y() )
                << ", " << fracStr( p.z() )
                << ")";
      REQUIRE( p.x() == 0.1 || p.x() == 0.4 || p.x() == 0.6 || p.x() == 0.9 );
    }
    std::cout << std::endl;
    //Rounded special positions are symmetrised to exactly 1/3, 2/3:
    const auto o2 = expand( "194", NC::Vector( 0.3333, 0.6667, 0.25 ) );
    REQUIRE( o2.positions.size() == 2 );
    REQUIRE( o2.positions.front().x() == 1.0 / 3.0 );
    REQUIRE( o2.positions.front().y() == 2.0 / 3.0 );
    REQUIRE( o2.positions.back().x() == 2.0 / 3.0 );
    REQUIRE( o2.positions.back().y() == 1.0 / 3.0 );
    std::cout << "  194 site (0.3333, 0.6667, 0.25): exactly (1/3, 2/3, 1/4)"
              << " and (2/3, 1/3, 3/4)" << std::endl;
    //Coordinates which are not simple rationals, or with a too large common
    //denominator (here about 1e16), are expanded in double arithmetic
    //instead:
    const double s2 = std::sqrt( 2.0 ) - 1.0;
    const auto o3 = expand( "62", NC::Vector( s2, 0.25, 0.2 ) );
    REQUIRE( o3.positions.size() == 4 && o3.positions.front().x() == s2 );
    const double big1 = 34000003.0 / 99999989.0;
    const double big2 = 46430011.0 / 99999971.0;
    const auto o4 = orbitOf( NC::SGSymmetry::get( NC::SpaceGroup( "191", NC::SpaceGroupParseMode::Relaxed ) ),
                             NC::Vector( big1, big2, 0.0873 ) );
    REQUIRE( o4.positions.size() == 24 );
    REQUIRE( o4.positions.front().x() == big1 );
    for ( auto& p : o4.positions ) {
      bool found = false;
      for ( auto& q : o4.positions )
        if ( pdist( p, NC::Vector( q.y(), q.x(), 1.0 - q.z() ) ) < 1e-15 )
          found = true;
      REQUIRE( found );//(y,x,-z) is an operation of 191
    }
    std::cout << "  Fallback to double arithmetic: OK" << std::endl;
  }
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Grouping of positions into orbits
////////////////////////////////////////////////////////////////////////////////

//Tests groupPositionsIntoOrbits: positions expanded from several sites (in a
//scrambled order, and shifted by lattice vectors) are grouped back into the
//orbits of those sites, and incomplete orbits are rejected.

namespace test_grouporbits {

namespace {
  struct Case {
    const char * sg;
    std::vector<NC::Vector> sites;
  };

  //Positions of the orbits of the sites, in a scrambled order and some
  //shifted by lattice vectors, with the site index of each:
  void expandAll( const NC::SGSymmetry& sym, const Case& c,
                  std::vector<NC::Vector>& pos, std::vector<unsigned>& site_of )
  {
    struct Entry { double key; NC::Vector p; unsigned site, idx; };
    std::vector<Entry> entries;
    for ( auto isite : NC::ncrange( c.sites.size() ) ) {
      const auto orbit = orbitOf( sym, c.sites.at( isite ) );
      for ( auto& p : orbit.positions ) {
        const double key = std::fmod( 13.7*p.x() + 7.3*p.y() + 3.1*p.z(), 1.0 );
        const NC::Vector shift( entries.size() % 2, 0.0,
                                entries.size() % 3 == 0 ? -1.0 : 0.0 );
        entries.push_back( { key, p + shift, static_cast<unsigned>( isite ),
                             static_cast<unsigned>( entries.size() ) } );
      }
    }
    //std::sort with the original order as tie-breaker (std::stable_sort
    //triggers deprecation warnings in some libstdc++ versions):
    std::sort( entries.begin(), entries.end(),
               []( const Entry& a, const Entry& b )
               { return a.key != b.key ? a.key < b.key : a.idx < b.idx; } );
    pos.clear();
    site_of.clear();
    for ( auto& e : entries ) {
      pos.push_back( e.p );
      site_of.push_back( e.site );
    }
  }
}

void run()
{
  const Case cases[] = {
    { "225", { NC::Vector( 0.0, 0.0, 0.0 ), NC::Vector( 0.25, 0.25, 0.25 ),
               NC::Vector( 0.34, 0.4643, 0.0873 ) } },
    { "227:2", { NC::Vector( 0.125, 0.125, 0.125 ),
                 NC::Vector( 0.5, 0.5, 0.5 ) } },
    { "194", { NC::Vector( 1.0/3, 2.0/3, 0.25 ),
               NC::Vector( 0.1, 0.2, 0.3 ) } },
    { "14:c1", { NC::Vector( 0.0, 0.0, 0.0 ), NC::Vector( 0.1, 0.2, 0.3 ) } },
  };
  for ( auto& c : cases ) {
    const auto& sym
      = NC::SGSymmetry::get( NC::SpaceGroup( c.sg, NC::SpaceGroupParseMode::Strict ) );
    std::vector<NC::Vector> pos;
    std::vector<unsigned> site_of;
    expandAll( sym, c, pos, site_of );
    const auto orbits = groupOrbits( sym, pos );
    REQUIRE( orbits.size() == c.sites.size() );
    std::size_t ntot = 0;
    std::cout << "  " << c.sg << ": " << pos.size() << " positions in "
              << orbits.size() << " orbits of sizes";
    for ( auto io : NC::ncrange( orbits.size() ) ) {
      const auto& orbit = orbits.at( io );
      REQUIRE( std::is_sorted( orbit.begin(), orbit.end() ) );
      if ( io > 0 )
        REQUIRE( orbits.at( io - 1 ).front() < orbit.front() );
      const unsigned site = site_of.at( orbit.front() );
      for ( auto i : orbit )
        REQUIRE( site_of.at( i ) == site );
      REQUIRE( orbit.size()
               == orbitOf( sym, c.sites.at( site ) ).positions.size() );
      ntot += orbit.size();
      std::cout << " " << orbit.size();
    }
    REQUIRE( ntot == pos.size() );
    std::cout << std::endl;
  }

  //Incomplete orbit:
  const auto& sym
    = NC::SGSymmetry::get( NC::SpaceGroup( "225", NC::SpaceGroupParseMode::Strict ) );
  std::vector<NC::Vector> pos;
  std::vector<unsigned> site_of;
  expandAll( sym, cases[0], pos, site_of );
  pos.pop_back();
  try {
    groupOrbits( sym, pos );
    REQUIRE( false );
  } catch ( NC::Error::BadInput& e ) {
    std::cout << "  225 with one position missing -> BadInput: " << e.what()
              << std::endl;
  }
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Structure symmetry check
////////////////////////////////////////////////////////////////////////////////

//Tests checkStructureSymmetry: structures expanded from a generic and a
//special site in each of the 530 settings must be matched to that setting
//with zero deviation, while perturbed or wrongly declared structures give the
//deviation. Also haveSameRotationParts, which tells settings apart only by
//their symmetry-equivalent reflections.

namespace test_structcheck {

namespace {
  //Positions by species: one generic orbit (species 0) and the orbit of the
  //origin (species 1):
  std::vector<std::vector<NC::Vector>> structure( const NC::SGSymmetry& sym )
  {
    std::vector<std::vector<NC::Vector>> res;
    res.push_back( orbitOf( sym, NC::Vector( 0.34, 0.4643, 0.0873 ) ).positions );
    res.push_back( orbitOf( sym, NC::Vector( 0.0, 0.0,
                                                           0.0 ) ).positions );
    return res;
  }
  bool sameRots( const char * a, const char * b )
  {
    return NC::haveSameRotationParts( NC::SGSymmetry::get( NC::SpaceGroup( a, NC::SpaceGroupParseMode::Relaxed ) ),
                                      NC::SGSymmetry::get( NC::SpaceGroup( b, NC::SpaceGroupParseMode::Relaxed ) ) );
  }
}

void run()
{
  unsigned nok = 0;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
    const bool rhombo = ( std::string( sg.choice() ) == "r" );
    const auto& sym = NC::SGSymmetry::get( sg );
    const auto res = NC::checkStructureSymmetry( sg.number(), rhombo,
                                                 structure( sym ) );
    REQUIRE( res.maxDeviation < 1e-12 );//(rounding only)
    if ( res.settingIdx.value().get() != hn ) {
      //A few settings have identical operations (e.g. 68:1 and 68:1ba-c),
      //in which case the first one is returned:
      const NC::SpaceGroup other{ res.settingIdx.value() };
      REQUIRE( res.settingIdx.value().get() < hn );
      auto ops1 = sym.operations();
      auto ops2 = NC::SGSymmetry::get( other ).operations();
      std::sort( ops1.begin(), ops1.end() );
      std::sort( ops2.begin(), ops2.end() );
      REQUIRE( ops1 == ops2 );
      std::cout << "  " << sg << " has the same operations as " << other
                << std::endl;
      continue;
    }
    ++nok;
  }
  std::cout << "Structures in all other " << nok << " settings matched to"
            << " their setting (deviation < 1e-12)" << std::endl;
  {
    auto pos = structure( NC::SGSymmetry::get( NC::SpaceGroup( "19", NC::SpaceGroupParseMode::Relaxed ) ) );
    pos.at(0).at(0) += NC::Vector( 0.01, 0.0, 0.0 );
    const auto res = NC::checkStructureSymmetry( 19, false, pos );
    const NC::SpaceGroup setting{ res.settingIdx.value() };
    std::cout << "SG-19 structure with one atom moved by 0.01: setting "
              << setting << ", deviation " << NC::fmt( res.maxDeviation, "%.3g" )
              << std::endl;
    REQUIRE( NC::floateq( res.maxDeviation, 0.01 ) );
  }
  for ( unsigned sgtrue : { 149u, 150u } ) {
    const unsigned sgclaim = ( sgtrue == 149 ? 150 : 149 );
    auto pos = structure( NC::SGSymmetry::get( NC::SpaceGroup( std::to_string( sgtrue ), NC::SpaceGroupParseMode::Relaxed ) ) );
    const auto res = NC::checkStructureSymmetry( sgclaim, false, pos );
    std::cout << "SG-" << sgtrue << " structure checked as SG-" << sgclaim
              << ": deviation " << NC::fmt( res.maxDeviation, "%.3g" ) << std::endl;
    REQUIRE( res.maxDeviation > 0.05 );
  }
  //Rhombohedral axes only match R settings:
  {
    auto pos = structure( NC::SGSymmetry::get( NC::SpaceGroup( "166:r", NC::SpaceGroupParseMode::Strict ) ) );
    REQUIRE( NC::checkStructureSymmetry( 166, true, pos ).settingIdx.value()
             == NC::SpaceGroup( "166:r", NC::SpaceGroupParseMode::Strict ).settingIdx() );
    REQUIRE( !NC::checkStructureSymmetry( 3, true, pos ).settingIdx.has_value() );
  }
  //Check against one given setting: the deviation is zero for the setting
  //the structure was expanded in, and not for another one:
  for ( auto pr : { std::make_pair( "227:2", "227:2" ),
                    std::make_pair( "227:2", "227:1" ),
                    std::make_pair( "14:c1", "14:c1" ),
                    std::make_pair( "14:c1", "14:b1" ) } ) {
    const NC::SpaceGroup sgtrue( pr.first, NC::SpaceGroupParseMode::Strict );
    const NC::SpaceGroup sgclaim( pr.second, NC::SpaceGroupParseMode::Strict );
    const auto res = NC::checkStructureSymmetry( sgclaim, structure( NC::SGSymmetry::get( sgtrue ) ) );
    std::cout << sgtrue << " structure checked against " << sgclaim
              << ": deviation ";
    if ( res.maxDeviation < 1e-12 )
      std::cout << "< 1e-12" << std::endl;
    else
      std::cout << NC::fmt( res.maxDeviation, "%.3g" ) << std::endl;
    REQUIRE( res.settingIdx.value() == sgclaim.settingIdx() );
    REQUIRE( sgtrue == sgclaim ? res.maxDeviation < 1e-12 : res.maxDeviation > 0.05 );
  }
  std::cout << "Same symmetry-equivalent reflections: 62 vs 62:cab "
            << sameRots( "62", "62:cab" ) << ", 227 vs 227:2 "
            << sameRots( "227", "227:2" ) << ", 14 vs 14:b2 "
            << sameRots( "14", "14:b2" ) << ", 14 vs 14:c1 "
            << sameRots( "14", "14:c1" ) << ", 166 vs 166:r "
            << sameRots( "166", "166:r" ) << std::endl;
  REQUIRE( sameRots( "62", "62:cab" ) && sameRots( "227", "227:2" )
           && sameRots( "14", "14:b2" ) && !sameRots( "14", "14:c1" )
           && !sameRots( "166", "166:r" ) );
}

}

//Settings in which a given unit cell can be expressed:
namespace test_cellsettings {

void run()
{
  struct Case { unsigned sg; NC::CellParameters cell; };
  const Case cases[] = { { 14, { 5, 6, 7, 90, 100, 90 } },
                         { 14, { 5, 6, 7, 90, 90, 100 } },
                         { 14, { 5, 6, 7, 90, 90, 90 } },
                         { 62, { 5, 6, 7, 90, 90, 90 } },
                         { 166, { 4, 4, 10, 90, 90, 120 } },
                         { 166, { 5, 5, 5, 60, 60, 60 } },
                         { 166, { 4, 0, 10, 0, 0, 0 } },
                         { 225, { 4, 0, 0, 0, 0, 0 } },
                         { 225, { 4, 5, 6, 90, 90, 90 } } };
  for ( auto& e : cases ) {
    const auto& c = e.cell;
    std::cout << "SG-" << e.sg << " with cell " << c.a << ' ' << c.b << ' '
              << c.c << " / " << c.alpha << ' ' << c.beta << ' ' << c.gamma
              << ':';
    const auto v = NC::settingsCompatibleWithCell( e.sg, c );
    for ( auto& sg : v )
      std::cout << ' ' << sg;
    if ( v.empty() )
      std::cout << " (none)";
    std::cout << std::endl;
  }
}

}

int main()
{
  std::cout << "==== Expansion of sites into orbits ====" << std::endl;
  test_orbit::run();
  std::cout << "==== Grouping of positions into orbits ====" << std::endl;
  test_grouporbits::run();
  std::cout << "==== Structure symmetry check ====" << std::endl;
  test_structcheck::run();
  std::cout << "==== Settings compatible with a cell ====" << std::endl;
  test_cellsettings::run();
  return 0;
}
