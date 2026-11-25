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


#include "NCrystal/internal/dyninfoutils/NCVDOSDebyeGroup.hh"
#include "NCrystal/internal/sabscatter/NCVDOSDebyeGroupScatter.hh"
#include "NCrystal/internal/sabscatter/NCSABScatter.hh"
#include "NCrystal/internal/vdos/NCVDOSDebyeXS.hh"
#include "NCrystal/internal/phys_utils/NCDebyeMSD.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCDataSources.hh"
#include <iostream>
#include <random>
#include <sstream>

//Tests of the grouped treatment of VDOSDebye entries (NCVDOSDebyeGroup.hh)
//and of VDOSDebyeGroupScatter.

namespace NC = NCrystal;
namespace NCDI = NCrystal::DynInfoUtils;
using GD = NCDI::VDOSDebyeGroupData;

#define REQUIRE(x) nc_assert_always(x)

namespace {
  using Entries = std::vector<NCDI::DebyeTempAndScale>;
  using Clusters = std::vector<std::vector<std::size_t>>;
  const NC::Temperature temp{ 300.0 };
  const NC::AtomMass mass{ 12.011 };

  //Pseudo-random entries with Debye temperatures in 250..450K:
  Entries testEntries()
  {
    Entries v;
    for ( auto i : NC::ncrange( 40 ) ) {
      const double x = std::fmod( i * 0.6180339887498949, 1.0 );
      const double y = std::fmod( i * 0.4142135623730951, 1.0 );
      v.push_back( { 0.001 + 0.03 * y,
                     NC::DebyeTemperature{ 250.0 + 200.0 * x } } );
    }
    return v;
  }

  double msdOf( const NCDI::DebyeTempAndScale& e )
  {
    return NC::debyeIsotropicMSD( e.debyeTemperature, temp, mass );
  }

  //Clusters as sets of Debye temperatures (independent of input order):
  std::vector<NC::VectD> asTDebye( const Entries& v, const Clusters& cl )
  {
    std::vector<NC::VectD> res;
    for ( auto& c : cl ) {
      res.emplace_back();
      for ( auto i : c )
        res.back().push_back( v.at( i ).debyeTemperature.dbl() );
    }
    return res;
  }

  template<class TErr>
  void expectError( std::function<void()> f )
  {
    bool threw = false;
    try { f(); } catch ( TErr& ) { threw = true; }
    REQUIRE( threw );
  }

  void testSettings()
  {
    for ( int raw : { 0, 1, 2, 3, 4, 5, 6, 100, 103, 106,
                      1000, 1003, 1005, 2000, 2005, 2006 } ) {
      const NC::VDOS::VDOSLux vl( raw );
      std::cout << "vdoslux=" << raw << ": lvl=" << vl.lvl() << " legacy="
                << vl.isLegacy() << " allowGrouping=" << vl.allowGrouping()
                << " tolerance(neff=1,3,4,10,100)=" << GD::tolerance( vl, 1.0 )
                << "," << GD::tolerance( vl, 3.0 ) << ","
                << GD::tolerance( vl, 4.0 ) << ","
                << GD::tolerance( vl, 10.0 ) << "," << GD::tolerance( vl, 100.0 )
                << " expansionLux="
                << vl.expansionLux().raw() << std::endl;
    }
  }

  void testClusters()
  {
    const Entries v = testEntries();
    Entries shuffled( v.rbegin(), v.rend() );
    std::swap( shuffled[3], shuffled[17] );
    std::size_t nprev = v.size() + 1;
    for ( double tol : { 0.0, 0.02, 0.05, 0.1, 0.2, 0.35, 10.0 } ) {
      const auto cl = GD::clusterEntries( v, temp, mass, tol );
      //Independent of entry order:
      REQUIRE( asTDebye( v, cl )
               == asTDebye( shuffled, GD::clusterEntries( shuffled, temp,
                                                         mass, tol ) ) );
      //A partition, with decreasing summed scales, contiguous msd ranges, and
      //entries by increasing msd:
      std::vector<int> seen( v.size(), 0 );
      double sprev = NC::kInfinity;
      std::vector<NC::PairDD> ranges;
      NC::StableSum sse, wtot;
      for ( auto& c : cl ) {
        REQUIRE( !c.empty() );
        NC::StableSum s, sx;
        for ( auto i : NC::ncrange( c.size() ) ) {
          ++seen.at( c[i] );
          if ( i )
            REQUIRE( msdOf( v[c[i-1]] ) <= msdOf( v[c[i]] ) );
          s.add( v[c[i]].scale );
          sx.add( v[c[i]].scale * std::log( msdOf( v[c[i]] ) ) );
        }
        REQUIRE( s.sum() <= sprev );
        sprev = s.sum();
        ranges.emplace_back( msdOf( v[c.front()] ), msdOf( v[c.back()] ) );
        const double mu = sx.sum() / s.sum();
        for ( auto i : c )
          sse.add( v[i].scale * NC::ncsquare( std::log( msdOf( v[i] ) ) - mu ) );
        wtot.add( s.sum() );
      }
      for ( auto n : seen )
        REQUIRE( n == 1 );
      std::sort( ranges.begin(), ranges.end() );
      for ( auto i : NC::ncrange( std::size_t( 1 ), ranges.size() ) )
        REQUIRE( ranges[i-1].second < ranges[i].first );
      REQUIRE( sse.sum() <= ( tol * tol + 1e-12 ) * wtot.sum() );
      REQUIRE( cl.size() <= nprev );
      nprev = cl.size();
      std::cout << "tolerance=" << tol << ": " << v.size() << " entries -> "
                << cl.size() << " clusters" << std::endl;
    }
    using BadInput = NC::Error::BadInput;
    expectError<BadInput>( [](){ GD::clusterEntries( {}, temp, mass, 0.1 ); } );
    expectError<BadInput>( [&v](){ GD::clusterEntries( v, temp, mass, -0.1 ); } );
  }

  //Max relative deviation of the cross sections from direct VDOSDebyeXS
  //evaluation of all entries (and checks consistency with the clusters):
  double checkXS( const GD& d, const Entries& v )
  {
    double wsum = 0.0;
    for ( auto& e : v )
      wsum += e.scale;
    double maxdev = 0.0;
    for ( double ekin : NC::logspace( -5.0, 1.0, 25 ) ) {
      const NC::NeutronEnergy en{ ekin };
      double ref = 0.0;
      for ( auto& e : v )
        ref += e.scale * NC::VDOS::VDOSDebyeXS( e.debyeTemperature, temp,
                                                 mass ).sigma( en ) / wsum;
      const double xs = d.xs().eval( ekin );
      maxdev = NC::ncmax( maxdev, NC::ncabs( xs / ref - 1.0 ) );
      if ( d.clusters().size() > 1 ) {
        double xssum = 0.0;
        for ( auto& c : d.clusters() )
          xssum += c.xs.value().eval( ekin );
        REQUIRE( NC::floateq( xssum, xs, 1e-12, 0.0 ) );
      }
    }
    return maxdev;
  }

  void testData()
  {
    const Entries v = testEntries();
    const NC::VDOS::VDOSLux lux3( 3 );
    auto d = NCDI::createVDOSDebyeGroupData( v, temp, mass, lux3, 3 );
    const double tol3 = GD::tolerance( lux3, GD::effectiveNumberOfEntries( v ) );
    const auto cl = GD::clusterEntries( v, temp, mass, tol3 );
    REQUIRE( d->nEntries() == v.size() );
    REQUIRE( d->clusters().size() == cl.size() && cl.size() > 1 );
    REQUIRE( d->temperature() == temp && d->mass() == mass );
    std::size_t ntot = 0;
    double stot = 0.0, sprev = NC::kInfinity;
    for ( auto& c : d->clusters() ) {
      REQUIRE( !c.sabext->supportsXS() && c.xs.has_value() );
      REQUIRE( c.scale <= sprev );
      sprev = c.scale;
      ntot += c.n;
      stot += c.scale;
    }
    REQUIRE( ntot == v.size() );
    REQUIRE( NC::floateq( stot, 1.0, 1e-12, 0.0 ) );
    REQUIRE( checkXS( *d, v ) < 3e-4 );
    std::cout << "VDOSDebyeGroupData: " << d->clusters().size()
              << " clusters, cross sections OK" << std::endl;

    //Cache: hit regardless of entry order, miss for other vdoslux:
    Entries shuffled( v.rbegin(), v.rend() );
    REQUIRE( &*d == &*NCDI::createVDOSDebyeGroupData( shuffled, temp, mass,
                                                      lux3, 3 ) );
    REQUIRE( &*d != &*NCDI::createVDOSDebyeGroupData( v, temp, mass,
                                                      NC::VDOS::VDOSLux( 2 ),
                                                      3 ) );

    //Clusters picked in proportion to their cross sections, checked via the
    //mean final energy:
    for ( double ekin : { 0.003, 0.025, 0.2 } ) {
      const NC::NeutronEnergy en{ ekin };
      const unsigned N = 100000;
      auto rng = NC::createBuiltinRNG( 12345 );
      NC::StableSum sum;
      for ( unsigned i = 0; i < N; ++i )
        sum.add( d->sampleScatter( *rng, en ).ekin.dbl() );
      NC::StableSum expected, var;
      const double xstot = d->xs().eval( ekin );
      for ( auto& c : d->clusters() ) {
        NC::StableSum sk, sk2;
        for ( unsigned i = 0; i < N; ++i ) {
          const double ef = c.sabext->sampleScatter( *rng, en ).ekin.dbl();
          sk.add( ef );
          sk2.add( ef * ef );
        }
        const double p = c.xs.value().eval( ekin ) / xstot;
        expected.add( p * sk.sum() / N );
        var.add( p * ( sk2.sum() / N ) );
      }
      const double sigma
        = std::sqrt( ( var.sum() - NC::ncsquare( expected.sum() ) ) / N );
      REQUIRE( NC::ncabs( sum.sum() / N - expected.sum() ) < 6.0 * sigma );
    }
    std::cout << "VDOSDebyeGroupData: cluster picking OK" << std::endl;

    //A single cluster: no separate table, and the same random numbers (and
    //samplings) as its expansion:
    const Entries one = { { 0.7, NC::DebyeTemperature{ 400.0 } } };
    auto d1 = NCDI::createVDOSDebyeGroupData( one, temp, mass, lux3, 3 );
    REQUIRE( d1->clusters().size() == 1 && !d1->clusters().front().xs.has_value() );
    REQUIRE( checkXS( *d1, one ) < 3e-4 );
    {
      auto rng1 = NC::createBuiltinRNG( 777 );
      auto rng2 = NC::createBuiltinRNG( 777 );
      for ( unsigned i = 0; i < 1000; ++i ) {
        const NC::NeutronEnergy en{ 0.001 * ( 1 + i % 50 ) };
        auto o1 = d1->sampleScatter( *rng1, en );
        auto o2 = d1->clusters().front().sabext->sampleScatter( *rng2, en );
        REQUIRE( o1.ekin.dbl() == o2.ekin.dbl() && o1.mu.dbl() == o2.mu.dbl() );
      }
    }
    std::cout << "VDOSDebyeGroupData: single cluster OK" << std::endl;

    using BadInput = NC::Error::BadInput;
    auto create = []( const Entries& e, int raw, int knllux )
    {
      NCDI::createVDOSDebyeGroupData( e, temp, mass, NC::VDOS::VDOSLux( raw ),
                                      knllux );
    };
    expectError<BadInput>( [&create](){ create( {}, 3, 3 ); } );
    expectError<BadInput>( [&create](){
      create( { { -1.0, NC::DebyeTemperature{ 300.0 } } }, 3, 3 ); } );
    for ( int raw : { 6, 103, 1003 } )
      expectError<BadInput>( [&create,&v,raw](){ create( v, raw, 3 ); } );
    expectError<BadInput>( [&create,&v](){ create( v, 3, -1 ); } );
  }

  std::string compNames( const char * cfgstr )
  {
    auto sc = NC::FactImpl::createScatter( NC::MatCfg( cfgstr ) );
    auto c = dynamic_cast<const NC::ProcImpl::ProcComposition*>( &*sc );
    std::ostringstream ss;
    if ( !c ) {
      ss << " " << sc->name();
      return ss.str();
    }
    for ( auto& comp : c->components() ) {
      ss << " " << comp.scale << "*" << comp.process->name();
      auto m = dynamic_cast<const NC::VDOSDebyeGroupScatter*>( &*comp.process );
      if ( m ) {
        NC::CachePtr cp;
        const NC::NeutronEnergy en{ 0.025 };
        const double xs = m->crossSectionIsotropic( cp, en ).dbl();
        REQUIRE( xs > 0.0 && std::isfinite( xs ) );
      }
    }
    return ss.str();
  }

  void testFactory()
  {
    //25 carbon labels with Debye temperatures 300..396K, hydrogen (one entry)
    //and oxygen (free gas, not affected):
    std::ostringstream ss;
    ss << "NCMAT v7\n@STATEOFMATTER\n  solid\n@DENSITY\n  1.0 g_per_cm3\n@ATOMDB\n";
    for ( auto i : NC::ncrange( 25 ) )
      ss << "  X" << i + 1 << " is C\n";
    for ( auto i : NC::ncrange( 25 ) )
      ss << "@DYNINFO\n  element X" << i + 1 << "\n  fraction 0.03\n"
         << "  type vdosdebye\n  debye_temp " << 300 + 4 * i << "\n";
    ss << "@DYNINFO\n  element H\n  fraction 0.15\n  type vdosdebye\n"
       << "  debye_temp 500\n";
    ss << "@DYNINFO\n  element O\n  fraction 0.1\n  type freegas\n";
    NC::DataSources::registerInMemoryFileData( "vdosdebyegroup.ncmat", ss.str() );
    for ( int lux : { 0, 3, 5, 103 } ) {
      std::ostringstream cfg;
      cfg << "vdosdebyegroup.ncmat;comp=inelas;vdoslux=" << lux;
      std::cout << "vdoslux=" << lux << ":" << compNames( cfg.str().c_str() )
                << std::endl;
    }
    //Mass outside the supported range (warning, separate expansions):
    NC::DataSources::registerInMemoryFileData(
      "vdosdebyegroup_heavy.ncmat",
      "NCMAT v7\n@STATEOFMATTER\n  solid\n@DENSITY\n  1.0 g_per_cm3\n"
      "@ATOMDB\n  U 300u 5fm 0b 0b\n"
      "@DYNINFO\n  element U\n  fraction 1\n  type vdosdebye\n  debye_temp 200\n" );
    //(created before printing, since it emits a warning):
    const std::string heavy
      = compNames( "vdosdebyegroup_heavy.ncmat;comp=inelas;vdoslux=1" );
    std::cout << "mass 300u:" << heavy << std::endl;
    //T>1000*TDebye (error message pointing to vdoslux=101):
    NC::DataSources::registerInMemoryFileData(
      "vdosdebyegroup_cold.ncmat",
      "NCMAT v7\n@STATEOFMATTER\n  solid\n@DENSITY\n  1.0 g_per_cm3\n"
      "@DYNINFO\n  element Al\n  fraction 1\n  type vdosdebye\n  debye_temp 0.2\n" );
    try {
      compNames( "vdosdebyegroup_cold.ncmat;comp=inelas;vdoslux=1" );
      REQUIRE( false );
    } catch ( NC::Error::CalcError& e ) {
      std::cout << "TDebye=0.2K: CalcError: " << e.what() << std::endl;
    }
  }
}

namespace {
  void testSmallGroups()
  {
    //Groups with n<=3 entries (random distinct Debye temperatures and
    //scales) get neither merging nor lowering at any vdoslux:
    std::mt19937 gen( 12345 );//fully specified, so portable
    auto rnd = [&gen]() { return gen() / 4294967296.0; };
    std::size_t ncheck = 0;
    for ( std::size_t n = 1; n <= 3; ++n ) {
      for ( int itry = 0; itry < 200; ++itry ) {
        Entries v;
        for ( std::size_t i = 0; i < n; ++i )
          v.push_back( { 1e-3 + rnd(),
                         NC::DebyeTemperature{ 100.0 + 50.0 * i
                                               + 40.0 * rnd() } } );
        const double neff = GD::effectiveNumberOfEntries( v );
        REQUIRE( neff <= n * ( 1.0 + 1e-12 ) );
        REQUIRE( !GD::lowersLevels( neff ) );
        for ( int raw = 0; raw <= 5; ++raw ) {
          const double tol = GD::tolerance( NC::VDOS::VDOSLux( raw ), neff );
          REQUIRE( tol == 0.0 );
          REQUIRE( GD::clusterEntries( v, temp, mass, tol ).size() == n );
          ++ncheck;
        }
      }
    }
    std::cout << "Small groups (n<=3): no merging or lowering in " << ncheck
              << " checks" << std::endl;
  }
}

int main()
{
  testSettings();
  testClusters();
  testData();
  testFactory();
  testSmallGroups();
  std::cout << "All ok" << std::endl;
  return 0;
}
