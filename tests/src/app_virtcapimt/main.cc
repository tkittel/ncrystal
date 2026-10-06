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

//Multi-threaded usage of the virtual C API (NCrystal/virtualapi/ncvirtapi.h):
//each thread works with its own clones (one per material), and must get
//exactly the same results as a single-threaded reference. For type 2, the
//threads also use the same info handles concurrently.

#include "NCrystal/ncrystal.h"
#include "NCrystal/virtualapi/ncvirtapi.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

namespace {

  const ncrystal_vapi_type1_v2_t * api = nullptr;
  const ncrystal_vapi_type2_v1_t * api2 = nullptr;

  void require( bool cond, const char * what )
  {
    if ( !cond ) {
      std::cout << "FAILED: " << what << std::endl;
      std::exit( 1 );
    }
  }

  struct RNG {
    std::uint64_t s;
    static double generate( void * p )
    {
      auto& r = *static_cast<RNG*>( p );
      r.s = r.s * 6364136223846793005ULL + 1442695040888963407ULL;
      return double( r.s >> 11 ) * ( 1.0 / 9007199254740992.0 );
    }
  };

  const std::vector<const char *> cfgs = {
    "stdlib::Al_sg225.ncmat",
    "stdlib::Polyethylene_CH2.ncmat;temp=250K",
    "stdlib::Ge_sg227.ncmat;mos=40arcmin;dir1=@crys_hkl:5,1,1@lab:0,0,1"
    ";dir2=@crys_hkl:0,-1,1@lab:0,1,0"
  };

  constexpr unsigned npoints = 400;
  constexpr unsigned nthreads = 8;

  //Neutron state for point i (energies from 0.001 to 0.5 eV, varying
  //directions):
  void neutron( unsigned i, double * n )
  {
    const double t = double( i ) / ( npoints - 1 );
    n[0] = 0.001 * std::pow( 500.0, t );
    const double phi = 7.0 * t, costh = 2.0 * t - 1.0;
    const double sinth = std::sqrt( 1.0 - costh * costh );
    n[1] = sinth * std::cos( phi );
    n[2] = sinth * std::sin( phi );
    n[3] = costh;
  }

  //Results (cross section and sampled outcome) for point i of material imat,
  //calculated with the given handle. The point is sampled with RNG seed i+1.
  struct Result {
    double xs;
    double out[4];
    bool operator==( const Result& o ) const
    {
      return xs == o.xs && out[0] == o.out[0] && out[1] == o.out[1]
        && out[2] == o.out[2] && out[3] == o.out[3];
    }
  };

  Result evaluate( ncrystal_vapi_t1v2_scatter_t * h, unsigned i )
  {
    ncrystal_vapi_error_t err;
    Result r;
    double n[4];
    neutron( i, n );
    require( api->cross_section( h, n, &r.xs, &err ) == 0, "cross_section" );
    RNG rng{ i + 1 };
    require( api->sample_scatter( h, RNG::generate, &rng, n, &err ) == 0,
             "sample_scatter" );
    std::copy( n, n + 4, r.out );
    return r;
  }

  //Type 2: scatter and absorption cross sections and sampling (the absorption
  //cross section is stored in out[3], after the sampled energy and the
  //first two components of the direction):
  Result evaluate2( ncrystal_vapi_t2v1_scatter_t * sc,
                    ncrystal_vapi_t2v1_absorption_t * ab, unsigned i )
  {
    ncrystal_vapi_error_t err;
    Result r;
    double n[4];
    neutron( i, n );
    require( api2->scatter_cross_section( sc, n, &r.xs, &err ) == 0,
             "scatter_cross_section" );
    require( api2->absorption_cross_section( ab, n, &r.out[3], &err ) == 0,
             "absorption_cross_section" );
    RNG rng{ i + 1 };
    require( api2->sample_scatter( sc, RNG::generate, &rng, n, &err ) == 0,
             "sample_scatter" );
    std::copy( n, n + 3, r.out );
    return r;
  }

  //Information from an info handle (which threads may share):
  std::vector<double> infoValues( const ncrystal_vapi_t2v1_info_t * h )
  {
    ncrystal_vapi_error_t err;
    std::vector<double> v = { api2->info_unique_id( h ),
                              api2->info_density( h ),
                              api2->info_number_density( h ) };
    double t;
    require( api2->info_temperature( h, &t, &err ) == 0, "temperature" );
    v.push_back( t );
    ncrystal_vapi_t2v1_component_t cmps[8];
    const std::size_t nc = api2->info_composition( h, 1, nullptr, nullptr,
                                                   cmps, 8, &err );
    require( nc > 0 && nc <= 8, "composition" );
    for ( std::size_t i = 0; i < nc; ++i ) {
      v.push_back( double( cmps[i].Z ) );
      v.push_back( double( cmps[i].A ) );
      v.push_back( cmps[i].fraction );
    }
    return v;
  }

  void testType2()
  {
    api2 = static_cast<const ncrystal_vapi_type2_v1_t*>
      ( ncrystal_access_virtual_c_api( 2001 ) );
    require( api2 != nullptr, "api2" );
    ncrystal_vapi_error_t err;

    //Base handles, and single-threaded reference results:
    std::vector<ncrystal_vapi_t2v1_info_t*> infos;
    std::vector<ncrystal_vapi_t2v1_scatter_t*> basesc;
    std::vector<ncrystal_vapi_t2v1_absorption_t*> baseab;
    std::vector<std::vector<Result>> ref( cfgs.size() );
    std::vector<std::vector<double>> refinfo;
    for ( std::size_t imat = 0; imat < cfgs.size(); ++imat ) {
      infos.push_back( api2->create_info( cfgs[imat], &err ) );
      basesc.push_back( api2->create_scatter( cfgs[imat], &err ) );
      baseab.push_back( api2->create_absorption( cfgs[imat], &err ) );
      require( infos.back() && basesc.back() && baseab.back(), "create" );
      refinfo.push_back( infoValues( infos.back() ) );
      for ( unsigned i = 0; i < npoints; ++i )
        ref[imat].push_back( evaluate2( basesc.back(), baseab.back(), i ) );
    }

    std::vector<std::vector<std::vector<Result>>> results( nthreads );
    std::vector<int> infook( nthreads, 0 );
    std::vector<std::thread> threads;
    for ( unsigned ithread = 0; ithread < nthreads; ++ithread ) {
      threads.emplace_back( [&,ithread]()
      {
        ncrystal_vapi_error_t terr;
        std::vector<ncrystal_vapi_t2v1_scatter_t*> sc;
        std::vector<ncrystal_vapi_t2v1_absorption_t*> ab;
        for ( std::size_t imat = 0; imat < cfgs.size(); ++imat ) {
          sc.push_back( api2->clone_scatter( basesc[imat], &terr ) );
          ab.push_back( api2->clone_absorption( baseab[imat], &terr ) );
          require( sc.back() && ab.back(), "clones in thread" );
        }
        auto& res = results[ithread];
        res.resize( cfgs.size(), std::vector<Result>( npoints ) );
        constexpr unsigned steps[nthreads] = { 1, 3, 7, 9, 11, 13, 17, 19 };
        bool ok = true;
        for ( unsigned k = 0; k < npoints; ++k ) {
          const unsigned i = ( k * steps[ithread] + ithread ) % npoints;
          for ( std::size_t imat = 0; imat < cfgs.size(); ++imat ) {
            res[imat][i] = evaluate2( sc[imat], ab[imat], i );
            if ( k % 50 == 0 && infoValues( infos[imat] ) != refinfo[imat] )
              ok = false;
          }
        }
        infook[ithread] = ok ? 1 : 0;
        for ( std::size_t imat = 0; imat < cfgs.size(); ++imat ) {
          api2->deallocate_scatter( sc[imat] );
          api2->deallocate_absorption( ab[imat] );
        }
      } );
    }
    for ( auto& t : threads )
      t.join();

    for ( unsigned ithread = 0; ithread < nthreads; ++ithread ) {
      require( infook[ithread] == 1, "same info values in all threads" );
      for ( std::size_t imat = 0; imat < cfgs.size(); ++imat )
        for ( unsigned i = 0; i < npoints; ++i )
          require( results[ithread][imat][i] == ref[imat][i],
                   "type 2: same results in threads as in the reference" );
    }
    for ( std::size_t imat = 0; imat < cfgs.size(); ++imat ) {
      api2->deallocate_info( infos[imat] );
      api2->deallocate_scatter( basesc[imat] );
      api2->deallocate_absorption( baseab[imat] );
    }
    std::cout << "Type 2: " << nthreads << " threads x " << cfgs.size()
              << " materials x " << npoints << " points (scattering and"
              << " absorption, and shared info handles): identical to the"
              << " single-threaded reference: OK" << std::endl;
  }

}

int main()
{
  api = static_cast<const ncrystal_vapi_type1_v2_t*>
    ( ncrystal_access_virtual_c_api( 1002 ) );
  require( api != nullptr, "api" );
  ncrystal_vapi_error_t err;

  //Base handles, and single-threaded reference results:
  std::vector<ncrystal_vapi_t1v2_scatter_t*> base;
  std::vector<std::vector<Result>> ref( cfgs.size() );
  for ( std::size_t imat = 0; imat < cfgs.size(); ++imat ) {
    base.push_back( api->create_scatter( cfgs[imat], &err ) );
    require( base.back() != nullptr, "create_scatter" );
    for ( unsigned i = 0; i < npoints; ++i )
      ref[imat].push_back( evaluate( base.back(), i ) );
  }

  //Each thread clones the base handles (concurrently), and processes all
  //points for all materials in a thread-specific order, interleaving the
  //materials. Thread 0 creates new handles instead of cloning:
  std::vector<std::vector<std::vector<Result>>> results( nthreads );
  std::vector<std::thread> threads;
  for ( unsigned ithread = 0; ithread < nthreads; ++ithread ) {
    threads.emplace_back( [&base,&results,ithread]()
    {
      ncrystal_vapi_error_t terr;
      std::vector<ncrystal_vapi_t1v2_scatter_t*> handles;
      for ( std::size_t imat = 0; imat < cfgs.size(); ++imat ) {
        handles.push_back( ithread == 0
                           ? api->create_scatter( cfgs[imat], &terr )
                           : api->clone_scatter( base[imat], &terr ) );
        require( handles.back() != nullptr, "handle in thread" );
      }
      auto& res = results[ithread];
      res.resize( cfgs.size(), std::vector<Result>( npoints ) );
      //Visit the points in different orders (the steps must be coprime with
      //npoints, to visit all points):
      constexpr unsigned steps[nthreads] = { 1, 3, 7, 9, 11, 13, 17, 19 };
      for ( unsigned k = 0; k < npoints; ++k ) {
        const unsigned i = ( k * steps[ithread] + ithread ) % npoints;
        for ( std::size_t imat = 0; imat < cfgs.size(); ++imat )
          res[imat][i] = evaluate( handles[imat], i );
      }
      for ( auto h : handles )
        api->deallocate_scatter( h );
    } );
  }
  for ( auto& t : threads )
    t.join();

  for ( unsigned ithread = 0; ithread < nthreads; ++ithread )
    for ( std::size_t imat = 0; imat < cfgs.size(); ++imat )
      for ( unsigned i = 0; i < npoints; ++i )
        require( results[ithread][imat][i] == ref[imat][i],
                 "same results in threads as in the reference" );
  for ( auto h : base )
    api->deallocate_scatter( h );
  std::cout << nthreads << " threads x " << cfgs.size() << " materials x "
            << npoints << " points: identical to the single-threaded"
            << " reference: OK" << std::endl;
  testType2();
  return 0;
}
