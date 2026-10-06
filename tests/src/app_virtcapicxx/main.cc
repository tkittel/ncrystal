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

//Tests of the NCrystalVAPIType2V1 C++ class from the example project in
//examples/virtualcapi_project (which applications copy into their own
//projects), comparing with the plain virtual C API. The NCrystalVAPIType2V1.hh
//and NCrystalVAPIType2V1.cc files next to this file are exact copies of those
//in the example project (as checked by "ncdevtool check copies").

#include "NCrystalVAPIType2V1.hh"
#include "NCrystal/ncrystal.h"
#include "NCrystal/virtualapi/ncvirtapi.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <thread>
#if defined (_WIN32) || defined (WIN32)
#  define NCTEST_WINDOWS
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

namespace {

  using NCV = NCrystalVAPIType2V1;

  void require( bool cond, const char * what )
  {
    if ( !cond ) {
      std::cout << "FAILED: " << what << std::endl;
      std::exit( 1 );
    }
  }

  //Let NCrystalVAPIType2V1::load() find the NCrystal library linked with this
  //test, via the NCRYSTAL_LIB environment variable (as in the Python API):
  void setupEnvironment()
  {
    const void * fct = reinterpret_cast<const void*>
      ( &ncrystal_access_virtual_c_api );
    std::string libpath;
#ifdef NCTEST_WINDOWS
    HMODULE hm = nullptr;
    require( GetModuleHandleExA( GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS
                                 | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                                 static_cast<LPCSTR>( fct ), &hm ) != 0,
             "GetModuleHandleExA" );
    char buf[4096];
    require( GetModuleFileNameA( hm, buf, sizeof(buf) ) > 0,
             "GetModuleFileNameA" );
    libpath = buf;
    _putenv_s( "NCRYSTAL_LIB", libpath.c_str() );
#else
    Dl_info dlinfo;
    require( dladdr( fct, &dlinfo ) != 0 && dlinfo.dli_fname, "dladdr" );
    libpath = dlinfo.dli_fname;
    setenv( "NCRYSTAL_LIB", libpath.c_str(), 1 );
#endif
#ifdef NCRYSTAL_NAMESPACE_PROTECTION
#  define NCTEST_STR2(x) #x
#  define NCTEST_STR(x) NCTEST_STR2(x)
    const char * ns = NCTEST_STR( NCRYSTAL_NAMESPACE_PROTECTION );
#else
    const char * ns = "";
#endif
#ifdef NCTEST_WINDOWS
    _putenv_s( "NCRYSTAL_LIB_NAMESPACE_PROTECTION", ns );
#else
    setenv( "NCRYSTAL_LIB_NAMESPACE_PROTECTION", ns, 1 );
#endif
  }

  const ncrystal_vapi_type2_v1_t * rawAPI()
  {
    return static_cast<const ncrystal_vapi_type2_v1_t *>
      ( ncrystal_access_virtual_c_api( 2001 ) );
  }

  NCV::Neutron neutron( double wl, double ux, double uy, double uz )
  {
    return NCV::Neutron{ ncrystal_wl2ekin( wl ), ux, uy, uz };
  }

  //Simple deterministic RNG:
  struct RNG {
    std::uint64_t s;
    double operator()()
    {
      s = s * 6364136223846793005ULL + 1442695040888963407ULL;
      return double( s >> 11 ) * ( 1.0 / 9007199254740992.0 );
    }
    static double generate( void * p ) { return ( *static_cast<RNG*>( p ) )(); }
  };

  const std::vector<std::string> cfgs = {
    "stdlib::Al_sg225.ncmat",
    "stdlib::Polyethylene_CH2.ncmat;temp=250K",
    "stdlib::Ge_sg227.ncmat;mos=40arcmin;dir1=@crys_hkl:5,1,1@lab:0,0,1"
    ";dir2=@crys_hkl:0,-1,1@lab:0,1,0"
  };

  const std::vector<double> wls = { 0.5, 1.0, 1.8, 3.0, 4.6, 20.0 };

  struct MyException : std::exception {
    const char * what() const noexcept override { return "MyException"; }
  };

  void testLoad()
  {
    auto api = NCV::load();
    require( api != nullptr, "load" );
    require( NCV::load() == api, "load returns the same object" );
    std::cout << "load: OK" << std::endl;
  }

  void testInfo()
  {
    auto api = NCV::load();
    auto raw = rawAPI();
    for ( auto& cfg : cfgs ) {
      auto info = api->createInfo( cfg );
      auto rh = raw->create_info( cfg.c_str(), nullptr );
      require( info.uniqueID() == std::uint64_t( raw->info_unique_id( rh ) ),
               "unique id" );
      require( info.density() == raw->info_density( rh ), "density" );
      require( info.numberDensity() == raw->info_number_density( rh ),
               "number density" );
      double t;
      require( raw->info_temperature( rh, &t, nullptr ) == 0, "raw temp" );
      require( info.temperature() == t, "temperature" );
      std::cout << cfg << ": density " << info.density() << " g/cm3, T="
                << info.temperature() << "K, composition:";
      for ( auto& c : info.composition() )
        std::cout << " (" << c.Z << "," << c.A << "," << c.fraction << ")";
      std::cout << std::endl;
      raw->deallocate_info( rh );
    }
    std::cout << "info: OK" << std::endl;
  }

  void testComposition()
  {
    auto api = NCV::load();
    auto info = api->createInfo( "stdlib::Polyethylene_CH2.ncmat" );
    unsigned ncalls = 0;
    auto natab = [&ncalls]( unsigned long Z )
      -> std::vector<std::pair<unsigned long,double>>
    {
      ++ncalls;
      if ( Z == 1 )
        return { { 1, 0.99985 }, { 2, 0.00015 } };
      if ( Z == 6 )
        return { { 12, 0.99 }, { 13, 0.01 } };
      return {};
    };
    auto cmp = info.composition( false, natab );
    require( ncalls == 2, "natab calls" );
    require( cmp.size() == 4, "four isotopes" );
    double sum = 0.0;
    std::cout << "isotopes:";
    for ( auto& c : cmp ) {
      std::cout << " (" << c.Z << "," << c.A << "," << c.fraction << ")";
      sum += c.fraction;
    }
    std::cout << std::endl;
    require( std::fabs( sum - 1.0 ) < 1e-12, "sum is 1" );
    //Natural elements are kept without natab:
    require( info.composition().size() == 2, "natural elements" );
    //Exceptions in natab are propagated unchanged:
    bool caught = false;
    try {
      info.composition( false, []( unsigned long )
                        -> std::vector<std::pair<unsigned long,double>>
                        { throw MyException(); } );
    } catch ( MyException& ) {
      caught = true;
    }
    require( caught, "natab exception propagated" );
    //Missing abundances give an NCrystal error:
    try {
      info.composition( false );
      require( false, "missing natab must fail" );
    } catch ( NCV::Error& e ) {
      std::cout << "expected error: [" << e.type() << "] " << e.what()
                << std::endl;
    }
    std::cout << "composition: OK" << std::endl;
  }

  void testProcesses()
  {
    auto api = NCV::load();
    auto raw = rawAPI();
    for ( auto& cfg : cfgs ) {
      auto sc = api->createScatter( cfg );
      auto ab = api->createAbsorption( cfg );
      auto rsc = raw->create_scatter( cfg.c_str(), nullptr );
      auto rab = raw->create_absorption( cfg.c_str(), nullptr );
      require( sc.isOriented() == ( raw->scatter_is_oriented( rsc ) != 0 ),
               "isOriented" );
      //Moving objects around, and cloning them:
      auto sc2 = sc.clone();
      NCV::Scatter sc3 = std::move( sc2 );
      auto ab2 = ab.clone();
      NCV::Absorption ab3 = ab.clone();
      ab3 = std::move( ab2 );
      RNG r1{ 17 }, r2{ 17 };
      std::function<double()> fr1 = [&r1]() { return r1(); };
      for ( double wl : wls ) {
        auto n = neutron( wl, 0.0, 0.6, 0.8 );
        double a[4] = { n.ekin, n.ux, n.uy, n.uz };
        double xs_sc, xs_ab;
        require( raw->scatter_cross_section( rsc, a, &xs_sc, nullptr ) == 0,
                 "raw xs" );
        require( raw->absorption_cross_section( rab, a, &xs_ab, nullptr ) == 0,
                 "raw xs" );
        require( sc.crossSection( n ) == xs_sc, "scatter xs" );
        require( sc3.crossSection( n ) == xs_sc, "scatter xs of clone" );
        require( ab.crossSection( n ) == xs_ab, "absorption xs" );
        require( ab3.crossSection( n ) == xs_ab, "absorption xs of clone" );
        //Sampling with std::function and with a function pointer, both
        //identical to the raw API with the same random numbers:
        for ( int k = 0; k < 20; ++k ) {
          auto n1 = neutron( wl, 0.0, 0.6, 0.8 );
          auto n2 = n1;
          double ra[4] = { n1.ekin, n1.ux, n1.uy, n1.uz };
          RNG r3 = r1;
          sc.sampleScatter( fr1, n1 );
          sc3.sampleScatter( RNG::generate, &r2, n2 );
          require( raw->sample_scatter( rsc, RNG::generate, &r3, ra,
                                        nullptr ) == 0, "raw sampling" );
          require( n1.ekin == ra[0] && n1.ux == ra[1] && n1.uy == ra[2]
                   && n1.uz == ra[3], "sampling as raw API" );
          require( n2.ekin == ra[0] && n2.ux == ra[1] && n2.uy == ra[2]
                   && n2.uz == ra[3], "sampling with function pointer" );
        }
      }
      std::cout << cfg << ": oriented=" << sc.isOriented()
                << ", xs(1.8Aa): scatter "
                << sc.crossSection( neutron( 1.8, 0.0, 0.6, 0.8 ) )
                << " barn, absorption "
                << ab.crossSection( neutron( 1.8, 0.0, 0.6, 0.8 ) )
                << " barn" << std::endl;
      raw->deallocate_scatter( rsc );
      raw->deallocate_absorption( rab );
    }
    std::cout << "scatter and absorption: OK" << std::endl;
  }

  void testErrors()
  {
    auto api = NCV::load();
    try {
      api->createScatter( "stdlib::doesnotexist.ncmat" );
      require( false, "missing file must fail" );
    } catch ( std::runtime_error& e ) {
      auto ncerr = dynamic_cast<NCV::Error*>( &e );
      require( ncerr && ncerr->type() == "FileNotFound", "FileNotFound" );
      std::cout << "expected error: [" << ncerr->type() << "] " << e.what()
                << std::endl;
    }
    auto info = api->createInfo( "phases<0.5*stdlib::Al_sg225.ncmat;temp=200K"
                                 "&0.5*stdlib::Cu_sg225.ncmat;temp=300K>" );
    try {
      info.temperature();
      require( false, "no temperature must fail" );
    } catch ( NCV::Error& e ) {
      require( e.type() == "MissingInfo", "MissingInfo" );
      std::cout << "expected error: [" << e.type() << "] " << e.what()
                << std::endl;
    }
    auto sc = api->createScatter( "stdlib::Al_sg225.ncmat" );
    auto n = neutron( 1.8, 0.0, 0.0, 0.0 );
    try {
      sc.crossSection( n );
      require( false, "invalid direction must fail" );
    } catch ( NCV::Error& e ) {
      std::cout << "expected error: [" << e.type() << "] " << e.what()
                << std::endl;
    }
    //Exceptions in the RNG are propagated unchanged, and the neutron is not
    //modified:
    n = neutron( 1.8, 0.0, 0.0, 1.0 );
    const auto n_orig = n;
    bool caught = false;
    try {
      sc.sampleScatter( []() -> double { throw MyException(); }, n );
    } catch ( MyException& ) {
      caught = true;
    }
    require( caught, "rng exception propagated" );
    require( n.ekin == n_orig.ekin && n.uz == n_orig.uz, "neutron unchanged" );
    std::cout << "errors: OK" << std::endl;
  }

  void testThreads()
  {
    auto api = NCV::load();
    auto info = api->createInfo( cfgs.at( 1 ) );
    auto sc = api->createScatter( cfgs.at( 1 ) );
    std::vector<double> ref;
    for ( int i = 1; i <= 200; ++i )
      ref.push_back( sc.crossSection( neutron( 0.05 * i, 0.0, 0.0, 1.0 ) ) );
    std::vector<int> ok( 6, 0 );
    std::vector<std::thread> threads;
    for ( std::size_t it = 0; it < ok.size(); ++it ) {
      threads.emplace_back( [&,it]()
      {
        auto mysc = sc.clone();
        bool good = ( info.density() == 0.92 );//shared Info object
        for ( int i = 200; i >= 1; --i )
          if ( mysc.crossSection( neutron( 0.05 * i, 0.0, 0.0, 1.0 ) )
               != ref.at( i - 1 ) )
            good = false;
        ok[it] = good ? 1 : 0;
      } );
    }
    for ( auto& t : threads )
      t.join();
    for ( auto v : ok )
      require( v == 1, "same results in threads" );
    std::cout << "threads: OK" << std::endl;
  }
}

int main()
{
  setupEnvironment();
  testLoad();
  testInfo();
  testComposition();
  testProcesses();
  testErrors();
  testThreads();
  return 0;
}
