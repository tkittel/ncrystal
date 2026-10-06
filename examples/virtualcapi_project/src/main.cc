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

// Example of using NCrystal through the NCrystalVAPIType2V1 class (from
// NCrystalVAPIType2V1.hh/.cc, which can be copied into any project), without
// any build-time dependency on NCrystal.

#include "NCrystalVAPIType2V1.hh"
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace {
  //A neutron with the given wavelength (in Aa), moving along the z-axis:
  NCrystalVAPIType2V1::Neutron neutronWithWavelength( double wl )
  {
    const double ekin = 0.081804209605330899 / ( wl * wl );//h^2/(2*m) [eV*Aa^2]
    return NCrystalVAPIType2V1::Neutron{ ekin, 0.0, 0.0, 1.0 };
  }
}

int main( int argc, char ** argv )
{
  const std::string cfgstr = ( argc > 1 ? argv[1]
                               : "stdlib::Polyethylene_CH2.ncmat;temp=250K" );

  //Load NCrystal (without NCrystal, an application could simply disable the
  //features using it):
  std::shared_ptr<const NCrystalVAPIType2V1> nc;
  try {
    nc = NCrystalVAPIType2V1::load();
  } catch ( std::runtime_error& e ) {
    std::cout << "NCrystal is not available: " << e.what() << std::endl;
    return 1;
  }

  //Material information:
  auto info = nc->createInfo( cfgstr );
  std::cout << "Material: " << cfgstr << "\n"
            << "  density: " << info.density() << " g/cm3\n"
            << "  number density: " << info.numberDensity() << " atoms/Aa^3\n"
            << "  temperature: " << info.temperature() << " K\n"
            << "  composition:";
  for ( auto& c : info.composition() )
    std::cout << " (Z=" << c.Z << ", A=" << c.A << ", " << c.fraction << ")";
  std::cout << std::endl;

  //The composition broken down into isotopes, using natural abundances from
  //the application (here a small table):
  auto natab = []( unsigned long Z )
    -> std::vector<std::pair<unsigned long,double>>
  {
    if ( Z == 1 )
      return { { 1, 0.99985 }, { 2, 0.00015 } };
    if ( Z == 6 )
      return { { 12, 0.9893 }, { 13, 0.0107 } };
    return {};//not known
  };
  try {
    std::cout << "  isotopes:";
    for ( auto& c : info.composition( false, natab ) )
      std::cout << " (Z=" << c.Z << ", A=" << c.A << ", " << c.fraction << ")";
    std::cout << std::endl;
  } catch ( NCrystalVAPIType2V1::Error& e ) {
    std::cout << " not available (" << e.what() << ")" << std::endl;
  }

  //Cross sections:
  auto scatter = nc->createScatter( cfgstr );
  auto absorption = nc->createAbsorption( cfgstr );
  std::cout << "Cross sections (oriented: "
            << ( scatter.isOriented() ? "yes" : "no" ) << "):" << std::endl;
  for ( double wl : { 0.5, 1.0, 1.8, 4.0, 10.0 } ) {
    auto n = neutronWithWavelength( wl );
    std::cout << "  " << wl << " Aa: scattering "
              << scatter.crossSection( n ) << " barn, absorption "
              << absorption.crossSection( n ) << " barn" << std::endl;
  }

  //Sampling scatterings, with random numbers from the application:
  std::mt19937_64 generator( 12345 );
  std::uniform_real_distribution<double> uniform( 0.0, 1.0 );
  std::function<double()> rng = [&]() { return uniform( generator ); };
  std::cout << "Scatterings of 1.8 Aa neutrons:" << std::endl;
  for ( int i = 0; i < 5; ++i ) {
    auto n = neutronWithWavelength( 1.8 );
    scatter.sampleScatter( rng, n );
    std::cout << "  ekin = " << n.ekin << " eV, direction = (" << n.ux
              << ", " << n.uy << ", " << n.uz << ")" << std::endl;
  }

  //Threads: each thread uses its own clone of the Scatter object (the Info
  //object could be used by all threads directly):
  std::vector<std::thread> threads;
  std::vector<double> results( 4 );
  for ( std::size_t ithread = 0; ithread < results.size(); ++ithread ) {
    threads.emplace_back( [&scatter,&results,ithread]()
    {
      auto sc = scatter.clone();
      double sum = 0.0;
      for ( int i = 1; i <= 1000; ++i )
        sum += sc.crossSection( neutronWithWavelength( 0.01 * i ) );
      results[ithread] = sum;
    } );
  }
  for ( auto& t : threads )
    t.join();
  std::cout << "Sum of cross sections in each of " << results.size()
            << " threads:";
  for ( double r : results )
    std::cout << " " << r;
  std::cout << std::endl;

  //Errors are thrown as exceptions:
  try {
    nc->createInfo( "stdlib::doesnotexist.ncmat" );
  } catch ( NCrystalVAPIType2V1::Error& e ) {
    std::cout << "Expected error of type " << e.type() << ": " << e.what()
              << std::endl;
  }
  return 0;
}
