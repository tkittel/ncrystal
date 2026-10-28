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
// Checks that all RandIdxByWeightPicker picks (NCRandUtils.hh) give exactly
// the same indices as the static oneShotPick on the same (prefix of the)
// cumulative weights, for identical random numbers, and consume the same
// number of random numbers.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>
#include <sstream>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  using NC::VectD;
  using NC::ncrange;
  using Picker = NC::RandIdxByWeightPicker;

  std::uint64_t g_nchecks = 0;

  VectD cumulFromWeights( const VectD& w )
  {
    VectD c;
    c.reserve( w.size() );
    double s = 0.0;
    for ( auto e : w )
      c.push_back( s += e );
    return c;
  }

  VectD cumulWithTies( NC::RNG& rng, std::size_t n, double zeroFrac )
  {
    VectD w(n,0.0);
    for ( auto& e : w )
      if ( rng.generate() >= zeroFrac )
        e = rng.generate();
    return cumulFromWeights( w );
  }

  void checkArray( const char* label, const VectD& cumul )
  {
    const Picker picker{ VectD(cumul) };
    REQUIRE( picker.size() == cumul.size() );
    auto rng = NC::getRNG();
    VectD uvals = { 0.0, std::numeric_limits<double>::denorm_min(), 1e-300,
                    0.5, std::nextafter(1.0,0.0), 1.0 };
    for ( unsigned i = 0; i < 2000; ++i )
      uvals.push_back( rng->generate() );
    const std::size_t nn = cumul.size();
    for ( auto n : { std::size_t(1), std::size_t(2), std::size_t(4),
                     std::size_t(5), std::size_t(7), nn/3, nn } ) {
      if ( n < 1 || n > nn )
        continue;
      const NC::Span<const double> pre( cumul.data(), cumul.data() + n );
      for ( auto u : uvals ) {
        const auto iref = Picker::oneShotPick( u, pre );
        REQUIRE( iref < n );
        REQUIRE( picker.pick( u, n ) == iref );
        REQUIRE( NC::pickRandIdxByWeight( u, pre ) == iref );
        if ( n == nn )
          REQUIRE( picker.pick( u ) == iref );
        ++g_nchecks;
      }
      //RNG versions, including the number of random numbers consumed:
      auto rng1 = NC::createBuiltinRNG( 12345 + n );
      auto rng2 = NC::createBuiltinRNG( 12345 + n );
      auto rng3 = NC::createBuiltinRNG( 12345 + n );
      for ( unsigned i = 0; i < 2000; ++i ) {
        const auto iref = Picker::oneShotPick( *rng1, pre );
        REQUIRE( ( n == nn ? picker.pick( *rng2 ) : picker.pick( *rng2, n ) ) == iref );
        REQUIRE( NC::pickRandIdxByWeight( *rng3, pre ) == iref );
        ++g_nchecks;
      }
      const double rnext = rng1->generate();
      REQUIRE( rng2->generate() == rnext );
      REQUIRE( rng3->generate() == rnext );
    }
    std::cout << "OK: " << label << " (n=" << nn << ")" << std::endl;
  }
}

int main()
{
  checkArray( "single entry", VectD{ 1.0 } );
  checkArray( "two entries", VectD{ 1.0, 3.0 } );
  checkArray( "small, with zero weights", VectD{ 0.0, 0.0, 1.0, 1.0, 2.5, 2.5 } );
  checkArray( "all weights zero", VectD( 10, 0.0 ) );
  auto rng = NC::getRNG();
  for ( double zeroFrac : { 0.0, 0.3, 0.7, 0.95 } ) {
    std::ostringstream label;
    label << "cumul-with-ties, zeroFrac=" << zeroFrac;
    checkArray( label.str().c_str(), cumulWithTies( *rng, 2000, zeroFrac ) );
  }
  checkArray( "cumul-with-ties, large", cumulWithTies( *rng, 150000, 0.3 ) );
  {
    //Extremely clustered weights (one huge, rest tiny), and leading zeros:
    VectD w(3000,1e-12);
    w[1500] = 1.0;
    w[0] = w[1] = 0.0;
    checkArray( "clustered", cumulFromWeights( w ) );
  }
  {
    //Values within a few ULP of each other:
    VectD c;
    double v = 1.0;
    for ( unsigned i = 0; i < 50; ++i ) {
      c.push_back( v );
      v = std::nextafter( v, 2.0 );
    }
    checkArray( "near-ULP-spaced", c );
  }
  std::cout << "All " << g_nchecks << " checks passed." << std::endl;
  return 0;
}
