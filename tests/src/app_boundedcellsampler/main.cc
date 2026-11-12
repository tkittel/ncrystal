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

#include "NCSABBoundedCellSample.hh"
#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/sab/NCSABCellInteg.hh"
#include "NCrystal/internal/sab/NCSABSurveyor.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include "NCTestUtils/NCTestStat.hh"
#include <iostream>

//Compares BoundedCellSampler with RefCellSampler in single S(alpha,beta)
//cells, using Kolmogorov-Smirnov tests of the sampled (alpha,beta) points.

#define REQUIRE(x) nc_assert_always(x)

namespace NC = NCrystal;
namespace NS = NCrystal::SAB;

namespace {

  struct TestCell {
    double e, a1, a2, b1, b2;
    std::array<double,4> s;
  };

  NS::CellData makeCell( const TestCell& t )
  {
    NS::CellData c;
    c.a1 = t.a1;
    c.a2 = t.a2;
    c.b1 = t.b1;
    c.b2 = t.b2;
    for ( auto i : NC::ncrange( 4 ) ) {
      c.S[i] = t.s[i];
      c.logS[i] = NS::safeLogOrElse( t.s[i], 0.0 );
    }
    return c;
  }

  //Probability of the edge at beta=b1, from the integrals within KB(E) of the
  //cell with S=0 on the other edge:
  double probEdge1( const NS::CellData& c, NC::EDivKT e )
  {
    const auto scheme = NC::SABCfg::IntegrationScheme::MaxPrec;
    NS::CellData c1 = c;
    c1.S[2] = c1.S[3] = c1.logS[2] = c1.logS[3] = 0.0;
    NS::CellData c2 = c;
    c2.S[0] = c2.S[1] = c2.logS[0] = c2.logS[1] = 0.0;
    NC::StableSumKahan sum;
    NS::StdLogLinCellIntegrator::integrateWithinKB( c1, e, scheme, sum );
    const double w1 = sum.sum();
    NS::StdLogLinCellIntegrator::integrateWithinKB( c2, e, scheme, sum );
    REQUIRE( sum.sum() > 0.0 );
    return w1 / sum.sum();
  }

  //Samples of BoundedCellSampler (after a pack/unpack round trip of its
  //preprocessed data) and RefCellSampler:
  struct Samples { NC::VectD a, b; };
  Samples sampleBounded( const TestCell& t, const NS::CellData& c,
                         NC::RNG& rng, unsigned n )
  {
    const NC::EDivKT e{ t.e };
    NS::SABCellSurvey cellsurv( t.a1, t.a2, t.b1, t.b2, e );
    using BCS = NS::BoundedCellSampler;
    auto bcsdata = BCS::prepareBCSData( probEdge1( c, e ), cellsurv, c, e );
    NC::MixedDataVector buf;
    buf.append( 123.0 );
    const auto idx = BCS::pack( buf, bcsdata );
    REQUIRE( idx > 0 && idx < buf.size_bytes() );
    buf.append( 456.0 );
    bcsdata = BCS::unpack( buf, idx );
    BCS sampler( c, bcsdata, e );
    Samples res;
    for ( unsigned i = 0; i < n; ++i ) {
      auto ab = sampler.sampleAlphaBeta( rng );
      REQUIRE( NC::ncsquare( ab.alpha - ab.beta ) <= 4.0 * t.e * ab.alpha );
      res.a.push_back( ab.alpha );
      res.b.push_back( ab.beta );
    }
    return res;
  }

  Samples sampleRef( const TestCell& t, const NS::CellData& c,
                     NC::RNG& rng, unsigned n )
  {
    NS::RefCellSampler sampler( c, NC::EDivKT{t.e} );
    Samples res;
    for ( unsigned i = 0; i < n; ++i ) {
      auto ab = sampler.sampleAlphaBeta( rng );
      res.a.push_back( ab.alpha );
      res.b.push_back( ab.beta );
    }
    return res;
  }

  //Worst KS p-value of the projections alpha, beta, alpha-beta, alpha+beta:
  double worstPValue( const Samples& x, const Samples& y )
  {
    auto proj = []( const Samples& s, int mode )
    {
      NC::VectD v;
      for ( auto i : NC::ncrange( s.a.size() ) ) {
        const double a = s.a[i];
        const double b = s.b[i];
        v.push_back( mode == 0 ? a : ( mode == 1 ? b
                                       : ( mode == 2 ? a - b : a + b ) ) );
      }
      return v;
    };
    double worst = 1.0;
    for ( int mode = 0; mode < 4; ++mode )
      worst = NC::ncmin( worst, nctest::ksPValue( proj( x, mode ),
                                                  proj( y, mode ) ) );
    return worst;
  }

}

int main()
{
  const std::array<double,4> one{ { 1.0, 1.0, 1.0, 1.0 } };
  //Same cells as tests/scripts/sabcellsample.py:
  const std::vector<TestCell> cells = {
    { 0.06, 6.6767857524277714, 6.75015702442328,
      5.4022664171285335, 5.5029404418379908, one },
    { 1.5, 1.0, 2.0, -2.0, -1.0, { { 1.0, 20.0, 3.0, 4.0 } } },
    { 0.01, 5.8, 6.0, 5.2, 5.35, { { 1e3, 1, 1e3, 1 } } },
    { 0.505, 0.0, 1.0, -1.0, -0.5, { { 1, 1e10, 1, 1e10 } } },
    { 0.505, 0.0, 1.0, -1.0, -0.5, { { 1000, 1, 3, 2000 } } },
    { 0.505, 0.0, 1.0, -1.0, -0.5, { { 0, 1, 0, 1 } } },
    { 0.505, 0.0, 1.0, -1.0, -0.5, { { 1, 0, 3, 0 } } },
    { 0.55, 0.0, 1.0, -1.0, -0.5, { { 1, 0, 3, 0 } } },
    { 0.7, 0.0, 1.0, -1.0, -0.5, { { 1, 0, 3, 0 } } },
    { 0.01, 5.8, 6.0, 5.2, 5.35, { { 1, 1, 0, 0 } } },
    { 0.01, 5.8, 6.0, 5.2, 5.35, { { 0, 0, 1, 1 } } },
    { 0.01, 5.8, 6.0, 5.2, 5.35, { { 1, 1, 1, 1 } } },
    { 0.01, 5.8, 6.0, 5.2, 5.35, { { 1e3, 1, 1, 1 } } },
    { 0.01, 5.8, 6.0, 5.2, 5.35, { { 1, 1e3, 1, 1 } } },
    { 0.01, 5.8, 6.0, 5.2, 5.35, { { 1, 1, 1e3, 1 } } },
    { 0.01, 5.8, 6.0, 5.2, 5.35, { { 1, 1, 1, 1e3 } } },
    { 0.0001, 5.5, 6.0, 5.2, 6.0, { { 1, 1, 0, 0 } } },
    { 0.0001, 5.5, 6.0, 5.2, 6.0, { { 0, 0, 1, 1 } } },
    { 0.0001, 5.5, 6.0, 5.2, 6.0, { { 1, 1, 1, 1 } } },
    { 0.0001, 5.5, 6.0, 5.2, 6.0, { { 1.0, 1e4, 0, 0 } } },
    { 0.0001, 5.5, 6.0, 5.2, 6.0, { { 0, 0, 1e4, 1.0 } } },
    { 0.0001, 5.5, 6.0, 5.2, 6.0, { { 1.0, 1e4, 1e4, 1.0 } } },
    { 0.000001, 5.5, 6.0, 5.2, 6.0, { { 1.0, 1e4, 1e4, 1.0 } } },
    { 2.0, 0.5, 12.0, -3.0, 4.0, { { 0, 1, 0, 1 } } },
    { 2.0, 0.5, 12.0, -3.0, 4.0, one },
    { 1.0, 0.5, 12.0, -3.0, 4.0, { { 0.1, 1000, 2000, 0.2 } } },
    { 2.0, 0.5, 12.0, -3.0, 4.0, { { 0.1, 1000, 1000, 0.1 } } },
    { 2.0, 0.5, 12.0, -3.0, 4.0, { { 0.1, 1000, 0, 0 } } },
    { 2.0, 0.5, 12.0, -3.0, 4.0, { { 0, 0, 1000, 0.1 } } },
    { 2.0, 0.5, 12.0, -3.0, 4.0, { { 1, 1000, 1, 1 } } },
    { 2.0, 0.5, 12.0, -3.0, 4.0, { { 0, 0, 1, 1 } } },
    { 1e-10, 5.5, 6.0, 5.2, 6.0, { { 1.0, 1e50, 0, 0 } } },
    { 1e-15, 5.5, 6.0, 5.2, 6.0, { { 1.0, 1e5, 0, 0 } } },
    { 1e-15, 5.5, 6.0, 5.2, 6.0, { { 1.0, 1e250, 0, 0 } } },
    { 1, 1e-50, 1e-10, -1, 1, one },
    { 0.2812650004, 5e-10, 8e-10, -3e-5, -2e-05, one },
    { 0.2812650005, 5e-10, 8e-10, -3e-5, -2e-05, one },
    { 0.000395856, 8.36433e-50, 8.36433e-10, -2.78659e-05, -2.60257e-08,
      { { 2.09324e-51, 2.09324e-11, 2.09321e-51, 2.09321e-11 } } },
    { 1.0, 0.5, 12.0, -3.0, 4.0, { { 0.1, 1, 20000, 0.2 } } },
    { 2.0, 0.5, 3.0, -10.0, 40.5, { { 1.0, 20.0, 3.0, 4.0 } } },
    { 1e-5, 0.5, 12.0, 11.0, 12.0, one },
    { 1e-5, 11.0, 13.0, 0.5, 12.0, one },
    { 1e-5, 0.5, 12.0, 0.5, 12.0, one },
    { 1e-5, 0.5, 12.0, 11.0, 12.0, { { 1.0, 20.0, 3.0, 4.0 } } },
    { 0.01, 11.35, 13, 10, 12.0, { { 0, 1, 1, 0 } } },
    { 1e-5, 11.0, 13.0, 0.5, 12.0, { { 1.0, 20.0, 3.0, 0 } } },
    { 0.01, 5.5-1e-7, 5.5, 5, 5.92, { { 1.0, 1e-200, 1.0, 1e-200 } } },
    { 0.01, 5.5-1e-10, 5.5, 5, 5.92, { { 1.0, 1e-20, 1.0, 1e-20 } } },
    { 0.01, 5.5-1e-11, 5.5, 5, 5.92, { { 1e-20, 1.0, 1.0, 1e-20 } } },
    { 1e-5, 11.0, 13.0, 0.5, 12.0, { { 1.0, 20.0, 3.0, 1e-300 } } },
    { 1e-5, 0.5, 12.0, 0.5, 12.0, { { 1.0, 20.0, 3.0e-4, 4.0 } } },
    { 0.01, 5.499, 5.5, 5, 5.92, { { 1.0, 99, 1.0, 99 } } },
    { 0.01, 5.499, 5.5, 5, 5.92, { { 1.0, 101, 1.0, 101 } } },
    { 0.01, 5.499, 5.5, 5, 5.92, { { 1.0, 0.011, 1.0, 0.011 } } },
    { 0.01, 5.5-1e-11, 5.5, 5, 5.92, { { 1.0, 0.009, 1.0, 0.009 } } },
    { 40, 158.28456094789, 160.60010716438572,
      -2.2701714000978403e-14, -2.1202501083299951e-17,
      { { 0.028331793397541188, 0.027532855700519738,
          0.027532855700519426, 0.026748257342144074 } } },
    { 9.4e-9, 0.0, 5.1, -0.0124, 3.95,
      { { 6.7e-11, 4.4e-12, 0.144, 2.3e-10 } } },
    { 0.0092, 7.6e-5, 1.15, 0.0107, 5.07, { { 7e-4, 2.3e-9, 0.087, 2.1e-8 } } },
    { 1.25e-6, 2.7e-6, 0.347, 0.0023, 0.75,
      { { 6.7e-6, 2.1e-3, 0.172, 6.9e-11 } } },
    { 1.0, 0.0505102572167269, 1.0, 0.5, 1.0, { { 1e-30, 1.0, 0.0, 1.0 } } },
  };

  constexpr unsigned nsample = 20000;
  auto rng = NC::createBuiltinRNG( 123 );
  std::vector<std::pair<unsigned,double>> worst;
  unsigned icell = 0;
  for ( auto& t : cells ) {
    ++icell;
    //Only cells partly within KB(E) (as in the ["sab","sglcell"] query):
    NS::SABSurveyor surv( NC::VectD{ t.a1, t.a2 }, NC::VectD{ t.b1, t.b2 } );
    if ( !( t.e > surv.data().front().e_touch ) ) {
      std::cout << "Cell #" << icell << ": outside KB(E), skipped"
                << std::endl;
      continue;
    }
    const auto c = makeCell( t );
    const auto sb = sampleBounded( t, c, *rng, nsample );
    const auto sr = sampleRef( t, c, *rng, nsample );
    worst.emplace_back( icell, worstPValue( sb, sr ) );
  }
  //Threshold as in tests/scripts/sabcellsample.py:
  const double threshold = NC::ncmax( 1e-3, 0.05 / ( 4.0 * worst.size() ) );
  bool ok = true;
  for ( auto& w : worst ) {
    const bool cell_ok = w.second > threshold;
    std::cout << "Cell #" << w.first << ": "
              << ( cell_ok ? "OK" : "FAILED" ) << std::endl;
    if ( !cell_ok ) {
      std::cout << "   worst p-value " << w.second << std::endl;
      ok = false;
    }
  }
  REQUIRE( ok );
  return 0;
}
