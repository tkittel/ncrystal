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

//Tests of the polynomial helpers in ncrystal_core/src/xstable/NCPolyBoundHelper.hh
//(through the functions for testing in its .cc file), and of the
//non-negativity check with Bernstein coefficients.

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include <iostream>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace XSTable {
    namespace detail {
      //Defined in ncrystal_core/src/xstable/NCPolyBoundHelper.cc (polynomials
      //as power-basis coefficients, lowest order first, at most 10):
      VectD polyMultiply( Span<const double> a, Span<const double> b );
      VectD polyComposeCubic( Span<const double> cubic, Span<const double> u );
      VectD polyBernsteinCoefficients( Span<const double> p );
      bool polyBernsteinNonNegative( Span<const double> bernstein );
    }
  }
}

namespace NCD = NCrystal::XSTable::detail;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  constexpr int nmax = 9;

  NC::VectD randomPoly( NC::RandXRSRImpl& rng, int deg )
  {
    NC::VectD p;
    for ( int k = 0; k <= deg; ++k )
      p.push_back( rng.generate() * 2.0 - 1.0 );
    return p;
  }

  double evalPoly( const NC::VectD& p, double z )
  {
    double r = 0.0;
    for ( std::size_t i = p.size(); i > 0; --i )
      r = r * z + p[i-1];
    return r;
  }

  double minOnGrid( const NC::VectD& p )
  {
    double m = NC::kInfinity;
    for ( int i = 0; i <= 20000; ++i )
      m = std::min( m, evalPoly( p, i / 20000.0 ) );
    return m;
  }

  //The polynomial of the Bernstein coefficients b, at z:
  double evalBernstein( const NC::VectD& b, double z )
  {
    NC::VectD tmp = b;
    const std::size_t n = b.size() - 1;
    for ( std::size_t r = 1; r <= n; ++r )
      for ( std::size_t i = 0; i <= n - r; ++i )
        tmp[i] = ( 1.0 - z ) * tmp[i] + z * tmp[i+1];
    return tmp[0];
  }

  bool nonNegative( const NC::VectD& p )
  {
    return NCD::polyBernsteinNonNegative( NCD::polyBernsteinCoefficients( p ) );
  }

}

int main()
{
  NC::RandXRSRImpl rng( 2468 );

  //Products, composition and Bernstein coefficients agree with direct
  //evaluation:
  for ( int k = 0; k < 2000; ++k ) {
    const auto p = randomPoly( rng, k % 4 );
    const auto q = randomPoly( rng, k % 3 );
    const auto u = randomPoly( rng, 2 );
    const auto c = randomPoly( rng, 3 );
    const auto pq = NCD::polyMultiply( p, q );
    const auto cu = NCD::polyComposeCubic( c, u );
    const auto b = NCD::polyBernsteinCoefficients( pq );
    REQUIRE( b.size() == nmax + 1 );
    for ( double z : { 0.0, 0.3, 0.77, 1.0 } ) {
      const double pz = evalPoly( p, z ), qz = evalPoly( q, z ), uz = evalPoly( u, z );
      REQUIRE( NC::floateq( evalPoly( pq, z ), pz * qz, 1e-13, 1e-13 ) );
      REQUIRE( NC::floateq( evalPoly( cu, z ),
                            c[0] + uz * ( c[1] + uz * ( c[2] + uz * c[3] ) ),
                            1e-12, 1e-12 ) );
      REQUIRE( NC::floateq( evalBernstein( b, z ), evalPoly( pq, z ), 1e-12, 1e-12 ) );
    }
  }
  std::cout << "Polynomial products, compositions and Bernstein coefficients: OK" << std::endl;

  //Non-negativity, compared with the minimum on a fine grid:
  unsigned nneg = 0, npos = 0;
  for ( int k = 0; k < 5000; ++k ) {
    auto p = randomPoly( rng, 1 + k % nmax );
    const double m = minOnGrid( p );
    const bool nonneg = nonNegative( p );
    if ( m < -1e-9 ) {
      REQUIRE( !nonneg );
      ++nneg;
    }
    if ( nonneg )
      REQUIRE( m >= -1e-12 );
    //Shifted to a positive minimum, it is always decided:
    p[0] += 1e-3 - m;
    REQUIRE( nonNegative( p ) );
    ++npos;
  }
  std::cout << "Non-negativity: " << nneg << " negative and " << npos
            << " positive random polynomials: OK" << std::endl;

  //Special cases: a narrow negative dip inside the interval, a double root
  //with a small positive margin, and a negative value only at an end point:
  {
    //100*(z-0.3)^2 - 1e-6 and (z-0.3)^2 + 1e-9:
    REQUIRE( !nonNegative( { 9.0 - 1e-6, -60.0, 100.0 } ) );
    REQUIRE( nonNegative( { 0.09 + 1e-9, -0.6, 1.0 } ) );
    REQUIRE( !nonNegative( { 0.999, -1.0 } ) );
  }
  std::cout << "Non-negativity special cases: OK" << std::endl;
  return 0;
}
