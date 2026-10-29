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

#include "NCrystal/internal/utils/NCFMA.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include <iostream>
#include <cstring>

namespace NC = NCrystal;

namespace {
  std::uint64_t dblBits( double x )
  {
    std::uint64_t b;
    std::memcpy( &b, &x, sizeof(b) );
    return b;
  }
  double fromBits( std::uint64_t b )
  {
    double x;
    std::memcpy( &x, &b, sizeof(x) );
    return x;
  }
}

int main()
{
  //softFMA must agree bit-for-bit with std::fma for random arguments in wide
  //exponent ranges, and for (near) cancellations of a*b and c:
  NC::RandXRSRImpl rng( 1234 );
  auto rdbl = [&rng]( int emin, int emax )
  {
    const std::uint64_t m = rng.genUInt64() & 0xFFFFFFFFFFFFFull;
    const auto e = static_cast<std::uint64_t>( emin + 1023
                                               + static_cast<int>( rng.genUInt64() % ( emax - emin + 1 ) ) );
    const std::uint64_t s = rng.genUInt64() & 1u;
    return fromBits( ( s << 63 ) | ( e << 52 ) | m );
  };
  std::uint64_t nchecks = 0, nbad = 0;
  auto check = [&]( double a, double b, double c )
  {
    ++nchecks;
    const double r1 = std::fma( a, b, c );
    const double r2 = NC::softFMA( a, b, c );
    if ( dblBits( r1 ) != dblBits( r2 ) && !( std::isnan(r1) && std::isnan(r2) ) ) {
      if ( nbad++ < 10 )
        std::cout << "Mismatch for a=" << a << " b=" << b << " c=" << c << std::endl;
    }
  };
  for ( int i = 0; i < 400000; ++i ) {
    check( rdbl(-30,30), rdbl(-30,30), rdbl(-60,60) );
    check( rdbl(-500,500), rdbl(-500,500), rdbl(-1000,1000) );
    const double a = rdbl(-10,10), b = rdbl(-10,10);
    check( a, b, -a * b );
    check( a, b, std::nextafter( -a * b, ( i % 2 ) ? 1e300 : -1e300 ) );
    check( a, b, -a * b * ( 1.0 + std::ldexp( double( i % 1024 ) - 512.0, -60 ) ) );
  }
  //Special values:
  const double sv[] = { 0.0, -0.0, 1.0, -1.0, 1e-310, -1e-310, 1e308, -1e308,
                        std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN() };
  for ( auto a : sv )
    for ( auto b : sv )
      for ( auto c : sv )
        check( a, b, c );
  nc_assert_always( nbad == 0 );
  std::cout << "softFMA agrees with std::fma in " << nchecks << " checks" << std::endl;
  return 0;
}
