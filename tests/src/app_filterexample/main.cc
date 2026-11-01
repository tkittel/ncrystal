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


//Test of the evaluation function of examples/ncrystal_example_filter.c, which
//is copied here verbatim (the copy is checked by "ncdevtool check copies"). It
//must give the same results as NCrystal's internal evaluation of the tables.

#include "NCrystal/internal/xstable/NCIsoMacroXSTable.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include <iostream>
#include <cstring>

namespace NC = NCrystal;
namespace NCX = NCrystal::XSTable;

/* Reference implementation for a function which checks the format of a       */
/* table, to protect against corrupted or truncated copies. It returns 1 if   */
/* the table is valid, and 0 otherwise. It should be called once, after the   */
/* table is obtained or copied (not for each evaluation).                     */
int filter_check_table( const double * table )
{
  const double magic = 4631638934963022.0;
  unsigned n;
  if ( !( table[0] == magic && table[1] == 1.0 ) )
    return 0;
  if ( !( table[2] >= 2.0 && table[2] <= 1e9 ) )
    return 0;
  n = (unsigned)table[2];
  if ( !( table[2] == (double)n ) )
    return 0;
  return table[3+3*n] == magic ? 1 : 0;
}

/* Reference implementation for a function which evaluates the table at the   */
/* given wavelength. The function only reads the table, and calls no other    */
/* functions, so it can also be used in code running on a GPU (for instance   */
/* with "#pragma acc routine seq" for OpenACC).                               */
double filter_macroxs( const double * table, double wavelength )
{
  /* Implements a cubic Hermite interpolation of the values and slopes at the */
  /* two neighbouring points, where at a discontinuity (two points with the   */
  /* same wavelength) the point above it is used. Beyond the last point, its  */
  /* value and slope are extrapolated linearly (clamped at 0).                */
  const unsigned n = (unsigned)table[2];
  const double * wl = table + 3;
  const double * macroxs = wl + n;
  const double * slope = macroxs + n;
  unsigned lo, hi, mid;
  double h, u, v, value;
  if ( !( wavelength > wl[0] ) )
    return macroxs[0];
  if ( wavelength >= wl[n-1] ) {
    if ( slope[n-1] == 0.0 )
      return macroxs[n-1];
    value = macroxs[n-1] + ( wavelength - wl[n-1] ) * slope[n-1];
    return value > 0.0 ? value : 0.0;
  }
  /* Binary search for the interval with wl[lo] <= wavelength < wl[hi]: */
  lo = 0;
  hi = n - 1;
  while ( hi - lo > 1 ) {
    mid = lo + ( hi - lo ) / 2;
    if ( wl[mid] <= wavelength )
      lo = mid;
    else
      hi = mid;
  }
  h = wl[hi] - wl[lo];
  u = ( wavelength - wl[lo] ) / h;
  v = 1.0 - u;
  return ( v * v * ( ( 1.0 + 2.0 * u ) * macroxs[lo] + u * h * slope[lo] )
           + u * u * ( ( 3.0 - 2.0 * u ) * macroxs[hi] - v * h * slope[hi] ) );
}

namespace {

  //A synthetic block in the format of NCX::IsoMacroXSTable::raw():
  NC::VectD makeBlock( const NC::VectD& wl, const NC::VectD& xs,
                       const NC::VectD& slope )
  {
    const double magic = NCX::IsoMacroXSTable::magic();
    const auto version = NCX::IsoMacroXSTable::formatVersion();
    NC::VectD res = { magic, static_cast<double>( version ),
                      static_cast<double>( wl.size() ) };
    for ( auto v : { &wl, &xs, &slope } )
      res.insert( res.end(), v->begin(), v->end() );
    res.push_back( magic );
    return res;
  }

  //Compare with the internal evaluation at the given wavelengths:
  void compare( NC::Span<const double> block, const NC::VectD& wls,
                const char * descr )
  {
    nc_assert_always( filter_check_table( block.data() ) == 1 );
    const auto n = static_cast<std::size_t>( block[2] );
    nc_assert_always( block.size() == 3 * n + 4 );
    const double * wl = block.data() + 3;
    for ( auto w : wls ) {
      const double a = filter_macroxs( block.data(), w );
      const double b = NCX::evalIsoMacroXSTable( wl, wl + n, wl + 2 * n, n, w );
      if ( !( a == b || NC::ncabs( a - b ) <= 1e-14 * NC::ncmax( NC::ncabs( b ), 1e-300 ) ) ) {
        std::cout << descr << ": FAILED at wavelength " << w << ": " << a
                  << " vs. " << b << std::endl;
        nc_assert_always( false );
      }
    }
    std::cout << descr << ": " << wls.size() << " wavelengths OK" << std::endl;
  }

}

int main()
{
  //Synthetic tables, with a discontinuity, and different slopes of the last
  //point (positive, negative and zero):
  const NC::VectD wl = { 0.0, 1.0, 2.0, 2.0, 3.0 };
  const NC::VectD xs = { 1.0, 2.0, 3.0, 1.0, 2.0 };
  const NC::VectD wls = { -1.0, 0.0, 0.5, 1.0, 1.5, 1.999, 2.0, 2.5, 3.0, 3.5, 4.0, 100.0,
                          std::numeric_limits<double>::infinity() };
  compare( makeBlock( wl, xs, { 1.0, 0.5, 2.0, -1.0, 1.0 } ), wls,
           "synthetic, rising" );
  compare( makeBlock( wl, xs, { 1.0, 0.5, 2.0, -1.0, -0.5 } ), wls,
           "synthetic, falling" );
  compare( makeBlock( wl, xs, { 1.0, 0.5, 2.0, -1.0, 0.0 } ), wls,
           "synthetic, flat" );

  //The format check rejects other blocks:
  {
    const NC::VectD good = makeBlock( wl, xs, { 1.0, 0.5, 2.0, -1.0, 1.0 } );
    nc_assert_always( filter_check_table( good.data() ) == 1 );
    nc_assert_always( good.front() == 4631638934963022.0 );
    nc_assert_always( good.back() == 4631638934963022.0 );
    const char magic_bytes[] = "NCryst0C";
    nc_assert_always( std::memcmp( &good.front(), magic_bytes, 8 ) == 0
                      || std::memcmp( &good.front(), "C0tsyrCN", 8 ) == 0 );
    auto modified = [&good]( std::size_t i, double v )
    {
      NC::VectD b = good;
      b.at(i) = v;
      return filter_check_table( b.data() );
    };
    nc_assert_always( modified( 0, 0.0 ) == 0 );//magic
    nc_assert_always( modified( 1, 2.0 ) == 0 );//version
    nc_assert_always( modified( 2, 4.0 ) == 0 );//n (trailing magic not found)
    nc_assert_always( modified( 2, 5.5 ) == 0 );//n not an integer
    nc_assert_always( modified( 2, 1.0 ) == 0 );//n too small
    nc_assert_always( modified( good.size() - 1, 0.0 ) == 0 );//trailing magic
    std::cout << "filter_check_table: OK" << std::endl;
  }

  //Tables of materials, at random wavelengths (also beyond the table), at
  //all table points, and next to them:
  NC::RandXRSRImpl rng( 12345 );
  for ( auto cfg : { "stdlib::Al_sg225.ncmat",
                     "stdlib::Be_sg194.ncmat;temp=80K",
                     "stdlib::Polyethylene_CH2.ncmat",
                     "phases<0.3*stdlib::Al_sg225.ncmat&0.7*stdlib::Cu_sg225.ncmat>" } ) {
    auto t = NCX::createIsoMacroXSTable( NC::MatCfg( cfg ),
                                         NCX::IsoMacroXSTableParams() );
    NC::VectD w;
    for ( unsigned i = 0; i < 100000; ++i )
      w.push_back( std::pow( 10.0, -8.0 + 12.0 * rng.generate() ) );
    for ( auto x : t.wavelengths() )
      for ( auto f : { 1.0, 1.0 - 1e-12, 1.0 + 1e-12 } )
        w.push_back( x * f );
    compare( t.raw(), w, cfg );
  }
  return 0;
}
