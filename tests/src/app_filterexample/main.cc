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

namespace NC = NCrystal;
namespace NCX = NCrystal::XSTable;

/* Reference implementation for a function which evaluates the table at the   */
/* given wavelength. The function only reads the arrays, and calls no other   */
/* functions, so it can also be used in code running on a GPU (for instance   */
/* with "#pragma acc routine seq" for OpenACC).                               */
double filter_macroxs( unsigned n, const double * wl, const double * macroxs,
                       const double * slope, double wavelength )
{
  /* Implements a cubic Hermite interpolation of the values and slopes at the */
  /* two neighbouring points, where at a discontinuity (two points with the   */
  /* same wavelength) the point above it is used. Beyond the last point, its  */
  /* value and slope are extrapolated linearly (clamped at 0).                */
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

  //Compare with the internal evaluation at the given wavelengths:
  void compare( const NC::VectD& wl, const NC::VectD& xs,
                const NC::VectD& slope, const NC::VectD& wls,
                const char * descr )
  {
    const auto n = static_cast<unsigned>( wl.size() );
    for ( auto w : wls ) {
      const double a = filter_macroxs( n, wl.data(), xs.data(), slope.data(), w );
      const double b = NCX::evalIsoMacroXSTable( wl.data(), xs.data(),
                                                 slope.data(), wl.size(), w );
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
  compare( wl, xs, { 1.0, 0.5, 2.0, -1.0, 1.0 }, wls, "synthetic, rising" );
  compare( wl, xs, { 1.0, 0.5, 2.0, -1.0, -0.5 }, wls, "synthetic, falling" );
  compare( wl, xs, { 1.0, 0.5, 2.0, -1.0, 0.0 }, wls, "synthetic, flat" );

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
    for ( auto x : t.wl )
      for ( auto f : { 1.0, 1.0 - 1e-12, 1.0 + 1e-12 } )
        w.push_back( x * f );
    compare( t.wl, t.xs, t.slope, w, cfg );
  }
  return 0;
}
