
/******************************************************************************/
/*                                                                            */
/*  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   */
/*                                                                            */
/*  Copyright 2015-2026 NCrystal developers                                   */
/*                                                                            */
/*  Licensed under the Apache License, Version 2.0 (the "License");           */
/*  you may not use this file except in compliance with the License.          */
/*  You may obtain a copy of the License at                                   */
/*                                                                            */
/*      http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                            */
/*  Unless required by applicable law or agreed to in writing, software       */
/*  distributed under the License is distributed on an "AS IS" BASIS,         */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/*  See the License for the specific language governing permissions and       */
/*  limitations under the License.                                            */
/*                                                                            */
/******************************************************************************/

/* Example showing how to use NCrystal for a filter or a window, where a      */
/* neutron beam is simply attenuated by a material: ncrystal_filtertable      */
/* provides a table of the macroscopic cross section, and filter_macroxs      */
/* below evaluates it. A neutron travelling a distance L (in cm) through the  */
/* material is transmitted with probability exp(-macroxs*L), so in a Monte    */
/* Carlo simulation, its weight can be multiplied by that factor.             */
/*                                                                            */
/* NOTE: Such usage is mainly intended for specialised use-cases who needs a  */
/* quick tabulated macroscopic cross section rather than the more detailed    */
/* cross sections and scattering sampling available in the usual NCrystal     */
/* interfaces.                                                                */

#include "NCrystal/cinterface/ncrystal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int main(void) {
  /* declarations first (to support ancient pedantic C compilers) */
  unsigned n, i;
  double * nc_wl;
  double * nc_macroxs;
  double * nc_slope;
  double * wl;
  double * macroxs;
  double * slope;
  double xs;
  const double wavelengths[] = { 1.0, 3.0, 3.9, 4.0, 6.0, 20.0 };

  /* A beryllium filter cooled to 80 K: */
  ncrystal_filtertable( "stdlib::Be_sg194.ncmat;temp=80K",
                        &n, &nc_wl, &nc_macroxs, &nc_slope, NULL );

  /* Copy the table to arrays owned by the application (e.g. allocated such   */
  /* that they are accessible on a GPU), and free the arrays from NCrystal:   */
  wl = (double*)malloc( n * sizeof(double) );
  macroxs = (double*)malloc( n * sizeof(double) );
  slope = (double*)malloc( n * sizeof(double) );
  if ( !wl || !macroxs || !slope ) {
    printf( "Memory allocation failed\n" );
    return 1;
  }
  memcpy( wl, nc_wl, n * sizeof(double) );
  memcpy( macroxs, nc_macroxs, n * sizeof(double) );
  memcpy( slope, nc_slope, n * sizeof(double) );
  ncrystal_dealloc_doubleptr( nc_wl );
  ncrystal_dealloc_doubleptr( nc_macroxs );
  ncrystal_dealloc_doubleptr( nc_slope );

  printf( "Be filter at 80K (table with %u points):\n", n );
  for ( i = 0; i < sizeof(wavelengths) / sizeof(wavelengths[0]); ++i ) {
    xs = filter_macroxs( n, wl, macroxs, slope, wavelengths[i] );
    printf( "  wavelength %g Aa: macroscopic cross section %g /cm"
            " (mean free path %g cm)\n", wavelengths[i], xs, 1.0 / xs );
  }

  free( wl );
  free( macroxs );
  free( slope );
  return 0;
}
