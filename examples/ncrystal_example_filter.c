
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
/* neutron beam is simply attenuated by a material:                           */
/* ncrystal_isomacroxstable provides a table of the macroscopic cross         */
/* section, and filter_macroxs below evaluates it. A neutron travelling a     */
/* distance L (in cm) through the material is transmitted with probability    */
/* exp(-macroxs*L), so in a Monte Carlo simulation, its weight can be         */
/* multiplied by that factor.                                                 */
/*                                                                            */
/* NOTE: Such usage is mainly intended for specialised use-cases who needs a  */
/* quick tabulated macroscopic cross section rather than the more detailed    */
/* cross sections and scattering sampling available in the usual NCrystal     */
/* interfaces.                                                                */
/*                                                                            */
/* The table returned by ncrystal_isomacroxstable is a single array of 3*n+4  */
/* doubles:                                                                   */
/*                                                                            */
/*   table[0]                 : magic number 4631638934963022.0               */
/*   table[1]                 : format version (1)                            */
/*   table[2]                 : number of points, n (at least 2)              */
/*   table[3] ... table[2+n]  : wavelengths (Aa), sorted, the first is 0      */
/*   table[3+n] ...           : macroscopic cross sections (1/cm)             */
/*   table[3+2*n] ...         : their slopes d(macroxs)/d(wavelength), in     */
/*                              1/(cm*Aa)                                     */
/*   table[3+3*n]             : magic number again                            */
/*                                                                            */
/* The magic number is an integer, whose bytes in memory read "NCryst0C" (on  */
/* little-endian machines). Between neighbouring points, the cross section is */
/* the cubic Hermite interpolation of their values and slopes, and it is      */
/* never negative. Two points with the same wavelength represent a            */
/* discontinuity (e.g. a Bragg edge) or a change of slope, and the second of  */
/* them is used at that wavelength. The first point is the limit for          */
/* wavelength -> 0, and beyond the last point its value and slope are         */
/* extrapolated linearly (clamped at 0).                                      */

#include "NCrystal/cinterface/ncrystal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int main(void) {
  /* declarations first (to support ancient pedantic C compilers) */
  unsigned i;
  size_t nbytes;
  double * nc_table;
  double * table;
  double xs;
  const double wavelengths[] = { 1.0, 3.0, 3.9, 4.0, 6.0, 20.0 };

  /* A beryllium filter cooled to 80 K: */
  nc_table = ncrystal_isomacroxstable( "stdlib::Be_sg194.ncmat;temp=80K",
                                       NULL );

  /* Copy the table to an array owned by the application (e.g. allocated such */
  /* that it is accessible on a GPU), and free the array from NCrystal:       */
  nbytes = ( 3 * (size_t)nc_table[2] + 4 ) * sizeof(double);
  table = (double*)malloc( nbytes );
  if ( !table ) {
    printf( "Memory allocation failed\n" );
    return 1;
  }
  memcpy( table, nc_table, nbytes );
  ncrystal_dealloc_doubleptr( nc_table );
  if ( !filter_check_table( table ) ) {
    printf( "Invalid table\n" );
    return 1;
  }

  printf( "Be filter at 80K (table with %u points):\n", (unsigned)table[2] );
  for ( i = 0; i < sizeof(wavelengths) / sizeof(wavelengths[0]); ++i ) {
    xs = filter_macroxs( table, wavelengths[i] );
    printf( "  wavelength %g Aa: macroscopic cross section %g /cm"
            " (mean free path %g cm)\n", wavelengths[i], xs, 1.0 / xs );
  }

  free( table );
  return 0;
}
