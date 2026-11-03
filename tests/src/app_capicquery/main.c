
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

#include "NCrystal/ncrystal.h"
#include <stdio.h>

/* Test of ncrystal_info_cquery (and that its availability can be detected */
/* with the preprocessor):                                                  */
#ifndef ncrystal_info_cquery
#  error "ncrystal_info_cquery not available"
#endif

static void query( ncrystal_info_t info, const char* key )
{
  ncrystal_cqueryresult_t res;
  int i;
  const int ok = ncrystal_info_cquery( info, key, &res );
  printf( "  %s:", key );
  if ( !ok ) {
    printf( " ERROR (is_error=%i): %s\n", res.is_error, res.val_str );
    return;
  }
  for ( i = 0; i < res.n_dbl; ++i )
    printf( " %.12g", res.val_dbl[i] );
  if ( res.has_str )
    printf( " \"%s\"", res.val_str );
  printf( "\n" );
}

static void queryall( const char* title, const char* cfgstr )
{
  const char* keys[] = { "spacegroup", "spacegroup_setting", "lattice_abc",
                         "lattice_angles", "cell_volume", "n_positions",
                         "n_atoms_average", "cell_numberdensity" };
  unsigned i;
  ncrystal_info_t info = ncrystal_create_info( cfgstr );
  printf( "%s (%s):\n", title, cfgstr );
  for ( i = 0; i < sizeof(keys)/sizeof(keys[0]); ++i )
    query( info, keys[i] );
  ncrystal_unref( &info );
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;
  queryall( "Crystal", "stdlib::Al_sg225.ncmat" );
  queryall( "Crystal, density override", "stdlib::Al_sg225.ncmat;density=0.5x" );
  ncrystal_register_in_mem_file_data(
    "nacl_occu.ncmat",
    "NCMAT v8\n@CELL\n  cubic 5.64\n@SPACEGROUP\n  225\n"
    "@ATOMSITES\n  Na 0 0 0\n  Cl 1/2 1/2 1/2\n"
    "@ATOMOCCUPANCIES\n  Na 0.9\n"
    "@DYNINFO\n  element Na\n  fraction 9/19\n  type vdosdebye\n"
    "  debye_temp 300\n"
    "@DYNINFO\n  element Cl\n  fraction 10/19\n  type vdosdebye\n"
    "  debye_temp 300\n" );
  queryall( "Partial occupancy", "nacl_occu.ncmat" );
  queryall( "Gas", "freegas::He/1e-3gcm3" );
  queryall( "Multiphase",
            "phases<0.5*stdlib::Al_sg225.ncmat&0.5*stdlib::Cu_sg225.ncmat>" );
  {
    ncrystal_info_t info = ncrystal_create_info( "stdlib::Al_sg225.ncmat" );
    printf( "List of keys:\n" );
    query( info, "list" );
    printf( "Bad keys:\n" );
    query( info, "foobar" );
    query( info, "" );
    ncrystal_unref( &info );
  }
  return 0;
}
