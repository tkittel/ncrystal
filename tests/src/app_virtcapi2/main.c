
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

/* Tests of type 2 version 1 of the virtual C API (NCrystal/virtualapi/       */
/* ncvirtapi.h), comparing with the normal C API and with type 1 version 2.  */

#include "NCrystal/ncrystal.h"
#include "NCrystal/virtualapi/ncvirtapi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static void require( int cond, const char * what )
{
  if ( !cond ) {
    printf( "FAILED: %s\n", what );
    exit( 1 );
  }
}

/* Simple deterministic RNG with numbers in [0,1): */
typedef struct { uint64_t s; } rng_t;

static double rng_fct( void * p )
{
  rng_t * r = (rng_t*)p;
  r->s = r->s * 6364136223846793005ULL + 1442695040888963407ULL;
  return (double)( r->s >> 11 ) * ( 1.0 / 9007199254740992.0 );
}

static const ncrystal_vapi_type2_v1_t * api = NULL;

#define NCFGS 4
static const char * cfgs[NCFGS] = {
  "stdlib::Al_sg225.ncmat",
  "stdlib::Polyethylene_CH2.ncmat;temp=250K",
  "stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat",
  "stdlib::Ge_sg227.ncmat;mos=40arcmin;dir1=@crys_hkl:5,1,1@lab:0,0,1"
  ";dir2=@crys_hkl:0,-1,1@lab:0,1,0"
};

#define NWL 6
static const double wls[NWL] = { 0.5, 1.0, 1.8, 3.0, 4.6, 20.0 };
static const double dir[3] = { 0.0, 0.6, 0.8 };

static void neutron( double * n, double wl )
{
  n[0] = ncrystal_wl2ekin( wl ); n[1] = dir[0]; n[2] = dir[1]; n[3] = dir[2];
}

static double ref_xs( ncrystal_process_t p, double wl )
{
  double xs = -1.0;
  double d[3];
  d[0] = dir[0]; d[1] = dir[1]; d[2] = dir[2];
  ncrystal_crosssection( p, ncrystal_wl2ekin( wl ),
                         (const double (*)[3])&d, &xs );
  return xs;
}

static void test_interface( void )
{
  api = (const ncrystal_vapi_type2_v1_t *)ncrystal_access_virtual_c_api( 2001 );
  require( api != NULL, "interface 2001 available" );
  require( api->interface_id == 2001, "interface_id" );
  require( api->struct_size == sizeof(ncrystal_vapi_type2_v1_t),
           "struct_size" );
  printf( "interface: OK\n" );
}

/* Info handles, compared with the normal C API: */
static void test_info( void )
{
  ncrystal_vapi_error_t err;
  unsigned i;
  for ( i = 0; i < NCFGS; ++i ) {
    ncrystal_vapi_t2v1_info_t * h;
    ncrystal_info_t ref;
    double t = -1.0;
    memset( &err, 0, sizeof(err) );
    h = api->create_info( cfgs[i], &err );
    require( h != NULL && err.code == 0, "create_info ok" );
    ref = ncrystal_create_info( cfgs[i] );
    require( api->info_density( h ) == ncrystal_info_getdensity( ref ),
             "density" );
    require( api->info_number_density( h )
             == ncrystal_info_getnumberdensity( ref ), "number density" );
    require( api->info_temperature( h, &t, &err ) == 0, "temperature ok" );
    require( t == ncrystal_info_gettemperature( ref ), "temperature" );
    printf( "info %s: density %g g/cm3, number density %g /Aa^3, T=%g K\n",
            cfgs[i], api->info_density( h ), api->info_number_density( h ),
            t );
    ncrystal_unref( &ref );
    api->deallocate_info( h );
  }
  api->deallocate_info( NULL );
  printf( "info: OK\n" );
}

/* Unique ids: the same for the same cfg-string (until caches are cleared): */
static void test_unique_id( void )
{
  ncrystal_vapi_error_t err;
  ncrystal_vapi_t2v1_info_t * a, * b, * c, * d;
  double ida, idb, idc, idd;
  a = api->create_info( "stdlib::Al_sg225.ncmat", &err );
  b = api->create_info( "stdlib::Al_sg225.ncmat", &err );
  c = api->create_info( "stdlib::Al_sg225.ncmat;temp=200K", &err );
  require( a && b && c, "create_info ok" );
  ida = api->info_unique_id( a );
  idb = api->info_unique_id( b );
  idc = api->info_unique_id( c );
  require( ida == idb, "same id for the same cfg-string" );
  require( ida != idc, "different id for different materials" );
  require( ida >= 0.0 && ida <= 9007199254740992.0
           && ida == (double)(uint64_t)ida,
           "id is an integral value in [0,2^53]" );
  /* Clearing the caches keeps existing handles valid, but new handles get */
  /* new ids:                                                              */
  api->clear_caches();
  require( api->info_unique_id( a ) == ida, "existing handle unchanged" );
  require( api->info_density( a ) > 2.0, "existing handle still works" );
  d = api->create_info( "stdlib::Al_sg225.ncmat", &err );
  require( d != NULL, "create_info after clear_caches ok" );
  idd = api->info_unique_id( d );
  require( idd != ida, "new id after clear_caches" );
  api->deallocate_info( a );
  api->deallocate_info( b );
  api->deallocate_info( c );
  api->deallocate_info( d );
  printf( "unique ids: OK\n" );
}

/* A natural abundance function, counting its calls via the state: */
static size_t natabund( void * state, unsigned long Z, unsigned long * A,
                        double * fraction, size_t capacity )
{
  require( capacity >= 128, "natabund capacity" );
  ++*(int*)state;
  if ( Z == 1 ) {
    A[0] = 1; fraction[0] = 0.99985;
    A[1] = 2; fraction[1] = 0.00015;
    return 2;
  }
  if ( Z == 5 ) {
    A[0] = 10; fraction[0] = 0.2;
    A[1] = 11; fraction[1] = 0.8;
    return 2;
  }
  if ( Z == 6 ) {
    A[0] = 12; fraction[0] = 0.99;
    A[1] = 13; fraction[1] = 0.01;
    return 2;
  }
  if ( Z == 8 ) {
    A[0] = 16; fraction[0] = 1.0;
    return 1;
  }
  return 0;
}

static void print_composition( const char * cfgstr, int prefer, int with_natab )
{
  ncrystal_vapi_error_t err;
  ncrystal_vapi_t2v1_component_t cmps[16];
  ncrystal_vapi_t2v1_info_t * h;
  size_t n, n2, i;
  double sum = 0.0;
  int ncalls = 0;
  h = api->create_info( cfgstr, &err );
  require( h != NULL, "create_info ok" );
  /* First the number of entries, then the entries: */
  memset( &err, 0, sizeof(err) );
  n = api->info_composition( h, prefer, with_natab ? natabund : NULL,
                             &ncalls, NULL, 0, &err );
  require( n > 0 && n <= 16 && err.code == 0, "composition count" );
  n2 = api->info_composition( h, prefer, with_natab ? natabund : NULL,
                              &ncalls, cmps, 16, &err );
  require( n2 == n, "composition entries" );
  printf( "composition of %s (prefer_natural_elements=%d%s):", cfgstr,
          prefer, with_natab ? ", with natabund" : "" );
  for ( i = 0; i < n; ++i ) {
    printf( " (%lu,%lu,%.6g)", cmps[i].Z, cmps[i].A, cmps[i].fraction );
    sum += cmps[i].fraction;
  }
  printf( "\n" );
  require( sum > 1.0 - 1e-12 && sum < 1.0 + 1e-12, "fractions sum to 1" );
  require( ( ncalls > 0 )
           == ( with_natab && ( !prefer || strstr( cfgstr, "phases" ) ) ),
           "natabund called when needed" );
  /* A too small capacity gives the total count, and only fills the first: */
  if ( n > 1 ) {
    ncrystal_vapi_t2v1_component_t one[1];
    require( api->info_composition( h, prefer, with_natab ? natabund : NULL,
                                    &ncalls, one, 1, &err ) == n,
             "small capacity" );
    require( one[0].Z == cmps[0].Z && one[0].A == cmps[0].A
             && one[0].fraction == cmps[0].fraction, "first entry" );
  }
  api->deallocate_info( h );
}

static void test_composition( void )
{
  /* Natural boron in one phase and enriched boron in the other, so the     */
  /* natural boron must be broken down into isotopes:                       */
  const char * mixed = "phases<0.5*stdlib::B4C_sg166_BoronCarbide.ncmat"
    "&0.5*solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11>";
  print_composition( "stdlib::Al_sg225.ncmat", 1, 0 );
  print_composition( "stdlib::Polyethylene_CH2.ncmat", 1, 0 );
  print_composition( "stdlib::Polyethylene_CH2.ncmat", 1, 1 );
  print_composition( "stdlib::Polyethylene_CH2.ncmat", 0, 1 );
  print_composition( "stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat", 1, 0 );
  print_composition( "stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat", 0, 1 );
  print_composition( "solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11", 1, 0 );
  print_composition( mixed, 1, 1 );
  printf( "composition: OK\n" );
}

/* Scatter and absorption handles, compared with the normal C API: */
static void test_processes( void )
{
  ncrystal_vapi_error_t err;
  unsigned i, j;
  for ( i = 0; i < NCFGS; ++i ) {
    ncrystal_vapi_t2v1_scatter_t * sc, * sc2;
    ncrystal_vapi_t2v1_absorption_t * ab, * ab2;
    ncrystal_scatter_t refsc = ncrystal_create_scatter( cfgs[i] );
    ncrystal_absorption_t refab = ncrystal_create_absorption( cfgs[i] );
    int oriented = !ncrystal_isnonoriented( ncrystal_cast_scat2proc( refsc ) );
    memset( &err, 0, sizeof(err) );
    sc = api->create_scatter( cfgs[i], &err );
    ab = api->create_absorption( cfgs[i], &err );
    require( sc && ab && err.code == 0, "create ok" );
    sc2 = api->clone_scatter( sc, &err );
    ab2 = api->clone_absorption( ab, &err );
    require( sc2 && ab2, "clone ok" );
    require( api->scatter_is_oriented( sc ) == oriented, "is_oriented" );
    require( api->scatter_is_oriented( sc2 ) == oriented, "is_oriented clone" );
    for ( j = 0; j < NWL; ++j ) {
      double n[4];
      double xs_sc = -1.0, xs_sc2 = -1.0, xs_ab = -1.0, xs_ab2 = -1.0;
      neutron( n, wls[j] );
      require( api->scatter_cross_section( sc, n, &xs_sc, &err ) == 0, "xs" );
      require( api->scatter_cross_section( sc2, n, &xs_sc2, &err ) == 0, "xs" );
      require( api->absorption_cross_section( ab, n, &xs_ab, &err ) == 0,
               "xs" );
      require( api->absorption_cross_section( ab2, n, &xs_ab2, &err ) == 0,
               "xs" );
      require( xs_sc == ref_xs( ncrystal_cast_scat2proc( refsc ), wls[j] ),
               "scatter xs as normal C API" );
      require( xs_ab == ref_xs( ncrystal_cast_abs2proc( refab ), wls[j] ),
               "absorption xs as normal C API" );
      require( xs_sc2 == xs_sc && xs_ab2 == xs_ab, "clones give the same xs" );
      if ( j == 2 )
        printf( "%s: oriented=%d, at 1.8Aa: scatter %.6g barn,"
                " absorption %.6g barn\n", cfgs[i], oriented, xs_sc, xs_ab );
    }
    api->deallocate_scatter( sc );
    api->deallocate_scatter( sc2 );
    api->deallocate_absorption( ab );
    api->deallocate_absorption( ab2 );
    ncrystal_unref( &refsc );
    ncrystal_unref( &refab );
  }
  api->deallocate_scatter( NULL );
  api->deallocate_absorption( NULL );
  printf( "scatter and absorption: OK\n" );
}

/* Sampling, compared with type 1 version 2 (same random numbers): */
static void test_sampling( void )
{
  const ncrystal_vapi_type1_v2_t * api1 =
    (const ncrystal_vapi_type1_v2_t *)ncrystal_access_virtual_c_api( 1002 );
  ncrystal_vapi_error_t err;
  unsigned i, k;
  require( api1 != NULL, "interface 1002 available" );
  for ( i = 0; i < NCFGS; ++i ) {
    ncrystal_vapi_t2v1_scatter_t * sc = api->create_scatter( cfgs[i], &err );
    ncrystal_vapi_t1v2_scatter_t * sc1 = api1->create_scatter( cfgs[i], &err );
    rng_t r, r1;
    require( sc && sc1, "create ok" );
    r.s = r1.s = 12345 + i;
    for ( k = 0; k < 100; ++k ) {
      double n[4], n1[4];
      neutron( n, wls[k % NWL] );
      neutron( n1, wls[k % NWL] );
      require( api->sample_scatter( sc, rng_fct, &r, n, &err ) == 0, "sample" );
      require( api1->sample_scatter( sc1, rng_fct, &r1, n1, &err ) == 0,
               "sample type 1" );
      require( memcmp( n, n1, sizeof(n) ) == 0, "same as type 1" );
    }
    api->deallocate_scatter( sc );
    api1->deallocate_scatter( sc1 );
  }
  printf( "sampling: OK\n" );
}

static void expect_error( int failed, const ncrystal_vapi_error_t * err,
                          const char * type, const char * what )
{
  require( failed && err->code != 0, what );
  require( strcmp( err->type, type ) == 0, what );
  printf( "  expected error (%s): [%s] %s\n", what, err->type, err->message );
}

static void test_errors( void )
{
  ncrystal_vapi_error_t err;
  ncrystal_vapi_t2v1_info_t * h;
  ncrystal_vapi_t2v1_absorption_t * ab;
  ncrystal_vapi_t2v1_component_t cmps[4];
  double t = -123.0, xs = -123.0;
  double n[4];
  printf( "errors:\n" );

  memset( &err, 0, sizeof(err) );
  expect_error( api->create_info( "stdlib::doesnotexist.ncmat", &err ) == NULL,
                &err, "FileNotFound", "missing file" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->create_absorption( NULL, &err ) == NULL, &err,
                "BadInput", "null cfgstr" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->clone_absorption( NULL, &err ) == NULL, &err,
                "BadInput", "clone of null handle" );

  /* Phases with different temperatures give no temperature: */
  h = api->create_info( "phases<0.5*stdlib::Al_sg225.ncmat;temp=200K"
                        "&0.5*stdlib::Cu_sg225.ncmat;temp=300K>", &err );
  require( h != NULL, "create_info ok" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->info_temperature( h, &t, &err ) != 0, &err,
                "MissingInfo", "no temperature" );
  require( t == -123.0, "temperature unchanged after error" );
  api->deallocate_info( h );

  /* Breaking down natural elements needs abundances: */
  h = api->create_info( "stdlib::Polyethylene_CH2.ncmat", &err );
  require( h != NULL, "create_info ok" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->info_composition( h, 0, NULL, NULL, cmps, 4, &err ) == 0,
                &err, "CalcError", "composition without natabund" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->info_composition( h, 1, NULL, NULL, NULL, 4, &err ) == 0,
                &err, "BadInput", "composition with null array" );
  api->deallocate_info( h );

  /* Invalid neutron for absorption: */
  ab = api->create_absorption( "stdlib::Al_sg225.ncmat", &err );
  require( ab != NULL, "create_absorption ok" );
  n[0] = -1.0; n[1] = 0.0; n[2] = 0.0; n[3] = 1.0;
  memset( &err, 0, sizeof(err) );
  expect_error( api->absorption_cross_section( ab, n, &xs, &err ) != 0, &err,
                "BadInput", "absorption with invalid neutron" );
  require( xs == -123.0, "xs unchanged after error" );
  api->deallocate_absorption( ab );
  printf( "errors: OK\n" );
}

int main( void )
{
  test_interface();
  test_info();
  test_unique_id();
  test_composition();
  test_processes();
  test_sampling();
  test_errors();
  return 0;
}
