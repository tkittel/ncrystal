
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

/* Tests of the virtual C API (NCrystal/virtualapi/ncvirtapi.h), comparing    */
/* with the results of the normal C API (which uses the same physics).        */

#include "NCrystal/ncrystal.h"
#include "NCrystal/virtualapi/ncvirtapi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#ifndef _WIN32
#  include <unistd.h>
#  include <sys/wait.h>
#endif

static void require( int cond, const char * what )
{
  if ( !cond ) {
    printf( "FAILED: %s\n", what );
    exit( 1 );
  }
}

/* Simple deterministic RNG with numbers in [0,1): */
typedef struct { uint64_t s; int nzero; } rng_t;

static double rng_fct( void * p )
{
  rng_t * r = (rng_t*)p;
  if ( r->nzero > 0 ) {
    --r->nzero;
    return 0.0;
  }
  r->s = r->s * 6364136223846793005ULL + 1442695040888963407ULL;
  return (double)( r->s >> 11 ) * ( 1.0 / 9007199254740992.0 );
}

static double rng_toolarge( void * p ) { (void)p; return 1.5; }
static double rng_nan( void * p ) { (void)p; return NAN; }

static const ncrystal_vapi_type1_v2_t * api = NULL;

#define NCFGS 5
static const char * cfgs[NCFGS] = {
  "stdlib::Al_sg225.ncmat",
  "stdlib::Polyethylene_CH2.ncmat;temp=250K",
  "stdlib::Cu_sg225.ncmat;dcutoff=0.5",
  "gasmix::He/2bar",
  "stdlib::Ge_sg227.ncmat;mos=40arcmin;dir1=@crys_hkl:5,1,1@lab:0,0,1"
  ";dir2=@crys_hkl:0,-1,1@lab:0,1,0"
};

#define NWL 9
static const double wls[NWL] = { 0.5, 1.0, 1.8, 2.3, 3.0, 4.0, 4.6, 7.0, 20.0 };
#define NDIRS 3
static const double dirs[NDIRS][3] = { { 0.0, 0.0, 1.0 },
                                       { 0.0, 0.6, 0.8 },
                                       { 0.48, 0.6, 0.64 } };

static double ref_xs( ncrystal_scatter_t sc, double ekin, const double * dir )
{
  double xs = -1.0;
  double d[3];
  d[0] = dir[0]; d[1] = dir[1]; d[2] = dir[2];
  ncrystal_crosssection( ncrystal_cast_scat2proc( sc ), ekin,
                         (const double (*)[3])&d, &xs );
  return xs;
}

static double api_xs( ncrystal_vapi_t1v2_scatter_t * h, double ekin,
                      const double * dir )
{
  ncrystal_vapi_error_t err;
  double n[4];
  double xs = -1.0;
  n[0] = ekin; n[1] = dir[0]; n[2] = dir[1]; n[3] = dir[2];
  require( api->cross_section( h, n, &xs, &err ) == 0, "cross_section ok" );
  return xs;
}

static ncrystal_vapi_t1v2_scatter_t * create( const char * cfgstr )
{
  ncrystal_vapi_error_t err;
  ncrystal_vapi_t1v2_scatter_t * h = api->create_scatter( cfgstr, &err );
  require( h != NULL, "create_scatter ok" );
  return h;
}

static void test_interface( void )
{
  const unsigned bad_ids[6] = { 0, 1, 1001, 1003, 2001, 4294967295u };
  unsigned i;
  api = (const ncrystal_vapi_type1_v2_t*)ncrystal_access_virtual_c_api( 1002 );
  require( api != NULL, "interface 1002 available" );
  require( api->interface_id == 1002, "interface_id" );
  require( api->struct_size == sizeof(ncrystal_vapi_type1_v2_t), "struct_size" );
  require( api->create_scatter && api->clone_scatter
           && api->deallocate_scatter && api->cross_section
           && api->sample_scatter, "all functions set" );
  require( ncrystal_access_virtual_c_api( 1002 ) == (const void*)api,
           "same struct every time" );
  for ( i = 0; i < 6; ++i )
    require( ncrystal_access_virtual_c_api( bad_ids[i] ) == NULL,
             "unknown interface id gives NULL" );
  printf( "interface: OK\n" );
}

/* Cross sections agree with the normal C API, also when calls for different */
/* materials are interleaved (which would give wrong results if caches were  */
/* shared between materials), and when repeated (cache hits):                */
static void test_xs( void )
{
  ncrystal_scatter_t refs[NCFGS];
  ncrystal_vapi_t1v2_scatter_t * hs[NCFGS];
  unsigned i, iwl, idir, irep;
  unsigned long ncalls = 0;
  for ( i = 0; i < NCFGS; ++i ) {
    refs[i] = ncrystal_create_scatter( cfgs[i] );
    hs[i] = create( cfgs[i] );
  }
  for ( irep = 0; irep < 2; ++irep )
    for ( iwl = 0; iwl < NWL; ++iwl )
      for ( idir = 0; idir < NDIRS; ++idir )
        for ( i = 0; i < NCFGS; ++i ) {
          double ekin = ncrystal_wl2ekin( wls[iwl] );
          double a = api_xs( hs[i], ekin, dirs[idir] );
          double b = ref_xs( refs[i], ekin, dirs[idir] );
          ++ncalls;
          if ( !( a == b ) ) {
            printf( "xs mismatch for %s at %g Aa: %.17g vs. %.17g\n",
                    cfgs[i], wls[iwl], a, b );
            require( 0, "xs agree" );
          }
        }
  for ( i = 0; i < NCFGS; ++i ) {
    printf( "xs(1.8 Aa, dir=(0,0,1)) = %.6g barn for %s\n",
            api_xs( hs[i], ncrystal_wl2ekin( 1.8 ), dirs[0] ), cfgs[i] );
    api->deallocate_scatter( hs[i] );
    ncrystal_unref( &refs[i] );
  }
  /* The oriented material depends on the direction: */
  {
    ncrystal_vapi_t1v2_scatter_t * h = create( cfgs[NCFGS-1] );
    double e = ncrystal_wl2ekin( 1.8 );
    require( api_xs( h, e, dirs[0] ) != api_xs( h, e, dirs[2] ),
             "oriented material depends on direction" );
    api->deallocate_scatter( h );
  }
  printf( "cross sections: %lu interleaved calls agree with the C API: OK\n",
          ncalls );
}

/* Clones share the physics, have their own caches, and outlive the original: */
static void test_clone( void )
{
  unsigned i, iwl;
  ncrystal_vapi_error_t err;
  for ( i = 0; i < NCFGS; ++i ) {
    ncrystal_vapi_t1v2_scatter_t * h = create( cfgs[i] );
    ncrystal_vapi_t1v2_scatter_t * c1 = api->clone_scatter( h, &err );
    ncrystal_vapi_t1v2_scatter_t * c2;
    require( c1 != NULL && c1 != h, "clone ok" );
    for ( iwl = 0; iwl < NWL; ++iwl ) {
      double e = ncrystal_wl2ekin( wls[iwl] );
      require( api_xs( h, e, dirs[1] ) == api_xs( c1, e, dirs[1] ),
               "clone gives same xs" );
    }
    api->deallocate_scatter( h );
    c2 = api->clone_scatter( c1, &err );
    require( c2 != NULL, "clone of clone ok" );
    for ( iwl = 0; iwl < NWL; ++iwl ) {
      double e = ncrystal_wl2ekin( wls[iwl] );
      require( api_xs( c1, e, dirs[2] ) == api_xs( c2, e, dirs[2] ),
               "clone of clone gives same xs" );
    }
    api->deallocate_scatter( c1 );
    api->deallocate_scatter( c2 );
  }
  api->deallocate_scatter( NULL );
  printf( "clones: OK\n" );
}

/* Sampling with the same random numbers gives exactly the same results as   */
/* the C API (ncrystal_samplescatter_rs), for the original and for clones:    */
static void test_sampling( void )
{
  unsigned i, iwl, isample;
  ncrystal_vapi_error_t err;
  for ( i = 0; i < NCFGS; ++i ) {
    ncrystal_scatter_t ref = ncrystal_create_scatter( cfgs[i] );
    ncrystal_vapi_t1v2_scatter_t * h = create( cfgs[i] );
    ncrystal_vapi_t1v2_scatter_t * c = api->clone_scatter( h, &err );
    double sum_ekin = 0.0;
    require( c != NULL, "clone ok" );
    for ( iwl = 0; iwl < NWL; ++iwl ) {
      for ( isample = 0; isample < 20; ++isample ) {
        rng_t r1, r2, r3;
        double n1[4], n2[4], ekin_ref;
        double dir_ref[3], dir_in[3];
        r1.s = r2.s = r3.s = 1000u * iwl + isample + 1;
        r1.nzero = r2.nzero = r3.nzero = 0;
        n1[0] = n2[0] = ncrystal_wl2ekin( wls[iwl] );
        n1[1] = n2[1] = dir_in[0] = dirs[2][0];
        n1[2] = n2[2] = dir_in[1] = dirs[2][1];
        n1[3] = n2[3] = dir_in[2] = dirs[2][2];
        require( api->sample_scatter( h, rng_fct, &r1, n1, &err ) == 0,
                 "sample_scatter ok" );
        require( api->sample_scatter( c, rng_fct, &r2, n2, &err ) == 0,
                 "sample_scatter with clone ok" );
        ncrystal_samplescatter_rs( rng_fct, &r3, ref, ncrystal_wl2ekin( wls[iwl] ),
                                   (const double (*)[3])&dir_in, &ekin_ref,
                                   (double (*)[3])&dir_ref );
        require( n1[0] == n2[0] && n1[1] == n2[1] && n1[2] == n2[2]
                 && n1[3] == n2[3], "clone samples identically" );
        require( n1[0] == ekin_ref && n1[1] == dir_ref[0]
                 && n1[2] == dir_ref[1] && n1[3] == dir_ref[2],
                 "same sampling as the C API" );
        require( n1[0] >= 0.0
                 && fabs( n1[1]*n1[1] + n1[2]*n1[2] + n1[3]*n1[3] - 1.0 ) < 1e-12,
                 "valid outcome" );
        sum_ekin += n1[0];
      }
    }
    printf( "sampling: mean final energy %.6g eV for %s\n",
            sum_ekin / ( NWL * 20 ), cfgs[i] );
    api->deallocate_scatter( h );
    api->deallocate_scatter( c );
    ncrystal_unref( &ref );
  }
  /* A random number of exactly 0 is accepted: */
  {
    ncrystal_vapi_t1v2_scatter_t * h = create( cfgs[1] );
    rng_t r;
    double n[4];
    r.s = 1; r.nzero = 3;
    n[0] = ncrystal_wl2ekin( 1.8 ); n[1] = 0.0; n[2] = 0.0; n[3] = 1.0;
    require( api->sample_scatter( h, rng_fct, &r, n, &err ) == 0,
             "random numbers of 0 accepted" );
    require( r.nzero == 0, "zeros were used" );
    api->deallocate_scatter( h );
  }
  printf( "sampling: OK\n" );
}

static void expect_error( int ret, const ncrystal_vapi_error_t * err,
                          const char * type, const char * what )
{
  require( ret != 0 && err->code != 0, what );
  if ( strcmp( err->type, type ) != 0 ) {
    printf( "Got error type %s, expected %s\n", err->type, type );
    require( 0, what );
  }
  printf( "  expected error (%s): [%s] %s\n", what, err->type, err->message );
}

static void test_errors( void )
{
  ncrystal_vapi_error_t err;
  ncrystal_vapi_t1v2_scatter_t * h = create( cfgs[0] );
  double n[4], n_orig[4], xs;
  rng_t r;
  char longcfg[3001];
  unsigned i;
  printf( "errors:\n" );

  memset( &err, 0, sizeof(err) );
  require( api->create_scatter( "stdlib::doesnotexist.ncmat", &err ) == NULL,
           "create_scatter fails" );
  expect_error( 1, &err, "FileNotFound", "missing file" );
  memset( &err, 0, sizeof(err) );
  require( api->create_scatter( NULL, &err ) == NULL, "null cfgstr" );
  expect_error( 1, &err, "BadInput", "null cfgstr" );
  memset( &err, 0, sizeof(err) );
  require( api->create_scatter( "stdlib::Al_sg225.ncmat;temp=-5K", &err )
           == NULL, "bad parameter" );
  require( err.code != 0, "bad parameter" );
  printf( "  expected error (bad parameter): [%s] %s\n", err.type, err.message );
  memset( &err, 0, sizeof(err) );
  require( api->clone_scatter( NULL, &err ) == NULL, "null handle" );
  expect_error( 1, &err, "BadInput", "clone of null handle" );

  /* Invalid neutron states (the outputs must be unchanged): */
  {
    const double bad[6][4] = { { -1.0, 0.0, 0.0, 1.0 },
                               { NAN, 0.0, 0.0, 1.0 },
                               { INFINITY, 0.0, 0.0, 1.0 },
                               { 0.025, 0.0, 0.0, 0.0 },
                               { 0.025, NAN, 0.0, 1.0 },
                               { 0.025, INFINITY, 0.0, 1.0 } };
    for ( i = 0; i < 6; ++i ) {
      memcpy( n, bad[i], sizeof(n) );
      xs = -123.0;
      memset( &err, 0, sizeof(err) );
      expect_error( api->cross_section( h, n, &xs, &err ), &err, "BadInput",
                    "cross_section with invalid neutron" );
      require( xs == -123.0, "xs unchanged after error" );
      memset( &err, 0, sizeof(err) );
      r.s = 1; r.nzero = 0;
      expect_error( api->sample_scatter( h, rng_fct, &r, n, &err ), &err,
                    "BadInput", "sample_scatter with invalid neutron" );
      require( memcmp( n, bad[i], sizeof(n) ) == 0, "neutron unchanged" );
    }
  }

  /* Null arguments and bad random numbers: */
  n[0] = n_orig[0] = 0.025; n[1] = n_orig[1] = 0.0;
  n[2] = n_orig[2] = 0.0; n[3] = n_orig[3] = 1.0;
  memset( &err, 0, sizeof(err) );
  expect_error( api->cross_section( NULL, n, &xs, &err ), &err, "BadInput",
                "cross_section with null handle" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->cross_section( h, NULL, &xs, &err ), &err, "BadInput",
                "cross_section with null neutron" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->cross_section( h, n, NULL, &err ), &err, "BadInput",
                "cross_section with null result" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->sample_scatter( h, NULL, NULL, n, &err ), &err,
                "BadInput", "sample_scatter with null rng" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->sample_scatter( h, rng_toolarge, NULL, n, &err ), &err,
                "BadInput", "rng returning 1.5" );
  memset( &err, 0, sizeof(err) );
  expect_error( api->sample_scatter( h, rng_nan, NULL, n, &err ), &err,
                "BadInput", "rng returning NaN" );
  require( memcmp( n, n_orig, sizeof(n) ) == 0, "neutron unchanged" );

  /* Long messages are truncated, and always null-terminated: */
  memset( longcfg, 'x', 3000 );
  longcfg[3000] = '\0';
  memset( &err, 'z', sizeof(err) );
  require( api->create_scatter( longcfg, &err ) == NULL, "long cfgstr" );
  require( err.code != 0, "long cfgstr error" );
  require( strlen( err.message ) == sizeof(err.message) - 1,
           "long message truncated" );
  require( strlen( err.type ) < sizeof(err.type), "type terminated" );
  printf( "  long message truncated to %u characters: OK\n",
          (unsigned)strlen( err.message ) );

  /* The handle still works after errors: */
  require( api->cross_section( h, n, &xs, &err ) == 0 && xs > 0.0,
           "handle works after errors" );
  api->deallocate_scatter( h );
  printf( "errors: OK\n" );
}

/* With a NULL error pointer, the error is printed and the program exits: */
static void test_null_err( void )
{
#ifndef _WIN32
  int fds[2];
  pid_t pid;
  int status = 0;
  char buf[512];
  ssize_t nread;
  require( pipe( fds ) == 0, "pipe" );
  fflush( stdout );
  pid = fork();
  require( pid >= 0, "fork" );
  if ( pid == 0 ) {
    close( fds[0] );
    dup2( fds[1], 2 );
    api->create_scatter( "stdlib::doesnotexist.ncmat", NULL );
    _exit( 0 );/* not reached */
  }
  close( fds[1] );
  nread = read( fds[0], buf, sizeof(buf) - 1 );
  close( fds[0] );
  require( nread > 0, "error message printed" );
  buf[nread] = '\0';
  if ( strchr( buf, '\n' ) )
    *strchr( buf, '\n' ) = '\0';
  require( waitpid( pid, &status, 0 ) == pid, "waitpid" );
  require( WIFEXITED( status ) && WEXITSTATUS( status ) == 1, "exit code 1" );
  printf( "null error pointer: printed \"%s\" and exited with code 1: OK\n",
          buf );
#else
  printf( "null error pointer: OK\n" );
#endif
}

int main( void )
{
  test_interface();
  test_xs();
  test_clone();
  test_sampling();
  test_errors();
  test_null_err();
  return 0;
}
