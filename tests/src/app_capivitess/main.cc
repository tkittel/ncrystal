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


// Small application exercising the C API in a manner similar to how the VITESS
// sample_ncrystal module does, so changes affecting such usage are noticed:
// own RNG and message/error handlers, material info (state of matter,
// structure, density, composition) and per-neutron cross sections and
// scatterings.

#include "NCrystal/ncrystal.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

namespace {

  constexpr double k_pi = 3.14159265358979323846;

  //Deterministic RNG standing in for the application's own generator:
  unsigned long long s_rngstate = 12345;
  double my_randgen()
  {
    s_rngstate = s_rngstate * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>( s_rngstate >> 11 ) * ( 1.0 / 9007199254740992.0 );
  }

  void my_msghandler( const char * msg, unsigned msgtype )
  {
    std::printf( "%s: %s\n", ( msgtype == 1 ? "Warning"
                               : ( msgtype == 0 ? "Info" : "Raw" ) ), msg );
  }

  //Errors end the application (here with status 0, since the test provokes an
  //error at the very end):
  void my_errhandler( char * errtype, char * errmsg )
  {
    std::printf( "ERROR [%s]: %s\n", errtype, errmsg );
    std::fflush( stdout );
    std::exit( 0 );
  }

  void show_info( const char * cfgstr )
  {
    ncrystal_info_t info = ncrystal_create_info( cfgstr );
    static const char * states[] = { "Unknown", "Solid", "Gas", "Liquid" };
    int state = ncrystal_info_getstateofmatter( info );
    std::printf( "  State of matter: %s\n", states[ state >= 0 && state <= 3 ? state : 0 ] );
    ncrystal_cqueryresult_t res;
    if ( ncrystal_info_cquery( info, "spacegroup", &res ) && res.n_dbl == 1 )
      std::printf( "  Crystal (space group: %g)\n", res.val_dbl[0] );
    else
      std::printf( "  No crystal structure\n" );
    std::printf( "  Number density: %.6g [1/Aa^3]\n",
                 ncrystal_info_getnumberdensity( info ) );
    std::printf( "  Density: %.6g [g/cm3]\n", ncrystal_info_getdensity( info ) );
    double sigma_coh = 0.0, sigma_inc = 0.0, sigma_abs = 0.0;
    const unsigned ncomp = ncrystal_info_ncomponents( info );
    for ( unsigned i = 0; i < ncomp; ++i ) {
      unsigned atomdata_index;
      double fraction;
      ncrystal_info_getcomponent( info, i, &atomdata_index, &fraction );
      ncrystal_atomdata_t atomdata = ncrystal_create_component_atomdata( info, i );
      const char * label;
      const char * description;
      double mass, incxs, cohsl, absxs;
      unsigned nsubcomp, z, a;
      ncrystal_atomdata_getfields( atomdata, &label, &description, &mass,
                                   &incxs, &cohsl, &absxs, &nsubcomp, &z, &a );
      std::printf( "  Component %u: %s (Z=%u, A=%u) fraction=%.6g mass=%.6g\n",
                   i, ( label && label[0] ? label : description ),
                   z, a, fraction, mass );
      const double cohxs = 0.04 * k_pi * cohsl * cohsl;//4*pi*b^2, fm^2 -> barn
      sigma_coh += fraction * cohxs;
      sigma_inc += fraction * incxs;
      sigma_abs += fraction * absxs;
      ncrystal_unref( &atomdata );
    }
    std::printf( "  Bound coh/inc/abs xs per atom: %.6g/%.6g/%.6g [barn]\n",
                 sigma_coh, sigma_inc, sigma_abs );
    ncrystal_unref( &info );
  }

  void scatter_neutrons( const char * cfgstr )
  {
    ncrystal_scatter_t scat = ncrystal_create_scatter( cfgstr );
    ncrystal_process_t proc = ncrystal_cast_scat2proc( scat );
    const double wavelengths[] = { 0.5, 1.8, 4.0, 10.0 };
    for ( double wl : wavelengths ) {
      const double ekin = ncrystal_wl2ekin( wl );
      const double dir[3] = { 0.0, 0.6, 0.8 };
      double xs;
      ncrystal_crosssection( proc, ekin, &dir, &xs );
      //Sample a number of scatterings, checking basic sanity:
      const int n = 1000;
      bool ok = true;
      for ( int i = 0; i < n; ++i ) {
        double ekin_fin;
        double dir_fin[3];
        ncrystal_samplescatter( scat, ekin, &dir, &ekin_fin, &dir_fin );
        const double dirmag2 = ( dir_fin[0]*dir_fin[0] + dir_fin[1]*dir_fin[1]
                                 + dir_fin[2]*dir_fin[2] );
        if ( !( ekin_fin > 0.0 && ncrystal_ekin2wl( ekin_fin ) > 0.0 )
             || std::fabs( dirmag2 - 1.0 ) > 1e-9 )
          ok = false;
      }
      std::printf( "  wl=%g Aa: xs=%.5g barn, %d scatterings %s\n",
                   wl, xs, n, ( ok ? "OK" : "FAILED" ) );
    }
    ncrystal_unref( &scat );
  }

}

int main()
{
  ncrystal_setrandgen( my_randgen );//before creating scatter objects
  ncrystal_setmsghandler( my_msghandler );
  ncrystal_setquietonerror( 1 );
  ncrystal_seterrhandler( my_errhandler );

  const char * cfg_sc = ( "stdlib::Ge_sg227.ncmat;dcutoff=0.5;mos=40.0arcsec"
                          ";dir1=@crys_hkl:5,1,1@lab:0,0,1"
                          ";dir2=@crys_hkl:0,-1,1@lab:0,1,0" );
  const char * cfgs[] = {
    cfg_sc,
    "stdlib::Al_sg225.ncmat;temp=200K",
    "gasmix::He/1bar",
    "stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat"
  };
  for ( const char * cfg : cfgs ) {
    std::printf( "Material %s:\n", cfg );
    show_info( cfg );
    scatter_neutrons( cfg );
  }

  std::printf( "Material with invalid cfg-string:\n" );
  show_info( "stdlib::Al_sg225.ncmat;temp=-5K" );
  std::printf( "Should not get here\n" );
  return 1;
}
