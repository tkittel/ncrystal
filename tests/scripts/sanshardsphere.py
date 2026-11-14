#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

# NEEDS: numpy mpmath

#Tests the hard-sphere SANS model (SANSSphereScatter): cross sections and
#sampled scattering angles against independent high precision evaluations, and
#the tabulation of the cross sections (IsoMacroXSTable).

#Before enable_fpe (mpmath 1.3 computes inf-inf when imported):
import mpmath # noqa F401
import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.xstable import IsoMacroXSTable
from NCrystalDev.misc import evaluate_query as ncquery
from NCTestUtils.mpmathctx import get_mpmath_context
from NCTestUtils.stat import kolmogorov_smirnov_pvalue
import numpy as np

R = 10.0#sphere radius (Aa)

NC.registerInMemoryFileData( 'Al_sans.ncmat', f"""NCMAT v8
@CELL
  cubic 4.04958
@SPACEGROUP
  225
@ATOMSITES
  Al 0 0 0
@OTHERPHASES
  0.1 void.ncmat
@CUSTOM_HARDSPHERESANS
  {R}
@DYNINFO
  element Al
  fraction 1
  type vdosdebye
  debye_temp 410
""" )

mp = get_mpmath_context( 40 )

def xs_shape( a ):
    #sigma(k) is proportional to this, with a=R*k:
    a = mp.mpf( a )
    v = ( 32*a**4 - 8*a**2 + 4*a*mp.sin(4*a) + mp.cos(4*a) - 1 ) / a**6
    return float( v )

def qr_cdf_unnormalised( x ):
    #Integral of [sin(t)-t*cos(t)]^2/t^5 from 0 to x (vectorised; the closed
    #form cancels badly for small x, where the series is used instead):
    x = np.asarray( x, dtype = float )
    xs = np.maximum( x, 1e-2 )
    closed = ( ( 2*xs**4 - 2*xs**2 + 2*xs*np.sin(2*xs) + np.cos(2*xs) - 1 )
               / ( 8*xs**4 ) )
    series = x**2 / 18 - x**4 / 180
    return np.where( x < 1e-2, series, closed )

def test_xs():
    sc = NC.createScatter( 'Al_sans.ncmat;comp=sans' )
    sc.dump()
    #Both the Taylor expansion (R*k < sqrt(0.9)) and the direct formula:
    wl = np.geomspace( 0.1, 1e4, 400 )
    xs = sc.xsect( wl = wl )
    ratio = np.array( [ x / xs_shape( R * 2*np.pi / w ) for x, w in zip( xs, wl ) ] )
    rel = np.abs( ratio / np.median( ratio ) - 1.0 )
    assert rel.max() < 1e-13, rel.max()
    #Limit for long wavelengths (shape -> 256/9):
    assert abs( sc.xsect( wl = 1e9 ) / ( np.median( ratio ) * 256/9 ) - 1.0 ) < 1e-13
    print('cross sections at 400 wavelengths OK (max rel. deviation < 1e-13)')

def test_sampling():
    sc = NC.createScatter( 'Al_sans.ncmat;comp=sans' )
    rng = np.random.default_rng( 123 )
    for wl in ( 2.0, 20.0, 100.0, 1000.0 ):
        ekin, mu = sc.sampleScatterIsotropic( ekin = NC.wl2ekin( wl ), repeat = 20000 )
        assert np.all( ekin == NC.wl2ekin( wl ) )#elastic
        k = 2*np.pi / wl
        x = R * k * np.sqrt( np.clip( 2.0 * ( 1.0 - mu ), 0.0, None ) )
        xmax = 2 * R * k
        u = qr_cdf_unnormalised( x ) / qr_cdf_unnormalised( xmax )
        p = kolmogorov_smirnov_pvalue( u, rng.uniform( size = 100000 ) )
        assert p > 1e-3, (wl, p)
        print(f'sampling at {wl:g} Aa OK (KS test on Q*R)')

def test_tabulation( cfg ):
    r = ncquery( ['isomacroxstable', cfg], huge_arrays = True )
    tol = r['tol']
    f = IsoMacroXSTable( cfg )
    wl = np.concatenate( [ np.geomspace( 1e-3, 1e5, 20001 ),
                           np.random.default_rng( 1 ).uniform( 1.0, 200.0, 5000 ) ] )
    exact = ( NC.createScatter(cfg).xsect( wl = wl ) + NC.createAbsorption(cfg).xsect( wl = wl )
              ) * NC.createInfo(cfg).numberdensity
    #The tolerance is guaranteed up to wlmax, beyond the table it is
    #extrapolated:
    m = wl <= r['wlmax']
    tab = f.xsect( wl = wl )
    reldev = np.abs( tab[m] - exact[m] ) / exact[m]
    assert reldev.max() <= tol, reldev.max()
    print(f'{cfg}: {len(f.table[0])} table points, within the tolerance {tol:g}'
          f' up to {r["wlmax"]:g} Aa')

def main():
    test_xs()
    test_sampling()
    test_tabulation( 'Al_sans.ncmat;comp=sans' )
    test_tabulation( 'Al_sans.ncmat' )
    test_tabulation( 'Al_sans.ncmat;knllux=-1' )

if __name__ == '__main__':
    main()
