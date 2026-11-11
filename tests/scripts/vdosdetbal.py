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

#Detailed balance, S(alpha,+beta)=exp(-beta)*S(alpha,-beta), of kernels
#built from a VDOS: the Keinert-Mattes H in ZrH VDOS, and a Debye model. Legacy
#beta grids are symmetric, allowing an exact check at all mirror nodes. Other
#beta grids are not, so there the mirror value is interpolated (linear in
#log(S)), except for the mirror nodes they do have (for the Debye model, also
#at its cutoff energy, the edge of the G1 range).

# NEEDS: numpy

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import numpy as np

VDOS = [0.00154, 0.00617, 0.01387, 0.02467, 0.03854, 0.05550, 0.07554, 0.09866, 0.12490, 0.15420,
        0.18650, 0.22200, 0.26050, 0.30220, 0.34690, 0.39470, 0.44550, 0.49950, 0.55650, 0.61660,
        *[0.0] * 46,
        0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00004, 0.00007, 0.00011, 0.00014, 0.00017,
        0.00020, 0.00023, 0.00026, 0.00067, 0.00144, 0.00221, 0.00298, 0.00376, 0.00453, 0.00530,
        0.00608, 0.00685, 0.01793, 0.02900, 0.04008, 0.05115, 0.06223, 0.07330, 0.08438, 0.09545,
        0.14760, 0.24070, 0.33390, 0.42700, 0.52020, 0.61330, 0.70650, 0.79960, 0.89280, 1.11600,
        1.38400, 1.70500, 2.08600, 2.53400, 3.05700, 3.66100, 4.35300, 5.14000, 6.02700, 7.01700,
        8.11200, 9.31100, 10.61000, 12.01000, 13.50000, 15.06000, 16.69000, 18.36000, 20.06000, 21.76000,
        23.44000, 25.07000, 26.62000, 28.07000, 29.39000, 30.55000, 31.54000, 32.33000, 32.91000, 33.26000,
        33.38000, 33.26000, 32.91000, 32.33000, 31.54000, 30.55000, 29.39000, 28.07000, 26.62000, 25.06000,
        23.44000, 21.76000, 20.06000, 18.36000, 16.69000, 15.06000, 13.50000, 12.01000, 10.61000, 9.31100,
        8.11200, 7.01700, 6.02700, 5.14000, 4.35300, 3.66100, 3.05700, 2.53400, 2.08600, 1.70500,
        1.38400, 1.11600, 0.89280, 0.72630, 0.55980, 0.45050, 0.34120, 0.27170, 0.20210, 0.15930,
        0.11640, 0.09082, 0.06519, 0.00000]          # 0.001 ... 0.180 eV, step 0.001 eV
NCMAT = ( 'NCMAT v7\n@DENSITY\n  1.0 g_per_cm3\n@DYNINFO\n  element H\n'
          '  fraction 1\n  type vdos\n  vdos_egrid 0.001 0.180\n'
          '  vdos_density ' + ' '.join(f'{v:g}' for v in VDOS) + '\n' )

_debye = """NCMAT v3
@CELL
 lengths 4.04958 4.04958 4.04958
 angles 90 90 90
@SPACEGROUP
 225
@ATOMPOSITIONS
 Al 0 0 0
 Al 0 1/2 1/2
 Al 1/2 0 1/2
 Al 1/2 1/2 0
@DEBYETEMPERATURE
 Al 410.4
"""

def load_kernel( temp, vdoslux, fn = 'h_zrhvdos.ncmat' ):
    info = NC.createInfo( f'{fn};temp={temp}K' )
    k = info.dyninfos[0].loadKernel( vdoslux = vdoslux )
    a, b = np.asarray( k['alpha'] ), np.asarray( k['beta'] )
    return a, b, np.asarray( k['sab'] ).reshape( len(b), len(a) ).T

def lnratio( x, y ):
    #log(x/y) where both are positive, nan elsewhere:
    ok = ( x > 0 ) & ( y > 0 )
    res = np.full( x.shape, np.nan )
    res[ok] = np.log( x[ok] ) - np.log( y[ok] )
    return res

def check_exact( temp, vdoslux = 1003, fn = 'h_zrhvdos.ncmat' ):
    _, b, S = load_kernel( temp, vdoslux, fn )
    nmirror, dmax = 0, 0.0
    for j in np.where( b > 0 )[0]:
        m = np.where( np.abs( b + b[j] ) <= 1e-12 * b[j] )[0]
        if not len(m):
            continue
        nmirror += 1
        d = lnratio( S[:,j], S[:,m[0]] * np.exp( -b[j] ) )
        d = d[ np.isfinite( d ) ]
        if len(d):
            dmax = max( dmax, float( np.abs( d ).max() ) )
    assert nmirror > 10, nmirror
    assert dmax < 1e-12, dmax
    print(f'{fn} vdoslux={vdoslux} T={temp}K: detailed balance exact at all'
          f' {nmirror} mirror nodes')

def check_interpolated( temp ):
    a, b, S = load_kernel( temp, 3 )
    bl = -b[ b < 0 ][::-1]
    Sl = S[ :, b < 0 ][ :, ::-1 ]
    sel = a >= 100
    ds = []
    for j in np.where( ( b >= 5 ) & ( b < bl[-1] ) )[0]:
        m = int( np.searchsorted( bl, b[j] ) )
        t = ( b[j] - bl[m-1] ) / ( bl[m] - bl[m-1] )
        s0, s1 = Sl[sel,m-1], Sl[sel,m]
        ok = ( s0 > 0 ) & ( s1 > 0 ) & ( S[sel,j] > 0 )
        lnmirror = ( 1 - t ) * np.log( s0[ok] ) + t * np.log( s1[ok] )
        ds.append( np.log( S[sel,j][ok] ) - ( lnmirror - b[j] ) )
    d = np.concatenate( ds )
    assert len(d) > 1000
    assert np.median( np.abs( d ) ) < 0.01, np.median( np.abs( d ) )
    print(f'vdoslux=3 T={temp}K: detailed balance holds at alpha>=100,'
          ' beta>=5 (median over interpolated mirrors)')

def main():
    NC.registerInMemoryFileData( 'h_zrhvdos.ncmat', NCMAT )
    NC.registerInMemoryFileData( 'al_debye.ncmat', _debye )
    for temp in ( 293.6, 1200.0 ):
        check_exact( temp )
        check_interpolated( temp )
    check_exact( 200.0, 3, 'al_debye.ncmat' )

if __name__ == '__main__':
    main()
