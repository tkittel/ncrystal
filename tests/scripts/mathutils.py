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

# NEEDS: numpy

# Test the "mathutils" query and the query input array mechanism
# (ncrystal_jsonquery_inputarray) used to pass numpy arrays to it.

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev.misc import evaluate_query as ncquery
from NCrystalDev.exceptions import NCBadInput
import numpy as np
import math

def expect_badinput( *query ):
    try:
        ncquery( list(query) )
    except NCBadInput as e:
        print(f'  Expected error: {e}')
        return
    raise SystemExit('Did not emit exception as required!')

def fourierint( fr, fi, w, *opts ):
    return ncquery( ['mathutils','fourierint',fr,
                     'none' if fi is None else fi, w] + list(opts),
                    huge_arrays = True )['F']

def test_gaussian():
    #f(t)=exp(-g t^2/2 + i w0 t) transforms to a Gaussian centred at w0:
    g, w0, h = 0.7, 1.3, 0.05
    t = np.arange( 0.0, 40.0, h )
    env = np.exp( -0.5 * g * t * t )
    fr, fi = env * np.cos( w0 * t ), env * np.sin( w0 * t )
    w = np.linspace( -5.0, 8.0, 131 )
    exact = np.exp( -(w-w0)**2 / ( 2 * g ) ) / math.sqrt( 2 * math.pi * g )
    sel = exact > 1e-8 * exact.max()
    print('Max relative errors for Gaussian test:')
    for scheme in ('trap','trapec','lin','cubic'):
        for ev in ('direct','fftlin','fftcubic'):
            F = fourierint( fr, fi, w, f'h={h}',
                            f'scheme={scheme}', f'eval={ev}' )
            assert F.shape == w.shape
            relerr = np.abs( F[sel]/exact[sel] - 1.0 ).max()
            print(f'  {scheme:6} {ev:8} : %s'%('<1e-5' if relerr < 1e-5
                                                else '<0.01' if relerr < 0.01
                                                else 'LARGE'))
    #Non-contiguous arrays are accepted:
    assert np.array_equal( fourierint( np.repeat(env,2)[::2], None, w, f'h={h}' ),
                           fourierint( env, None, w, f'h={h}' ) )
    #Real even input with fi=None is the same as fi=0:
    a = fourierint( env, None, w, f'h={h}' )
    b = fourierint( env, np.zeros_like(env), w, f'h={h}' )
    assert np.array_equal( a, b )

def test_errors():
    print('Errors:')
    arr = np.linspace( 1.0, 0.0, 5 )
    w = np.linspace( 0.0, 1.0, 3 )
    q = ['mathutils','fourierint']
    expect_badinput( *q, arr, 'none', 'junk', 'h=1' )
    expect_badinput( *q )
    expect_badinput( *q, arr, 'none', w )
    expect_badinput( *q, arr, 'none', w, 'h=1', 'scheme=bla' )
    expect_badinput( *q, arr, arr[:-1], w, 'h=1' )
    expect_badinput( 'mathutils', 'bla' )
    expect_badinput( *q, np.ones((2,2)), 'none', w, 'h=1' )
    expect_badinput( *q, 117.0, 'none', w, 'h=1' )
    expect_badinput( *q, [1.0,2.0], 'none', w, 'h=1' )
    expect_badinput( *q, (1.0,2.0), 'none', w, 'h=1' )
    expect_badinput( *q, np.arange(5,dtype=np.int32), 'none', w, 'h=1' )
    expect_badinput( *q, arr.astype(np.float32), 'none', w, 'h=1' )

def test_placeholders():
    print('Placeholders:')
    from NCrystalDev._chooks import _get_raw_cfcts
    reg = lambda a : _get_raw_cfcts()['jsonquery_inputarrays']([a])[0]
    arr = np.linspace( 1.0, 0.0, 5 )
    w = np.linspace( 0.0, 1.0, 3 )
    q = ['mathutils','fourierint']
    p = reg( arr )
    assert p.startswith('__ncrystal__dblinarray::') and p.endswith('::5')
    #Input arrays can only be used once, and are discarded after each query
    #(even if unused or in case of errors):
    pw = reg( w )
    ncquery( q + [p, 'none', pw, 'h=1'] )
    expect_badinput( *q, p, 'none', w, 'h=1' )
    expect_badinput( *q, arr, 'none', pw, 'h=1' )
    p, pw = reg( arr ), reg( w )
    expect_badinput( *q, p, 'none', pw )
    expect_badinput( *q, p, 'none', w, 'h=1' )
    #Tampered size:
    p = reg( arr )
    expect_badinput( *q, p[:-1]+'4', 'none', w, 'h=1' )
    #Empty arrays are allowed:
    p = reg( np.zeros(0) )
    assert p.endswith('::0')
    F = ncquery( q + [arr, 'none', p, 'h=1'], huge_arrays=True )['F']
    assert F.shape == (0,)

def main():
    test_gaussian()
    test_errors()
    test_placeholders()
    print('All OK')

if __name__ == '__main__':
    main()
