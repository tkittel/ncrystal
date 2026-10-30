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

def test_tail():
    #Geometric tail on truncated exp(-g t + i w0 t) reproduces the long-range
    #trapezoid result (identical up to rounding):
    g, w0, h = 0.05, 2.0, 0.1
    w = np.linspace( -5.0, 9.0, 141 )
    def f( tmax ):
        t = np.arange( 0.0, tmax + 0.5 * h, h )
        env = np.exp( -g * t )
        return env * np.cos( w0 * t ), env * np.sin( w0 * t )
    a = fourierint( *f(14.0), w, f'h={h}', 'scheme=trap', 'tail=geometric' )
    b = fourierint( *f(1000.0), w, f'h={h}', 'scheme=trap' )
    assert np.abs( a - b ).max() < 1e-12 * np.abs( b ).max()
    #Fully decayed (zero last value) means no tail:
    z = np.zeros(3)
    fr0, fi0 = f(14.0)
    fr0, fi0 = np.concatenate([fr0,z]), np.concatenate([fi0,z])
    assert np.array_equal( fourierint( fr0, fi0, w, f'h={h}', 'scheme=trap',
                                       'tail=geometric' ),
                           fourierint( fr0, fi0, w, f'h={h}', 'scheme=trap' ) )
    c = fourierint( *f(14.0), w, f'h={h}', 'scheme=trap' )
    assert np.abs( c - b ).max() > 1e-2 * np.abs( b ).max()
    print('Geometric tail OK')

def convolve( xf, f, xg, g, z, *opts ):
    return ncquery( ['mathutils','convolve',xf,f,xg,g,z] + list(opts),
                    huge_arrays = True )['h']

def test_convolve():
    #Exact convolution of piecewise linear functions (with steps at the ends)
    #vs. a fine brute-force trapezoid evaluation (error ~1e-5 from the steps):
    xf = np.array([-1.0,0.0,0.5,2.0])
    f = np.array([0.3,1.0,-0.2,0.7])
    xg = np.array([-0.5,0.1,0.4])
    g = np.array([0.5,2.0,0.1])
    z = np.linspace( -3.0, 4.0, 71 )
    h = convolve( xf, f, xg, g, z, 'scheme=exactlin' )
    assert np.array_equal( h, convolve( xg, g, xf, f, z, 'scheme=exactlin' ) )
    xx = np.linspace( -5.0, 5.0, 1000001 )
    fi = np.interp( xx, xf, f, left = 0.0, right = 0.0 )
    ref = np.array( [ np.trapezoid( fi * np.interp( v - xx, xg, g, left = 0.0,
                                                    right = 0.0 ), xx )
                      for v in z ] )
    assert np.abs( h - ref ).max() < 3e-5
    assert h[0] == 0.0 and h[-1] == 0.0
    #Gaussians (sigma 1 and 0.5), all schemes close to the exact result:
    x1 = np.linspace( -8.0, 8.0, 321 )
    x2 = np.linspace( -4.0, 4.0, 161 )
    def gauss( x, s ):
        return np.exp( -0.5 * (x/s)**2 ) / ( s * math.sqrt( 2 * math.pi ) )
    z = np.linspace( -3.0, 3.0, 61 ) + 0.0123#not aligned with input nodes
    exact = gauss( z, math.sqrt(1.25) )
    print('Max relative errors for Gaussian convolution:')
    for scheme in ('direct','exactlin','fft'):
        for interp in ('lin','lagrange4','spline','pchip'):
            if scheme == 'exactlin' and interp != 'lin':
                continue
            h = convolve( x1, gauss(x1,1.0), x2, gauss(x2,0.5), z,
                          f'scheme={scheme}', f'interp={interp}' )
            relerr = np.abs( h / exact - 1.0 ).max()
            print(f'  {scheme:8} {interp:9} : %s'%(
                '<1e-6' if relerr < 1e-6 else '<1e-3' if relerr < 1e-3 else
                '<1e-2' if relerr < 1e-2 else 'LARGE'))

def test_overshoot():
    #x^2 on [-1,1] with sharp cutoffs between the nodes of a coarse grid,
    #convolved with a narrow Gaussian: spline interpolation rings into
    #negative values, pchip (and lin) do not:
    xf = np.linspace( -0.3, 0.3, 301 )
    f = np.exp( -0.5 * ( xf / 0.03 )**2 )
    xg = np.arange( -2.0, 2.05, 0.1 ) + 0.037
    g = np.where( np.abs( xg ) <= 1.0, xg * xg, 0.0 )
    z = np.linspace( -1.5, 1.5, 301 )
    print('Most negative value relative to max for sharp cutoffs:')
    for scheme in ('direct','fft'):
        for interp in ('lin','spline','pchip'):
            h = convolve( xf, f, xg, g, z, f'scheme={scheme}',
                          f'interp={interp}', 'step=0.006' )
            r = h.min() / h.max()
            print(f'  {scheme:6} {interp:6} : %s'%( '>=-1e-15' if r >= -1e-15
                                                    else '<-1e-3' if r < -1e-3
                                                    else 'small' ))

def test_lines():
    #Oscillator lines: exact sum rules and detailed balance, and weights vs.
    #the closed form for small A (w_{+-1} ~ A e^{+-betaE/2}/(2 sinh)):
    print('Oscillator lines:')
    for A, bE in ( (1e-6,1.0), (0.3,2.0), (5.0,0.5), (50.0,30.0), (2.0,317.0) ):
        r = ncquery( ['mathutils','osclines',f'A={A}',f'betaE={bE}'],
                     huge_arrays = True )
        s, w = r['shifts'], r['weights']
        assert np.all( np.diff( s ) > 0 ) and np.all( w > 0 )
        assert abs( w.sum() - 1.0 ) < 1e-12
        assert abs( ( w * s ).sum() / ( -A * bE ) - 1.0 ) < 1e-10
        k = np.rint( -s / bE ).astype( int )
        d = dict( zip( k, w ) )
        for kk in k:
            if kk > 0 and -kk in d and d[kk] > 1e-12 * w.max():
                assert abs( d[-kk] / d[kk] * math.exp( kk * bE ) - 1 ) < 1e-12
        print(f'  A={A} betaE={bE}: lines k={k.min()}..{k.max()} OK')
    r = ncquery( ['mathutils','osclines','A=1e-6','betaE=1'], huge_arrays = True )
    d = dict( zip( np.rint( -r['shifts'] ).astype( int ), r['weights'] ) )
    assert abs( d[1] / ( 1e-6 * math.exp( 0.5 ) / ( 2 * math.sinh( 0.5 ) ) )
                - 1 ) < 1e-5
    #Shift-and-add of a Gaussian: exact for shifted copies:
    x = np.linspace( -10.0, 10.0, 2001 )
    g = np.exp( -0.5 * x * x )
    sh = np.array( [ -2.0, 0.5, 3.0 ] )
    wt = np.array( [ 0.2, 0.5, 0.3 ] )
    z = np.linspace( -4.0, 4.0, 81 ) + 0.0123
    exact = sum( wk * np.exp( -0.5 * ( z - sk )**2 ) for sk, wk in zip( sh, wt ) )
    for interp in ( 'lin', 'spline', 'pchip' ):
        h = ncquery( ['mathutils','shiftadd',x,g,sh,wt,z,f'interp={interp}'],
                     huge_arrays = True )['h']
        relerr = np.abs( h / exact - 1 ).max()
        print(f'  shiftadd {interp:6}: %s'%( '<1e-7' if relerr < 1e-7 else
                                             '<1e-4' if relerr < 1e-4 else
                                             'LARGE' ))

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
    expect_badinput( *q, arr, 'none', w, 'h=1', 'scheme=lin', 'tail=geometric' )
    expect_badinput( *q, arr[::-1].copy(), 'none', w, 'h=1', 'scheme=trap',
                     'tail=geometric' )
    expect_badinput( 'mathutils','osclines','A=1' )
    expect_badinput( 'mathutils','osclines','A=-1','betaE=1' )
    expect_badinput( 'mathutils','shiftadd',w,w,w,w[:-1],w )
    q2 = ['mathutils','convolve']
    expect_badinput( *q2, arr, arr, w, w, w, 'scheme=bla' )
    expect_badinput( *q2, w, w, w, w, w, 'interp=bla' )
    expect_badinput( *q2, w, w, w, w, w, 'scheme=exactlin', 'interp=spline' )
    expect_badinput( *q2, w, w, w, w[:-1], w )
    expect_badinput( *q2, arr, arr, w, w, w )
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
    test_tail()
    test_convolve()
    test_overshoot()
    test_lines()
    test_errors()
    test_placeholders()
    print('All OK')

if __name__ == '__main__':
    main()
