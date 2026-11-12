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


# Tests of the cross-platform reproducible arithmetic in _stablemath.

import NCTestUtils.enable_fpe # noqa F401
import math
from fractions import Fraction
import NCrystalDev._stablemath as sm

def test_cos_sin():
    #Exact values:
    for a, c in ( (0,1.0), (60,0.5), (90,0.0), (120,-0.5), (180,-1.0),
                  (240,-0.5), (270,0.0), (300,0.5), (360,1.0), (-60,0.5),
                  (420,0.5), (90.0,0.0) ):
        assert sm.cos_deg( a ) == c, a
    for a, sv in ( (0,0.0), (30,0.5), (90,1.0), (150,0.5), (180,0.0),
                   (270,-1.0), (-30,-0.5) ):
        assert sm.sin_deg( a ) == sv, a
    #Other values: close to math.cos/sin, and printed as reference (the
    #reference log thus also checks cross-platform reproducibility):
    angles = ( 1e-9, 0.5, 30, 45, 84, 89.999999, 90.000001, 97, 105, 109.47,
               119.9954, 135, 179.5, 200, 359.9, -45 )
    for a in angles:
        c, s = sm.cos_deg( a ), sm.sin_deg( a )
        assert abs( c - math.cos( math.radians( a ) ) ) < 1e-15
        assert abs( s - math.sin( math.radians( a ) ) ) < 1e-15
        print( f'cos_deg({a!r}) = {c!r}, sin_deg({a!r}) = {s!r}' )
    #Symmetries hold exactly (when the shifted angle is exact as a float):
    def exact( x, y ):
        return Fraction( x ) == Fraction( y )
    for a in angles:
        c, s = sm.cos_deg( a ), sm.sin_deg( a )
        assert sm.cos_deg( -a ) == c and sm.sin_deg( -a ) == -s
        if exact( 180 - a, 180 - Fraction( a ) ):
            assert sm.cos_deg( 180 - a ) == -c and sm.sin_deg( 180 - a ) == s
        if exact( a + 360, Fraction( a ) + 360 ):
            assert sm.cos_deg( a + 360 ) == c and sm.sin_deg( a + 360 ) == s

def test_snap_rational():
    for x, expect in ( ( 0.5, Fraction(1,2) ),
                       ( 1/3, Fraction(1,3) ),
                       ( -2/3 + 1e-12, Fraction(-2,3) ),
                       ( 0.1234567, Fraction(0.1234567) ),
                       ( 1.0, Fraction(1) ) ):
        assert sm.snap_rational( x ) == expect, x
    print( 'snap_rational(0.1234567, maxden=1000) is exact:',
           sm.snap_rational( 0.1234567 ) == Fraction( 0.1234567 ) )

def test_transform_sym():
    u = [ [ 0.0123, 0.0021, -0.0032 ],
          [ 0.0021, 0.0087, 0.0011 ],
          [ -0.0032, 0.0011, 0.0154 ] ]
    #Identity: unchanged.
    ident = [ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ]
    assert sm.transform_sym( ident, u ) == u
    #Hexagonal 3-fold rotation (integer matrix) and a basis change with
    #fractions: exactly the correctly rounded result, and exactly symmetric.
    for name, m in ( ( '3-fold', [ [ 0, -1, 0 ], [ 1, -1, 0 ], [ 0, 0, 1 ] ] ),
                     ( 'C-centring', [ [ Fraction(1,2), Fraction(1,2), 0 ],
                                       [ Fraction(-1,2), Fraction(1,2), 0 ],
                                       [ 0, 0, 1 ] ] ) ):
        r = sm.transform_sym( m, u )
        fm = [ [ Fraction( v ) for v in row ] for row in m ]
        fu = [ [ Fraction( v ) for v in row ] for row in u ]
        exact = [ [ sum( fm[i][k] * fu[k][n] * fm[j][n]
                         for k in range(3) for n in range(3) )
                    for j in range(3) ] for i in range(3) ]
        assert r == [ [ float( v ) for v in row ] for row in exact ]
        assert all( r[i][j] == r[j][i] for i in range(3) for j in range(3) )
        print( f'{name}:', [ [ repr( v ) for v in row ] for row in r ] )

def test_scale():
    vals = [ 0.0123, 0.0087, -0.0032, 0.0 ]
    r = sm.scale( vals, 0.0134366, 0.013429928642 )
    assert r == [ float( Fraction( v ) * Fraction( 0.0134366 )
                         / Fraction( 0.013429928642 ) ) for v in vals ]
    assert sm.scale( vals, 2, 2 ) == vals
    print( 'scale:', [ repr( v ) for v in r ] )

def main():
    test_cos_sin()
    test_snap_rational()
    test_transform_sym()
    test_scale()

if __name__ == '__main__':
    main()
