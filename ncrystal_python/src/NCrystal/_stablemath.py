
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

"""

Internal helpers for cross-platform reproducible arithmetic: results are
bitwise identical on all platforms. The sine and cosine of angles in degrees
come from the C++ library (stable_sincos_deg), the rest uses exact rational
arithmetic (fractions) and correctly rounded conversions to float. Unlike
e.g. math.cos (libm) or numpy matrix products (SIMD, FMA, summation order),
they are thus free of platform dependent ulp noise.

"""

__all__ = []

import functools as _functools
from fractions import Fraction as _Fraction


def _frac( x ):
    return x if isinstance( x, _Fraction ) else _Fraction( x )

def sincos_deg( x ):
    """(sin,cos) of the angle x in degrees, from stable_sincos_deg of the C++
    library (exact for multiples of 30 and 45 degrees)."""
    return _sincos_deg_impl( float( x ) )

@_functools.lru_cache( maxsize = 256 )#(callers repeat the few cell angles)
def _sincos_deg_impl( x ):
    from .misc import evaluate_query
    ( s, c ), = evaluate_query( [ 'util', 'sincosdeg', '@'+repr( x ) ] )
    return float( s ), float( c )

def cos_deg( x ):
    """Cosine of the angle x in degrees (cf. sincos_deg)."""
    return sincos_deg( x )[1]

def sin_deg( x ):
    """Sine of the angle x in degrees (cf. sincos_deg)."""
    return sincos_deg( x )[0]

def snap_rational( x, maxden = 1000, tol = 1e-9 ):
    """x as a Fraction: the nearest one with denominator at most maxden if
    within tol, otherwise the exact value of x."""
    f = _frac( x )
    s = f.limit_denominator( maxden )
    return s if abs( s - f ) <= tol else f

def matmul( a, b ):
    """Exact product of matrices (nested sequences of numbers), as nested
    lists of Fractions."""
    a = [ [ _frac( v ) for v in row ] for row in a ]
    b = [ [ _frac( v ) for v in row ] for row in b ]
    return [ [ sum( ( a[i][k] * b[k][j] for k in range( len( b ) ) ),
                    _Fraction( 0 ) )
               for j in range( len( b[0] ) ) ] for i in range( len( a ) ) ]

def transform_sym( m, u ):
    """M.U.M^T of square matrices (nested sequences of numbers), computed
    exactly and returned as nested lists of floats (each correctly
    rounded)."""
    mt = [ list( col ) for col in zip( *m ) ]
    res = matmul( matmul( m, u ), mt )
    return [ [ float( v ) for v in row ] for row in res ]

def scale( values, num, den ):
    """The values (sequence of numbers) times num/den, each computed exactly
    and returned as a correctly rounded float."""
    f = _frac( num ) / _frac( den )
    return [ float( _frac( v ) * f ) for v in values ]
