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

# Tests the fmt_value_rtol and fmt_vector_rtol number formatting helpers.

import NCTestUtils.enable_fpe # noqa F401
import random
from NCrystalDev._fmt import fmt_value_rtol, fmt_vector_rtol, _fmt_ndigits

rtols = [ 0.0, 5e-14, 1e-9, 5e-6, 1e-3 ]

def check_value( x, rtol ):
    s = fmt_value_rtol( x, rtol )
    y = float( s )
    if x == 0.0:
        assert s == '0'
    elif rtol == 0.0:
        assert y == x, (x,s)
    else:
        assert abs( y - x ) <= rtol * abs( x ), (x,rtol,s)
        #No shorter number of significant digits suffices:
        digits = s.lstrip('-.0').split('e')[0].replace('.','').rstrip('0')
        nd = len( digits )
        if nd > 1:
            t = _fmt_ndigits( x, nd - 1, False )
            assert abs( float( t ) - x ) > rtol * abs( x ), (x,rtol,s,t)
    return s

def check_grid( grid, rtol ):
    res = fmt_vector_rtol( grid, rtol, is_grid = True )
    vals = [ float( s ) for s in res ]
    for x, s, y in zip( grid, res, vals ):
        if x == 0.0:
            assert s == '0'
        else:
            assert abs( y - x ) <= rtol * abs( x ), (x,rtol,s)
    for i in range( len( grid ) - 1 ):
        dx = grid[i+1] - grid[i]
        assert abs( ( vals[i+1] - vals[i] ) - dx ) <= rtol * abs( dx )
        assert vals[i+1] > vals[i]
    return res

def main():
    rng = random.Random( 123 )
    for x in ( 0.0, 300.0, 0.123456789, -0.123456789, 1e-300, 2.5e-7,
               0.06000000000000001, 1234567.89, 1.0/3.0 ):
        print( f'{x!r}: ' + ', '.join( f'rtol={r:g} -> {check_value(x,r)}'
                                        for r in rtols ) )
    n = 0
    for _ in range( 20000 ):
        x = rng.choice( (1,-1) ) * 10**rng.uniform( -30, 30 )
        for r in rtols:
            check_value( x, r )
            n += 1
    print(f'Checked {n} random values')
    #Grids, including neighbours much closer than the tolerance:
    grid = [ 1.0, 1.00001, 1.00002, 2.0, 3.0 ]
    for r in rtols:
        print(f'grid {grid} with rtol={r:g}: {check_grid( grid, r )}')
    ngrids = 0
    for _ in range( 200 ):
        g, x = [], rng.uniform( 0.0, 1.0 )
        for i in range( rng.randint( 2, 100 ) ):
            g.append( x )
            x += 10**rng.uniform( -12, 0 ) * x
        for r in rtols:
            check_grid( g, r )
            ngrids += 1
    print(f'Checked {ngrids} random grids')

if __name__ == '__main__':
    main()
