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
from NCrystalDev._fmt import fmt_value_rtol, fmt_vector_rtol

rtols = [ 0.0, 5e-14, 1e-9, 5e-6, 1e-3 ]

def check_values( xs, rtol ):
    #Within rtol/2 (exact for rtol below 1e-14, which needs more than 15
    #digits), and idempotent (all values in two batched calls):
    texts = fmt_vector_rtol( xs, rtol, strip_leading_zero = False )
    ys = [ float( s ) for s in texts ]
    for x, s, y in zip( xs, texts, ys ):
        if x == 0.0:
            assert s == '0'
        elif rtol < 1e-14:
            assert y == x, (x,s)
        else:
            assert abs( y - x ) <= 0.5 * rtol * abs( x ), (x,rtol,s)
    assert fmt_vector_rtol( ys, rtol, strip_leading_zero = False ) == texts
    return texts

def check_value( x, rtol ):
    s = fmt_value_rtol( x, rtol )
    assert check_values( [ x ], rtol ) == [ s ]
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
    xs = [ rng.choice( (1,-1) ) * 10**rng.uniform( -30, 30 )
           for _ in range( 20000 ) ]
    for r in rtols:
        check_values( xs, r )
    print(f'Checked {len(xs)*len(rtols)} random values')
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
