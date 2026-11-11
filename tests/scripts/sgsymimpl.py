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

# Unit tests of the tolerance helpers in _sgsymimpl.

import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.common import ensure_error
import NCrystalDev._sgsymimpl as sgs
from NCrystalDev.exceptions import NCBadInput

def show( **kw ):
    r = sgs._resolve_tolerances( **kw )
    print( f'{kw!r:50} -> ' + ', '.join( f'{k}={v:g}'
                                         for k, v in sorted( r.items() ) ) )
    return r

def test_resolve():
    #Default position tolerance is 0.01Aa, scaling linearly with tol:
    d = show()
    assert d['symprec'] == 0.01
    assert d['uij'] == d['occu'] == sgs._TOL_DEFAULT
    t = 5 * sgs._TOL_DEFAULT
    d = show( tol = t )
    assert abs( d['symprec'] - 0.05 ) < 1e-15
    assert d['uij'] == d['occu'] == t
    #Individual overrides:
    d = show( tol = t, tol_symprec = 0.003 )
    assert d['symprec'] == 0.003 and d['uij'] == t
    d = show( tol_uij = 0.02, tol_occu = 0.03 )
    assert d['symprec'] == 0.01 and d['uij'] == 0.02 and d['occu'] == 0.03

def test_anchor():
    #The 0.01Aa anchor does not depend on the value of _TOL_DEFAULT:
    orig = sgs._TOL_DEFAULT
    try:
        sgs._TOL_DEFAULT = 0.02
        d = show()
        assert d['symprec'] == 0.01 and d['uij'] == 0.02
        d = show( tol = 0.04 )
        assert abs( d['symprec'] - 0.02 ) < 1e-15
    finally:
        sgs._TOL_DEFAULT = orig

def test_errors():
    for kw in ( { 'tol': 0 }, { 'tol': -1e-3 }, { 'tol': 1.0 },
                { 'tol': 0.0501 }, { 'tol_uij': 0.06 }, { 'tol_symprec': 0.51 },
                { 'tol': float('nan') }, { 'tol': 'abc' },
                { 'tol_symprec': 0.0 }, { 'tol_symprec': float('inf') },
                { 'tol_uij': 2.0 }, { 'tol_occu': -0.1 } ):
        with ensure_error( NCBadInput, None ):
            sgs._resolve_tolerances( **kw )
    #The largest values:
    d = show( tol = 0.05 )
    assert abs( d['symprec'] - 0.5 ) < 1e-12
    show( tol = 1e-3, tol_symprec = 0.5, tol_uij = 0.05, tol_occu = 0.05 )

def test_compare():
    rc = sgs._relclose
    assert rc( 0.0, 0.0, 1e-3 )
    assert rc( 1.0, 1.001, 1e-3 ) and not rc( 1.0, 1.0011, 1e-3 )
    assert rc( -2.0, -2.002, 1e-3 ) and not rc( 0.0, 1e-12, 1e-3 )
    ec = sgs._eigclose
    e = ( 0.001, 0.005, 0.010 )
    #Tolerance is relative to the largest eigenvalue:
    assert ec( e, ( 0.001005, 0.005, 0.010 ), 1e-3 )
    assert not ec( e, ( 0.001, 0.005, 0.010015 ), 1e-3 )
    print( 'compare OK' )

def test_group():
    g = sgs._group_by_tol
    def close( a, b ):
        return sgs._relclose( a, b, 0.01 )
    #A chain 1.0, 1.008, 1.016: the third is compared with the
    #representative (1.0) only, so it starts a new group:
    print( g( [ 1.0, 1.008, 1.016, 2.0, 1.005, 1.02 ], close ) )
    print( g( [], close ) )
    print( g( [ 'a', 'b', 'a' ], lambda a, b : a == b ) )

test_resolve()
test_anchor()
test_errors()
test_compare()
test_group()
