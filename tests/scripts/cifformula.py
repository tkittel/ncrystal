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


# Tests of the chemical formula helpers used to check loaded CIF data against
# the formula stated in it (_cifcomposerimpl).

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev._cifcomposerimpl import ( _decode_formula, _formula_to_dict,
                                          _formulas_incompatible )

def test_decode():
    for s, expect in ( ( 'Al2 O3', { 'Al': 2, 'O': 3 } ),
                       ( 'Fe2O3', { 'Fe': 2, 'O': 3 } ),
                       ( 'C12 H19 N O2', { 'C': 12, 'H': 19, 'N': 1, 'O': 2 } ),
                       ( 'H2 O H2', { 'H': 4, 'O': 1 } ),
                       ( 'D2 O', { 'D': 2, 'O': 1 } ),
                       ( 'Si O2 ', { 'Si': 1, 'O': 2 } ),
                       ( 'Xq2', None ),
                       ( '', None ) ):
        r = _decode_formula( s )
        print( f'{s!r:>16} -> {r}' )
        assert r == expect, ( s, r )
    assert _formula_to_dict( [ ( 'O', 2 ), ( 'O', 1 ), ( 'Na', 0 ) ] ) == { 'O': 3 }
    assert _formula_to_dict( None ) is None and _formula_to_dict( [] ) is None

def test_incompatible():
    for f1, f2, expect in ( ( 'Al2 O3', 'Al4 O6', False ),#(same proportions)
                            ( 'Al2 O3', [ ( 'Al', 2 ), ( 'O', 3.015 ) ], False ),#(within 1%)
                            ( 'Al2 O3', [ ( 'Al', 2 ), ( 'O', 3.1 ) ], True ),
                            ( 'Al2 O3', 'Al2 O4', True ),
                            ( 'Al2 O3', 'Al2 Cr O3', True ),
                            ( 'H2 O', 'D2 O', False ),#(isotopes as H)
                            ( 'H2 O', 'Xq', False ),#(not decodable: no check)
                            ( [ ( 'Na', 4.0 ), ( 'Cl', 4.0 ) ], 'Na Cl', False ) ):
        r = _formulas_incompatible( f1, f2 )
        print( f'{f1!r} versus {f2!r}: incompatible={r}' )
        assert r == expect, ( f1, f2 )

test_decode()
test_incompatible()
