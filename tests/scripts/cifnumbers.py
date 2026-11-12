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


# Tests of _cifreadimpl.parse_cif_number: values of CIF numeric strings and their
# rounding errors (half the step of the last printed digit).

import NCTestUtils.enable_fpe # noqa F401
import math
from NCrystalDev._cifreadimpl import parse_cif_number as vu

def check( s, value, rounding_error ):
    r = vu( s )
    assert r is not None, s
    v, e = r
    assert v == value and math.copysign( 1.0, v ) == math.copysign( 1.0, value ), ( s, v, value )
    assert e == rounding_error, ( s, e, rounding_error )
    print( f'{s!r:>16} -> value {v!r}, rounding error {e!r}' )

def test_decimal():
    check( '0.0283', 0.0283, 5e-05 )
    check( '-0.0032', -0.0032, 5e-05 )
    check( '+0.5', 0.5, 0.05 )
    check( '.34', 0.34, 0.005 )
    check( '-.5', -0.5, 0.05 )
    check( '1.', 1.0, 0.5 )
    check( '12.75', 12.75, 0.005 )
    check( '0.343400', 0.3434, 5e-07 )#(trailing zeros are printed digits)
    check( '0.000', 0.0, 0.0005 )#(a rounded zero, not exact)
    check( '-0.000', 0.0, 0.0005 )
    check( ' 0.012 ', 0.012, 0.0005 )
    check( '007.50', 7.5, 0.005 )

def test_uncertainty():
    #(the standard uncertainty in parentheses is ignored):
    check( '0.345(34)', 0.345, 0.0005 )
    check( '0.0283(4)', 0.0283, 5e-05 )
    check( '1.23e-4(2)', 1.23e-4, 5e-07 )
    check( '12(3)', 12.0, 0.5 )#(an integer with uncertainty is measured)
    check( '0(1)', 0.0, 0.5 )

def test_integer():
    check( '2344', 2344.0, 0.5 )
    check( '-7', -7.0, 0.5 )
    check( '100', 100.0, 0.5 )
    #Zero without decimals: exact (a constraint in practice):
    for s in ( '0', '-0', '+0', '00', '0e-3', '0E5' ):
        check( s, 0.0, 0.0 )

def test_exponent():
    check( '1234.2e-6', 1234.2e-6, 5e-08 )
    check( '0.0001E2', 0.01, 0.005 )
    check( '3400e-6', 3400e-6, 5e-07 )#(trailing zeros of the mantissa count)
    check( '1e-3', 1e-3, 0.0005 )
    check( '5.e-3', 5e-3, 0.0005 )
    check( '2.83E+2', 283.0, 0.5 )
    check( '-6.02e23', -6.02e23, 5e+20 )
    check( '.5e1', 5.0, 0.5 )

def test_invalid():
    for s in ( '?', '.', '', ' ', '-', '+', 'e3', '.e3', '1,5', 'inf', '-inf',
               'nan', '0x10', '1.2.3', '1e', '1e+', '1.2(3)e-4', '(3)',
               '1.2(3', '1.2(a)', '1 2', 'abc', "'0.1'", '1.0d3',
               '1e400', '-1e400', '0.0e9999', '1.5e999999' ):
        assert vu( s ) is None, s
    print( 'invalid strings rejected' )

def test_extreme_exponents():
    #(tiny values: below the smallest double, so no rounding error):
    check( '1e-500', 0.0, 0.0 )
    check( '1e-1000000000', 0.0, 0.0 )#(must not hang)
    check( '1.5e-320', 1.5e-320, 5e-322 )#(subnormal)
    check( '1e300', 1e300, 5e299 )

def main():
    test_decimal()
    test_uncertainty()
    test_integer()
    test_exponent()
    test_invalid()
    test_extreme_exponents()

if __name__ == '__main__':
    main()
