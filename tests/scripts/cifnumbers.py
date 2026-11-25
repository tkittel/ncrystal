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


# Tests of _cifimpl.parse_cif_number: values of CIF numeric strings and their
# rounding errors (half the step of the last printed digit).

import NCTestUtils.enable_fpe # noqa F401
import math
from fractions import Fraction
from NCrystalDev._cifimpl import parse_cif_number as vu

def up( x ):
    #Smallest double not below the decimal value of x (e.g. 5e-05):
    exact = Fraction( repr( x ) )
    r = float( exact )
    return math.nextafter( r, math.inf ) if Fraction( r ) < exact else r

def check( s, value, rounding_error ):
    #(rounding_error given as the intended decimal value, e.g. 5e-05):
    r = vu( s )
    assert r is not None, s
    v, e = r
    assert v == value and math.copysign( 1.0, v ) == math.copysign( 1.0, value ), ( s, v, value )
    assert e == up( rounding_error ), ( s, e, rounding_error )
    assert Fraction( e ) >= Fraction( repr( rounding_error ) )
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
    #Zeros without uncertainty are exact (e.g. 0 printed with %.3f):
    for s in ( '0.000', '-0.000', '0.', '.0', '0.0e3' ):
        check( s, 0.0, 0.0 )
    check( ' 0.012 ', 0.012, 0.0005 )
    check( '007.50', 7.5, 0.005 )

def test_uncertainty():
    #(the standard uncertainty in parentheses is ignored):
    check( '0.345(34)', 0.345, 0.0005 )
    check( '0.0283(4)', 0.0283, 5e-05 )
    check( '1.23e-4(2)', 1.23e-4, 5e-07 )
    check( '12(3)', 12.0, 0.5 )#(an integer with uncertainty is measured)
    check( '0(1)', 0.0, 0.5 )
    check( '0.000(2)', 0.0, 0.0005 )#(a refined zero is not exact)

def test_integer():
    check( '2344', 2344.0, 0.5 )
    check( '-7', -7.0, 0.5 )
    check( '100', 100.0, 0.5 )
    #Zero: exact (a constraint in practice):
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
               '1e400', '-1e400', '1.5e999999', '1e309' ):
        assert vu( s ) is None, s
    print( 'invalid strings rejected' )

def test_extreme_exponents():
    #(tiny values: below the smallest double, so no rounding error):
    check( '1e-500', 0.0, 0.0 )
    check( '1e-1000000000', 0.0, 0.0 )#(must not hang)
    assert vu( '1e-' + '9' * 5000 ) == ( 0.0, 0.0 )#(beyond int() digit limit)
    assert vu( '0.0e' + '9' * 5000 ) == ( 0.0, 0.0 )
    for s in ( '0.0e9999', '0.0e310', '0.00e400' ):
        check( s, 0.0, 0.0 )
    assert vu( '0.0e310(1)' ) is None#(rounding error not a finite double)
    #(leading zeros in the exponent do not count as digits):
    assert vu( '1e+000000000005' ) == ( 1e5, 5e4 )
    assert vu( '1e-00000000000000001' ) == ( 0.1, up( 0.05 ) )
    check( '1.5e-320', 1.5e-320, 5e-322 )#(subnormal)
    check( '1e300', 1e300, 5e299 )
    check( '0.0e308(1)', 0.0, 5e306 )
    check( '1e308', 1e308, 5e307 )

def test_floor():
    #Rounding errors are at least 1e-15 times the value (the conversion of the
    #printed number to a double), for numbers with very many digits:
    for s in ( '0.12345678901234567890', '-123456.789012345678' ):
        v, e = vu( s )
        assert e == abs( v ) * 1e-15, s
        print( f'{s!r:>24} -> rounding error {e!r} (floor)' )
    check( '0.1234567890123', 0.1234567890123, 5e-14 )#(above the floor)

def test_round_up():
    #Rounding errors are never below the exact value (e.g. the nearest double
    #to 5e-05 is below 5e-05), and exact ones are unchanged:
    for s, exact in ( ( '0.0283', '5e-5' ), ( '0.1', '0.05' ),
                      ( '1.23e-4', '5e-7' ), ( '2344', '0.5' ),
                      ( '7.5', '0.05' ), ( '3400e-6', '5e-7' ) ):
        e, exact = vu( s )[1], Fraction( exact )
        assert Fraction( e ) >= exact > Fraction( math.nextafter( e, -math.inf ) ), s
    assert vu( '2344' )[1] == 0.5 and vu( '1.' )[1] == 0.5
    print( 'rounding errors are the smallest doubles not below the exact values' )

def main():
    test_decimal()
    test_uncertainty()
    test_integer()
    test_exponent()
    test_invalid()
    test_extreme_exponents()
    test_floor()
    test_round_up()

if __name__ == '__main__':
    main()
