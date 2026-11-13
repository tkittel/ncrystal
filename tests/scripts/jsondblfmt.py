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

#Exact JSON output of doubles from the C++ layer (streamJSON and the arrays of
#streamJSONHugeDblVect via streamJSONFastNotPretty), checked against
#a platform independent reference computed with Python's correctly rounded
#formatting, and decoded back to the identical Python floats.

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev.misc import evaluate_query
import ctypes
import json
import locale
import math
import os
import random
import struct
import sys

def bits( x ):
    return struct.unpack( '<Q', struct.pack( '<d', x ) )[0]

def ref_shortstr( x ):
    #Reference for NCrystal's dbl2shortstr: %.15g if it round-trips, otherwise
    #%.17g (both with fewer digits for denormals).
    nstrip, t = 0, x
    while abs( t ) < sys.float_info.min:
        nstrip += 1
        t *= 10
    if nstrip:
        nstrip -= 1
    s = f'{x:.{max(0,14-nstrip)+1}g}'
    if float( s ) != x:
        s = f'{x:.{max(1,16-nstrip)+1}g}'
    return s

def ref_streamjson( x ):
    if math.isinf( x ):
        return '1.0e99999' if x > 0 else '-1.0e99999'
    if x == 0.0:
        return '0.0'
    s = ref_shortstr( x )
    if not any( c in s for c in '.eE' ):
        return f'{x:.1f}'
    return s

def ref_fastnotpretty( x ):
    if math.isinf( x ) or x == 0.0:
        return ref_streamjson( x )
    s = f'{x:.17g}'
    if abs( x ) < 1e17 and x == math.floor( x ):
        s += '.0'
    return s

def run_query( values ):
    raw = evaluate_query( ['util','jsontest']
                          + [ str( bits(x) ) for x in values ], unpack = False )
    def tokens( key ):
        return raw.split( f'"{key}":[', 1 )[1].split( ']', 1 )[0].split( ',' )
    return json.loads( raw ), tokens( 'scalar' ), tokens( 'array' )

def edge_values():
    return [ 1e16, 10461373122004696.0, 2.0**53, 2.0**53 + 2, -1e22, 1e17,
             1e-17, 0.1, 5.0, -5.0, 1234.5, 3.0e-5, 1e-300, -2.5e-300,
             5e-324, -5e-324, 2.2250738585072014e-308,
             2.225073858507201e-308, 1.7976931348623157e308,
             -1.7976931348623157e308, 1.0, 0.0, -0.0, math.inf, -math.inf ]

def random_values( n ):
    rng = random.Random( 123456 )
    res = []
    while len( res ) < n:
        b = rng.getrandbits( 64 )
        if ( b >> 52 ) & 0x7ff == 0x7ff:
            continue#inf or nan (signalling nans would trigger the fpe trap)
        res.append( struct.unpack( '<d', struct.pack( '<Q', b ) )[0] )
    return res

def check( values, label ):
    d, tok_scalar, tok_array = run_query( values )
    for key, toks, reffct in ( ( 'scalar', tok_scalar, ref_streamjson ),
                               ( 'array', tok_array, ref_fastnotpretty ) ):
        assert len( toks ) == len( values ) == len( d[key] )
        nbad = 0
        for x, tok, dec in zip( values, toks, d[key] ):
            ok = ( tok == reffct( x ) and isinstance( dec, float )
                   and dec == x )
            if not ok:
                nbad += 1
                if nbad <= 5:
                    print( f'  MISMATCH ({key}): {x!r}: got {tok},'
                           f' expected {reffct(x)}' )
        print( f'{label}: {len(values)} values, {key} output as expected'
               f' for {len(values)-nbad}' )
        assert nbad == 0

def check_locales( values ):
    #Output must not depend on the decimal point of the C locale. Uses the C
    #library's setlocale directly (POSIX only), since it returns NULL for
    #unavailable locales:
    if os.name != 'posix':
        return
    setlocale = ctypes.CDLL( None ).setlocale
    setlocale.restype = ctypes.c_char_p
    setlocale.argtypes = [ ctypes.c_int, ctypes.c_char_p ]
    args = ['util','jsontest'] + [ str( bits(x) ) for x in values ]
    ref = evaluate_query( args, unpack = False )
    for name in ( 'de_DE.UTF-8', 'de_DE.utf8', 'fr_FR.UTF-8', 'fr_FR.utf8',
                  'ps_AF.UTF-8', 'ps_AF.utf8', 'fa_IR.UTF-8', 'fa_IR.utf8' ):
        if setlocale( locale.LC_NUMERIC, name.encode() ) is None:
            continue
        raw = evaluate_query( args, unpack = False )
        setlocale( locale.LC_NUMERIC, b'C' )
        assert raw == ref, f'output differs in locale {name}'

def main():
    check( edge_values(), 'edge values' )
    for x in edge_values()[:8]:
        print( f'  {x!r} -> {ref_streamjson(x)} / {ref_fastnotpretty(x)}' )
    check( random_values( 20000 ), 'random bit patterns' )
    check_locales( edge_values() + random_values( 2000 ) )

if __name__ == '__main__':
    main()
