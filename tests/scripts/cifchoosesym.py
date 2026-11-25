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

# Tests of _cifsymimpl.determine_setting: the space group setting of CIF data,
# from the listed symmetry operations (cross-checked against the symbol), or
# else from the symbol (several candidates resolved by the cell), or from an
# explicit override, in default and permissive mode.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev._cifreadimpl import CIFContent
from NCrystalDev._cifsymimpl import determine_setting
from NCrystalDev.misc import evaluate_query
from NCrystalDev._common import WarningSpy

cubic = { 'a': '5.43', 'b': '5.43', 'c': '5.43',
          'alpha': '90', 'beta': '90', 'gamma': '90' }
hexagonal = { 'a': '3.0', 'b': '3.0', 'c': '10.0',
              'alpha': '90', 'beta': '90', 'gamma': '120' }
rhombohedral = { 'a': '4.0', 'b': '4.0', 'c': '4.0',
                 'alpha': '70', 'beta': '70', 'gamma': '70' }

def ops_of( setting ):
    return evaluate_query( [ 'sgsym', 'info', 'strict', setting ] )['ops']

def content( cell, ops = None, hall = None, hm = None, itnumber = None,
             text = '' ):
    return CIFContent( text = text, cif_raw = None, cell = cell,
                       listed_ops = ops, hall = hall, hm = hm,
                       itnumber = itnumber, formula = None, sites = [],
                       atom_type_number_in_cell = [], temperatures = {} )

def test( descr, c, permissive = False, override = None ):
    with WarningSpy( block = True ) as warnings:
        try:
            s = determine_setting( c, permissive = permissive,
                                   override_spacegroup = override )
        except NC.NCBadInput as e:
            print( f'{descr}: rejected: {e}' )
            s = None
        else:
            print( f'{descr}: {s}' )
    for w, _ in warnings:
        print( f'   warning: {w}' )
    return s

def main():
    ops227 = ops_of( '227:2' )
    #A: listed operations (any order, with duplicates), cross-checked against
    #the symbol:
    assert test( 'listed operations', content( cubic, ops227 ) ) == '227:2'
    test( 'listed operations, reversed with duplicates',
          content( cubic, ops227[::-1] + ops227[:3] ) )
    test( 'listed operations and matching symbol',
          content( cubic, ops227, hm = 'F d -3 m' ) )
    test( 'listed operations and unrecognised symbol (ignored)',
          content( cubic, ops227, hm = 'F d -3 m (bla)' ) )
    test( 'shared operations of number 68',
          content( { 'a': '5', 'b': '6', 'c': '7', 'alpha': '90',
                     'beta': '90', 'gamma': '90' }, ops_of( '68:1' ) ) )
    #A.2: unusable operations are an error, in permissive mode the symbol is
    #used:
    for descr, ops in ( ( 'missing operation', ops227[:-1] ),
                        ( 'unparsable operation', ops227[:-1] + [ 'x+y+z' ] ),
                        ( 'origin shift', [ 'x,y,z', '-x+1/2,-y+1/2,-z+1/2' ] ) ):
        for permissive in ( False, True ):
            test( f'{descr}, permissive={permissive}',
                  content( cubic, ops, hm = 'F d -3 m:2' ),
                  permissive = permissive )
    test( 'missing operation, without symbol, permissive',
          content( cubic, ops227[:-1] ), permissive = True )
    test( 'missing operation, pymatgen file',
          content( cubic, ops227[:-1], text = '# generated using pymatgen' ) )
    #A.3: operations and symbol of different space groups:
    for permissive in ( False, True ):
        test( f'symbol mismatch, permissive={permissive}',
              content( cubic, ops227, hm = 'F m -3 m' ),
              permissive = permissive )
    test( 'Hall symbol mismatch', content( cubic, ops227, hall = '-F 4 2 3',
                                           hm = 'F d -3 m' ) )
    #B: no listed operations, the symbol (Hall, H-M or number):
    test( 'Hall symbol', content( cubic, hall = '-F 4vw 2vw 3' ) )
    test( 'H-M symbol with origin choice', content( cubic, hm = 'F d -3 m:1' ) )
    test( 'H-M symbol, one setting', content( cubic, hm = 'Fm-3m' ) )
    test( 'unrecognised Hall symbol, H-M used',
          content( cubic, hall = 'bla', hm = 'Fm-3m' ) )
    test( 'IT number, one setting', content( cubic, itnumber = '225' ) )
    test( 'R -3 m, hexagonal axes', content( hexagonal, hm = 'R -3 m' ) )
    test( 'R -3 m, rhombohedral axes', content( rhombohedral, hm = 'R -3 m' ) )
    test( 'R -3 m, unfitting cell', content( cubic | { 'gamma': '100' },
                                             hm = 'R -3 m' ) )
    test( 'ambiguous origin choice', content( cubic, hm = 'F d -3 m' ) )
    test( 'ambiguous IT number', content( cubic, itnumber = '227' ) )
    test( 'no symmetry information', content( cubic ) )
    test( 'unrecognised symbol', content( cubic, hm = 'bla' ) )
    #C: explicit override (ignoring the symmetry information of the data):
    test( 'override, H-M symbol',
          content( cubic, ops227, hm = 'F d -3 m' ), override = 'F d -3 m:1' )
    test( 'override, setting', content( cubic, hm = 'F d -3 m' ),
          override = '227:1' )
    test( 'override, ambiguous', content( cubic ), override = 'F d -3 m' )
    test( 'override, unknown', content( cubic ), override = 'bla' )

if __name__ == '__main__':
    main()
