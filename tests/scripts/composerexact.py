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

# Exact NCMAT output of the NCMATComposer for crystals (without spglib):
# special positions written as fractions however they were given, values
# within the number format, idempotent output (writing again what was loaded
# gives the same text), the dynamic info fractions, and the verification of
# the space group (including a higher symmetry in the same coordinates).

import NCTestUtils.reprint_escaped_warnings # noqa F401
from NCTestUtils.common import ensure_error

import NCrystalDev as NC

def lines( c, *starts ):
    return [ ll.strip() for ll in c.create_ncmat().splitlines()
             if ll.strip().startswith( starts ) ]

def hexc( sites, *, use_sites = True, sg = 194 ):
    c = NC.NCMATComposer()
    c.set_cellsg( a = 3.0, b = 3.0, c = 5.0, alpha = 90, beta = 90,
                  gamma = 120, spacegroup = sg )
    ( c.set_atomsites if use_sites else c.set_atompos )( sites )
    c.set_dyninfo_debyetemp( 'C', 300 )
    return c

def test_special_positions():
    #The site 2c of P6_3/mmc, however precisely given, is written as exact
    #fractions (the sites are snapped as when the data is loaded):
    for x, y in ( ( 0.3333, 0.6667 ), ( 0.33333, 0.66667 ),
                  ( 0.333333333, 0.666666667 ), ( 1/3, 2/3 ) ):
        print( f'site ({x}, {y}, 1/4):', lines( hexc( [ ('C',x,y,0.25) ] ),
                                                 'C ' ) )
    #Also when given as all positions:
    c = hexc( [ ('C',0.3333,0.6667,0.25), ('C',0.6667,0.3333,0.75) ],
              use_sites = False )
    print( 'positions (0.3333, 0.6667, 1/4) etc.:', lines( c, 'C ' ) )
    #A general position is written as given (within the number format, here
    #also with rounding noise in a small coordinate removed):
    print( 'general position:',
           lines( hexc( [ ('C',0.1234,0.2345,0.008389999999999898) ],
                        sg = 1 ), 'C ' ) )

def test_idempotent():
    #Composing a material from the loaded data and writing it again gives the
    #same text, also with calculated values (noise below the number format):
    c = NC.NCMATComposer()
    c.set_cellsg( a = 4.0 * ( 1 + 1e-15 ), b = 4.0 * ( 1 + 1e-15 ),
                  c = 5.0 / 3 * 3, alpha = 90, beta = 90, gamma = 90,
                  spacegroup = 123 )
    c.set_atomsites( [ ('Al',0,0,0), ('Cu',0.5,0.5,0.5),
                       ('O',0.1+0.2,0.3,0.0) ] )
    for e in ( 'Al', 'Cu', 'O' ):
        c.set_dyninfo_debyetemp( e, 1000.0 / 3 )
    def data( text ):
        #(without the comments, which say where data was transferred from):
        return [ ll for ll in text.splitlines()
                 if not ll.strip().startswith('#') ]
    t1 = data( c.create_ncmat() )
    t2 = data( NC.NCMATComposer( data = '\n'.join( t1 ) ).create_ncmat() )
    print( 'idempotent:', t1 == t2 )
    assert t1 == t2
    print( '\n'.join( '  ' + ll for ll in t1 ) )

def test_fractions():
    #Dynamic info fractions: exact from the counts, exact with occupancies
    #written as fractions, and decimals with measured occupancies:
    def fcc( occu ):
        c = NC.NCMATComposer()
        c.set_cellsg_cubic( 4.0, spacegroup = 225 )
        c.set_atomsites( [ ('Al',0,0,0,occu), ('O',0.5,0.5,0.5) ] )
        c.set_dyninfo_debyetemp( 'Al', 300 )
        c.set_dyninfo_debyetemp( 'O', 300 )
        return c
    for occu in ( 1.0, 0.9, 22/25, 0.883 ):
        print( f'occupancy {occu}:',
               lines( fcc( occu ), 'fraction', 'Al ' ) )

def test_verification():
    #Positions which do not have the symmetry of the space group, cell
    #parameters breaking its constraints, and a structure with a higher
    #symmetry in the same coordinates are rejected:
    c = hexc( [ ('C',0.3333,0.6667,0.25), ('C',0.1,0.2,0.75) ],
              use_sites = False )
    with ensure_error( NC.NCBadInput ):
        c.create_ncmat()
    c = NC.NCMATComposer()
    c.set_cellsg( a = 4.0, b = 4.1, c = 4.0, alpha = 90, beta = 90,
                  gamma = 90, spacegroup = 225 )
    c.set_atompos( [ ('Al',0,0,0), ('Al',0,0.5,0.5), ('Al',0.5,0,0.5),
                     ('Al',0.5,0.5,0) ] )
    c.set_dyninfo_debyetemp( 'Al', 300 )
    with ensure_error( NC.NCBadInput, 'Lattice parameters a and b must be'
                       ' equal for this space group (got a=4 and b=4.1)' ):
        c.create_ncmat()
    #(also tiny deviations):
    c.set_cellsg( a = 4.0, b = 4.0 * ( 1 + 1e-12 ), c = 4.0, alpha = 90,
                  beta = 90, gamma = 90, spacegroup = 225 )
    with ensure_error( NC.NCBadInput, 'Lattice parameters a and b must be'
                       ' equal for this space group (got a=4 and'
                       ' b=4.000000000004)' ):
        c.create_ncmat()
    c.set_cellsg_cubic( 4.0, spacegroup = 221 )
    with ensure_error( NC.NCBadInput, 'The crystal structure has the higher'
                       ' symmetry of space group setting 225 (in the same'
                       ' coordinates) than the given 221. Use'
                       ' refine_crystal_structure() to change to it, or'
                       ' verify_crystal_structure=False to keep the given'
                       ' one.' ):
        c.create_ncmat()
    t = c.create_ncmat( verify_crystal_structure = False )
    print( 'with verify_crystal_structure=False:',
           [ ll.strip() for ll in t.splitlines() if ll.strip() == '221' ] )

test_special_positions()
test_idempotent()
test_fractions()
test_verification()
