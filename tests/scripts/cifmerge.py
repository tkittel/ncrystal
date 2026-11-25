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


# Tests of _cifsitesimpl.merge_atoms: CIF sites at the same positions are
# merged into atoms with a composition, and inconsistencies are rejected.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev._cifsitesimpl import merge_atoms

def atom( lbl, elem, coords, occ = 1.0, uiso = 0.01, aniso = None ):
    return { 'cif_label': lbl, 'element_marker': elem,
             'expanded_coords': [ tuple( c ) for c in coords ],
             'rotations': [ None ] * len( coords ), 'occupancy': occ,
             'uiso': uiso, 'aniso': aniso }

orbit_a = [ ( 0.0, 0.0, 0.0 ), ( 0.5, 0.5, 0.5 ) ]
orbit_a2 = [ ( 0.5, 0.5, 0.5 ), ( 0.0, 0.0, 0.0 ) ]#(same, other order)
orbit_b = [ ( 0.25, 0.25, 0.25 ), ( 0.75, 0.75, 0.75 ) ]

def show( descr, atoms ):
    print( f'{descr}:' )
    for a in atoms:
        print( '  ', a['cif_labels'], a['composition'], 'occupancy',
               a['occupancy'], 'positions', a['positions'] )

def expect_error( descr, *args ):
    try:
        merge_atoms( *args )
    except NC.NCBadInput as e:
        print( f'{descr}: rejected: {e}' )
    else:
        raise SystemExit( f'{descr}: not rejected' )

def main():
    a = merge_atoms( [ atom( 'Na1', 'Na', orbit_a ),
                       atom( 'Cl1', 'Cl', orbit_b ) ], 4 )
    show( 'separate sites', a )
    assert len( a ) == 2 and a[0]['occupancy'] == 1.0
    #Mixed occupancy on one site (in any order of the positions):
    a = merge_atoms( [ atom( 'Al1', 'Al', orbit_a, occ = 0.75 ),
                       atom( 'Cr1', 'Cr', orbit_a2, occ = 0.25 ) ], None )
    show( 'mixed site', a )
    assert len( a ) == 1 and a[0]['composition'] == ( ( 0.25, 'Cr' ), ( 0.75, 'Al' ) )
    #Partial occupancy (composition normalised, occupancy kept):
    a = merge_atoms( [ atom( 'O1', 'O', orbit_a, occ = 0.4 ),
                       atom( 'F1', 'F', orbit_a, occ = 0.4 ) ], None )
    show( 'partially occupied mixed site', a )
    assert a[0]['occupancy'] == 0.8 and a[0]['composition'] == ( ( 0.5, 'F' ), ( 0.5, 'O' ) )
    #Rounding of occupancies slightly above 1 is accepted:
    a = merge_atoms( [ atom( 'A', 'Al', orbit_a, occ = 0.5000004 ),
                       atom( 'B', 'Cr', orbit_a, occ = 0.5000004 ) ], None )
    assert a[0]['occupancy'] == 1.0
    expect_error( 'occupancy above 1',
                  [ atom( 'A', 'Al', orbit_a, occ = 0.7 ),
                    atom( 'B', 'Cr', orbit_a, occ = 0.7 ) ], None )
    expect_error( 'same site given twice (occupancy 2)',
                  [ atom( 'A', 'Al', orbit_a ), atom( 'B', 'Al', orbit_a ) ], None )
    expect_error( 'different displacements at the same positions',
                  [ atom( 'A', 'Al', orbit_a, occ = 0.5 ),
                    atom( 'B', 'Cr', orbit_a, occ = 0.5, uiso = 0.02 ) ], None )
    expect_error( 'different tensors at the same positions',
                  [ atom( 'A', 'Al', orbit_a, occ = 0.5, aniso = { 'u11': 0.01 } ),
                    atom( 'B', 'Cr', orbit_a, occ = 0.5, aniso = { 'u11': 0.02 } ) ],
                  None )
    #Orbits of different sizes are never merged:
    a = merge_atoms( [ atom( 'A', 'Al', orbit_a ),
                       atom( 'B', 'Cr', orbit_a[:1] ) ], None )
    assert len( a ) == 2
    #Atom count check (1% tolerance):
    merge_atoms( [ atom( 'Na1', 'Na', orbit_a ) ], 2 )
    merge_atoms( [ atom( 'Na1', 'Na', orbit_a ) ], 2.01 )
    expect_error( 'wrong atom count', [ atom( 'Na1', 'Na', orbit_a ) ], 3 )
    print( 'atom count checks done' )

if __name__ == '__main__':
    main()
