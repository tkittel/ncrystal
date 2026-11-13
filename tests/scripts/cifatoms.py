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


# NEEDS: gemmi spglib

# Tests of the atoms of loaded CIF data: isotopes (D, T), elements deduced from
# labels, mixed sites and their remapping, and the stated atom counts.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev as NC
from NCrystalDev.cifutils import CIFLoader

def cif( sites, extra = '' ):
    return '\n'.join( [ 'data_test', '_cell_length_a 4.0',
                        '_cell_length_b 4.0', '_cell_length_c 4.0',
                        '_cell_angle_alpha 90', '_cell_angle_beta 90',
                        '_cell_angle_gamma 90', "_space_group_name_H-M_alt 'P 1'",
                        '_space_group_IT_number 1', 'loop_', '_atom_site_label',
                        '_atom_site_type_symbol', '_atom_site_fract_x',
                        '_atom_site_fract_y', '_atom_site_fract_z',
                        '_atom_site_occupancy' ] + sites ) + '\n' + extra

def show( descr, lc ):
    print( f'{descr}:' )
    for a in lc.atoms:
        print( '  ', a['cif_labels'], a['composition'], 'occupancy',
               a['occupancy'], 'positions', len( a['equivalent_positions'] ) )
    return lc

def main():
    lc = show( 'isotopes', CIFLoader( cif( [ 'D1 D 0 0 0 1',
                                             'T1 T 0.5 0.5 0.5 1',
                                             'O1 O 0.25 0.25 0.25 1' ] ),
                                      quiet = True, refine_with_spglib = False ) )
    comps = sorted( a['composition'][0][1] for a in lc.atoms )
    assert comps == [ 'H2', 'H3', 'O' ], comps
    #Tritium from the type symbol alone (gemmi does not know it), with charge:
    lc = show( 'tritium type symbol', CIFLoader( cif( [ 'H1 T 0 0 0 1',
                                                        'H2 T+ 0.5 0.5 0.5 1' ] ),
                                                 quiet = True,
                                                 refine_with_spglib = False ) )
    assert [ a['composition'] for a in lc.atoms ] == [ ( ( 1.0, 'H3' ), ) ]
    #Element deduced from the label when the type symbol is missing:
    lc = show( 'element from label', CIFLoader( cif( [ 'CaM1 . 0 0 0 1',
                                                       'OA1 . 0.5 0.5 0.5 1' ] ),
                                                quiet = True,
                                                refine_with_spglib = False ) )
    assert sorted( a['composition'][0][1] for a in lc.atoms ) == [ 'Ca', 'O' ]
    #A stated atom count which is not a number is not used:
    extra = 'loop_\n_atom_type_symbol\n_atom_type_number_in_cell\nNa ?\nCl 1\n'
    lc = show( 'atom count "?"', CIFLoader( cif( [ 'Na1 Na 0 0 0 1',
                                                   'Cl1 Cl 0.5 0.5 0.5 1' ], extra ),
                                            quiet = True, refine_with_spglib = False ) )
    #A stated atom count which is wrong is an error:
    extra = 'loop_\n_atom_type_symbol\n_atom_type_number_in_cell\nNa 1\nCl 2\n'
    try:
        CIFLoader( cif( [ 'Na1 Na 0 0 0 1', 'Cl1 Cl 0.5 0.5 0.5 1' ], extra ),
                   quiet = True, refine_with_spglib = False )
    except NC.NCBadInput as e:
        print( 'wrong atom count: rejected:', e )
    else:
        raise SystemExit( 'wrong atom count not rejected' )
    #Mixed site, with one of its elements remapped: still one atom per site.
    lc = show( 'mixed site', CIFLoader( cif( [ 'Al1 Al 0 0 0 0.75',
                                               'Cr1 Cr 0 0 0 0.25',
                                               'O1 O 0.5 0.5 0.5 1' ],
                                             "_publ_section_title 'Test'\n" ),
                                        quiet = True ) )
    descr0 = list( lc.extracted_description )
    c = lc.create_ncmat_composer( quiet = True, remap = [ ( 'Cr', 'Fe' ) ] )
    lc.create_ncmat_composer( quiet = True, remap = [ ( 'Cr', 'Fe' ) ] )
    assert list( lc.extracted_description ) == descr0#(loader not modified)
    txt = c.create_ncmat()
    print( "\n".join( ln for ln in txt.splitlines()
                      if ln.startswith( ( '@', '  ' ) ) ) )
    info = NC.directLoad( txt ).info
    assert abs( sum( ai.averageNumberOfAtoms for ai in info.atominfos ) - 2.0 ) < 1e-12
    print( 'atoms per cell after remap:', sum( ai.averageNumberOfAtoms for ai in info.atominfos ) )
    #Remap with a stated formula which can not be decoded (not checked):
    lc = CIFLoader( cif( [ 'Al1 Al 0 0 0 0.75', 'Cr1 Cr 0 0 0 0.25' ],
                         "_chemical_formula_sum 'Xq7 Zz'\n" ), quiet = True )
    c = lc.create_ncmat_composer( quiet = True, remap = [ ( 'Cr', 'Fe' ) ] )
    txt = c.create_ncmat()
    NC.directLoad( txt )
    print( 'undecodable formula with remap:',
           [ ln.strip() for ln in txt.splitlines()
             if ' is ' in ln and not ln.startswith( '#' ) ] )

if __name__ == '__main__':
    main()
