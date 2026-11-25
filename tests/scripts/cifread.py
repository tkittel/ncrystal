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


# NEEDS: gemmi

# Tests of _cifreadimpl.read_cif: CIF text read with gemmi into a CIFContent
# (printed strings) used by the rest of the CIF loading.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev as NC
from NCrystalDev._cifreadimpl import read_cif

def cif( *, hm = 'P 1 21/c 1', itno = '14', ops = None, extra = '',
         sites = None, aniso = None, cell = '5.1 6.3 7.2 90 105 90',
         opstag = '_space_group_symop_operation_xyz' ):
    a, b, c, al, be, ga = cell.split()
    lines = [ 'data_test', f'_cell_length_a {a}', f'_cell_length_b {b}',
              f'_cell_length_c {c}', f'_cell_angle_alpha {al}',
              f'_cell_angle_beta {be}', f'_cell_angle_gamma {ga}' ]
    if hm is not None:
        lines.append( f"_space_group_name_H-M_alt '{hm}'" )
    if itno is not None:
        lines.append( f'_space_group_IT_number {itno}' )
    if ops:
        lines += [ 'loop_', opstag ] + [ f"'{o}'" for o in ops ]
    lines += [ 'loop_', '_atom_site_label', '_atom_site_type_symbol',
               '_atom_site_fract_x', '_atom_site_fract_y',
               '_atom_site_fract_z', '_atom_site_occupancy',
               '_atom_site_U_iso_or_equiv',
               '_atom_site_symmetry_multiplicity' ]
    lines += sites or [ 'O1 O 0.79 0.87 0.63 1 0.0123(4) 4',
                        'Si1 Si 0.11 0.22 0.33 0.5 0.0071 4',
                        'Fe1 Fe 0 0 0 1 ? 2' ]
    if aniso:
        lines += [ 'loop_', '_atom_site_aniso_label' ] + [
            f'_atom_site_aniso_U_{k}' for k in ( '11','22','33','12','13','23' ) ]
        lines += aniso
    return '\n'.join( lines ) + '\n' + extra

def read( txt ):
    return read_cif( txt, quiet = True )

def show( descr, c ):
    print( f'{descr}:' )
    print( '   cell:', c.cell )
    for s in c.sites:
        print( '   site:', s )
    for k in ( 'listed_ops', 'hall', 'hm', 'itnumber', 'formula',
               'atom_type_number_in_cell' ):
        print( f'   {k}:', getattr( c, k ) )

def expect_error( descr, txt, detail = True ):
    try:
        read( txt )
    except NC.NCBadInput as e:
        #(without the details from gemmi, which might change between versions):
        print( f'{descr}: rejected: {e if detail else str(e).split("(")[0].rstrip()}' )
    else:
        raise SystemExit( f'{descr}: not rejected' )

def main():
    ops14 = [ 'x, y, z', '-x, y+1/2, -z+1/2', '-x, -y, -z', 'x, -y+1/2, z+1/2' ]
    c = read( cif( ops = ops14, extra = "_chemical_formula_sum 'Fe O Si'\n",
                   aniso = [ 'O1 0.012(2) 0.009 0.015 0.002 -0.003 0.001',
                             'Si1 0.01 ? 0.01 0 0 0' ] ) )
    show( 'monoclinic, listed operations, tensors', c )
    o1, si1, fe1 = c.sites
    #(printed strings, with uncertainties):
    assert o1['uij'] == [ '0.012(2)', '0.009', '0.015', '0.002', '-0.003', '0.001' ]
    assert o1['uiso'] == '0.0123(4)' and si1['occ'] == '0.5'
    #(tensors with missing values are ignored; null values are None):
    assert si1['uij'] is None and fe1['uij'] is None and fe1['uiso'] is None
    assert c.listed_ops == ops14 and c.formula == 'Fe O Si'
    assert c.cell == { 'a': '5.1', 'b': '6.3', 'c': '7.2', 'alpha': '90',
                          'beta': '105', 'gamma': '90' }
    #Deprecated _symmetry_* tags, no H-M symbol, only the IT number:
    c = read( cif( ops = ops14, hm = None, itno = None,
                   opstag = '_symmetry_equiv_pos_as_xyz',
                   extra = "_symmetry_Int_Tables_number 14\n"
                   "_symmetry_space_group_name_Hall '-P 2ybc'\n" ) )
    show( 'old tags, no H-M symbol', c )
    assert c.listed_ops == ops14 and c.itnumber == '14'
    assert c.hall == '-P 2ybc' and c.hm is None
    #No listed operations:
    assert read( cif() ).listed_ops is None
    #Atom counts, a missing type symbol, B instead of U:
    c = read( cif( extra = 'loop_\n_atom_type_symbol\n_atom_type_number_in_cell\n'
                   'O 4\nSi 2\nFe ?\n',
                   sites = [ 'O1 O 0.79 0.87 0.63 1 0.0123 4',
                             'CaM1 . 0.11 0.22 0.33 1 0.0071 4' ] ).replace(
                                 '_atom_site_U_iso_or_equiv',
                                 '_atom_site_B_iso_or_equiv' ) )
    show( 'atom counts, missing type symbol, Biso', c )
    assert c.atom_type_number_in_cell == [ '4', '2', None ]
    assert c.sites[1]['type_symbol'] is None
    assert c.sites[0]['biso'] == '0.0123' and c.sites[0]['uiso'] is None
    #Several blocks: the one with sites is used:
    c = read( 'data_first\n_cell_length_a 1\n' + cif() )
    assert len( c.sites ) == 3
    print( 'block with the sites selected among two' )
    #Missing cell angles are 90 degrees, but missing lengths are an error:
    c = read( '\n'.join( ln for ln in cif().splitlines()
                         if not ln.startswith( '_cell_angle' ) ) )
    assert all( c.cell[k] == '90' for k in ( 'alpha', 'beta', 'gamma' ) )
    print( 'missing cell angles: 90 degrees' )
    expect_error( 'missing cell length',
                  '\n'.join( ln for ln in cif().splitlines()
                             if not ln.startswith( '_cell_length_b' ) ) )
    expect_error( 'two blocks with sites', cif() + cif().replace( 'data_test', 'data_b' ) )
    expect_error( 'syntax error', 'data_x\n_a 1 2\n', detail = False )
    expect_error( 'no blocks', '# nothing\n' )
    #Aniso rows go to the sites with their label in order (labels reused):
    c = read( cif( sites = [ 'O1 O 0.1 0.2 0.3 1 ? 4', 'O1 O 0.4 0.2 0.3 1 ? 4',
                             'Si1 Si 0.6 0.2 0.3 1 ? 4' ],
                   aniso = [ 'O1 0.011 0.012 0.013 0 0 0',
                             'Si1 0.021 0.022 0.023 0 0 0',
                             'O1 0.031 0.032 0.033 0 0 0' ] ) )
    print( 'reused labels, tensors:', [ s['uij'] and s['uij'][0]
                                        for s in c.sites ] )
    assert [ s['uij'][0] for s in c.sites ] == [ '0.011', '0.031', '0.021' ]
    from NCrystalDev.cifutils import CIFLoader
    #All-zero tensors are placeholders for no tensor (the Uiso is used), and
    #errors about tensors name the CIF site:
    for descr, aniso in ( ( 'all-zero tensor', '0 0 0 0 0 0' ),
                          ( 'tensor not positive definite',
                            '0.01 -0.01 0.01 0 0 0' ) ):
        try:
            lc = CIFLoader( cif( sites = [ 'Fe1 Fe 0 0 0 1 0.012 2' ],
                                 aniso = [ f'Fe1 {aniso}' ], ops = ops14 ),
                            quiet = True, refine_with_spglib = False )
        except NC.NCBadInput as e:
            print( f'{descr}: rejected: {e}' )
        else:
            s, = lc.sites
            print( f'{descr}: uiso {s["uiso"]}, tensor {s["aniso"]}' )
    #Missing coordinates are rejected when loading (but unoccupied sites are
    #ignored):
    from NCrystalDev.cifutils import CIFLoader
    for sites in ( [ 'O1 O ? ? ? 1 ? 4' ],
                   [ 'O1 O 0.1 0.2 0.3 1 ? 4', 'O2 O ? ? ? 0 ? 4' ] ):
        try:
            lc = CIFLoader( cif( sites = sites, ops = ops14 ), quiet = True,
                            refine_with_spglib = False )
        except NC.NCBadInput as e:
            print( f'sites {sites}: rejected: {e}' )
        else:
            print( f'sites {sites}: loaded {len(lc.sites)} site(s)' )

if __name__ == '__main__':
    main()
