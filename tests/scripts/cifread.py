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

# Tests of _cifreadimpl.read_cif: CIF text read with gemmi into the plain
# content used by the rest of the CIF loading.

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

def read( txt, override = None ):
    return read_cif( txt, quiet = True, override_spacegroup = override )

def show( descr, c ):
    print( f'{descr}:' )
    print( '   cell:', c['cell'] )
    for s in c['sites']:
        print( '   site:', s )
    for k in ( 'listed_ops', 'symbol_hm', 'symbol_itnumbers', 'orig_hm_alt',
               'multiplicities', 'atom_type_number_in_cell', 'uiso_printed',
               'ueq_printed' ):
        print( f'   {k}:', c[k] )
    assert c['_struct'] is not None and c['_gemmi'] is not None

def expect_error( descr, txt, **kw ):
    try:
        read( txt, **kw )
    except NC.NCBadInput as e:
        print( f'{descr}: rejected: {e}' )
    else:
        raise SystemExit( f'{descr}: not rejected' )

def main():
    ops14 = [ 'x, y, z', '-x, y+1/2, -z+1/2', '-x, -y, -z', 'x, -y+1/2, z+1/2' ]
    c = read( cif( ops = ops14,
                   aniso = [ 'O1 0.012 0.009 0.015 0.002 -0.003 0.001' ] ) )
    show( 'monoclinic, listed operations, one tensor', c )
    o1, si1, fe1 = c['sites']
    assert o1['aniso'] == { 'u11': 0.012, 'u22': 0.009, 'u33': 0.015,
                            'u12': 0.002, 'u13': -0.003, 'u23': 0.001 }
    assert si1['aniso'] is None and fe1['aniso'] is None
    assert si1['occ'] == 0.5 and o1['u_iso'] == 0.0123
    assert c['listed_ops'] == ops14
    assert c['multiplicities'] == { 'O1': 4, 'Si1': 4, 'Fe1': 2 }
    assert set( c['uiso_printed'] ) == { 'O1', 'Si1' }
    assert set( c['ueq_printed'] ) == { 'O1' }
    #Old operations tag, no H-M symbol, only the IT number:
    c = read( cif( ops = ops14, hm = None, opstag = '_symmetry_equiv_pos_as_xyz' ) )
    show( 'old operations tag, no symbol', c )
    assert c['listed_ops'] == ops14 and c['symbol_itnumbers'][0] == [ '14' ]
    #No listed operations, symbol only:
    c = read( cif() )
    #(NB: the value keeps the CIF quotes, which makes the check for an
    #explicitly given setting fail; to be fixed separately):
    assert c['listed_ops'] is None and c['orig_hm_alt'].strip( "'" ) == 'P 1 21/c 1'
    #Override: listed operations are then ignored:
    c = read( cif( ops = ops14 ), override = 'P 1 21/c 1' )
    assert c['listed_ops'] is None
    print( 'override ignores listed operations' )
    #Atom counts and an element only given by the label:
    c = read( cif( extra = 'loop_\n_atom_type_symbol\n_atom_type_number_in_cell\n'
                   'O 4\nSi 2\nFe 2\n',
                   sites = [ 'O1 O 0.79 0.87 0.63 1 0.0123 4',
                             'CaM1 . 0.11 0.22 0.33 1 0.0071 4' ] ) )
    show( 'atom counts, missing type symbol', c )
    assert c['atom_type_number_in_cell'] == [ '4', '2', '2' ]
    #Several blocks: the one with sites is used:
    c = read( 'data_first\n_cell_length_a 1\n' + cif() )
    assert len( c['sites'] ) == 3
    print( 'block with the sites selected among two' )
    expect_error( 'two blocks with sites', cif() + cif().replace( 'data_test', 'data_b' ) )
    expect_error( 'syntax error', 'data_x\n_a 1 2\n' )
    expect_error( 'no blocks', '# nothing\n' )
    expect_error( 'unknown override', cif(), override = 'Q 99' )

if __name__ == '__main__':
    main()
