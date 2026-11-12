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

# Anisotropic displacements from CIF data (cifutils): the CIF sites and their
# uij values are kept. Uiso versus Ueq is resolved when the data is read: a
# Uiso differing from the Ueq by more than tol_uij is rejected (or used, in
# permissive mode), and the tensor is always rescaled to the Uiso, so each
# site carries a single magnitude.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev as NC
from NCrystalDev.cifutils import CIFLoader
from NCrystalDev._sgsymimpl import _cif_uij_ueq, _relclose
from NCrystalDev._cifreadimpl import uiso_from_cif, ueq_from_cif

cell = { 'a': 5.1, 'b': 6.3, 'c': 7.2,
         'alpha': 90.0, 'beta': 105.0, 'gamma': 90.0 }
#O site not at the representative NCrystal would choose (0.21,0.13,0.37):
uij_o = { 'u11': 0.012, 'u22': 0.009, 'u33': 0.015,
          'u12': 0.002, 'u13': -0.003, 'u23': 0.001 }
uij_si = { 'u11': 0.006, 'u22': 0.007, 'u33': 0.008,
           'u12': 0.0, 'u13': 0.001, 'u23': 0.0 }
ueq_o = _cif_uij_ueq( uij_o, cell )
ueq_si = _cif_uij_ueq( uij_si, cell )

def aniso_line( lbl, u ):
    return ' '.join( [ lbl ] + [ f'{u[k]:.4f}' for k in ( 'u11','u22','u33',
                                                          'u12','u13','u23' ) ] )

cif = f"""data_test
_cell_length_a {cell['a']}
_cell_length_b {cell['b']}
_cell_length_c {cell['c']}
_cell_angle_alpha {cell['alpha']}
_cell_angle_beta {cell['beta']}
_cell_angle_gamma {cell['gamma']}
_space_group_name_H-M_alt 'P 1 21/c 1'
_space_group_IT_number 14
loop_
_space_group_symop_operation_xyz
'x, y, z'
'-x, y+1/2, -z+1/2'
'-x, -y, -z'
'x, -y+1/2, z+1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
O1 O 0.79 0.87 0.63 {ueq_o:.6f}
Si1 Si 0.11 0.22 0.33 {1.2*ueq_si:.6f}
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
{aniso_line('O1',uij_o)}
{aniso_line('Si1',uij_si)}
"""

def main():
    #Si1 has Uiso=1.2*Ueq: rejected, unless permissive (using the Uiso):
    try:
        CIFLoader( cif )
    except NC.NCBadInput as e:
        print( 'Rejected:', e )
    else:
        raise SystemExit('Uiso/Ueq mismatch not rejected')
    loader = CIFLoader( cif, permissive = True )
    composer = loader.create_ncmat_composer( uiso_temperature = 300 )
    txt = composer.create_ncmat()
    i = txt.index('@ATOMSITES')
    print( txt[i:txt.index('@',i+1)], end = '' )
    info = NC.directLoad( txt, cfg_params = 'inelas=0;temp=300K' ).info
    for ai in info.atominfos:
        lbl = ai.atomData.displayLabel()
        s, = ai.sites
        u = uij_o if lbl == 'O' else uij_si
        ucif = [ u[k] for k in ( 'u11','u22','u33','u12','u13','u23' ) ]
        print( f'{lbl} loaded at {tuple(s["position"])}, uij:',
               ', '.join( f'{x:.5g}' for x in s['uij'] ) )
        if lbl == 'O':
            #(the inversion maps the CIF site to the representative, keeping
            #the tensor, and msd=Uiso=Ueq at 300K):
            assert all( abs( a - b ) < 1e-6 for a, b in zip( s['uij'], ucif ) )

ukeys = ( 'u11','u22','u33','u12','u13','u23' )
cells = { 'orthorhombic': { 'a': 5.1, 'b': 6.3, 'c': 7.2,
                            'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0 },
          'triclinic': { 'a': 5.1, 'b': 6.3, 'c': 7.2,
                         'alpha': 84.0, 'beta': 105.0, 'gamma': 97.0 } }
#Tensor of a single site in P-1 (any tensor is allowed on a general position):
uij_gen = { 'u11': 0.0123, 'u22': 0.0087, 'u33': 0.0154,
            'u12': 0.0021, 'u13': -0.0032, 'u23': 0.0011 }

def p1bar_cif( c, uiso ):
    uiso_str = '?' if uiso is None else repr( uiso )
    return f"""data_test
_cell_length_a {c['a']}
_cell_length_b {c['b']}
_cell_length_c {c['c']}
_cell_angle_alpha {c['alpha']}
_cell_angle_beta {c['beta']}
_cell_angle_gamma {c['gamma']}
_space_group_name_H-M_alt 'P -1'
_space_group_IT_number 2
loop_
_space_group_symop_operation_xyz
'x, y, z'
'-x, -y, -z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
Al1 Al 0.11 0.22 0.33 {uiso_str}
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
{aniso_line('Al1',uij_gen)}
"""

def scenario( cellname, descr, uiso, expect_ok, **kw ):
    c = cells[cellname]
    if uiso is not None:
        import gemmi
        blk = gemmi.cif.read_string( p1bar_cif( c, uiso ) ).sole_block()
        ( ui, ui_err ), = uiso_from_cif( gemmi, blk ).values()
        ( ue, ue_err ), = ueq_from_cif( gemmi, blk ).values()
        print( f'  {descr}: printed Uiso {ui:.6g} (rounding error {ui_err:.3g}),'
               f' Ueq {ue:.6g} (rounding error {ue_err:.3g})' )
    try:
        lc = CIFLoader( p1bar_cif( c, uiso ), quiet = True, **kw )
    except NC.NCBadInput as e:
        print( f'  {descr}: rejected: {e}' )
        assert not expect_ok
        return
    assert expect_ok
    s, = lc.sites
    u = s['aniso']
    #A single magnitude: the Ueq of the (rescaled) tensor, in the final cell,
    #is the uiso:
    assert _relclose( _cif_uij_ueq( u, lc.cellsg ), s['uiso'], 1e-9 )
    if uiso is not None:
        assert s['uiso'] == uiso
    print( f'  {descr}: uiso={s["uiso"]!r} uij:',
           ' '.join( repr( u[k] ) for k in ukeys ) )

def test_uiso_vs_ueq():
    for cellname, c in cells.items():
        ueq = _cif_uij_ueq( uij_gen, c )
        print( f'{cellname} cell, Ueq of the CIF tensor: {ueq:.8g}' )
        def off( rel, ueq = ueq ):
            return float( f'{ueq*(1.0+rel):.6g}' )
        scenario( cellname, 'no Uiso', None, True )
        #(differences explained by rounding of the printed values are fine):
        scenario( cellname, 'Uiso rounded to 4 decimals',
                  float( f'{ueq:.4f}' ), True )
        scenario( cellname, 'Uiso rounded to 3 decimals',
                  float( f'{ueq:.3f}' ), True )
        scenario( cellname, 'Uiso 0.05% off', off( 5e-4 ), True )
        #(the tensor components are printed with 4 decimals, so rounding
        #explains differences of a few permille):
        scenario( cellname, 'Uiso 0.5% off', off( 5e-3 ), True )
        scenario( cellname, 'Uiso 2% off', off( 0.02 ), False )
        scenario( cellname, 'Uiso 2% off, permissive', off( 0.02 ), True,
                  permissive = True )
        scenario( cellname, 'Uiso 2% off, tol_uij=0.03', off( 0.02 ), True,
                  tol_uij = 0.03 )
        scenario( cellname, 'Uiso 20% off', off( 0.2 ), False )
        scenario( cellname, 'Uiso 20% off, tol_uij=0.03', off( 0.2 ), False,
                  tol_uij = 0.03 )

if __name__ == '__main__':
    main()
    test_uiso_vs_ueq()
