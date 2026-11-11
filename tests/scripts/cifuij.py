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
# uij values are kept, and Uiso values inconsistent with Ueq give a warning.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev as NC
from NCrystalDev.cifutils import CIFLoader
from NCrystalDev._sgsymimpl import _cif_uij_ueq

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
    loader = CIFLoader( cif )
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

if __name__ == '__main__':
    main()
