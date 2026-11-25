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

# NEEDS: numpy spglib

# Tests the temperature associated with Uiso values when converting CIF data to
# NCMAT: explicitly chosen, from the CIF temperature fields, or the fallback,
# and the comments, warnings and errors in default and permissive mode.

import contextlib
import io

import NCrystalDev as NC
import NCrystalDev.cifutils as nc_cifutils
from NCrystalDev.msg import MsgCtrl


def cifdata( fields, with_uiso = True ):
    lines = [ 'data_test', "_symmetry_space_group_name_H-M 'F m -3 m'",
              '_cell_length_a 5.64', '_cell_length_b 5.64', '_cell_length_c 5.64',
              '_cell_angle_alpha 90', '_cell_angle_beta 90', '_cell_angle_gamma 90' ]
    lines += [ f'{k} {v}' for k, v in fields.items() ]
    lines += [ 'loop_', '_atom_site_label', '_atom_site_type_symbol',
               '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z' ]
    if with_uiso:
        lines.append( '_atom_site_U_iso_or_equiv' )
    lines += [ 'Na1 Na 0 0 0' + ( ' 0.0123' if with_uiso else '' ),
               'Cl1 Cl 0.5 0.5 0.5' + ( ' 0.0105' if with_uiso else '' ) ]
    return '\n'.join( lines ) + '\n'

def test( descr, fields, *, chosen = None, permissive = False, with_uiso = True ):
    print( f'==> {descr} (fields: {fields}, chosen: {chosen},'
           f' permissive: {permissive}{"" if with_uiso else ", no Uiso"})' )
    with MsgCtrl( spy = ['warning'], block = ['warning'] ) as mc:
        try:
            with contextlib.redirect_stdout( io.StringIO() ):
                lc = nc_cifutils.CIFLoader( cifdata( fields, with_uiso ),
                                            permissive = permissive )
                txt = lc.create_ncmat( uiso_temperature = chosen )
        except NC.NCBadInput as e:
            print( f'   ERROR: {e}' )
            txt = None
    for msg, _ in mc.warnings:
        print( f'   warning: {msg}' )
    if txt is None:
        return
    lines = txt.splitlines()
    i0 = next( ( i for i, e in enumerate( lines ) if 'Uiso (msd)' in e ), None )
    if i0 is None:
        print( '   (no comment about the Uiso temperature)' )
    else:
        for line in lines[i0:]:
            if not line.startswith( '#' ) or line.strip() == '#':
                break
            print( f'   {line}' )
    for line in lines:
        if 'Fallback (dummy) Debye' in line:
            print( f'   {line}' )
    for line in lines:
        if line.strip().startswith( ( 'msd ', 'debye_temp ' ) ):
            print( f'   {line.strip()}' )

def main():
    a, c = '_diffrn_ambient_temperature', '_cell_measurement_temperature'
    test( 'from the CIF', { a: '100(2)', c: '96(2)' } )
    test( 'from the CIF, only cell field', { c: '150' } )
    test( 'no temperature', {} )
    test( 'absent values', { a: '?', c: '.' } )
    test( 'inconsistent fields', { a: '100(2)', c: '293(2)' } )
    test( 'inconsistent fields', { a: '100(2)', c: '293(2)' }, permissive = True )
    test( 'invalid field', { a: '100(2)', c: 'room' } )
    test( 'invalid field', { a: '100(2)', c: 'room' }, permissive = True )
    test( 'only invalid field', { a: '5000' }, permissive = True )
    test( 'zero kelvin field', { a: '0' } )
    test( 'chosen', { a: '100(2)' }, chosen = 300 )
    test( 'chosen, no fields', {}, chosen = 300 )
    test( 'chosen, inconsistent fields', { a: '100(2)', c: '293(2)' }, chosen = 300 )
    test( 'chosen, invalid', {}, chosen = 0.05 )
    test( 'no Uiso values', { a: '100(2)', c: '293(2)' }, with_uiso = False )

if __name__ == '__main__':
    main()
