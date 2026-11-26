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

# Tests of cifutils.spacegroup_settings: each listed setting name,
# Hermann-Mauguin and Hall symbol selects its setting as override_spacegroup,
# except Hall symbols which also read as Hermann-Mauguin symbols of other
# settings (rejected as ambiguous unless prefixed) and those shared by two
# settings (resolved by the cell, which here fits both).

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.cifutils import spacegroup_settings
from NCrystalDev.misc import evaluate_query

#(an atom at a general position, unoccupied enough to be no split site in
#any setting; the cell need not fit):
cif = ( "data_x\n_cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n"
        "_atom_site_label Fe1\n_atom_site_fract_x 0.1234\n"
        "_atom_site_fract_y 0.2345\n_atom_site_fract_z 0.3456\n"
        "_atom_site_occupancy 0.01\n" )

def setting_of( override ):
    try:
        return evaluate_query( [ 'cif', 'structure', 'permissive', 'RAWDATA',
                                 f'override_spacegroup={override}',
                                 'allow_cell_misfit=1', 'RAWDATA', 'x.cif',
                                 'RAWDATA', cif ] )['spacegroup_setting']
    except NC.NCBadInput as e:
        return f'rejected: {e}'

def main():
    rows = spacegroup_settings()
    assert len( rows ) == 530 and len( set( rows ) ) == 530
    for r in rows[:3] + rows[-2:]:
        print( r )
    nexc = 0
    for s, hm, hall in rows:
        for kind, o in ( ( 'setting', s ), ( 'H-M', hm ), ( 'Hall', hall ) ):
            r = setting_of( o )
            if r != s:
                nexc += 1
                print( f'{s} {kind} "{o}": {r}' )
                #(the prefixed symbol is always accepted, with the cell
                #resolving the shared Hall symbols of number 68):
                rp = setting_of( f'{kind}:{o}' )
                assert kind == 'Hall' and ( rp == s or '68:' in rp ), rp
    print( f'{len( rows )} settings, {nexc} exceptions' )
    #The same table from the command line tool:
    import contextlib
    import io

    import NCrystalDev.cli as nc_cli
    buf = io.StringIO()
    with contextlib.redirect_stdout( buf ):
        nc_cli.run( 'cif2ncmat', '--list-spacegroups' )
    lines = buf.getvalue().splitlines()
    assert len( lines ) == 531
    print( 'ncrystal_cif2ncmat --list-spacegroups:' )
    print( '\n'.join( lines[:4] + [ '...' ] + lines[-2:] ) )

if __name__ == '__main__':
    main()
