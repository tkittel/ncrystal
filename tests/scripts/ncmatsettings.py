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

# NEEDS: numpy

# NCMAT v8: @SPACEGROUP with a setting choice ("227:2", "14:c1", ...),
# checked against the atom positions.

import NCTestUtils.enable_fpe  # noqa F401
import NCrystalDev as NC
from NCTestUtils.common import ensure_error

#A monoclinic structure with its 2-fold axis along c, which needs the
#setting 3:c (the default 3:b has the 2-fold axis along b):
p2c = """NCMAT v8
@CELL
  lengths 5 6 7
  angles 90 90 90
@SPACEGROUP
  {sg}
@ATOMPOSITIONS
  Al 0.1 0.2 0.3
  Al 0.9 0.8 0.3
@DYNINFO
  element Al
  fraction 1
  type vdosdebye
  debye_temp 400
"""

def load( sg = '3:c', version = 8 ):
    t = p2c.format( sg = sg ).replace( 'NCMAT v8', f'NCMAT v{version}' )
    NC.registerInMemoryFileData( 'p2c.ncmat', t )
    return NC.createInfo( 'p2c.ncmat;inelas=0;dcutoff=1.0' )

def main():
    info = load()
    print( 'P2 with 2-fold along c (3:c):', info.structure_info['n_atoms'],
           'atoms, setting', info.structure_info['spacegroup_setting'] )
    assert info.structure_info['spacegroup_setting'] == '3:c'
    fam = { e.hkl_label: sorted( zip( map(int,e.h), map(int,e.k), map(int,e.l) ) )
            for e in info.hklObjects() }
    for lbl, eqv in fam.items():
        #Laue class 2/m with the 2-fold along c: (h,k,l) ~ (-h,-k,l) (and Friedel):
        for (h,k,ll) in eqv:
            assert ( -h, -k, ll ) in eqv or ( h, k, -ll ) in eqv, ( lbl, eqv )
    print( f'  {len(fam)} hkl families, all with the 2/m (unique axis c) equivalences' )
    with ensure_error( NC.NCBadInput ):
        load( version = 7 )#settings need v8
    with ensure_error( NC.NCBadInput ):
        load( sg = '3:x' )#no such setting
    with ensure_error( NC.NCBadInput ):
        load( sg = '3:b' )#positions lack the symmetry of the declared setting

if __name__ == '__main__':
    main()
