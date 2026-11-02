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

# Test of inelas=vdosdebye (selected by inelas=auto for crystalline materials
# with Debye temperatures but no @DYNINFO sections), which must give the same
# cross sections as the corresponding @DYNINFO entries of type vdosdebye, both
# with the next-gen and the legacy (knllux=-1) kernel processing.

import NCrystalDev as NC

_cell = """NCMAT v3
@CELL
 lengths {a} {a} {a}
 angles 90 90 90
@SPACEGROUP
 225
@ATOMPOSITIONS
{pos}@DEBYETEMPERATURE
{tdebye}"""

_fcc = ( (0,0,0), (0,'1/2','1/2'), ('1/2',0,'1/2'), ('1/2','1/2',0) )

def register( name, a, elems, tdebyes ):
    #elems: element on the fcc sites, and optionally on the octahedral sites:
    pos = ''.join( f' {elems[0]} {x} {y} {z}\n' for x,y,z in _fcc )
    if len(elems) == 2:
        pos += ''.join( f' {elems[1]} {x} {y} {z}\n'
                        for x,y,z in ( ('1/2','1/2','1/2'), ('1/2',0,0),
                                       (0,'1/2',0), (0,0,'1/2') ) )
    data = _cell.format( a = a, pos = pos,
                         tdebye = ''.join( f' {e} {t}\n'
                                           for e,t in zip(elems,tdebyes) ) )
    NC.registerInMemoryFileData( f'{name}_debyeonly.ncmat', data )
    dyninfo = ''.join( f'@DYNINFO\n element {e}\n fraction {1/len(elems)}\n'
                       ' type vdosdebye\n' for e in elems )
    NC.registerInMemoryFileData( f'{name}_dyninfo.ncmat', data + dyninfo )

def main():
    register( 'Al', 4.04958, ['Al'], [410.4] )
    register( 'NaCl', 5.64, ['Na','Cl'], [290.0, 250.0] )
    energies = [ 1e-5, 1e-3, 0.025, 0.1, 1.0, 10.0, 100.0 ]
    for name in ( 'Al', 'NaCl' ):
        for knllux in ( 3, 0, -1 ):
            res = []
            for f in ( 'debyeonly', 'dyninfo' ):
                cfg = f'{name}_{f}.ncmat;temp=200K;knllux={knllux};comp=inelas'
                sc = NC.createScatter( cfg )
                res.append( [ sc.xsect( ekin = e ) for e in energies ] )
            assert res[0] == res[1], (name, knllux)
            print(f'{name} knllux={knllux}: '
                  + ' '.join( f'{x:.6g}' for x in res[0] ) )

if __name__ == '__main__':
    main()
