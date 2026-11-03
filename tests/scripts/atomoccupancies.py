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


# Parsing of @ATOMOCCUPANCIES sections in NCMAT v8 files, and their effect on
# the composition and density.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import NCTestUtils.reprint_escaped_warnings # noqa F401

def nacl( occu_lines, *, fracs = ('1/2','1/2'), version = 8,
          sites = True, extra = '' ):
    if sites:
        atoms = '@ATOMSITES\n  Na 0 0 0\n  Cl 1/2 1/2 1/2\n'
    else:
        atoms = ( '@ATOMPOSITIONS\n'
                  + ''.join( f'  Na {p}\n' for p in ( '0 0 0', '0 1/2 1/2',
                                                      '1/2 0 1/2', '1/2 1/2 0' ) )
                  + ''.join( f'  Cl {p}\n' for p in ( '1/2 1/2 1/2', '1/2 0 0',
                                                      '0 1/2 0', '0 0 1/2' ) ) )
    occu = ( '@ATOMOCCUPANCIES\n' + ''.join( f'  {ll}\n' for ll in occu_lines )
             if occu_lines is not None else '' )
    return f"""NCMAT v{version}
@CELL
  cubic 5.64
@SPACEGROUP
  225
{atoms}{occu}{extra}@DYNINFO
  element Na
  fraction {fracs[0]}
  type vdosdebye
  debye_temp 300
@DYNINFO
  element Cl
  fraction {fracs[1]}
  type vdosdebye
  debye_temp 300
"""

def load( data ):
    return NC.directLoad( data, cfg_params = 'inelas=0' ).info

def show( title, info ):
    comp = ' '.join( f'{a.displayLabel()}={f:.12g}' for f,a in info.composition )
    print( f'{title}: n_positions={info.structure_info["n_positions"]}'
           f' numberdensity={info.numberdensity:.12g} {comp}' )

def expect_error( title, data ):
    try:
        load( data )
    except NC.NCBadInput as e:
        print( f'{title}: NCBadInput: {e}' )
        return
    raise SystemExit(f'{title}: expected NCBadInput')

ref = load( nacl( None ) )
show( 'no occupancies', ref )
for sites in ( True, False ):
    kind = 'sites' if sites else 'positions'
    full = load( nacl( ['Na 1','Cl 1'], sites = sites ) )
    assert full.numberdensity == ref.numberdensity
    assert full.composition[0][0] == ref.composition[0][0]
    part = load( nacl( ['Na 0.9'], fracs = ('9/19','10/19'), sites = sites ) )
    show( f'Na 0.9 ({kind})', part )
    assert abs( part.numberdensity / ref.numberdensity - 7.6/8 ) < 1e-14
    part2 = load( nacl( ['Na 9/10','Cl 2/3'], fracs = ('27/47','20/47'),
                        sites = sites ) )
    show( f'Na 9/10, Cl 2/3 ({kind})', part2 )

#The same chemical element in two roles with different occupancies, separated
#with a generic label (Na on both rock salt positions):
twona = load( """NCMAT v8
@CELL
  cubic 5.64
@SPACEGROUP
  225
@ATOMDB
  X1 is Na
@ATOMSITES
  Na 0 0 0
  X1 1/2 1/2 1/2
@ATOMOCCUPANCIES
  Na 0.9
  X1 1/2
@DYNINFO
  element Na
  fraction 9/14
  type vdosdebye
  debye_temp 300
@DYNINFO
  element X1
  fraction 5/14
  type vdosdebye
  debye_temp 250
""" )
show( 'Na in two roles (0.9 and 1/2)', twona )
assert abs( twona.numberdensity * twona.structure_info['volume'] - 5.6 ) < 1e-12

#Exposure via queries, Python AtomInfo and the dump:
def show_exposure( title, info ):
    si = info.structure_info
    print( f'{title}: n_positions={si["n_positions"]}'
           f' n_atoms_average={si["n_atoms_average"]:.12g}'
           f' occupancies={info._query("occupancies")}' )
    for ai in info.atominfos:
        print( f'  {ai} numberOfPositions={ai.numberOfPositions}'
               f' averageNumberOfAtoms={ai.averageNumberOfAtoms:.12g}'
               f' occupancy={ai.occupancy:.12g}' )
    for line in info.dump_str().splitlines():
        if 'Atoms / unit cell' in line or ' position' in line and '[' in line:
            print( f'  dump: {line.strip()}' )

show_exposure( 'exposure, no occupancies', ref )
show_exposure( 'exposure, Na 9/10, Cl 2/3', part2 )
show_exposure( 'exposure, Na in two roles', twona )
assert ref.structure_info['n_atoms_average'] == 8.0
assert twona.structure_info['n_atoms_average'] == 5.6

#The obsolete AtomInfo.count alias still works, but warns (once):
for _ in range(3):
    _c = part2.atominfos[0].count
    assert _c == part2.atominfos[0].numberOfPositions
    print( 'AtomInfo.count:', _c )

expect_error( 'zero', nacl( ['Na 0'] ) )
expect_error( 'negative', nacl( ['Na -0.1'] ) )
expect_error( 'above one', nacl( ['Na 1.01'] ) )
expect_error( 'twice', nacl( ['Na 0.9','Na 0.9'], fracs = ('9/19','10/19') ) )
expect_error( 'unknown label', nacl( ['K 0.9'] ) )
expect_error( 'missing value', nacl( ['Na'] ) )
expect_error( 'extra value', nacl( ['Na 0.9 1'] ) )
expect_error( 'not a number', nacl( ['Na abc'] ) )
expect_error( 'empty section', nacl( [] ) )
expect_error( 'unweighted fractions (sites)', nacl( ['Na 0.9'] ) )
expect_error( 'unweighted fractions (positions)', nacl( ['Na 0.9'], sites = False ) )
expect_error( 'NCMAT v7', nacl( ['Na 0.9'], version = 7, sites = False,
                                fracs = ('9/19','10/19') ) )
