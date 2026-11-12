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


# NEEDS: gemmi numpy

# Tests of _cifsymimpl.decide_symmetry and of the expansion of positions with
# the chosen symmetry (expand_position, image_rotations).

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import warnings
import NCrystalDev as NC
import NCrystalDev._cifsymimpl as sym
from NCrystalDev._cifreadimpl import read_cif
from NCrystalDev._sgsymimpl import _resolve_tolerances

def cif( cell, hm_alt, itno, sites, hm_old = None ):
    a, b, c, al, be, ga = cell
    lines = [ 'data_test', f'_cell_length_a {a}', f'_cell_length_b {b}',
              f'_cell_length_c {c}', f'_cell_angle_alpha {al}',
              f'_cell_angle_beta {be}', f'_cell_angle_gamma {ga}',
              f'_space_group_IT_number {itno}' ]
    if hm_alt:
        lines.append( f'_space_group_name_H-M_alt {hm_alt}' )
    if hm_old:
        lines.append( f'_symmetry_space_group_name_H-M {hm_old}' )
    lines += [ 'loop_', '_atom_site_label', '_atom_site_type_symbol',
               '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z',
               '_atom_site_occupancy' ] + sites
    return '\n'.join( lines ) + '\n'

def decide( txt, refine = True, permissive = False ):
    content = read_cif( txt, quiet = True, override_spacegroup = None )
    with warnings.catch_warnings( record = True ) as wl:
        warnings.simplefilter( 'always' )
        choice = sym.decide_symmetry( content, permissive = permissive,
                                      tols = _resolve_tolerances( None ),
                                      quiet = True, refine_with_spglib = refine,
                                      override_spacegroup = None,
                                      cifdata = txt )
    return content, choice, [ str( w.message ) for w in wl ]

def show( descr, choice, ww ):
    print( f'{descr}: space group {choice["sgnumber"]} ({choice["hm"]}),'
           f' {len( choice["images"] )} images, p1_sites:'
           f' {None if choice["p1_sites"] is None else len( choice["p1_sites"] )},'
           f' expect_sgno: {choice["expect_sgno"]}' )
    for w in ww:
        print( '     warning:', w[:120] )

hexcell = ( 4.9, 4.9, 13.4, 90, 90, 120 )
site_hex = [ 'Fe1 Fe 0 0 0.1 1' ]

def test_setting_ambiguity():
    #R-3 (148) has hexagonal and rhombohedral settings:
    for descr, hm_alt in ( ( 'setting given (quoted)', "'R -3 :H'" ),
                           ( 'setting given (other spacing)', "'R -3:H'" ),
                           ( 'setting given (unquoted)', "R-3:H" ),
                           ( 'setting not given', None ) ):
        _, c, ww = decide( cif( hexcell, hm_alt, 148, site_hex,
                                hm_old = None if hm_alt else "'R -3'" ) )
        show( descr, c, ww )
        assert c['sgnumber'] == 148
        assert any( 'multiple choices' in w for w in ww ) == ( hm_alt is None )

def test_incompatible_cell():
    #Tetragonal space group, but a != b by 0.1%:
    cell = ( 4.0, 4.004, 6.0, 90, 90, 90 )
    txt = cif( cell, "'P 4/m m m'", 123, [ 'Al1 Al 0 0 0 1' ] )
    try:
        decide( txt, refine = False )
    except NC.NCBadInput as e:
        print( 'incompatible cell, no refinement: rejected:', e )
    else:
        raise SystemExit( 'not rejected' )
    _, c, ww = decide( txt )
    show( 'incompatible cell', c, ww )
    assert c['sgnumber'] == 1 and c['expect_sgno'] == 123
    assert c['images'] == [] and len( c['p1_sites'] ) == 1

def maps( op, R, x0, x ):
    #Whether the gemmi operation op has rotation R and maps x0 to x (mod 1):
    if [ [ e / op.DEN for e in row ] for row in op.rot ] != [ list( r ) for r in R ]:
        return False
    y = op.apply_to_xyz( list( x0 ) )
    return all( abs( d - round( d ) ) < 1e-9 for d in ( a - b for a, b in zip( y, x ) ) )

def test_expand_and_rotations():
    content, c, _ = decide( cif( hexcell, "'R -3 :H'", 148,
                                 [ 'Fe1 Fe 0.3333 0.6667 0.1 1',
                                   'O1 O 0.1 0.2 0.3 1' ] ) )
    cell = content['cell']
    f1a = 1.0 / ( ( cell['a'] + cell['b'] + cell['c'] ) / 3.0 )
    for site in content['sites']:
        pos0 = sym.wrap_position( c, site['fract'] )
        with warnings.catch_warnings( record = True ) as wl:
            warnings.simplefilter( 'always' )
            coords = sym.expand_position( c, pos0, f1a )
        print( f'{site["label"]}: {len( coords )} positions, first'
               f' {tuple( round( x, 12 ) for x in coords[0] )}' )
        for w in wl:
            print( '     warning:', str( w.message ) )
        rots = sym.image_rotations( c, coords )
        assert len( rots ) == len( coords )
        assert rots[0] == ( ( 1, 0, 0 ), ( 0, 1, 0 ), ( 0, 0, 1 ) )
        #Each rotation is that of an operation mapping the first to its image:
        for R, x in zip( rots, coords ):
            assert any( maps( op, R, coords[0], x ) for op in c['sg'].operations() )
    #Site 6c of R-3:H (on the 3-fold axis): 6 positions, snapped to 1/3,2/3:
    assert len( sym.expand_position( c, ( 1/3, 2/3, 0.1 ), f1a ) ) == 6
    assert len( sym.expand_position( c, ( 0.1, 0.2, 0.3 ), f1a ) ) == 18

test_setting_ambiguity()
test_incompatible_cell()
test_expand_and_rotations()
