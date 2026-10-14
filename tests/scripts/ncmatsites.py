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

# NCMAT v8: @SPACEGROUP with a setting choice ("227:2", "14:c1", ...) and the
# @ATOMSITES section (one position per orbit, expanded at load time).

import re

import NCrystalDev as NC
import NCTestUtils.enable_fpe  # noqa F401
import numpy as np
from NCTestUtils.common import ensure_error


def as_sites_text( name ):
    """The stdlib file rewritten as NCMAT v8 with @ATOMSITES (from the sites of
    the loaded Info) and the detected setting in @SPACEGROUP."""
    info = NC.createInfo( f'{name};inelas=0' )
    setting = info.structure_info['spacegroup_setting']
    lines = []
    for ai in info.atominfos:
        lbl = ai.atomData.displayLabel()
        for s in ai.sites:
            p = s['position']
            lines.append( f'  {lbl} {p[0]!r} {p[1]!r} {p[2]!r}' )
    text = NC.createTextData( name ).rawData
    text = re.sub( r'^NCMAT v\d+', 'NCMAT v8', text, count = 1, flags = re.MULTILINE )
    text = re.sub( r'^@SPACEGROUP\n.*?\n(?=@)', f'@SPACEGROUP\n  {setting}\n', text, flags = re.MULTILINE|re.DOTALL )
    text = re.sub( r'^@ATOMPOSITIONS\n.*?\n(?=@)', '@ATOMSITES\n' + '\n'.join(lines) + '\n', text, flags = re.MULTILINE|re.DOTALL )
    return text, info

def same_positions( info1, info2 ):
    def pset( info ):
        res = set()
        for ai in info.atominfos:
            for p in ai.positions:
                res.add( ( ai.atomData.displayLabel(), ) + tuple( round( x % 1.0, 6 ) % 1.0 for x in p ) )
        return res
    return pset( info1 ) == pset( info2 )

def main():
    egrid = np.geomspace( 1e-3, 5.0, 300 )
    for name in ( 'Al_sg225.ncmat', 'Si_sg227.ncmat', 'Al2O3_sg167_Corundum.ncmat',
                  'PbO-alpha_sg129_Litharge.ncmat', 'Mg_sg194.ncmat' ):
        text, info_orig = as_sites_text( name )
        nsites = sum( len( ai.sites ) for ai in info_orig.atominfos )
        NC.registerInMemoryFileData( 'sites.ncmat', text )
        info = NC.createInfo( 'sites.ncmat;inelas=0' )
        assert same_positions( info, info_orig ), name
        assert info.structure_info['spacegroup_setting'] == info_orig.structure_info['spacegroup_setting']
        def rounded( info ):
            return [ [ ( tuple( round(x,9) for x in s['position'] ), s['multiplicity'] )
                       for s in ai.sites ] for ai in info.atominfos ]
        assert rounded( info ) == rounded( info_orig ), name
        xs = NC.createScatter( 'sites.ncmat;inelas=0' ).xsect( egrid )
        xs_orig = NC.createScatter( f'{name};inelas=0' ).xsect( egrid )
        assert np.allclose( xs, xs_orig, rtol = 1e-12, atol = 0.0 ), name
        print( f'{name}: {nsites} sites -> {info.structure_info["n_atoms"]} atoms'
               f' (setting {info.structure_info["spacegroup_setting"]}), same positions and xs' )

    #A monoclinic structure with its 2-fold axis along c, previously rejected
    #(only unique axis b was supported), loads with the declared setting 3:c:
    p2c = """NCMAT v8
@CELL
  lengths 5 6 7
  angles 90 90 90
@SPACEGROUP
  3:c
@ATOMSITES
  Al 0.1 0.2 0.3
@DYNINFO
  element Al
  fraction 1
  type vdosdebye
  debye_temp 400
"""
    NC.registerInMemoryFileData( 'p2c.ncmat', p2c )
    info = NC.createInfo( 'p2c.ncmat;inelas=0;dcutoff=1.0' )
    print( 'P2 with 2-fold along c (3:c):', info.structure_info['n_atoms'], 'atoms,',
           info.atominfos[0].sites, 'setting', info.structure_info['spacegroup_setting'] )
    assert info.structure_info['n_atoms'] == 2
    fam = { e.hkl_label: ( e.mult, sorted( zip( map(int,e.h), map(int,e.k), map(int,e.l) ) ) )
            for e in info.hklObjects() }
    for lbl, (mult, eqv) in fam.items():
        #Laue class 2/m with the 2-fold along c: (h,k,l) ~ (-h,-k,l) (and Friedel):
        for (h,k,ll) in eqv:
            assert ( -h, -k, ll ) in eqv or ( h, k, -ll ) in eqv, ( lbl, eqv )
    print( f'  {len(fam)} hkl families, all with the 2/m (unique axis c) equivalences' )

    #Errors:
    base = p2c.replace( '3:c', '{sg}' ).replace( '@ATOMSITES', '{sec}' )
    def load( sg = '3:c', sec = '@ATOMSITES', version = 8, extra = '' ):
        t = base.format( sg = sg, sec = sec ).replace( 'NCMAT v8', f'NCMAT v{version}' ) + extra
        NC.registerInMemoryFileData( 'bad.ncmat', t )
        return NC.createInfo( 'bad.ncmat;inelas=0' )
    with ensure_error( NC.NCBadInput ):
        load( version = 7 )#settings and @ATOMSITES need v8
    with ensure_error( NC.NCBadInput ):
        load( sg = '3:x' )#no such setting
    with ensure_error( NC.NCBadInput ):
        load( extra = '@ATOMPOSITIONS\n  Al 0.1 0.2 0.3\n' )#both sections
    with ensure_error( NC.NCBadInput ):
        #sites need a space group:
        NC.registerInMemoryFileData( 'bad.ncmat', p2c.replace( '@SPACEGROUP\n  3:c\n', '' ) )
        NC.createInfo( 'bad.ncmat;inelas=0' )
    #With sites, any setting generates a consistent structure, so the setting
    #must be right by itself; with positions, it is checked against them:
    with ensure_error( NC.NCBadInput ):
        NC.registerInMemoryFileData( 'bad.ncmat', p2c.replace( '3:c', '3' ).replace(
            '@ATOMSITES\n  Al 0.1 0.2 0.3\n', '@ATOMPOSITIONS\n  Al 0.1 0.2 0.3\n  Al 0.9 0.8 0.3\n' ) )
        NC.createInfo( 'bad.ncmat;inelas=0' )
    print( 'Invalid inputs rejected' )

    #NCMAT v8 settings policy: a bare number is only allowed for the 140
    #numbers with a single setting; the 90 others must declare the setting
    #explicitly (incl. the orthorhombic standard axis order, 17:abc).
    #Pre-v8 files deduce the setting from the positions. Error messages list
    #the settings:
    def errmsg( f ):
        try:
            f()
        except NC.NCBadInput as e:
            return str(e)
        raise RuntimeError('did not fail')
    diamond_o2 = """NCMAT v{v}
@CELL
  cubic 5.43
@SPACEGROUP
  {sg}
@ATOMPOSITIONS
  Si 1/8 1/8 1/8
  Si 1/8 5/8 5/8
  Si 5/8 1/8 5/8
  Si 5/8 5/8 1/8
  Si 7/8 7/8 7/8
  Si 7/8 3/8 3/8
  Si 3/8 7/8 3/8
  Si 3/8 3/8 7/8
@DYNINFO
  element Si
  fraction 1
  type vdosdebye
  debye_temp 400
"""
    def load_diamond( v, sg ):
        NC.registerInMemoryFileData( 'd.ncmat', diamond_o2.format( v = v, sg = sg ) )
        return NC.createInfo( 'd.ncmat;inelas=0' )
    info = load_diamond( 7, '227' )
    print( 'v7, bare 227, positions in origin choice 2: setting deduced as',
           info.structure_info['spacegroup_setting'] )
    assert info.structure_info['spacegroup_setting'] == '227:2'
    msg = errmsg( lambda : load_diamond( 8, '227' ) )
    print( 'v8, bare 227:', msg )
    assert '227:1' in msg and '227:2' in msg
    msg = errmsg( lambda : load_diamond( 8, '227:3' ) )
    print( 'v8, 227:3:', msg )
    assert '227:1 (default), 227:2' in msg
    msg = errmsg( lambda : load_diamond( 8, '227:1' ) )
    print( 'v8, 227:1 with positions in origin choice 2:', msg )
    assert 'consistent with the setting 227:2' in msg
    info = load_diamond( 8, '227:2' )
    print( 'v8, 227:2: loads with setting', info.structure_info['spacegroup_setting'] )
    assert info.structure_info['spacegroup_setting'] == '227:2'

    #Orthorhombic: the standard axis order is declared as 17:abc (17:cab etc.
    #are the permutations); a declared 17:abc with positions actually in
    #17:cab is detected from the positions (the cell cannot tell them apart):
    from NCrystalDev.misc import evaluate_query
    def sg17( sgspec, declared = '17:abc' ):
        pos = evaluate_query( [ 'util', 'atompositions', 'strict', sgspec, 'Al,0.1,0.2,0.3' ] )
        text = """NCMAT v8
@CELL
  lengths 5 6 7
  angles 90 90 90
@SPACEGROUP
  """ + declared + """
@ATOMPOSITIONS
""" + ''.join( f'  Al {x!r} {y!r} {z!r}\n' for lbl,(x,y,z) in pos ) + """@DYNINFO
  element Al
  fraction 1
  type vdosdebye
  debye_temp 400
"""
        NC.registerInMemoryFileData( 'sg17.ncmat', text )
        return NC.createInfo( 'sg17.ncmat;inelas=0' )
    info = sg17( '17:abc' )
    print( 'v8, 17:abc with positions in the standard axis order: setting',
           info.structure_info['spacegroup_setting'] )
    assert info.structure_info['spacegroup_setting'] == '17:abc'
    msg = errmsg( lambda : sg17( '17:cab' ) )
    print( 'v8, 17:abc declared with positions in 17:cab:', msg )
    assert 'declared space group setting is 17:abc' in msg
    assert 'consistent with the setting 17:cab' in msg
    assert '17:cab (consistent with the positions)' in msg
    msg = errmsg( lambda : sg17( '17', declared = '17' ) )
    print( 'v8, bare 17:', msg )
    assert '17:abc (default), 17:cab, 17:bca' in msg

    #Rhombohedral numbers must declare the axes in v8:
    hexcell = """NCMAT v8
@CELL
  lengths 4.76 4.76 12.99
  angles 90 90 120
@SPACEGROUP
  {sg}
@ATOMSITES
  Al 0 0 0.35
@DYNINFO
  element Al
  fraction 1
  type vdosdebye
  debye_temp 400
"""
    def load_r( sg ):
        NC.registerInMemoryFileData( 'r.ncmat', hexcell.format( sg = sg ) )
        return NC.createInfo( 'r.ncmat;inelas=0' )
    msg = errmsg( lambda : load_r( '167' ) )
    print( 'v8, bare 167:', msg )
    assert '167:h (default), 167:r' in msg
    msg = errmsg( lambda : load_r( '167:r' ) )
    print( 'v8, 167:r with hexagonal cell:', msg )
    info = load_r( '167:h' )
    print( 'v8, 167:h:', info.structure_info['n_atoms'], 'atoms, setting',
           info.structure_info['spacegroup_setting'] )
    assert info.structure_info['n_atoms'] == 12

if __name__ == '__main__':
    main()
