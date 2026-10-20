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

def site_errors():
    #Sites given twice (also as a symmetry-equivalent site) or close to a
    #special position, and fractions not matching the generated atoms:
    def load( sites, sg = '225', fracs = None ):
        els = sorted( { line.split()[0] for line in sites } )
        fracs = fracs or { e : f'1/{len(els)}' for e in els }
        t = ( f'NCMAT v8\n@CELL\n  cubic 4.05\n@SPACEGROUP\n  {sg}\n@ATOMSITES\n'
              + ''.join( f'  {line}\n' for line in sites ) )
        for e in els:
            t += ( f'@DYNINFO\n  element {e}\n  fraction {fracs[e]}\n'
                   '  type vdosdebye\n  debye_temp 400\n' )
        NC.registerInMemoryFileData( 'sites.ncmat', t )
        return NC.createInfo( 'sites.ncmat;inelas=0;dcutoff=0.5' )
    for sites, sg in ( ( [ 'Al 1/4 1/4 1/4', 'Al 1/4 1/4 1/4' ], '225' ),
                       ( [ 'Al 0 0 1', 'Al 0 0 0' ], '225' ),
                       ( [ 'Al 0 0 0', 'Cr 0 0 -1' ], '225' ),
                       ( [ 'Al 0 0 0', 'Al 1/2 1/2 0' ], '225' ),
                       ( [ 'Si 1/8 1/8 1/8', 'Si 7/8 3/8 3/8' ], '227:2' ),
                       ( [ 'Al 0 0 0', 'Cr 1/2 1/2 0' ], '225' ),
                       ( [ 'Al 0.003 0 0' ], '225' ) ):
        print( f'Sites in {sg}:', ', '.join( sites ) )
        with ensure_error( NC.NCBadInput ):
            load( sites, sg )
    #Within the tolerance a site is moved onto the special position:
    info = load( [ 'Al 0.0001 0 0' ] )
    print( 'Site Al (0.0001, 0, 0) in 225:', info.structure_info['n_atoms'], 'atoms' )
    assert info.structure_info['n_atoms'] == 4
    info = load( [ 'Al 0 0 0', 'Cr 1/2 1/2 1/2' ] )
    assert info.structure_info['n_atoms'] == 8
    with ensure_error( NC.NCBadInput ):
        load( [ 'Al 0 0 0', 'Cr 1/2 1/2 1/2' ], fracs = { 'Al' : 0.3, 'Cr' : 0.7 } )

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

    #A monoclinic structure with its 2-fold axis along c loads with the
    #declared setting 3:c:
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
    print( 'Invalid inputs rejected' )
    site_errors()
    #With sites, any setting generates a consistent structure, so the setting
    #must be right by itself; with positions and a bare number, the setting is
    #deduced from them:
    NC.registerInMemoryFileData( 'pos.ncmat', p2c.replace( '3:c', '3' ).replace(
        '@ATOMSITES\n  Al 0.1 0.2 0.3\n', '@ATOMPOSITIONS\n  Al 0.1 0.2 0.3\n  Al 0.9 0.8 0.3\n' ) )
    info = NC.createInfo( 'pos.ncmat;inelas=0' )
    print( 'Bare 3 with positions in 3:c: setting',
           info.structure_info['spacegroup_setting'] )
    assert info.structure_info['spacegroup_setting'] == '3:c'

if __name__ == '__main__':
    main()
