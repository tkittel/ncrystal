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

# Tests related to the generation of lists of HKL planes.

import math
import re
import time

import NCrystalDev as NC
import NCTestUtils.enable_fpe  # noqa: F401
from NCTestUtils.env import ncsetenv


# Tests selection of a single HKL plane with the FILLHKL_SELECTHKL environment
# variable (for validation plots), with and without space group. The selected
# plane must get the values of its family in the full list of planes.
def test_select():
    def strip_sg( txt ):
        return re.sub( r'@SPACEGROUP\s*\n\s*\d+\s*\n', '', txt )

    def planes( txt, cfg ):
        NC.clearCaches()
        info = NC.directLoad( txt, cfg + ';comp=bragg', doScatter = False,
                              doAbsorption = False ).info
        return [ ( [ (int(a),int(b),int(c)) for a,b,c in zip(e.h,e.k,e.l) ],
                   e.mult, float(e.d), float(e.f2) ) for e in info.hklObjects() ]

    def select( txt, cfg, hkl ):
        ncsetenv( 'FILLHKL_SELECTHKL', hkl )
        try:
            return planes( txt, cfg )
        finally:
            ncsetenv( 'FILLHKL_SELECTHKL', None )

    def errmsg( e ):
        #Remove the build dependent env var prefix:
        return re.sub( r'NCRYSTAL[A-Z]*_FILLHKL', 'NCRYSTAL_FILLHKL', str(e) )

    def test( fn, cfg, hkls ):
        txt = NC.createTextData( f'stdlib::{fn}' ).rawData
        for label, t in ( ( 'SG', txt ), ( 'no SG', strip_sg( txt ) ) ):
            full = planes( t, cfg )
            for hkl in hkls:
                sel = tuple( int(x) for x in hkl.split(',') )
                msel = tuple( -x for x in sel )
                fam = [ p for p in full if sel in p[0] or msel in p[0] ]
                assert len(fam) <= 1
                try:
                    res = select( t, cfg, hkl )
                except NC.NCCalcError as e:
                    assert not fam
                    assert 'is not present' in str(e)
                    print(f'  {fn} [{label}] {hkl}: CalcError ({errmsg(e)[:40]}...)')
                    continue
                if not res:
                    #Plane outside the d-spacing range:
                    assert not fam
                    info = NC.directLoad( t, cfg + ';comp=bragg',
                                          doScatter = False,
                                          doAbsorption = False ).info
                    d = info.dspacingFromHKL( *sel )
                    assert not ( info.hklDLower() <= d <= info.hklDUpper() )
                    print(f'  {fn} [{label}] {hkl}: outside d-spacing range')
                    continue
                assert fam, 'selected plane not in full list'
                assert len(res) == 1
                hkllist, mult, d, f2 = res[0]
                assert hkllist == [ sel ] and mult == 2
                _, _, dref, f2ref = fam[0]
                if label == 'SG':
                    #Same code path and representative, so identical values:
                    assert d == dref and f2 == f2ref
                else:
                    #Family averages over (possibly) fewer members:
                    assert abs( d - dref ) < 1e-12 * dref
                    assert abs( f2 - f2ref ) < 1e-6 * f2ref
                print(f'  {fn} [{label}] {hkl}: d={d:.6g} F2={f2:.6g}')

    test( 'Al_sg225.ncmat', 'dcutoff=0.3',
          [ '1,1,1', '-1,-1,-1', '2,0,0', '0,0,-2', '0,-2,0', '-3,1,1',
            '1,-3,1', '1,1,0', '5,5,5', '10,10,10' ] )
    test( 'GaSe_sg194_GalliumSelenide.ncmat', 'dcutoff=0.4;dcutoffup=6',
          [ '1,0,0', '0,1,0', '-1,1,0', '1,-2,0', '2,-1,3', '0,0,2', '0,0,1',
            '0,0,4', '0,0,-4' ] )
    test( 'Al2O3_sg167_Corundum.ncmat', 'dcutoff=0.5',
          [ '1,0,4', '0,1,-4', '-1,0,-4', '1,1,3', '1,1,-3', '0,0,6' ] )

    #Invalid values:
    t = NC.createTextData( 'stdlib::Al_sg225.ncmat' ).rawData
    for bad in ( '1,1', '1,1,1,1', 'a,b,c', '1.5,1,1', '0,0,0', ' ' ):
        for label, tt in ( ( 'SG', t ), ( 'no SG', strip_sg( t ) ) ):
            try:
                select( tt, 'dcutoff=0.3', bad )
            except NC.NCBadInput as e:
                assert 'FILLHKL_SELECTHKL' in str(e)
                print(f'  invalid "{bad}" [{label}]: BadInput')
            else:
                raise RuntimeError(f'invalid value "{bad}" accepted')

# Tests that physics models work with a single plane selected via the
# FILLHKL_SELECTHKL environment variable (in case it is ever exposed as a cfg
# parameter): single crystals and layered crystals must scatter exactly on the
# selected (h,k,l) and (-h,-k,-l) planes, and powders must scale as d*|F|^2.
def test_selectphys():
    def strip_sg( txt ):
        return re.sub( r'@SPACEGROUP\s*\n\s*\d+\s*\n', '', txt )

    def load( txt, cfg, sel ):
        ncsetenv( 'FILLHKL_SELECTHKL', sel )
        try:
            NC.clearCaches()
            return NC.directLoad( txt, 'comp=bragg;' + cfg, doAbsorption = False )
        finally:
            ncsetenv( 'FILLHKL_SELECTHKL', None )

    def norm( v ):
        n = math.sqrt( sum( x*x for x in v ) )
        return tuple( x / n for x in v )

    def bragg_direction( g, d, wl ):
        """Incident direction satisfying the Bragg condition for the plane with
    normal g (in the lab frame) and d-spacing d."""
        g = norm( g )
        a = ( 0.0, 1.0, 0.0 ) if abs( g[0] ) > 0.9 else ( 1.0, 0.0, 0.0 )
        perp = norm( ( g[1]*a[2] - g[2]*a[1], g[2]*a[0] - g[0]*a[2],
                       g[0]*a[1] - g[1]*a[0] ) )
        s = wl / ( 2 * d )
        c = math.sqrt( 1 - s*s )
        return tuple( -s*a + c*b for a, b in zip( g, perp ) )

    def check_single_crystal( fn, cfg, gdir_lab, wl, selections, layered = False ):
        """For each selected plane, gdir_lab(hkl) gives its normal in the lab
    frame. Scattering at the exact Bragg geometry must have momentum transfer
    along the normal, with |q|=2pi/d. For single crystals, the planes of the
    full list must not add anything (they are not in the Bragg condition), but
    for layered crystals they might (due to rotation around the c-axis)."""
        txt = NC.createTextData( f'stdlib::{fn}' ).rawData
        full = NC.directLoad( txt, 'comp=bragg;' + cfg, doAbsorption = False )
        ekin = NC.wl2ekin( wl )
        k = 2 * math.pi / wl
        for sel, other in selections:
            hkl = tuple( int(x) for x in sel.split(',') )
            xs = []
            for label, t in ( ( 'SG', txt ), ( 'no SG', strip_sg( txt ) ) ):
                m = load( t, cfg, sel )
                assert m.info.nHKL() == 1
                e = next( m.info.hklObjects() )
                assert e.mult == 2 and ( int(e.h[0]), int(e.k[0]), int(e.l[0]) ) == hkl
                d = float( e.d )
                g = norm( gdir_lab( hkl ) )
                u = bragg_direction( g, d, wl )
                x = float( m.scatter.xsect( wl = wl, direction = u ) )
                assert x > 0.0
                xs.append( x )
                xfull = float( full.scatter.xsect( wl = wl, direction = u ) )
                if layered:
                    assert x <= xfull * ( 1 + 1e-9 )
                else:
                    assert abs( x - xfull ) < 1e-6 * x
                #No scattering at the Bragg geometry of another plane of the family:
                ho = tuple( int(v) for v in other.split(',') )
                uo = bragg_direction( norm( gdir_lab( ho ) ), d, wl )
                assert float( m.scatter.xsect( wl = wl, direction = uo ) ) == 0.0
                #Scattered neutrons have q along +-g, and |q|=2pi/d:
                ef, uf = m.scatter.sampleScatter( ekin, u, repeat = 100 )
                for ef_i, uf_i in zip( ef, zip( *uf ) ):
                    assert abs( ef_i - ekin ) < 1e-12 * ekin
                    q = tuple( k * ( a - b ) for a, b in zip( uf_i, u ) )
                    qmag = math.sqrt( sum( x*x for x in q ) )
                    assert abs( qmag - 2 * math.pi / d ) < 1e-6 * qmag
                    cosang = abs( sum( a*b for a, b in zip( q, g ) ) ) / qmag
                    assert cosang > math.cos( math.radians( 5.0 ) )
            assert abs( xs[0] - xs[1] ) < 1e-9 * xs[0]
            print(f'  {fn} {sel}: d={d:.5g}Aa, xsect at Bragg geometry {xs[0]:.4g} barn,'
                  f' zero for ({other}) (same without space group)')

    #Al (cubic), crystal axes aligned with the lab axes:
    check_single_crystal( 'Al_sg225.ncmat',
                          'mos=0.3deg;dir1=@crys_hkl:0,0,1@lab:0,0,1;'
                          'dir2=@crys_hkl:1,0,0@lab:1,0,0',
                          lambda hkl : hkl, 2.0,
                          [ ( '1,1,1', '1,1,-1' ), ( '0,0,2', '2,0,0' ),
                            ( '-3,1,1', '1,3,1' ) ] )

    #Pyrolytic graphite (layered crystal model), (0,0,l) planes along lab z:
    check_single_crystal( 'C_sg194_pyrolytic_graphite.ncmat',
                          'mos=2deg;lcaxis=0,0,1;dir1=@crys_hkl:0,0,1@lab:0,0,1;'
                          'dir2=@crys_hkl:1,0,0@lab:1,0,0',
                          lambda hkl : hkl, 3.0,
                          [ ( '0,0,2', '1,0,0' ), ( '0,0,-4', '1,0,0' ) ],
                          layered = True )

    #Powder: cross sections below the Bragg thresholds of two selected planes
    #must be proportional to d*|F|^2 (multiplicity 2 for both):
    al = NC.createTextData( 'stdlib::Al_sg225.ncmat' ).rawData
    for t in ( al, strip_sg( al ) ):
        res = []
        for sel in ( '1,1,1', '2,0,0', '-3,1,1' ):
            m = load( t, '', sel )
            e = next( m.info.hklObjects() )
            assert m.info.braggthreshold == 2 * float( e.d )
            assert float( m.scatter.xsect( wl = 2.01 * float( e.d ) ) ) == 0.0
            res.append( ( float( m.scatter.xsect( wl = 1.0 ) ),
                          float( e.d ) * float( e.f2 ) ) )
        for x, df2 in res[1:]:
            assert abs( ( x / res[0][0] ) / ( df2 / res[0][1] ) - 1 ) < 1e-9
    print('  Powder cross sections of selected planes scale as d*|F|^2')

def main():
    test_select()
    test_selectphys()

if __name__ == '__main__':
    main()
