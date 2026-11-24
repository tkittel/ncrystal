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

# Test the various ways of specifying the orientation of single crystals, by
# verifying that equivalent specifications give identical cross sections.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev.core as nccore
from NCrystalDev.constants import wl2ekin
import math
import random

def xs_values( cfgstr ):
    #Cross sections at a fixed set of neutron energies and directions.
    sc = nccore.createScatter( cfgstr )
    rng = random.Random( 12345 )
    res = []
    for i in range( 1500 ):
        ekin = wl2ekin( rng.uniform( 0.5, 6.0 ) )
        direction = tuple( rng.gauss( 0.0, 1.0 ) for i in range(3) )
        res.append( sc.crossSection( ekin, direction ) )
    assert sum( 1 for e in res if e > 0.0 ) > 50
    return res

def same_xs( xs1, xs2 ):
    return all( abs( a - b ) <= 1e-9 * max( a, b ) for a, b in zip(xs1,xs2) )

def require_equivalent( base, *orientations ):
    xs = [ xs_values( base + ';' + o ) for o in orientations ]
    for e in xs[1:]:
        assert same_xs( xs[0], e )
    return xs[0]

def require_different( xs, base, *orientations ):
    for o in orientations:
        assert not same_xs( xs, xs_values( base + ';' + o ) )

def test_crysdirs():
    #Cubic: Lattice directions and plane normals coincide (also with the
    #Cartesian axes).
    require_equivalent(
        'Ge_sg227.ncmat;dcutoff=0.5;mos=2deg;comp=bragg',
        'dir1=@crys_hkl:1,1,0@lab:0,0,1;dir2=@crys_hkl:0,0,1@lab:0,1,0',
        'dir1=@crys_uvw:1,1,0@lab:0,0,1;dir2=@crys_uvw:0,0,1@lab:0,1,0',
        'dir1=@crys_uvw:2,2,0@lab:0,0,5;dir2=@crys_hkl:0,0,1@lab:0,1,0',
        'dir1=@crys:1,1,0@lab:0,0,1;dir2=@crys_uvw:0,0,1@lab:0,1,0',
    )
    #Hexagonal: b is 120 degrees from a (which is along x), with b* along y
    #and c along c* (and z).
    be = 'Be_sg194.ncmat;dcutoff=0.5;mos=2deg;comp=bragg'
    y120 = math.sqrt(0.75)
    xs = require_equivalent(
        be,
        'dir1=@crys_uvw:0,1,0@lab:0,0,1;dir2=@crys_uvw:0,0,1@lab:0,1,0',
        f'dir1=@crys:-0.5,{y120!r},0@lab:0,0,1;dir2=@crys:0,0,1@lab:0,1,0',
        f'dir1=@crys:-0.5,{y120!r},0@lab:0,0,1;dir2=@crys_hkl:0,0,1@lab:0,1,0',
        'dir1=@crys_uvw:0,1,0@lab:0,0,1;dir2=@crys_hkl:0,0,3@lab:0,1,0',
    )
    require_different(
        xs, be,
        'dir1=@crys_hkl:0,1,0@lab:0,0,1;dir2=@crys_uvw:0,0,1@lab:0,1,0',
    )
    require_equivalent(
        be,
        'dir1=@crys_hkl:0,1,0@lab:0,0,1;dir2=@crys_uvw:1,0,0@lab:0,1,0',
        'dir1=@crys:0,1,0@lab:0,0,1;dir2=@crys:1,0,0@lab:0,1,0',
    )
    #Triclinic. The direction [uvw] is in the plane (hkl) if hu+kv+lw=0, so
    #the following are consistent (even with a tiny dirtol). Also note that
    #the part of dir2 along dir1 is always ignored.
    tri = 'CaSiO3_sg2_Wollastonite.ncmat;dcutoff=1.0;mos=2deg;comp=bragg'
    xs = require_equivalent(
        tri + ';dirtol=1e-9',
        'dir1=@crys_hkl:1,0,1@lab:0,0,1;dir2=@crys_uvw:1,2,-1@lab:0,1,0',
        'dir1=@crys_hkl:2,0,2@lab:0,0,1;dir2=@crys_uvw:-1,-2,1@lab:0,-1,0',
        'dir1=@crys_hkl:1,0,1@lab:0,0,1;dir2=@crys_uvw:1,2,-1@lab:0,1,0',
        'dir1=@crys_uvw:1,2,-1@lab:0,1,0;dir2=@crys_hkl:1,0,1@lab:0,0,1',
    )
    require_different(
        xs, tri,
        'dir1=@crys_hkl:1,0,1@lab:0,0,1;dir2=@crys_hkl:1,2,-1@lab:0,1,0'
        ';dirtol=180deg',
    )

def rotate( v, axis, angle ):
    #Right-handed rotation of v around axis (Rodrigues' formula).
    n = math.sqrt( sum( e*e for e in axis ) )
    k = [ e / n for e in axis ]
    kxv = ( k[1]*v[2] - k[2]*v[1], k[2]*v[0] - k[0]*v[2], k[0]*v[1] - k[1]*v[0] )
    kdotv = sum( a*b for a, b in zip(k,v) )
    c, s = math.cos( angle ), math.sin( angle )
    return tuple( v[i]*c + kxv[i]*s + k[i]*kdotv*(1.0-c) for i in range(3) )

def vstr( v ):
    return f'{v[0]!r},{v[1]!r},{v[2]!r}'

def test_rot_impl( base, c1, c2 ):
    lab1, lab2 = ( 0.1, 0.2, 1.0 ), ( 0.0, 1.0, -0.2 )
    def orient( axis = None, angle = 0.0 ):
        #dir1 and dir2, with lab directions optionally rotated.
        a, b = ( ( rotate( v, axis, angle ) for v in ( lab1, lab2 ) )
                 if axis else ( lab1, lab2 ) )
        return f'dir1={c1}@lab:{vstr(a)};dir2={c2}@lab:{vstr(b)}'
    o = orient()
    xs = require_equivalent(
        base, o,
        o + ';rot=0,0,1:0',
        o + ';rot=0,1,0:360deg',
        o + ';rot=1,-2,3:-720deg',
    )
    deg = math.pi / 180
    for axis, angle, anglestrs in [
            ( (0,1,0), 30*deg,
              ( '30deg', '1800arcmin', '108000 arcsec', '-330deg',
                repr(30*deg), repr(30*deg)+'rad' ) ),
            ( (1,-2,3), -77.7*deg, ( '-77.7deg', ) ),
            ( (0.3,0,-1), 180*deg, ( '180deg', '-180deg' ) ),
    ]:
        #Opposite axis and angle is the same rotation:
        neg = vstr( [ -e for e in axis ] ) + f':{-angle!r}'
        xs_rot = require_equivalent(
            base,
            orient( axis, angle ),
            o + ';rot=' + neg,
            *( o + f';rot={vstr(axis)}:{a}' for a in anglestrs ),
        )
        assert not same_xs( xs, xs_rot )
    #The rotation is applied after dir1 and dir2, also if these are not
    #consistent (in which case the lab part of dir2 along dir1 is ignored):
    require_equivalent(
        base + ';dirtol=180deg',
        orient( (0,1,0), 0.3 ),
        o + ';rot=0,1,0:0.3',
        'rot=0,1,0:0.3;' + o,
    )

def test_rot():
    test_rot_impl( 'Ge_sg227.ncmat;dcutoff=0.5;mos=2deg;comp=bragg',
                   '@crys_hkl:5,1,1', '@crys_hkl:0,-1,1' )
    test_rot_impl( 'C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.5;mos=2deg'
                   ';comp=bragg', '@crys_hkl:0,0,1', '@crys_uvw:1,0,0' )

def test_arbitrary():
    #The cross section of a neutron along the lab direction of dir1 does not
    #depend on the rotation of the crystal around that direction.
    def xs_along( cfgstr, direction ):
        sc = nccore.createScatter( cfgstr )
        res = [ sc.crossSection( wl2ekin( 0.5 + 0.002 * i ), direction )
                for i in range( 3000 ) ]
        assert sum( 1 for e in res if e > 0.0 ) > 50
        return res
    base = 'Ge_sg227.ncmat;dcutoff=0.5;mos=2deg;comp=bragg;'
    lab1 = ( 0.1, 0.2, 1.0 )
    #NB: Using directions which are not normals of any reflecting plane, since
    #the cross section of a neutron exactly along such a normal is not
    #numerically stable at back-scattering.
    for c1 in ( '@crys_hkl:3,5,10', '@crys_uvw:2,1,0.1', '@crys:0.1,-0.3,1',
                '@crys:1,0.01,0.02' ):
        dir1 = f'dir1={c1}@lab:{vstr(lab1)}'
        xs = [ xs_along( base + dir1 + ';' + dir2, lab1 )
               for dir2 in ( 'dir2=arbitrary',
                             'dir2=arbitrary;dirtol=1e-9',
                             'dir2=@crys:1,2,3@lab:0,1,0;dirtol=180deg',
                             'dir2=@crys_hkl:-3,0,1@lab:1,1,0;dirtol=180deg' ) ]
        assert all( same_xs( xs[0], e ) for e in xs[1:] )
        #Also when subsequently rotated:
        axis, angle = ( 1.0, -2.0, 0.5 ), 0.7
        rot = f';rot={vstr(axis)}:{angle!r}'
        xs_rot = xs_along( base + dir1 + ';dir2=arbitrary' + rot,
                           rotate( lab1, axis, angle ) )
        assert same_xs( xs[0], xs_rot )
        assert not same_xs( xs[0], xs_along( base + dir1 + ';dir2=arbitrary'
                                             + rot, lab1 ) )
    #For layered crystals with dir1 along lcaxis, the rotation around dir1 is
    #of no importance at all:
    pg = 'C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.5;mos=2deg;comp=bragg'
    require_equivalent(
        pg,
        'dir1=@crys_hkl:0,0,1@lab:0.1,0.2,1;dir2=arbitrary',
        'dir1=@crys_hkl:0,0,1@lab:0.1,0.2,1;dir2=@crys:1,0,0@lab:0,1,0'
        ';dirtol=180deg',
        'dir1=@crys_uvw:0,0,1@lab:0.1,0.2,1;dir2=@crys_hkl:1,5,0@lab:1,0,0'
        ';dirtol=180deg',
    )

def test_lcaxis():
    #In hexagonal crystals c is along c* and z, b is 120 degrees from a
    #(which is along x), and b* is along y.
    y120 = math.sqrt(0.75)
    o = 'dir1=@crys_hkl:1,0,3@lab:0.1,0.2,1;dir2=@crys_uvw:1,0,0@lab:0,1,-0.2'
    for mat in ( 'C_sg194_pyrolytic_graphite.ncmat', 'Be_sg194.ncmat' ):
        base = mat + ';dcutoff=0.5;mos=2deg;comp=bragg;dirtol=180deg;' + o
        xs = require_equivalent(
            base,
            'lcaxis=0,0,1',
            'lcaxis=@crys:0,0,2',
            'lcaxis=@crys_hkl:0,0,1',
            'lcaxis=@crys_hkl:0,0,-7',
            'lcaxis=@crys_uvw:0,0,0.5',
        )
        xs_b = require_equivalent(
            base,
            'lcaxis=@crys_uvw:0,1,0',
            f'lcaxis=-0.5,{y120!r},0',
            'lcaxis=@crys_uvw:1,2,0;lcaxis=@crys_uvw:0,-3,0',
        )
        xs_bstar = require_equivalent(
            base,
            'lcaxis=@crys_hkl:0,1,0',
            'lcaxis=0,1,0',
            'lcaxis=@crys_uvw:1,2,0',
        )
        assert not same_xs( xs, xs_b )
        assert not same_xs( xs, xs_bstar )
        assert not same_xs( xs_b, xs_bstar )
    #An axis in the lab frame. Here the orientation takes the Cartesian
    #crystal directions (x,y,z) to the lab directions (z,x,y):
    o = 'dir1=@crys_hkl:0,0,1@lab:0,1,0;dir2=@crys_uvw:1,0,0@lab:0,0,1'
    for mat in ( 'C_sg194_pyrolytic_graphite.ncmat', 'Ge_sg227.ncmat' ):
        base = mat + ';dcutoff=0.5;mos=2deg;comp=bragg;'
        xs = require_equivalent(
            base + o,
            'lcaxis=@lab:0.2,1,0.5',
            'lcaxis=@crys:0.5,0.2,1',
            'lcaxis=@lab:-0.4,-2,-1',
        )
        require_different( xs, base + o, 'lcaxis=@lab:0.5,0.2,1' )
        #An axis in the lab frame is not affected by a rotation of the
        #crystal:
        axis, angle = ( 1.0, -2.0, 0.5 ), 0.7
        lab1, lab2 = ( vstr( rotate( v, axis, angle ) )
                       for v in ( (0,1,0), (0,0,1) ) )
        o_rot = ( f'dir1=@crys_hkl:0,0,1@lab:{lab1}'
                  f';dir2=@crys_uvw:1,0,0@lab:{lab2}' )
        xs_rot = require_equivalent(
            base + 'lcaxis=@lab:0.2,1,0.5',
            o + f';rot={vstr(axis)}:{angle!r}',
            o_rot,
        )
        assert not same_xs( xs, xs_rot )
        #...unlike an axis in the crystal:
        lc_rot = vstr( rotate( (0.2,1,0.5), axis, angle ) )
        require_equivalent(
            base,
            'lcaxis=@crys:0.5,0.2,1;' + o_rot,
            'lcaxis=@crys:0.5,0.2,1;' + o + f';rot={vstr(axis)}:{angle!r}',
            f'lcaxis=@lab:{lc_rot};' + o_rot,
        )
    #dir2 can be left out when lcaxis is along dir1 (the PG file itself has
    #lcaxis=0,0,1):
    pg = 'C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.5;mos=2deg;comp=bragg'
    require_equivalent(
        pg,
        'dir1=@crys_hkl:0,0,1@lab:0.1,0.2,1',
        'dir1=@crys_uvw:0,0,-3@lab:-0.1,-0.2,-1',
        'dir1=@crys:0,0,1@lab:0.1,0.2,1',
        'dir1=@crys_hkl:0,0,1@lab:0.1,0.2,1;lcaxis=@crys_hkl:0,0,2',
        'dir1=@crys_hkl:0,0,1@lab:0.1,0.2,1;dir2=arbitrary',
        'dir1=@crys_hkl:0,0,1@lab:0.1,0.2,1;dir2=@crys:1,2,3@lab:1,0,0'
        ';dirtol=180deg',
    )
    require_equivalent(
        'Be_sg194.ncmat;dcutoff=0.5;mos=2deg;comp=bragg',
        'lcaxis=@crys_uvw:0,1,0;dir1=@crys_uvw:0,1,0@lab:0.1,0.2,1',
        f'lcaxis=@crys_uvw:0,1,0;dir1=@crys:-0.5,{y120!r},0@lab:0.1,0.2,1',
        f'lcaxis=-0.5,{y120!r},0;dir1=@crys_uvw:0,2,0@lab:0.1,0.2,1',
        'lcaxis=@crys_uvw:0,1,0;dir1=@crys_uvw:0,1,0@lab:0.1,0.2,1'
        ';dir2=@crys_hkl:1,0,3@lab:1,0,0;dirtol=180deg',
    )

def test_bad():
    from NCrystalDev.exceptions import NCBadInput
    ge = 'Ge_sg227.ncmat;dcutoff=0.5'
    o = ';mos=2deg;dir1=@crys_hkl:0,0,1@lab:0,0,1;dir2=@crys_hkl:1,0,0@lab:1,0,0'
    nccore.createScatter( ge + o + ';rot=0,1,0:1deg' )
    for bad in [ ge + ';rot=0,1,0:1deg',
                 ge + o + ';rot=0,0,0:1deg',
                 ge + o + ';rot=0,1,0',
                 ge + o + ';rot=0,1,0:1degree',
                 ge + o + ';rot=@lab:0,1,0:1deg',
                 ge + ';mos=2deg;dir1=arbitrary;dir2=arbitrary',
                 ge + ';mos=2deg;dir2=arbitrary',
                 ge + ';dir1=@crys_hkl:0,0,1@lab:0,0,1;dir2=arbitrary',
                 ge + o.replace('dir2=','dir2=arbitrary;dir3='),
                 #dir2 can not be left out:
                 ge + ';mos=2deg;dir1=@crys_hkl:0,0,1@lab:0,0,1',
                 ge + ';mos=2deg;dir1=@crys_hkl:0,0,1@lab:0,0,1'
                 ';lcaxis=@lab:0,0,1',
                 ge + ';mos=2deg;dir1=@crys_hkl:0,0,1@lab:0,0,1'
                 ';lcaxis=@crys_hkl:0,1,1',
                 ge + ';mos=2deg;dir1=@crys_hkl:0,0,1@lab:0,0,1'
                 ';lcaxis=0,1e-2,1',
                 ( 'Be_sg194.ncmat;mos=2deg;dir1=@crys_uvw:0,1,0@lab:0,0,1'
                   ';lcaxis=@crys_hkl:0,1,0' ),
                 #Bad syntax:
                 ge + o + ';lcaxis=@crys_xyz:0,0,1',
                 ge + o + ';lcaxis=@lab:0,0,0',
                 ge + o + ';lcaxis=@lab:0,0',
                ]:
        try:
            nccore.createScatter( bad )
        except NCBadInput:
            continue
        raise SystemExit( 'Did not fail: ' + bad )

def main():
    test_crysdirs()
    test_rot()
    test_arbitrary()
    test_lcaxis()
    test_bad()

if __name__ == '__main__':
    main()
