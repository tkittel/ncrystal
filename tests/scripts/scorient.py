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

def main():
    test_crysdirs()

if __name__ == '__main__':
    main()
