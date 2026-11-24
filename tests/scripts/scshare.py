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

# Verify that the Bragg diffraction models of single crystals are shared between
# materials differing only in the orientation of the crystal.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev.core as nccore
from NCrystalDev.constants import wl2ekin
import random

def bragg_uids( sc ):
    #Returns ( name and uid of the model wrapped in a RotatedScatter,
    #          uid of any PowderBragg component )
    s = sc.getSummary(short=False)
    comps = ( [ c for _, c in s['specific']['components'] ]
              if s['name'] == 'ProcComposition' else [ s ] )
    rot = [ c for c in comps if c['name'] == 'RotatedScatter' ]
    powder = [ c['uid'] for c in comps if c['name'] == 'PowderBragg' ]
    assert len(rot) == 1 and len(powder) <= 1
    assert rot[0]['isOriented']
    (scale, model), = rot[0]['specific']['components']
    assert scale == 1.0
    return ( model['name'], model['uid'] ), ( powder[0] if powder else None )

#A lab vector (x,y,z) with orientation o1 corresponds to (z,x,y) with o2:
o1 = ( ';dir1=@crys_hkl:5,1,1@lab:0,0,1;dir2=@crys_hkl:0,-1,1@lab:0,1,0'
       ';dirtol=180deg' )
o2 = ( ';dir1=@crys_hkl:5,1,1@lab:1,0,0;dir2=@crys_hkl:0,-1,1@lab:0,0,1'
       ';dirtol=180deg' )
o3 = ';dir1=@crys:1,0,0@lab:0,1,1;dir2=@crys:0,1,0@lab:1,0,0'

def check_sharing( base, name, has_powder, variations ):
    keepalive = []
    def load( cfgstr ):
        keepalive.append( nccore.createScatter( cfgstr ) )
        return bragg_uids( keepalive[-1] )
    ref = load( base + o1 )
    assert ref[0][0] == name
    assert ( ref[1] is not None ) == has_powder
    #Same model no matter the orientation, or which other components are
    #present:
    for o in ( o1, o2, o3 ):
        for extra in ( '', ';inelas=0', ';comp=bragg', ';vdoslux=1' ):
            assert load( base + o + extra ) == ref
    #But not if parameters affecting the model itself are changed:
    seen = set( ref ) - { None }
    for v in variations:
        for o in ( o1, o3 ):
            ids = load( base + o + v )
            assert ids[0][0] == name
            assert not seen.intersection( ids )
            assert ids == load( base + o2 + v )
        seen.update( e for e in ids if e is not None )
    #The cache does not survive a clearance:
    assert load( base + o2 ) == ref
    nccore.clearCaches()
    assert not seen.intersection( load( base + o2 ) )

def check_xs( base ):
    #Cross sections must be the same if both the crystal and the neutron are
    #rotated.
    sc1 = nccore.createScatter( base + o1 )
    sc2 = nccore.createScatter( base + o2 )
    assert bragg_uids( sc1 ) == bragg_uids( sc2 )
    rng = random.Random( 12345 )
    n = 0
    for i in range( 4000 ):
        ekin = wl2ekin( rng.uniform( 0.5, 6.0 ) )
        x, y, z = ( rng.gauss( 0.0, 1.0 ) for i in range(3) )
        xs1 = sc1.crossSection( ekin, ( x, y, z ) )
        xs2 = sc2.crossSection( ekin, ( z, x, y ) )
        assert abs( xs1 - xs2 ) <= 1e-9 * max( xs1, xs2 )
        if xs1 > 0.0:
            n += 1
    assert n > 100

def main():
    ge = 'Ge_sg227.ncmat;dcutoff=0.3;mos=2deg'
    check_sharing( ge, 'SCBragg', True,
                   [ ';mos=1deg', ';mosprec=1e-4', ';sccutoff=0.5',
                     ';sccutoff=0', ';temp=200K', ';dcutoff=0.35' ] )
    check_sharing( ge + ';sccutoff=0', 'SCBragg', False, [] )
    pg = 'C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.3;mos=2deg'
    check_sharing( pg, 'LCBragg', True,
                   [ ';lcmode=3', ';lcaxis=0,1,1', ';mos=1deg' ] )
    check_xs( ge + ';comp=bragg' )
    check_xs( pg + ';comp=bragg' )

if __name__ == '__main__':
    main()
