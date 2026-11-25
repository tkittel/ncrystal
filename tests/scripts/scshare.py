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
import math
import random

def crystal_component( sc ):
    #The SCBragg or LCBragg component.
    s = sc.getSummary(short=False)
    comps = ( [ c for _, c in s['specific']['components'] ]
              if s['name'] == 'ProcComposition' else [ s ] )
    cryst = [ c for c in comps if c['name'] in ('SCBragg','LCBragg') ]
    assert len(cryst) == 1
    assert cryst[0]['isOriented']
    return cryst[0], comps

def bragg_uids( sc ):
    #Returns ( name and model_uid of the SCBragg or LCBragg component,
    #          uid of any PowderBragg component )
    cryst, comps = crystal_component( sc )
    powder = [ c['uid'] for c in comps if c['name'] == 'PowderBragg' ]
    assert len(powder) <= 1
    return ( ( cryst['name'], cryst['specific']['model_uid'] ),
             ( powder[0] if powder else None ) )

#A lab vector (x,y,z) with orientation o1 corresponds to (z,x,y) with o2:
o1 = ( ';dir1=@crys_hkl:5,1,1@lab:0,0,1;dir2=@crys_hkl:0,-1,1@lab:0,1,0'
       ';dirtol=180deg' )
o2 = ( ';dir1=@crys_hkl:5,1,1@lab:1,0,0;dir2=@crys_hkl:0,-1,1@lab:0,0,1'
       ';dirtol=180deg' )
o3 = ';dir1=@crys:1,0,0@lab:0,1,1;dir2=@crys:0,1,0@lab:1,0,0'
o4 = o3 + ';rot=1,2,3:17deg'
o5 = ';dir1=@crys_hkl:5,1,1@lab:0,0,1;dir2=arbitrary'

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
    for o in ( o1, o2, o3, o4, o5 ):
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

def check_descriptions():
    #The model specific descriptions of SCBragg and LCBragg.
    sc = nccore.createScatter( 'Ge_sg227.ncmat;dcutoff=0.5;mos=2deg'
                               ';sccutoff=0;comp=bragg' + o1 )
    s = sc.getSummary(short=False)
    assert s['name'] == 'SCBragg'
    d = s['specific']
    assert sorted( d.keys() ) == [ 'crystal2lab', 'dmax', 'dmin', 'model_uid',
                                   'mos', 'nfamilies', 'summarystr' ]
    assert abs( d['mos'] - math.radians(2.0) ) < 1e-14
    assert 0.5 <= d['dmin'] < d['dmax'] < 4.0
    assert d['nfamilies'] > 3
    assert d['summarystr'].startswith('nfamilies=')
    pg = nccore.createScatter( 'C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.5'
                               ';mos=2deg;sccutoff=0;comp=bragg'
                               ';dir1=@crys_hkl:0,0,1@lab:0,1,0' )
    s = pg.getSummary(short=False)
    assert s['name'] == 'LCBragg'
    d = s['specific']
    assert sorted( d.keys() ) == [ 'crystal2lab', 'dmax', 'dmin', 'lcaxis',
                                   'lcmode', 'model_uid', 'mos',
                                   'nplanesets', 'summarystr' ]
    assert d['lcaxis'] == [ 0.0, 0.0, 1.0 ] and d['lcmode'] == 0
    assert 0.5 <= d['dmin'] < d['dmax'] < 4.0 and d['nplanesets'] > 3
    #The lab direction of the layer axis (crystal z) is y:
    col_z = [ row[2] for row in d['crystal2lab'] ]
    assert max( abs(a-b) for a,b in zip( col_z, [0.0,1.0,0.0] ) ) < 1e-12
    assert 'lcaxis=(0,0,1)' in d['summarystr']
    assert 'lcmode' not in d['summarystr']

def main():
    check_descriptions()
    ge = 'Ge_sg227.ncmat;dcutoff=0.3;mos=2deg'
    check_sharing( ge, 'SCBragg', True,
                   [ ';mos=1deg', ';mosprec=1e-4', ';sccutoff=0.5',
                     ';sccutoff=0', ';temp=200K', ';dcutoff=0.35' ] )
    check_sharing( ge + ';sccutoff=0', 'SCBragg', False, [] )
    pg = 'C_sg194_pyrolytic_graphite.ncmat;dcutoff=0.3;mos=2deg'
    check_sharing( pg, 'LCBragg', True,
                   [ ';lcmode=3', ';lcaxis=0,1,1', ';mos=1deg',
                     ';lcaxis=@crys_uvw:1,0,2' ] )
    #Also when leaving out dir2:
    sc1 = nccore.createScatter( pg + ';dir1=@crys_hkl:0,0,1@lab:0,0,1' )
    sc2 = nccore.createScatter( pg + ';dir1=@crys_hkl:0,0,1@lab:0,1,1' )
    sc3 = nccore.createScatter( pg + o1 )
    assert bragg_uids( sc1 ) == bragg_uids( sc2 ) == bragg_uids( sc3 )
    check_xs( ge + ';comp=bragg' )
    check_xs( pg + ';comp=bragg' )

if __name__ == '__main__':
    main()
