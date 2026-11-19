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


# Tests NCMATComposer.set_number_format.

import NCTestUtils.enable_fpe # noqa F401
import random
import NCrystalDev as NC

rng = random.Random( 42 )
n = 40
#Irregular grid with some very close neighbours, and a density ending with 0:
egrid = [ 0.01 ]
for i in range( n - 1 ):
    egrid.append( egrid[-1] * ( 1.0 + ( 1e-7 if i % 10 == 3 else 0.05 )
                                * rng.uniform( 0.5, 1.5 ) ) )
density = [ rng.uniform( 0.1, 1.0 ) for i in range( n - 1 ) ] + [ 0.0 ]
temperature = 293.15 + rng.uniform( 0.0, 1.0 )
matdens = 2.6989999 + 1e-9 * rng.uniform( 0.0, 1.0 )

def compose( **kwargs ):
    c = NC.NCMATComposer()
    c.set_composition( 'Al', 'Al' )
    c.set_density( matdens, 'g/cm3' )
    c.set_state_of_matter( 'solid' )
    c.set_default_temperature( temperature )
    c.set_dyninfo_vdos( 'Al', vdos_egrid = egrid, vdos = density )
    if kwargs:
        c.set_number_format( **kwargs )
    return c

def read_back( c ):
    info = NC.directLoad( c.create_ncmat(), doScatter = False,
                          doAbsorption = False ).info
    di = info.dyninfos[0]
    return ( list( di.vdosOrigEgrid() ), list( di.vdosOrigDensity() ),
             info.getTemperature().dbl()
             if hasattr( info.getTemperature(), 'dbl' )
             else info.getTemperature() )

def check( lbl, c, srtol, vrtol ):
    eg, dens, temp = read_back( c )
    assert len( eg ) == n and len( dens ) == n
    for x, y in zip( density, dens ):
        assert abs( y - x ) <= vrtol * abs( x ), (lbl,x,y)
    for x, y in zip( egrid, eg ):
        assert abs( y - x ) <= vrtol * abs( x ), (lbl,x,y)
    for i in range( n - 1 ):
        dx, dy = egrid[i+1] - egrid[i], eg[i+1] - eg[i]
        assert abs( dy - dx ) <= vrtol * dx, (lbl,i)
    assert abs( temp - temperature ) <= srtol * temperature, lbl
    print(f'{lbl}: {len(c.create_ncmat())} bytes, values within tolerance')

check( 'defaults', compose(), 1e-14, 1e-6 )
check( 'vector_rtol=1e-3', compose( vector_rtol = 1e-3 ), 5e-14, 1e-3 )
check( 'exact', compose( scalar_rtol = 0, vector_rtol = 0 ), 0.0, 0.0 )
assert read_back( compose( scalar_rtol = 0, vector_rtol = 0 ) )[:2] \
    == ( egrid, density )
print('exact: identical values read back')

#Each setting leaves the other one unchanged:
c = compose( vector_rtol = 1e-3 )
c.set_number_format( scalar_rtol = 0 )
d = compose( scalar_rtol = 0, vector_rtol = 1e-3 )
assert c.create_ncmat() == d.create_ncmat()
c.set_number_format( vector_rtol = 0 )
assert c.create_ncmat() == compose( scalar_rtol = 0,
                                    vector_rtol = 0 ).create_ncmat()
print('set_number_format only changes the given settings')

for kwargs in ( { 'scalar_rtol' : 1e-5 }, { 'vector_rtol' : 2e-3 },
                { 'scalar_rtol' : -1.0 }, { 'vector_rtol' : 1e-15 },
                { 'scalar_rtol' : None } ):
    try:
        compose( **kwargs )
    except NC.NCBadInput as e:
        print(f'set_number_format({kwargs}): BadInput: {e}')
    else:
        raise RuntimeError('not rejected')
try:
    NC.NCMATComposer().set_number_format( 1e-6 )
except TypeError:
    print('set_number_format(1e-6): TypeError (keyword arguments only)')
else:
    raise RuntimeError('not rejected')
