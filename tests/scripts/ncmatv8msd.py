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

# The NCMAT v8 msd keyword of vdosdebye @DYNINFO sections ("msd <value>
# at_temperature <T>").

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.vdos import debyeTempFromIsotropicMSD
from NCTestUtils.common import ensure_error

def ncmat( dyninfo, *, cell = True, som = None ):
    s = 'NCMAT v8\n'
    if cell:
        s += ( '@CELL\n  cubic 4.04958\n@SPACEGROUP\n  225\n@ATOMPOSITIONS\n'
               '  Al 0 0 0\n  Al 0 1/2 1/2\n  Al 1/2 0 1/2\n  Al 1/2 1/2 0\n' )
    else:
        s += '@DENSITY\n  2.7 g_per_cm3\n'
    if som:
        s += f'@STATEOFMATTER\n  {som}\n'
    return s + '@DYNINFO\n  element Al\n  fraction 1\n' + dyninfo

def load( content, cfg = ';dcutoff=0.5' ):
    NC.registerInMemoryFileData( 'v8msd.ncmat', content )
    return NC.createInfo( 'v8msd.ncmat' + cfg )

def main():
    msd, T, mass = 0.0085, 293.15, 26.9815385
    dt_expected = debyeTempFromIsotropicMSD( msd = msd, mass = mass,
                                             temperature = T )
    info = load( ncmat( f'  type vdosdebye\n  msd {msd} at_temperature {T}\n' ),
                 cfg = f';dcutoff=0.5;temp={T}' )
    di = info.dyninfos[0]
    print( f'vdosdebye with msd {msd} at {T} K -> Debye temperature'
           f' {di.debyeTemperature():.6g} K (expected {dt_expected:.6g} K)' )
    assert abs( di.debyeTemperature() - dt_expected ) < 1e-6 * dt_expected
    ai_msd = info.atominfos[0].msd
    print( f'  msd at {T} K from the Info object: {ai_msd:.6g} Aa^2' )
    assert abs( ai_msd - msd ) < 1e-6 * msd
    #Same with an explicit debye_temp gives the same physics:
    info2 = load( ncmat( f'  type vdosdebye\n  debye_temp {dt_expected!r}\n' ),
                  cfg = f';dcutoff=0.5;temp={T}' )
    assert abs( info2.atominfos[0].msd - ai_msd ) < 1e-9
    print( '  identical to an explicit debye_temp' )
    #Non-crystalline solids may use it too, other states not:
    load( ncmat( f'  type vdosdebye\n  msd {msd} at_temperature {T}\n',
                 cell = False, som = 'solid' ), cfg = '' )
    print( 'non-crystalline solid with msd: loaded' )
    with ensure_error( NC.NCBadInput ):
        load( ncmat( f'  type vdosdebye\n  msd {msd} at_temperature {T}\n',
                     cell = False, som = 'liquid' ), cfg = '' )
    with ensure_error( NC.NCBadInput ):
        load( ncmat( f'  type vdosdebye\n  msd {msd}\n' ) )
    with ensure_error( NC.NCBadInput ):
        load( ncmat( f'  type vdosdebye\n  msd {msd} at_temperature {T}\n'
                     f'  debye_temp 400\n' ) )

if __name__ == '__main__':
    main()
