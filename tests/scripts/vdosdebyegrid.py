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

#Checks the idealised Debye VDOS grid (createVDOSDebye): layout, msd vs the
#analytic Debye model, and that the Python and C++ versions agree (via cross
#sections of "vdosdebye" vs "vdos" with the Python grid).

import NCrystalDev as NC
import NCrystalDev.vdos as ncv

kB = NC.constant_boltzmann

def check( debye_temp, temp, element ):
    (emin, emax), dens = ncv.createVDOSDebye( debye_temp )
    h = kB * debye_temp / 1000.5
    assert len(dens) == 11
    assert dens[-1] == 0.0 and dens[-2] > 0.0
    assert abs( emin / h - 991.0 ) < 1e-9 and abs( emax / h - 1001.0 ) < 1e-9
    mass = NC.atomDB(element).averageMassAMU()
    msd = ncv.analyseVDOS( emin, emax, dens, temp, mass )['msd']
    msd_ref = ncv.debyeIsotropicMSD( debye_temperature = debye_temp,
                                     temperature = temp, mass = mass )
    #(VDOSEval itself is good to ~5e-6 at low T, the old grid was off by 5e-5)
    assert abs( msd / msd_ref - 1.0 ) < 1e-5
    def mat( dyn ):
        return '\n'.join( [ 'NCMAT v7', '@STATEOFMATTER', '  solid',
                            '@DENSITY', '  1.0 g_per_cm3', '@DYNINFO',
                            f'  element {element}', '  fraction 1' ] + dyn )
    fn_d = f'debyegrid_d_{element}_{debye_temp}.ncmat'
    fn_v = f'debyegrid_v_{element}_{debye_temp}.ncmat'
    NC.registerInMemoryFileData( fn_d, mat( [ '  type vdosdebye',
                                              f'  debye_temp {debye_temp}' ] ) )
    NC.registerInMemoryFileData( fn_v, mat( [ '  type vdos',
                                              f'  vdos_egrid {emin:.17g} {emax:.17g}',
                                              '  vdos_density '
                                              + ' '.join(f'{d:.17g}' for d in dens) ] ) )
    #(vdoslux=103: expansions, as for "vdos", not the grouped treatment)
    cfg = f';temp={temp}K;vdoslux=103;comp=inelas'
    sc_d = NC.createScatter( fn_d + cfg )
    sc_v = NC.createScatter( fn_v + cfg )
    for e in ( 1e-5, 0.0253, 0.5 ):
        xd, xv = sc_d.xsect( wl = NC.ekin2wl(e) ), sc_v.xsect( wl = NC.ekin2wl(e) )
        #(not exact: vdosdebye rounds its cache key, e.g. the mass, to 1e-7)
        assert abs( xd / xv - 1.0 ) < 1e-7
    print(f'TDebye={debye_temp}K T={temp}K {element}: ok')

check( 100, 293.15, 'Pb' )
check( 400, 293.15, 'O' )
check( 600, 10, 'H' )
check( 2000, 1000, 'C' )
