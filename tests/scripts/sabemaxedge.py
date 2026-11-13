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

#Emax of SAB kernels where a few cells dominate all others by many orders of
#magnitude, and of the smallest allowed kernels (5x5).

import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.dirs import get_named_test_data_dir
from NCrystalDev.ncmat import NCMATComposer
from NCrystalDev.misc import evaluate_query
import NCrystalDev as NC
import numpy as np

def emax_div_kT( cfgstr ):
    d = evaluate_query( ['sab','proctries',cfgstr,'','0','1'] )
    return d['sabproc']['egrid'][-1]

def register( name, alpha, beta, sab ):
    c = NCMATComposer()
    c.set_composition( 'H', 'H' )
    c.set_density( 1.0, 'g/cm3' )
    c.set_state_of_matter( 'liquid' )
    c.set_dyninfo_scatknl( 'H', alphagrid = alpha, betagrid = beta,
                           temperature = 293.6, sab = sab )
    NC.registerInMemoryFileData( name, c.create_ncmat( meta_data = False ) )

def main():
    f = get_named_test_data_dir('synth') / 'sabsynth_dominant_cells.ncmat'
    info = NC.createInfo( f )
    k = info.dyninfos[0].loadKernel()
    a, b, s = k['alpha'], k['beta'], np.array( k['sab'] ).copy()
    #Same kernel with the dominant cells only 1e10 above the rest:
    s[ s > 1.0 ] *= 1e-20
    register( 'dominant_1e10.ncmat', a, b, s )
    e30, e10 = emax_div_kT( str(f) ), emax_div_kT( 'dominant_1e10.ncmat' )
    print( f'Emax/kT with dominant cells x1e30: {e30:.4g}, x1e10: {e10:.4g}' )
    assert abs( e30 / e10 - 1.0 ) < 1e-6
    kT = NC.constants.constant_boltzmann * 293.6
    sc = NC.createScatter( f'{f};comp=inelas' )
    print( 'xs at E/kT=0.5,1,2:',
           ' '.join( f'{sc.xsect( ekin = t*kT ):.4g}' for t in (0.5,1.0,2.0) ) )

    #Smallest allowed grid, with a Gaussian S:
    a5, b5 = np.linspace( 0.1, 20.0, 5 ), np.linspace( -10.0, 10.0, 5 )
    s5 = np.exp( -b5[:,None]**2 / 20 - a5[None,:] / 20 )
    register( 'gauss5x5.ncmat', a5, b5, s5 )
    print( f'5x5 kernel: Emax/kT {emax_div_kT("gauss5x5.ncmat"):.4g}' )
    sc = NC.createScatter( 'gauss5x5.ncmat;comp=inelas' )
    print( f'5x5 kernel: xs at E/kT=1: {sc.xsect( ekin = kT ):.4g}' )

if __name__ == '__main__':
    main()
