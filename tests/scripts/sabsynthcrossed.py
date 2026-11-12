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

#The integral of S within the kinematic boundary must increase with energy, also
#while a dominant cell with extreme S values is being crossed. This requires
#its crossed parts to be integrated precisely enough at every knllux level.

import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.dirs import get_named_test_data_dir
import NCrystalDev as NC

f = get_named_test_data_dir('synth') / 'sabsynth_crossed_dominant_cell.ncmat'
kT = NC.constants.constant_boltzmann * 293.6
for knllux in range(7):
    sc = NC.createScatter( f'{f};knllux={knllux};comp=inelas' )
    xs = [ sc.xsect( ekin = t * kT ) for t in ( 15.0, 18.0, 20.5, 21.0, 25.0 ) ]
    print( f'knllux={knllux}: xs at E/kT=15,18,20.5,21,25:',
           ' '.join( f'{x:.4g}' for x in xs ) )
