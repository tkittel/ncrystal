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

#NCMAT v8 requires vdos_density to end with a value of 0, older versions not.

import NCrystalDev as NC


def mat( version, densities ):
    return '\n'.join( [ f'NCMAT v{version}', '@STATEOFMATTER', '  solid',
                        '@DENSITY', '  1.0 g_per_cm3', '@DYNINFO',
                        '  element Al', '  fraction 1', '  type vdos',
                        '  vdos_egrid 0.01 0.05',
                        '  vdos_density ' + ' '.join(densities) ] )

for version in ( 7, 8 ):
    for dens in ( ['1','2','3','2','1'], ['1','2','3','2','1','0'] ):
        fn = f'vdoszero_v{version}_{len(dens)}.ncmat'
        NC.registerInMemoryFileData( fn, mat( version, dens ) )
        try:
            NC.createInfo( fn )
            res = 'loaded'
        except NC.NCBadInput as e:
            res = f'BadInput: {e}'
        print(f'NCMAT v{version}, last density {dens[-1]}: {res}')
