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

#At very low temperatures, VDOS kernel beta grids must end where
#S(alpha,beta)=exp(-beta)*S(alpha,-beta) underflows, keep the requested number
#of beta points, and not trigger the "Discarding N edges" warning (which would
#show up in the output).

import NCrystalDev as NC


def check( cfgstr ):
    info = NC.createInfo( cfgstr )
    for di in info.dyninfos:
        if not isinstance( di, NC.Info.DI_VDOS ):
            continue
        k = di.loadKernel( vdoslux = 1 )
        b = k['beta']
        s = k['sab'].reshape( len(b), len(k['alpha']) )
        assert len(b) == 240#from VDOS::gridDimFromLux
        assert s[-1].any() and s[0].any()
        print(f'{cfgstr} {di.atomData.displayLabel()}: nbeta={len(b)}'
              f' betamax={b[-1]:.4g}')


for c in ( 'stdlib::ThO2_sg225_ThoriumDioxide.ncmat;temp=0.1',
           'stdlib::CeO2_sg225_CeriumOxide.ncmat;temp=0.1',
           'stdlib::ThO2_sg225_ThoriumDioxide.ncmat;temp=1' ):
    check( c )
