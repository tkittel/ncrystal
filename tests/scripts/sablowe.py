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

# Next-gen kernel sampling and cross sections at E=0 and tiny energies:
# finite results, E'>=0, mu in [-1,1], and isotropic mu below the energy grid.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import numpy as np

cfgs = [ 'stdlib::Al_sg225.ncmat;temp=20K;comp=inelas;knllux=3',
         'stdlib::Polyethylene_CH2.ncmat;comp=inelas;knllux=3' ]
energies = [ 0.0, 5e-324, 1e-300, 1e-20, 1e-12, 1e-8 ]

for cfg in cfgs:
    print(f'==> {cfg}')
    sc = NC.createScatter( cfg )
    xs0 = sc.crossSectionIsotropic( 0.0 )
    assert np.isfinite( xs0 ) and xs0 > 0.0
    for e in energies:
        xs = sc.crossSectionIsotropic( e )
        assert np.isfinite( xs ) and xs > 0.0
        #sigma(E)*sqrt(E) is constant below the energy grid:
        if e >= 1e-300:
            assert abs( xs * np.sqrt( e ) / ( xs0 * np.sqrt( 5e-324 ) ) - 1 ) < 1e-2
        ekin, mu = sc.sampleScatterIsotropic( e, repeat = 10000 )
        assert np.isfinite( ekin ).all() and ( ekin >= 0.0 ).all()
        assert np.isfinite( mu ).all() and ( np.abs( mu ) <= 1.0 ).all()
        assert abs( mu.mean() ) < 0.05
        #(no values for E=0 or subnormal E, where some platforms flush to 0)
        val = f'xs*sqrt(E/eV)={xs*np.sqrt(e):.4g} barn, ' if e >= 1e-300 else ''
        print(f'  E={e:g} eV: {val}checks OK')
