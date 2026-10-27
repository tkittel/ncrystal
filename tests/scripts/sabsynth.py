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

# Kernel sampling with some of the synthetic kernels in tests/data/synth
# (long_sabsynth.py covers all): SABProcessor against the reference sampler,
# and its choice of sampling method.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCTestUtils.dirs import get_named_test_data_dir
from NCTestUtils.sabchecks import check_method_choice, compare_with_ref

synthdir = get_named_test_data_dir('synth')
for name in [ 'minimal', 'diagonal_sliver', 'alpha_far', 'tiny_values' ]:
    cfg = f'{synthdir}/sabsynth_{name}.ncmat;knllux=3'
    for ekin in [ 1e-5, 0.05 ]:
        assert compare_with_ref( cfg, '', ekin, 5000, 1000 ) > 1e-3
    check_method_choice( cfg, '', nperpt = 100, name = name )

#Cross sections below the energy grid, sigma(E)*sqrt(E) being constant there
#(even for extreme S values, where intermediate results could under- or
#overflow):
for name in [ 'tiny_values', 'huge_values' ]:
    sc = NC.createScatter( f'{synthdir}/sabsynth_{name}.ncmat;knllux=3'
                           ';comp=inelas' )
    ref = sc.crossSectionIsotropic( 1e-100 ) * 1e-50
    for ekin in [ 1e-300, 1e-200, 1e-20 ]:
        assert abs( sc.crossSectionIsotropic( ekin ) * ekin**0.5 / ref - 1 ) < 1e-6
    print(f'{name}: sigma(E)*sqrt(E) constant below the energy grid')
xs0 = NC.createScatter( f'{synthdir}/sabsynth_tiny_values.ncmat;knllux=3'
                        ';comp=inelas' ).crossSectionIsotropic( 0.0 )
assert 0.0 < xs0 < float('inf')
