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

# Kernel sampling with all the synthetic kernels in tests/data/synth:
# SABProcessor against the reference sampler, and its choice of sampling
# method.

import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.dirs import get_named_test_data_dir
from NCTestUtils.sabchecks import check_method_choice, compare_with_ref

synthdir = get_named_test_data_dir('synth')
names = [ 'alpha_far', 'diagonal_sliver', 'fine_beta_coarse_alpha',
          'huge_values', 'minimal', 'noise', 'spike_alphamin',
          'steep_alpha_coarse_beta', 'subnormal_values', 'tiny_values',
          'zero_holes' ]
for name in names:
    cfg = f'{synthdir}/sabsynth_{name}.ncmat;knllux=3'
    #subnormal_values: samplers disagree there (to be investigated):
    if name != 'subnormal_values':
        for ekin in [ 1e-5, 1e-3, 0.05 ]:
            assert compare_with_ref( cfg, '', ekin, 20000, 2000 ) > 1e-3
    #noise: the strip ladder leaves full-cell sampling at a few grid points
    #where a dedicated strip sampler would be cheaper (to be improved):
    check_method_choice( cfg, '', require_strips = False,
                         check_fc = ( name != 'noise' ), name = name )
