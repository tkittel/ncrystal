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

# Kernel sampling with the synthetic "noise" kernel in tests/data/synth: the
# choice of sampling method of SABProcessor.

import NCTestUtils.enable_fpe # noqa F401
from NCTestUtils.sabchecks import check_synth_kernel

#FIXME: The strip ladder leaves full-cell sampling at a few grid points where
#a dedicated strip sampler would be cheaper (to be improved), hence
#check_fc=False:
check_synth_kernel( 'noise', require_strips = False, check_fc = False )
