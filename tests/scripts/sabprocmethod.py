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

# Checks the sampling method chosen by SABProcessor at each energy grid point:
# beta-strip sampling only where its expected tries are below those of
# full-cell sampling, full-cell sampling elsewhere not much worse than a
# dedicated strip sampler, and measured strip tries as expected.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
from NCTestUtils.sabchecks import check_method_choice as check

check( 'stdlib::Al_sg225.ncmat;temp=20K;knllux=3', '' )
check( 'stdlib::Polyethylene_CH2.ncmat;knllux=3', 'H' )
check( 'stdlib::CaH2_sg62_CalciumHydride.ncmat;temp=20K;knllux=3', 'H' )
#Hard case (quasi-elastic peak on the beta=0 row):
check( 'H_from_ethanol_liquid_200K_sabsmall.ncmat;knllux=3', '' )
