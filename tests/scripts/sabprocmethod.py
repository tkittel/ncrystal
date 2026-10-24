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


# Checks the choice of sampling method at each energy grid point of
# SABProcessor (via the ["sab","proctries"] query): beta-strip sampling only
# where its expected tries are below the expected full-cell tries, and measured
# mean tries consistent with the expected ones of the method used.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
from NCrystalDev.misc import evaluate_query as ncquery
import math

cases = [
    ( 'stdlib::Al_sg225.ncmat;vdoslux=1', '' ),
    ( 'stdlib::Al_sg225.ncmat;vdoslux=2003', '' ),
    ( 'stdlib::Polyethylene_CH2.ncmat;vdoslux=1', 'H' ),
    ( 'stdlib::Polyethylene_CH2.ncmat;vdoslux=2003', 'H' ),
    ( 'benzene_solid_100K_sabsmall.ncmat', 'H' ),
    ( 'C_from_benzene_solid_20K_sabsmall.ncmat', '' ),
    ( 'Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat', 'Li' ),
    ( 'H_from_ethanol_liquid_200K_sabsmall.ncmat', '' ),
]

nperpoint = 400

def check( cfgstr, lbl ):
    r = ncquery( [ 'sab', 'proctries', cfgstr, lbl, str(nperpoint), '1234' ] )
    sp = r['sabproc']
    egrid = sp['egrid']
    fc_tries = [ 1.0 / ar for ar in sp['diagnostics']['full_cell_AR'] ]
    si = sp['sample_info']
    strip = dict( zip( si['epts_with_beta_strip_sampling'],
                       si['beta_strip_expected_tries'] ) )
    problems = []
    if si['epts_with_bounded_cell_sampling']:
        problems.append('bounded-cell sampling used')
    for i, tstrip in strip.items():
        if not tstrip < fc_tries[i]:
            problems.append(f'strip sampling at point {i} although'
                            f' {tstrip:g} >= {fc_tries[i]:g} tries')
    #Tries are geometrically distributed (variance t*(t-1) for mean t):
    for i, measured in enumerate( r['mean_tries'] ):
        t = strip.get( i, fc_tries[i] )
        sigma = math.sqrt( t * ( t - 1.0 ) / nperpoint )
        if not abs( measured - t ) < 6.0 * sigma + 1e-3 * t:
            problems.append(f'point {i}: mean tries {measured:g}, expected'
                            f' {t:g}')
    print(f'{cfgstr} [{lbl}]: {len(egrid)} grid points,'
          f' beta-strip sampling used: {"yes" if strip else "no"}')
    for p in problems:
        print(f'  PROBLEM: {p}')
    return not problems

def main():
    nfail = sum( 0 if check(*c) else 1 for c in cases )
    if nfail:
        raise SystemExit(f'ERROR: {nfail} cases failed')
    print('All cases OK')

if __name__ == '__main__':
    main()
