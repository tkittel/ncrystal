
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

# Shared checks of S(alpha,beta) sampling, used by several test scripts.

import math

from NCrystalDev.misc import evaluate_query as ncquery

from NCTestUtils.sabsampleutils import onedim_projections
from NCTestUtils.stat import kolmogorov_smirnov_pvalue


def pval_from_samples( samples1, samples2 ):
    ksp = kolmogorov_smirnov_pvalue
    pvals = {}
    pvals_worst = 2.0
    for key, fct in onedim_projections:
        tmp = fct(samples1.alpha,samples1.beta)
        pv = ksp(fct(samples2.alpha,samples2.beta),tmp)
        pvals_worst = min( pv, pvals_worst )
        pvals[key] = pv
    return pvals, pvals_worst

nperpt = 400

def check_method_choice( cfg, lbl ):
    d = ncquery( ['sab','proctries',cfg,lbl,str(nperpt),'1234'] )
    sp = d['sabproc']
    ar = sp['diagnostics']['full_cell_AR']
    dedicated = sp['diagnostics']['dedicated_strip_expected_tries']
    si = sp['sample_info']
    strip = dict( zip( si['epts_with_beta_strip_sampling'],
                       si['beta_strip_expected_tries'] ) )
    assert strip, 'beta-strip sampling not used anywhere'
    assert not si['epts_with_bounded_cell_sampling']
    for k, (a, ded, meas) in enumerate( zip( ar, dedicated, d['mean_tries'] ) ):
        if k in strip:
            t = strip[k]
            assert t * a < 1.0, f'strips chosen at grid point {k} but more expensive'
            #Mean of nperpt geometric-like try counts (6 sigma):
            tol = 6 * math.sqrt( max( t - 1.0, 0.0 ) / t / nperpt ) + 0.01
            assert abs( meas / t - 1 ) < tol, f'grid point {k}: {meas} vs {t} tries'
        elif ded > 0.0 and a > 0.0:
            assert 1.0 / a < 2.0 * ded, f'full-cell sampling much worse at grid point {k}'
    print(f'{cfg} [{lbl}]: sampling methods and tries OK')
