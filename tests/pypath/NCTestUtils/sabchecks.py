
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

import numpy as np
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

def check_method_choice( cfg, lbl, nperpt = 400, require_strips = True,
                         check_fc = True, name = None ):
    d = ncquery( ['sab','proctries',cfg,lbl,str(nperpt),'1234'] )
    sp = d['sabproc']
    ar = sp['diagnostics']['full_cell_AR']
    dedicated = sp['diagnostics']['dedicated_strip_expected_tries']
    si = sp['sample_info']
    strip = dict( zip( si['epts_with_beta_strip_sampling'],
                       si['beta_strip_expected_tries'] ) )
    assert strip or not require_strips, 'beta-strip sampling not used anywhere'
    for k, (a, ded, meas) in enumerate( zip( ar, dedicated, d['mean_tries'] ) ):
        if k == 0:
            continue#sampled just below the grid, by the below-grid sampler
        if k in strip:
            t = strip[k]
            assert t * a < 1.0, f'strips chosen at grid point {k} but more expensive'
            #Within maxRelTries=1.2 of a sampler built at this grid point:
            assert t <= 1.2 * ded, f'strips at grid point {k}: {t} vs {ded} tries'
            #Mean of nperpt geometric-like try counts (6 sigma):
            tol = 6 * math.sqrt( max( t - 1.0, 0.0 ) / t / nperpt ) + 0.01
            assert abs( meas / t - 1 ) < tol, f'grid point {k}: {meas} vs {t} tries'
        elif check_fc and ded > 0.0 and a > 0.0:
            assert 1.0 / a < 2.0 * ded, f'full-cell sampling much worse at grid point {k}'
    print(f'{name or cfg} [{lbl}]: sampling methods and tries OK')

class _AB:
    def __init__( self, d ):
        self.alpha, self.beta = np.asarray(d['alpha']), np.asarray(d['beta'])

def compare_with_ref( cfg, lbl, ekins, nsample, seed ):
    #Worst KS p-values of the 1D projections of (alpha,beta) samples from
    #SABProcessor and the reference sampler, at each of the energies (one
    #query each, so the samplers are only initialised once):
    n, e = str(nsample), ','.join( str(x) for x in ekins )
    proc = ncquery( ['sab','proc',cfg,lbl,n,str(seed),e],
                    huge_arrays = True )['sample']
    ref = ncquery( ['sab','refsample',cfg,lbl,n,str(seed+1),e,'ref'],
                   huge_arrays = True )['refsample']
    if len( ekins ) == 1:
        proc, ref = [ proc ], [ ref ]
    return [ pval_from_samples( _AB(p), _AB(r) )[1]
             for p, r in zip( proc, ref ) ]
