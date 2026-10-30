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

# Test the ["vdos","gn",...] query, which gives access to the G_n functions as
# built by VDOSGn with a selectable VDOSGn::Cfg (default, legacy or maxlux).

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query as ncquery
from NCrystalDev.vdos import extractGn
from NCrystalDev.exceptions import NCBadInput
import numpy as np

cfg = 'stdlib::Al_sg225.ncmat;temp=293.15K'

def gn_query( nmax, *extra ):
    r = ncquery( ['vdos','gn',cfg,'',str(nmax)] + list(extra),
                 huge_arrays = True )
    res = []
    for e in r['Gn']:
        emin, emax = e['energy_range']
        res.append( ( np.linspace( emin, emax, len(e['values']) ),
                      e['values'] ) )
    return r, res

def main():
    print( ncquery(['vdos','list']) )
    info = NC.createInfo( cfg )
    di = info.dyninfos[0]
    allgn = {}
    for c in ( None, 'default', 'legacy', 'maxlux' ):
        r, gn = gn_query( 6, *( [f'cfg={c}'] if c else [] ) )
        allgn[c] = gn
        print(f'cfg={c}: input={r["input"]} n={[e["n"] for e in r["Gn"]]}')
        for n, ( e, g ) in enumerate( gn, start = 1 ):
            assert abs( np.trapezoid( g, e ) - 1.0 ) < 1e-6
            assert g.min() >= 0.0
    #No cfg means default:
    for a, b in zip( allgn[None], allgn['default'] ):
        assert np.array_equal( a[1], b[1] )
    #Same as extractGn (which goes through the raw C API, up to rounding):
    for n in ( 1, 4, 6 ):
        e, g = extractGn( di, n, di.atomData.averageMassAMU(),
                          di.temperature, 1.0 )
        g = np.asarray( g )
        v = allgn['default'][n-1][1]
        sel = g > 1e-6 * g.max()
        assert v.size == g.size
        assert np.abs( v[sel] / g[sel] - 1.0 ).max() < 1e-8
    #G1 is identical for all cfgs, and legacy thins later than default:
    for c in ( 'legacy', 'maxlux' ):
        assert np.array_equal( allgn[c][0][1], allgn['default'][0][1] )
    print('legacy G4 has more points than default G4:',
          allgn['legacy'][3][1].size > allgn['default'][3][1].size )
    for bad in ( ['0'], ['3','cfg=bla'], ['3','cfg=legacy','x'], [] ):
        try:
            ncquery( ['vdos','gn',cfg,''] + bad )
        except NCBadInput as e:
            print('Expected error: ...'+str(e).split('] (',1)[-1][:60]+'...')
        else:
            raise SystemExit('Did not emit exception as required!')
    print('All OK')

if __name__ == '__main__':
    main()
