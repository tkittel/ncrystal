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

#Per-row checks of the SABAnalyser (NCSABAnalyser): every row it accepts for
#Teff or msd must agree with the VDOSEval-derived truth, and kernels cut at the
#beta grid ends must give either a refusal or a correct value. Only counts
#are printed, keeping the reference log robust against platform differences.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query as q
from NCrystalDev.ncmat import NCMATComposer
import numpy as np

_cases = [ ('stdlib::Al_sg225.ncmat','Al'),
           ('stdlib::Polyethylene_CH2.ncmat','H'),
           ('stdlib::UO2_sg225_UraniumDioxide.ncmat;temp=20','U') ]
_lux = 3

def analyse( cfgstr, lbl, diag = False ):
    return q(['sab','analyse',cfgstr,lbl]+(['diag'] if diag else []))['sabanalyse']

def relerr( a, b ):
    return abs( a / b - 1.0 )

def check_rows( cfgstr, lbl ):
    r = analyse( f'{cfgstr};vdoslux={_lux}', lbl, diag = True )
    d = r['diagnostics']
    tt, mt = r['teff_vdos'], r['msd_vdos']
    nt = nm = 0
    for t, st in zip( d['teff_row'], d['teff_row_status'] ):
        if st == 0:
            nt += 1
            assert relerr( t, tt ) < 0.01
    for m, st in zip( d['msd_row'], d['msd_row_status'] ):
        if st == 0:
            nm += 1
            assert relerr( m, mt ) < 0.01
    assert relerr( r['auto']['teff'], tt ) < 0.003
    assert relerr( r['auto']['msd'], mt ) < 0.003
    print(f'{cfgstr} {lbl}: all accepted rows within 1% of the truth'
          f' (teff rows: {"many" if nt>=50 else nt}, msd rows:'
          f' {"many" if nm>=50 else nm})')
    return tt, mt

def check_cut_kernels( cfgstr, lbl, tt, mt ):
    info = NC.createInfo( cfgstr )
    di = next( d for d in info.dyninfos if d.atomData.displayLabel() == lbl )
    k = di.loadKernel( vdoslux = _lux )
    a, b = np.asarray( k['alpha'] ), np.asarray( k['beta'] )
    sab = np.asarray( k['sab'] ).reshape( len(b), len(a) )
    cuts = [ ('beta<=2', b <= 2.0), ('beta<=0', b <= 0.0),
             ('beta>=0.2*betamin', b >= 0.2 * b[0]),
             ('beta>=0.05*betamin', b >= 0.05 * b[0]) ]
    for name, mask in cuts:
        c = NCMATComposer()
        c.set_composition( 'X', di.atomData.displayLabel() )
        c.set_density( 1.0 )
        c.set_state_of_matter( 'solid' )
        c.set_dyninfo_scatknl( 'X', alphagrid = a, betagrid = b[mask],
                               temperature = info.getTemperature(),
                               sab = sab[mask,:].flatten() )
        fn = f'sabanalyse_rows_{lbl}_{len(b[mask])}.ncmat'
        c.register_as( fn )
        r = q(['sab','analyse',fn,''])['sabanalyse']['auto']
        te, ms = r['teff'], r['msd']
        assert te is None or relerr( te, tt ) < 0.01
        assert ms is None or relerr( ms, mt ) < 0.01
        #NB: not printing which (could differ between platforms near the
        #thresholds):
        print(f'   cut kernel {name}: teff and msd refused or correct')

def main():
    for cfgstr, lbl in _cases:
        tt, mt = check_rows( cfgstr, lbl )
        check_cut_kernels( cfgstr, lbl, tt, mt )

if __name__ == '__main__':
    main()
