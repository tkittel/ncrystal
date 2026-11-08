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

#Dedicated test of the legacy inelastic physics (knllux=-1 and vdoslux
#1000..1005), so its coverage does not depend on other tests (and can be
#removed together with the legacy code): kernel processing without a suggested
#Emax, sampling above Emax and at very low energies, the legacy VDOS expansion
#error, and the sharing of SABScatterHelper objects.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
import NCrystalDev.core as nccore
from NCrystalDev.exceptions import NCCalcError
from NCTestUtils.common import ensure_error
from NCTestUtils.sab import extract_scathelper_uids

def xs_and_sampling( cfgstr, ekins ):
    sc = nccore.createScatter( cfgstr + ';comp=inelas' )
    assert len( extract_scathelper_uids( sc ) ) > 0
    print(f'{cfgstr}:')
    for e in ekins:
        xs = sc.crossSectionIsotropic( e )
        ef, mu = sc.sampleScatterIsotropic( e, repeat = 1000 )
        assert xs > 0.0
        assert ( ef > 0.0 ).all() and ( abs( mu ) <= 1.0 ).all()
        print(f'   E={e:g}eV: xs={xs:.3g}b, 1000 scatterings sampled')

def main():
    #Direct kernel without suggested Emax (determined by the legacy code), and
    #sampling far above that Emax and at very low energies:
    xs_and_sampling( 'benzene_solid_100K_sabsmall.ncmat;knllux=-1',
                     [ 1e-7, 0.025, 2.0, 50.0 ] )
    #VDOS-based kernel:
    xs_and_sampling( 'stdlib::Al_sg225.ncmat;vdoslux=1003;knllux=-1',
                     [ 1e-7, 0.025, 20.0 ] )
    #Legacy VDOS expansion error (it can not reach the needed energy):
    with ensure_error( NCCalcError, None ):
        nccore.createScatter( 'stdlib::CaO2H2_sg164_CalciumHydroxide.ncmat;'
                              'temp=1000K;vdoslux=1000;comp=inelas' )
    print('Legacy VDOS expansion error: OK')
    #Sharing of SABScatterHelper objects:
    su = extract_scathelper_uids
    al = 'stdlib::Al_sg225.ncmat;knllux=-1'
    u1 = su( nccore.load( f'{al};temp=300K' ) )
    u2 = su( nccore.load( f'{al}  ;  temp=300K' ) )
    u3 = su( nccore.load( f'{al};temp=200K' ) )
    assert len(u1) == 1 and u1 == u2 and u1 != u3
    nccore.clearCaches()
    assert su( nccore.load( f'{al};temp=300K' ) ) != u1
    print('Sharing of SABScatterHelper objects: OK')

if __name__ == '__main__':
    main()
