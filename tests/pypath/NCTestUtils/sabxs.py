
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

# adaptation of .xs utils for inelastic sab
import re

import NCrystalDev.core as nccore
from NCrystalDev._numpy import _np_geomspace, _np_linspace
from NCrystalDev.constants import constant_boltzmann, wl2ekin


def _rdtol_for_cfgstr( cfgstr ):
    #The acceptable cross-platform reldiff depends on the vdoslux setting:
    #vdoslux=2000 is not a priority (loosest tolerance), vdoslux=2001 is
    #the historical default tolerance, and vdoslux>=2002 (more refined next-
    #gen settings, where a still-open residual is being chased) get the
    #tightest tolerance:
    m = re.search( r';vdoslux=(\d+)', cfgstr )
    vdoslux = int( m.group(1) ) if m else None
    if vdoslux == 2000:
        return 1e-4
    if vdoslux == 2001:
        return 1e-5
    if vdoslux is not None and vdoslux > 2001:
        return 1e-6
    return 1e-5  #legacy/unspecified vdoslux: unchanged historical default


def run( testgroup ):
    from .xs import XSMonitor
    def testlist_filtered():
        yield from test_list_gen( testgroup )
    mon = XSMonitor( refdatadir = f'sabxs_{testgroup}',
                     matloadfct = _load_fct,
                     egridgenfct = _egrid_fct,
                     testlistgenfct = testlist_filtered,
                     test_rdtol = _rdtol_for_cfgstr )
    mon.run()

_test_focus = ( 'Al_sg225.ncmat',
                'CaH2_sg62_CalciumHydride.ncmat',
                )
#                'Li2O_sg225_LithiumOxide.ncmat' )



def test_list_gen( testgroup ):
    assert testgroup in ('A','B','C','D')

    if testgroup=='B':
        for t in [293.15,10,1000]:
            yield from [
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2000;knllux=0;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2000;knllux=1;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2000;knllux=2;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2000;knllux=3;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2000;knllux=4;temp={t:g}',

                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2001;knllux=0;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2001;knllux=1;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2001;knllux=2;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2001;knllux=3;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2001;knllux=4;temp={t:g}',

                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2002;knllux=0;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2002;knllux=1;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2002;knllux=2;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2002;knllux=3;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2002;knllux=4;temp={t:g}',

                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2003;knllux=0;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2003;knllux=1;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2003;knllux=2;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2003;knllux=3;temp={t:g}',
                f'stdlib::Li2O_sg225_LithiumOxide.ncmat;vdoslux=2003;knllux=4;temp={t:g}',
            ]

    if testgroup=='A':
        import NCTestUtils.enable_testdatapath # noqa F401
        for t in [293.15,10,1000]:
            yield from [
                f'Li_from_Li2O.ncmat;vdoslux=2000;knllux=0;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2000;knllux=1;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2000;knllux=2;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2000;knllux=3;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2000;knllux=4;temp={t:g}',

            f'O_from_Li2O.ncmat;vdoslux=2000;knllux=0;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2000;knllux=1;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2000;knllux=2;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2000;knllux=3;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2000;knllux=4;temp={t:g}',

            f'Li_from_Li2O.ncmat;vdoslux=2001;knllux=0;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2001;knllux=1;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2001;knllux=2;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2001;knllux=3;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2001;knllux=4;temp={t:g}',

            f'O_from_Li2O.ncmat;vdoslux=2001;knllux=0;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2001;knllux=1;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2001;knllux=2;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2001;knllux=3;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2001;knllux=4;temp={t:g}',
                #Keep the "next-gen but with free-gas extender"
                #comparison mode (knllux 200-206) exercised (a single
                #cfg, to not add measurably to the suite runtime):
                f'O_from_Li2O.ncmat;vdoslux=2001;knllux=203;temp={t:g}',

            f'Li_from_Li2O.ncmat;vdoslux=2002;knllux=0;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2002;knllux=1;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2002;knllux=2;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2002;knllux=3;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2002;knllux=4;temp={t:g}',

            f'O_from_Li2O.ncmat;vdoslux=2002;knllux=0;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2002;knllux=1;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2002;knllux=2;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2002;knllux=3;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2002;knllux=4;temp={t:g}',

            f'Li_from_Li2O.ncmat;vdoslux=2003;knllux=4;temp={t:g}',
                f'Li_from_Li2O.ncmat;vdoslux=2004;knllux=4;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2003;knllux=4;temp={t:g}',
                f'O_from_Li2O.ncmat;vdoslux=2004;knllux=4;temp={t:g}',
            ]
        #Direct-kernel material (auto-detected Teff drives the SCT
        #extension; NB no knllux 20x comparison entry: the total xs
        #above Emax is continuity-anchored and hence insensitive to the
        #extender model, giving a byte-identical reference file):
        yield 'stdlib::LiquidWaterH2O_T293.6K.ncmat;vdoslux=2001;knllux=1'
        #Trimmed JENDL-5 solid fixtures (mixed per-element teff/msd
        #outcomes resp. zero-row teff refusal; cf. the file headers):
        yield 'benzene_solid_100K_sabsmall.ncmat;vdoslux=2001;knllux=1'
        yield 'C_from_benzene_solid_20K_sabsmall.ncmat;vdoslux=2001;knllux=1'


def _load_fct( cfgstr ):
    if 'knllux=' not in (''.join(cfgstr.split())):
        #NB: add knllux=3 if knllux not specified, to get through the migration
        #period where knllux=-1 is the default.
        cfgstr += ';knllux=3'
    #in any case, only load inelastic:
    cfgstr += ';comp=inelas'
    return nccore.load(cfgstr)

def _egrid_fct( mat ):
    kT = mat.info.getTemperature()*constant_boltzmann
    n = 10
    e = set(_np_linspace( kT*1e-6, kT*1000, n ))
    e |= set(_np_geomspace( kT*1e-6, kT*1000, n ))
    e |= set(_np_geomspace(1e-10, 1e3, n ) )
    e |= {wl2ekin(wl) for wl in _np_linspace(0.1,15.0,n)}
    return e
