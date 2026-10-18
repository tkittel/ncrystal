
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

import re

import NCrystalDev.core as nccore
from NCrystalDev._numpy import _np_geomspace, _np_linspace
from NCrystalDev.constants import constant_boltzmann, wl2ekin

# Monitoring of inelastic cross sections from the next-gen VDOS expansion and
# kernel integration (vdoslux=2000..2006), on explicit lists of cfg-strings
# chosen to cover the physics (VDOS shape, mass, Debye temperature,
# temperature) and the numerical settings, as well as known troublemakers.
#
# Each test is a short script calling run(testname); the lists below define
# the tests. Reference data lives in tests/data/sabxs/<testname>.txt.

def run( testname ):
    import NCTestUtils.enable_testdatapath # noqa F401

    from .xs import XSMonitor
    _register_synthetic_materials()
    mon = XSMonitor( testname = testname,
                     cfglist = _tests[testname](),
                     matloadfct = _load_fct,
                     egridgenfct = _egrid_fct,
                     rdtolfct = _rdtol_for_cfgstr,
                     datadir = 'sabxs' )
    mon.run()

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
    return 1e-5

def _load_fct( cfgstr ):
    #only inelastic scattering is monitored:
    return nccore.load( cfgstr + ';comp=inelas' )

def _egrid_fct( mat ):
    kT = mat.info.getTemperature()*constant_boltzmann
    n = 10
    e = set(_np_linspace( kT*1e-6, kT*1000, n ))
    e |= set(_np_geomspace( kT*1e-6, kT*1000, n ))
    e |= set(_np_geomspace(1e-10, 1e3, n ) )
    e |= {wl2ekin(wl) for wl in _np_linspace(0.1,15.0,n)}
    return e

################################################################################
# Helpers for defining the lists

def _rotate( entries, *, vdoslux = None, knllux = None ):
    #Expand (material,temperatures) entries into cfg-strings, assigning the
    #next (vdoslux,knllux) pair to each (material,temperature) from a cyclic
    #sequence in which every vdoslux=2000..2005 meets every knllux=0..5
    #once per 36 cfg-strings. Fixed values can be given instead:
    i = 0
    for mat, temps in entries:
        for t in temps:
            v = vdoslux if vdoslux is not None else 2000 + i % 6
            k = knllux if knllux is not None else ( i % 6 + i // 6 ) % 6
            i += 1
            tstr = '' if t is None else f';temp={t:g}'
            yield f'{mat}{tstr};vdoslux={v};knllux={k}'

def _sl( name ):
    return f'stdlib::{name}.ncmat'

#Typical temperatures: 0.1K is the lowest supported value.
_Tcold = ( 5, 20 )
_Twide = ( 5, 20, 293.15, 1000, 3000 )

################################################################################
# Synthetic Debye-model materials (in-memory, no data files): single-element
# materials spanning Debye temperature and mass.

_synthetic_debye = [ ( 'H', 100.0 ), ( 'H', 1000.0 ), ( 'C', 400.0 ),
                     ( 'Al', 400.0 ), ( 'Pb', 90.0 ), ( 'U', 200.0 ) ]

def _synth_name( elem, tdebye ):
    return f'sabxs_debye_{elem}_{tdebye:g}K.ncmat'

_synth_done = [ False ]
def _register_synthetic_materials():
    if _synth_done[0]:
        return
    _synth_done[0] = True
    from NCrystalDev.ncmat import NCMATComposer
    for elem, td in _synthetic_debye:
        c = NCMATComposer()
        c.set_composition( 'X', elem )
        c.set_density( 1.0, 'g/cm3' )
        c.set_dyninfo_debyetemp( 'X', td )
        c.register_as( _synth_name( elem, td ) )

################################################################################
# The tests

def _list_trouble():
    #Known troublemakers from the next-gen development (and Li2O, the
    #original worst offender of the legacy algorithms):
    yield from _rotate( [ ( 'Li_from_Li2O.ncmat', ( 10, 293.15, 1000 ) ),
                          ( 'O_from_Li2O.ncmat', ( 10, 293.15, 1000 ) ),
                          ( _sl('UO2_sg225_UraniumDioxide'), ( 5, ) ),
                          ( _sl('ThO2_sg225_ThoriumDioxide'), ( 293.15, ) ),
                          ( _sl('Al_sg225'), ( 293.15, ) ),
                          ( _sl('CaH2_sg62_CalciumHydride'), ( 10, 1000 ) ),
                          ( _sl('Li2O_sg225_LithiumOxide'), ( 10, ) ) ] )
    #The 1e-15 -> 1e-5 amplification case behind the reducer noise fix:
    yield 'Li_from_Li2O.ncmat;temp=10;vdoslux=2004;knllux=4'
    #Next-gen kernel but with the free-gas extender:
    yield 'O_from_Li2O.ncmat;temp=293.15;vdoslux=2001;knllux=203'
    #Direct kernels (fixed temperature), incl. trimmed JENDL-5 solids with
    #per-element teff/msd outcomes resp. zero-row teff refusal:
    yield from _rotate( [ ( _sl('LiquidWaterH2O_T293.6K'), ( None, None ) ),
                          ( 'benzene_solid_100K_sabsmall.ncmat', ( None, None ) ),
                          ( 'C_from_benzene_solid_20K_sabsmall.ncmat',
                            ( None, None ) ) ] )

def _list_testdata():
    #VDOS fixtures in tests/data at and around their design temperatures:
    yield from _rotate( [ ( 'vdos_paraH2.ncmat', ( 5, 20, 50 ) ),
                          ( 'vdos_orthoD2.ncmat', ( 5, 20, 50 ) ),
                          ( 'vdos_liquidH2.ncmat', ( 14, 20, 40, 293.15 ) ),
                          ( 'vdos_hydrogen_osc.ncmat', ( 5, 14, 30 ) ),
                          ( 'vdos_clathrate_H.ncmat', ( 5, 100, 293.15 ) ),
                          ( 'vdos_clathrate_O2.ncmat', ( 5, 100, 293.15 ) ),
                          ( 'vdos_graphiteoxide_D.ncmat', ( 20, 293.15, 1000 ) ),
                          ( 'Al2O3_Alvdos_from_Al225.ncmat', ( 20, 293.15, 1000 ) ),
                          ( 'Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat',
                            ( None, None ) ) ] )
    #Expected refusal (excessive number of phonon orders needed):
    yield 'vdos_hydrogen_osc.ncmat;temp=293.15;vdoslux=2001;knllux=1'

def _list_light():
    yield from _rotate( [
        #Issue #349 (jumps in T dependence), incl. a close temperature pair:
        ( _sl('Polyethylene_CH2'), ( 20, 250, 293.15, 306.265, 306.27, 350,
                                     1000, 3000 ) ),
        ( _sl('LiH_sg225_LithiumHydride'), ( 5, 293.15, 1000, 3000 ) ),
        ( _sl('CaH2_sg62_CalciumHydride'), ( 20, 293.15, 3000 ) ),
        ( _sl('MgH2_sg136_MagnesiumHydride'), ( 1, 293.15, 1000 ) ),
        ( _sl('MgD2_sg136_MagnesiumDeuteride'), ( 5, 293.15, 1000 ) ),
        #Issue #262:
        ( _sl('CaO2H2_sg164_CalciumHydroxide'), ( 293.15, 700, 3000 ) ),
        ( _sl('LiquidWaterH2O_T293.6K'), ( None, None, None ) ),
        ( _sl('LiquidHeavyWaterD2O_T293.6K'), ( None, None, None ) ),
        ( 'solid::H2O/1gcm3', ( 5, 293.15, 1000 ) ) ] )

def _list_heavy():
    yield from _rotate( [ ( _sl('UO2_sg225_UraniumDioxide'), _Twide ),
                          ( _sl('ThO2_sg225_ThoriumDioxide'), _Twide ),
                          ( _sl('PbS_sg225_LeadSulfide'), _Twide ),
                          ( _sl('KBr_sg225_PotassiumBromide'), _Twide ),
                          ( 'solid::Pb/11.3gcm3', ( 5, 293.15, 1000 ) ) ] )

def _list_highdebye():
    _T = ( 5, 20, 293.15, 1000, 3000, 5000 )
    yield from _rotate( [ ( _sl('C_sg227_Diamond'), _T ),
                          ( _sl('Be_sg194'), _T ),
                          ( _sl('BeO_sg186'), _T ) ] )

def _list_simple():
    yield from _rotate( [ ( _sl('Al_sg225'), _Twide ),
                          ( _sl('Cu_sg225'), _Twide ),
                          ( _sl('MgO_sg225_Periclase'), _Twide ),
                          ( _sl('Li3N_sg191_LithiumNitride'), _Twide ),
                          ( _sl('Li2O_sg225_LithiumOxide'), ( 293.15, ) ) ] )

def _list_debye():
    entries = [ ( _sl('Ca_sg229_Calcium-gamma'), ( 5, 293.15, 1000 ) ),
                ( _sl('CaSiO3_sg2_Wollastonite'), ( 20, 293.15, 3000 ) ),
                ( _sl('B4C_sg166_BoronCarbide'), ( 20, 293.15, 3000 ) ),
                ( _sl('PEEK_C19H12O3'), ( 20, 293.15 ) ),
                ( 'solid::B4C/2.52gcm3', ( 5, 293.15, 1000 ) ),
                ( 'solid::CH2/0.95gcm3', ( 5, 293.15 ) ) ]
    for elem, td in _synthetic_debye:
        #T/Tdebye from ~1e-3 to ~10:
        temps = sorted({ float(f'{td*f:.3g}') for f in (0.01,0.3,3.0) })
        entries.append( ( _synth_name(elem,td), temps ) )
    yield from _rotate( entries )

def _list_lowtemp():
    yield from _rotate( [
        ( _sl('Al_sg225'), ( 0.1, 1 ) ),        #issue #109
        ( _sl('MgH2_sg136_MagnesiumHydride'), ( 0.1, 1 ) ),  #issue #109
        ( _sl('Be_sg194'), ( 0.1, 1 ) ),        #issue #109
        ( _sl('Polyethylene_CH2'), ( 0.1, 1 ) ),
        ( _sl('C_sg227_Diamond'), ( 0.1, ) ),
        ( _sl('UO2_sg225_UraniumDioxide'), ( 0.1, 1 ) ),
        ( 'Li_from_Li2O.ncmat', ( 0.1, 1 ) ),
        ( _sl('Ca_sg229_Calcium-gamma'), ( 0.1, 1 ) ),
        ( 'solid::CH2/0.95gcm3', ( 0.1, ) ),
        ( 'vdos_paraH2.ncmat', ( 0.1, 1 ) ),
        ( _synth_name('H',100.0), ( 0.1, ) ) ] )

def _list_vdoslux6():
    for i, ( mat, t ) in enumerate( [
            ( _sl('Al_sg225'), 293.15 ), ( _sl('Polyethylene_CH2'), 293.15 ),
            ( _sl('UO2_sg225_UraniumDioxide'), 1000 ),
            ( 'Li_from_Li2O.ncmat', 10 ), ( _sl('Be_sg194'), 20 ),
            ( _sl('CaH2_sg62_CalciumHydride'), 293.15 ) ] ):
        yield f'{mat};temp={t:g};vdoslux=2006;knllux={i%4}'

def _list_knllux6():
    for i, ( mat, t ) in enumerate( [
            ( _sl('Al_sg225'), 293.15 ), ( _sl('Polyethylene_CH2'), 20 ),
            ( 'O_from_Li2O.ncmat', 1000 ), ( _sl('UO2_sg225_UraniumDioxide'),
                                            293.15 ) ] ):
        yield f'{mat};temp={t:g};vdoslux={2000+i%3};knllux=6'

def _list_maxlux():
    #Both kernel luxury settings at their maximum (1-2 minutes for most
    #materials, ~30s for diamond):
    yield f'{_sl("C_sg227_Diamond")};temp=293.15;vdoslux=2006;knllux=6'

_tests = {
    'sabxs_trouble' : lambda : list(_list_trouble()),
    'sabxs_testdata' : lambda : list(_list_testdata()),
    'long_sabxs_light' : lambda : list(_list_light()),
    'long_sabxs_heavy' : lambda : list(_list_heavy()),
    'long_sabxs_highdebye' : lambda : list(_list_highdebye()),
    'long_sabxs_simple' : lambda : list(_list_simple()),
    'long_sabxs_debye' : lambda : list(_list_debye()),
    'long_sabxs_lowtemp' : lambda : list(_list_lowtemp()),
    'long_sabxs_vdoslux6' : lambda : list(_list_vdoslux6()),
    'long_sabxs_knllux6' : lambda : list(_list_knllux6()),
    'long_sabxs_maxlux' : lambda : list(_list_maxlux()),
}
