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


# Checks that the stdlib files of NCrystal 4.4.6 (copies in tests/data/refnc4d4,
# without the two big liquid water kernel files) still give exactly the same
# materials as the current stdlib files with the same names.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCTestUtils.dirs import get_named_test_data_dir
from NCrystalDev.datasrc import browseFiles

#Stdlib files deliberately changed since 4.4.6 must be added here (with the
#reason) to be excluded from the comparison:
excluded_files = {
    'SiO2-beta_sg180_BetaQuartz.ncmat':
    'O z-coordinates snapped to exact fractions for release 5.0.0',
    'void.ncmat':
    'state of matter (gas) declared for release 5.0.0',
}

#Stdlib files removed or renamed since 4.4.6 must be added here:
removed_files = []

#Files for which inelastic cross sections are also compared:
inelastic_files = [ 'Al_sg225.ncmat',
                    'Ca_sg229_Calcium-gamma.ncmat',
                    'MgD2_sg136_MagnesiumDeuteride.ncmat',
                    'Polyethylene_CH2.ncmat',
                    'UO2_sg225_UraniumDioxide.ncmat' ]

wls = [ 0.3, 1.0, 1.8, 2.5, 4.0, 8.0 ]

def material_values( cfgstr, inelastic ):
    info = NC.createInfo( cfgstr )
    v = [ info.isCrystalline(), str(info.stateOfMatter()),
          info.density, info.numberdensity ]
    if info.hasStructureInfo():
        si = info.structure_info
        v += [ (k,si[k]) for k in sorted(si) ]
    for ai in info.atominfos:
        v += [ ai.atomData.description(False), ai.numberOfPositions,
               ai.msd, ai.debyeTemperature,
               sorted( tuple(p) for p in ai.positions ) ]
    for di in info.dyninfos:
        v += [ type(di).__name__, di.atomData.description(False),
               di.fraction, di.temperature ]
        if isinstance( di, NC.Info.DI_VDOS ):
            v += [ tuple(di.vdos_egrid), tuple(di.vdos_density) ]
    v += list( NC.createScatter( cfgstr + ';inelas=0' ).xsect( wl = wls ) )
    v += list( NC.createAbsorption( cfgstr ).xsect( wl = wls ) )
    if inelastic:
        sc = NC.createScatter( cfgstr + ';elas=0;vdoslux=1' )
        v += list( sc.xsect( wl = wls ) )
    return repr(v)

def main():
    refdir = get_named_test_data_dir('refnc4d4')
    reffiles = sorted( f.name for f in refdir.glob('*.ncmat') )
    stdlib = { f.name for f in browseFiles( factory = 'stdlib' ) }
    for fn in ( list(excluded_files) + removed_files + inelastic_files ):
        assert fn in reffiles, f'{fn} not in tests/data/refnc4d4'
    ncompared = 0
    for fn in reffiles:
        if fn in removed_files:
            assert fn not in stdlib, f'{fn} is in removed_files but exists'
            print(f'{fn}: removed or renamed since 4.4.6')
            continue
        if fn not in stdlib:
            raise SystemExit( f'No stdlib file {fn} (removed or renamed?):'
                              ' add it to removed_files at the top of'
                              ' this test.' )
        refname = f'refnc4d4_{fn}'
        NC.registerInMemoryFileData( refname, (refdir/fn).read_text() )
        inel = fn in inelastic_files
        same = ( material_values( f'virtual::{refname}', inel )
                 == material_values( f'stdlib::{fn}', inel ) )
        if fn in excluded_files:
            if same:
                raise SystemExit( f'{fn} is in excluded_files but is no'
                                  ' longer different: remove it there.' )
            print(f'{fn}: excluded ({excluded_files[fn]})')
            continue
        if not same:
            raise SystemExit( f'{fn} differs from the 4.4.6 version. If'
                              ' deliberate, add it to excluded_files at the'
                              ' top of this test.' )
        ncompared += 1
        print(f'{fn}: identical{" (incl. inelastic)" if inel else ""}')
    print(f'All {ncompared} compared files identical.')

if __name__ == '__main__':
    main()
