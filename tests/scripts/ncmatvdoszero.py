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

#NCMAT v8 requires vdos_density to end with a value of 0, older versions not.
#Excess trailing zeros are ignored everywhere VDOS curves enter NCrystal.

import NCrystalDev as NC


def mat( version, densities, emax = '0.05', emin = '0.01' ):
    return '\n'.join( [ f'NCMAT v{version}', '@STATEOFMATTER', '  solid',
                        '@DENSITY', '  1.0 g_per_cm3', '@DYNINFO',
                        '  element Al', '  fraction 1', '  type vdos',
                        f'  vdos_egrid {emin} {emax}',
                        '  vdos_density ' + ' '.join(densities) ] )

for version in ( 7, 8 ):
    for dens in ( ['1','2','3','2','1'], ['1','2','3','2','1','0'] ):
        fn = f'vdoszero_v{version}_{len(dens)}.ncmat'
        NC.registerInMemoryFileData( fn, mat( version, dens ) )
        try:
            NC.createInfo( fn )
            res = 'loaded'
        except NC.NCBadInput as e:
            res = f'BadInput: {e}'
        print(f'NCMAT v{version}, last density {dens[-1]}: {res}')

#The same curve with no final zero (where allowed), one, or three trailing
#zeros (with emax extended accordingly) must give identical results via NCMAT
#v7 or v8 data, the NCMATComposer, and the C API (analyseVDOS). All energies
#are multiples of 2^-7 eV, so the three are exactly the same curve:
emin = '0.0078125'
curves = [ ( '0.0390625', [1,2,3,2,1] ),
           ( '0.046875', [1,2,3,2,1,0] ),
           ( '0.0625', [1,2,3,2,1,0,0,0] ) ]
wls = [ 0.5, 1.0, 2.0, 4.0, 8.0 ]
def material_values( data ):
    NC.registerInMemoryFileData( 'vdoszero_cmp.ncmat', data )
    NC.clearCaches()
    info = NC.createInfo( 'vdoszero_cmp.ncmat' )
    di = info.dyninfos[0]
    sc = NC.createScatter( 'vdoszero_cmp.ncmat;vdoslux=1;comp=inelas' )
    return repr( [ di.vdos_egrid, list(di.vdos_density),
                   list( sc.xsect( wl = wls ) ) ] )
results = []
for emax, dens in curves:
    sdens = [ str(d) for d in dens ]
    for version in ( 7, 8 ):
        if version == 8 and dens[-1] != 0:
            continue
        results.append( ( f'NCMAT v{version}, {len(dens)} values',
                          material_values( mat( version, sdens, emax, emin ) ) ) )
    c = NC.NCMATComposer()
    c.set_composition( 'Al', 'Al' )
    c.set_density( 1.0, 'g/cm3' )
    c.set_state_of_matter( 'solid' )
    c.set_dyninfo_vdos( 'Al', vdos_egrid = [ float(emin), float(emax) ],
                        vdos = dens )
    results.append( ( f'NCMATComposer, {len(dens)} values',
                      material_values( c.create_ncmat() ) ) )
    from NCrystalDev.vdos import analyseVDOS
    results.append( ( f'analyseVDOS, {len(dens)} values',
                      repr( analyseVDOS( float(emin), float(emax), dens,
                                           300.0, 27.0 ) ) ) )
refs = {}
for lbl, r in results:
    kind = 'analyseVDOS' if lbl.startswith('analyseVDOS') else 'material'
    assert r == refs.setdefault( kind, r ), lbl
    print(f'{lbl}: identical results')

#Fewer than 5 values after ignoring excess trailing zeros is an error:
bad = [ '1', '2', '3', '0', '0' ]
for version in ( 7, 8 ):
    NC.registerInMemoryFileData( f'vdoszero_bad_v{version}.ncmat',
                                 mat( version, bad ) )
    try:
        NC.createInfo( f'vdoszero_bad_v{version}.ncmat' )
    except NC.NCBadInput as e:
        print(f'NCMAT v{version}, {" ".join(bad)}: BadInput: {e}')
    else:
        raise RuntimeError('not rejected')
try:
    NC.NCMATComposer().set_dyninfo_vdos( 'Al', vdos_egrid = [ 0.01, 0.05 ],
                                         vdos = [ 1, 2, 3, 0, 0 ] )
except NC.NCBadInput as e:
    print(f'NCMATComposer, {" ".join(bad)}: BadInput: {e}')
else:
    raise RuntimeError('not rejected')
