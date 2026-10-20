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

# The NCMAT v8 msd keyword of vdosdebye and scatknl @DYNINFO sections ("msd
# <value> [at_temperature <T>]"), in files and as written by
# NCMATComposer.set_dyninfo_msd, and the debye_temp keyword of scatknl sections.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
import NCrystalDev as NC
from NCrystalDev.vdos import debyeTempFromIsotropicMSD, debyeIsotropicMSD
from NCTestUtils.common import ensure_error

def ncmat( dyninfo, *, cell = True, som = None ):
    s = 'NCMAT v8\n'
    if cell:
        s += ( '@CELL\n  cubic 4.04958\n@SPACEGROUP\n  225\n@ATOMPOSITIONS\n'
               '  Al 0 0 0\n  Al 0 1/2 1/2\n  Al 1/2 0 1/2\n  Al 1/2 1/2 0\n' )
    else:
        s += '@DENSITY\n  2.7 g_per_cm3\n'
    if som:
        s += f'@STATEOFMATTER\n  {som}\n'
    return s + '@DYNINFO\n  element Al\n  fraction 1\n' + dyninfo

def load( content, cfg = ';dcutoff=0.5' ):
    NC.registerInMemoryFileData( 'v8msd.ncmat', content )
    return NC.createInfo( 'v8msd.ncmat' + cfg )

def main():
    vdosdebye()
    scatknl()

def vdosdebye():
    msd, T, mass = 0.0085, 293.15, 26.9815385
    dt_expected = debyeTempFromIsotropicMSD( msd = msd, mass = mass,
                                             temperature = T )
    info = load( ncmat( f'  type vdosdebye\n  msd {msd} at_temperature {T}\n' ),
                 cfg = f';dcutoff=0.5;temp={T}' )
    di = info.dyninfos[0]
    print( f'vdosdebye with msd {msd} at {T} K -> Debye temperature'
           f' {di.debyeTemperature():.6g} K (expected {dt_expected:.6g} K)' )
    assert abs( di.debyeTemperature() - dt_expected ) < 1e-6 * dt_expected
    ai_msd = info.atominfos[0].msd
    print( f'  msd at {T} K from the Info object: {ai_msd:.6g} Aa^2' )
    assert abs( ai_msd - msd ) < 1e-6 * msd
    #Same with an explicit debye_temp gives the same physics:
    info2 = load( ncmat( f'  type vdosdebye\n  debye_temp {dt_expected!r}\n' ),
                  cfg = f';dcutoff=0.5;temp={T}' )
    assert abs( info2.atominfos[0].msd - ai_msd ) < 1e-9
    print( '  identical to an explicit debye_temp' )
    #Non-crystalline solids may use it too, other states not:
    load( ncmat( f'  type vdosdebye\n  msd {msd} at_temperature {T}\n',
                 cell = False, som = 'solid' ), cfg = '' )
    print( 'non-crystalline solid with msd: loaded' )
    with ensure_error( NC.NCBadInput ):
        load( ncmat( f'  type vdosdebye\n  msd {msd} at_temperature {T}\n',
                     cell = False, som = 'liquid' ), cfg = '' )
    with ensure_error( NC.NCBadInput ):
        load( ncmat( f'  type vdosdebye\n  msd {msd}\n' ) )
    with ensure_error( NC.NCBadInput ):
        load( ncmat( f'  type vdosdebye\n  msd {msd} at_temperature {T}\n'
                     f'  debye_temp 400\n' ) )
    #NCMATComposer.set_dyninfo_msd writes the keyword (and hence v8):
    c = NC.NCMATComposer( 'Al_sg225.ncmat' )
    c.set_dyninfo_msd( 'Al', msd = msd, temperature = T )
    out = c.create_ncmat()
    print( 'Composer output:' )
    print( '\n'.join( line for line in out.splitlines()
                      if not line.startswith('#') ) )
    assert out.startswith( 'NCMAT v8\n' )
    assert f'msd {msd!r} at_temperature {T!r}'.replace('0.0085','.0085') in out
    info3 = load( out, cfg = f';dcutoff=0.5;temp={T}' )
    assert abs( info3.atominfos[0].msd - msd ) < 1e-6 * msd
    print( '  round trip preserves the msd' )
    #Non-crystalline: msd implies a solid (written unless set), and a
    #contradicting state of matter is an error:
    c = NC.NCMATComposer()
    c.set_dyninfo_msd( 'Al', msd = msd, temperature = T, fraction = 1.0 )
    c.set_density( 2.7, 'g/cm3' )
    import re
    assert re.search( r'@STATEOFMATTER\s+solid\s', c.create_ncmat() )
    print( 'non-crystalline Composer material with msd is written as a solid' )
    c.set_state_of_matter( 'liquid' )
    with ensure_error( NC.NCBadInput,
                       'The msd/debye_temp dyninfo parameters describe atomic'
                       ' displacements in a solid, which is inconsistent with'
                       ' the declared state of matter "liquid"' ):
        c.create_ncmat()

def li2o_kernels( li, o ):
    #The two 10K scatknl sections of the Li2O test file, with extra lines:
    t = NC.createTextData( 'Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat' )
    t = t.rawData
    sec_li, sec_o = t[t.index('@DYNINFO'):].split('@DYNINFO')[1:]
    assert 'element     Li' in sec_li and 'element     O' in sec_o
    return f'@DYNINFO{sec_li}  {li}\n@DYNINFO{sec_o}  {o}\n'

def li2o_structure():
    lines, keep = [], False
    for line in NC.createTextData( 'stdlib::Li2O_sg225_LithiumOxide.ncmat' ):
        if line.startswith('@'):
            keep = line.startswith( ('@CELL','@SPACEGROUP','@ATOMPOSITIONS') )
        if keep:
            lines.append( line + '\n' )
    return ''.join( lines )

def scatknl():
    #Kernels at 10K, so the material temperature is 10K. The msd is
    #transferred to 10K via the Debye model when given at another
    #temperature, using the same functions as the loader (hence the tight
    #tolerances):
    T = 10.0
    def msd_at_T( msd, at_T, mass ):
        dt = debyeTempFromIsotropicMSD( msd = msd, temperature = at_T,
                                        mass = mass )
        return debyeIsotropicMSD( debye_temperature = dt, temperature = T,
                                  mass = mass )

    def load_cryst( li, o = 'debye_temp 600' ):
        info = load( 'NCMAT v8\n' + li2o_structure() + li2o_kernels( li, o ),
                     cfg = ';temp=10K;dcutoff=0.5' )
        ai = { a.atomData.displayLabel() : a for a in info.atominfos }
        return ai['Li'], ai['O']

    mass_li = load_cryst( 'msd 0.02' )[0].atomData.averageMassAMU()
    mass_o = load_cryst( 'msd 0.02' )[1].atomData.averageMassAMU()
    dt_li = debyeTempFromIsotropicMSD( msd = 0.03, temperature = 300.0,
                                       mass = mass_li )
    print( 'Crystalline Li2O with scatknl sections at 10K:' )
    for li, expected_msd in ( ( 'msd 0.02', 0.02 ),
                              ( 'msd 0.02 at_temperature 10', 0.02 ),
                              ( 'msd 0.03 at_temperature 300',
                                msd_at_T( 0.03, 300.0, mass_li ) ),
                              ( f'debye_temp {dt_li!r}',
                                msd_at_T( 0.03, 300.0, mass_li ) ),
                              ( 'debye_temp 400',
                                debyeIsotropicMSD( debye_temperature = 400.0,
                                                   temperature = T,
                                                   mass = mass_li ) ) ):
        ai_li, ai_o = load_cryst( li )
        #The exact dt_li value is used as input, but printing it in full would
        #expose platform-dependent last digits:
        label = li.replace( repr(dt_li), f'{dt_li:.6g}' )
        print( f'  Li "{label}": msd={ai_li.msd:.6g} Aa^2'
               f' Debye temperature={ai_li.debyeTemperature:.6g} K' )
        assert abs( ai_li.msd - expected_msd ) < 1e-12 * expected_msd
        dt_expected = debyeTempFromIsotropicMSD( msd = expected_msd,
                                                 temperature = T,
                                                 mass = mass_li )
        assert abs( ai_li.debyeTemperature - dt_expected ) < 1e-9 * dt_expected
        o_msd = debyeIsotropicMSD( debye_temperature = 600.0, temperature = T,
                                   mass = mass_o )
        assert abs( ai_o.msd - o_msd ) < 1e-12 * o_msd
    #Without keywords for O, its msd comes from the kernel analysis (only
    #needing to be sane here), not affecting the explicit Li value:
    ai_li, ai_o = load_cryst( 'msd 0.03 at_temperature 300', o = '' )
    print( f'  O without keywords: msd={ai_o.msd:.3g} Aa^2 (from kernel)' )
    assert 1e-4 < ai_o.msd < 0.1
    assert ai_li.msd == load_cryst( 'msd 0.03 at_temperature 300' )[0].msd

    #Non-crystalline solid: the msd enters only the incoherent elastic
    #physics, so compare cross sections of equivalent descriptions:
    def elas_xs( li, o = 'debye_temp 600' ):
        content = ( 'NCMAT v8\n@STATEOFMATTER\n  solid\n'
                    '@DENSITY\n  2.013 g_per_cm3\n' + li2o_kernels( li, o ) )
        NC.registerInMemoryFileData( 'v8msd.ncmat', content )
        sc = NC.createScatter( 'v8msd.ncmat;temp=10K;comp=elas' )
        return tuple( sc.xsect( wl = wl ) for wl in ( 0.5, 1.0, 1.8, 4.0 ) )
    ref = elas_xs( 'msd 0.03 at_temperature 300' )
    print( 'Non-crystalline Li2O, elastic xs at 0.5, 1, 1.8, 4 Aa:',
           ', '.join( f'{x:.5g}' for x in ref ), 'barn' )
    def same( a, b ):
        return all( abs( x - y ) <= 1e-12 * abs( y ) for x, y in zip( a, b ) )
    assert same( elas_xs( f'debye_temp {dt_li!r}' ), ref )
    msd10 = msd_at_T( 0.03, 300.0, mass_li )
    assert same( elas_xs( f'msd {msd10!r}' ), ref )
    assert same( elas_xs( f'msd {msd10!r} at_temperature 10' ), ref )
    #(and the value matters):
    assert not same( elas_xs( 'msd 0.03' ), ref )
    print( '  identical for the equivalent msd and debye_temp entries' )

if __name__ == '__main__':
    main()
