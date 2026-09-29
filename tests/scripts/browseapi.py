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

# Test the NCrystal.browse Python API.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.stabilise_ncpprint # noqa F401
import NCrystalDev as NC
import NCrystalDev.browse as nb
from NCrystalDev.misc import evaluate_query
from NCTestUtils.common import ensure_error

_crystal = """NCMAT v7
#
#   A small Al crystal.
#
#   Mentions Togo.
@CELL
 cubic 4.04958
@SPACEGROUP
 225
@ATOMPOSITIONS
 Al 0 1/2 1/2
 Al 0 0 0
 Al 1/2 1/2 0
 Al 1/2 0 1/2
@DEBYETEMPERATURE
 Al 400
"""

_gas = """NCMAT v7
# ------------------
#  Helium gas.
#  Second line.
#
#  Other paragraph.
@STATEOFMATTER
  gas
@DENSITY
  0.5 kg_per_m3
@DYNINFO
  element He
  fraction 1
  type freegas
"""

def nthreads():
    return evaluate_query(['util','factorythreads'])['nthreads']

def first( b ):
    return next( iter( b ) )

def show( title, b ):
    print(f'{title}: {b!r} {b.names()}')

def lazlau( fmt ):
    #Small .laz/.lau data (high dcutoff):
    from NCrystalDev.mcstasutils import cfgstr_2_hkl
    return '\n'.join( cfgstr_2_hkl( cfgstr = ( 'stdlib::Al_sg225.ncmat;'
                                               'dcutoff=1.0' ),
                                    tgtformat = fmt, verbose = False,
                                    fp_format = '%.8g' ) ) + '\n'

def test_entries():
    b = nb.DataBrowser( factory = 'virtual' )
    print('==> Entries:')
    for e in b:
        print(f'  {e} name={e.name!r} factory={e.factory!r}'
              f' comments={e.comments!r}')
        print(f'     props: {e.props}' if e.props
              else f'     ERROR {e.error!r}')
    p = first(b.match('crystal')).props
    assert p.sg == 225 and p.natoms == 4
    assert p.composition == ((13,((0,1.0),)),)
    assert p.crystalsystem == 'cubic' and p.crystal is True
    assert abs( p.braggthreshold - 2*4.04958/3**0.5 ) < 1e-9
    assert abs( p.cohsl - 3.449 ) < 1e-3
    assert abs( p.mass - 26.9815 ) < 1e-3 and p.incohxs > 0.0
    assert p.vdosdebye and not p.scatknl
    with ensure_error(AttributeError,'PhysicsProps has no attribute "foo"'):
        _ = p.foo
    g = first(b.match('gas.ncmat')).props
    assert g.sg is None and g.crystalsystem is None
    assert g.braggthreshold is None and g.crystal is False
    e = first(b.match('broken'))
    assert e.props is None and e.error
    e = first(b.match('crystal'))
    m = e.load(';temp=100K')
    assert m.info.getTemperature() == 100.0 and m.info.hasStructureInfo()
    print('Entries OK')

def test_selection():
    b = nb.DataBrowser()
    al = b.match('Al_sg225')
    show('Al_sg225',al)
    assert al.names() == ['virtual::Al_sg225.ncmat','stdlib::Al_sg225.ncmat']
    v = nb.DataBrowser( factory = 'virtual' )
    show('virtual',v)
    assert len(v) == 4 and bool(v) and not v.match('nonexistent')
    show('.NCMAT',v.match('.NCMAT'))
    show('stdlib::Be',b.match('stdlib::Be'))
    show('two patterns',v.match('crystal','gas'))
    show('search togo',v.search('togo'))
    show('search togo small',v.search('togo','small'))
    show('chained search',v.search('togo').search('mentions'))
    show('where formula',v.where('formula=he'))#case-insensitive
    show('where not He',v.where('formula!=He'))
    show('where several',v.where('crystal=yes','sg==225','natoms>3'))
    show('where None',v.where('sg > 200'))#missing values never match
    show('where None !=',v.where('sg != 225'))
    show('where text',v.where('state=GAS',"formula='He'"))
    show('where float',v.where('density>=1.5e0','cohsl>.5'))
    show('chain',v.match('.ncmat').search('gas').where('absxs<1'))
    #Selections never modify the original:
    assert len(v) == 4 and v.names()[0] == 'virtual::Al_sg225.ncmat'
    #Constructor arguments:
    show('constructor',nb.DataBrowser('gas','crystal',factory='virtual',
                                      search='small',where='crystal=yes'))
    print('Selection OK')

def test_output():
    import json

    import NCrystalDev._common as nc_common
    v = nb.DataBrowser(factory='virtual')
    print('==> dump():')
    v.dump()
    print('==> dump() of search:')
    v.search('togo').dump()
    d = json.loads(v.to_json())
    assert d[2]['props']['sg'] == 225 and d[2]['path'] is None
    print('to_json keys:',list(d[2]))
    print('==> to_json() of crystal:')
    nc_common.ncpprint(json.loads(v.match('crystal').to_json()))
    assert repr(v.match('crystal')) == 'DataBrowser(1 entry)'
    print('Output OK')

def test_errors():
    def bad( msg, fct ):
        with ensure_error(NC.NCBadInput,msg):
            fct()
    v = nb.DataBrowser(factory='virtual')
    bad('Unknown TextData factory: "nonexistent"',
        lambda : nb.DataBrowser(factory='nonexistent'))
    propnames = [ n for n,d in nb.physics_props_doc( conditions = True ) ]
    bad('Invalid property "foo" in condition "foo > 1" (must be one'
        ' of: '+', '.join(propnames)+')',lambda : v.where('foo > 1'))
    bad('Invalid property "SG" in condition "SG=225" (names are'
        ' case-sensitive, did you mean "sg"?)',lambda : v.where('SG=225'))
    bad('Invalid condition "crystal" (must be like "absxs>100" or sg==225)',
        lambda : v.where('crystal'))
    bad('Invalid condition "absxs >" (missing value)',
        lambda : v.where('absxs >'))
    bad('Invalid condition "absxs>x" (value "x" is not a number)',
        lambda : v.where('absxs>x'))
    bad('Invalid condition "elements=B" ("elements" can not be used in'
        ' conditions)',lambda : v.where('elements=B'))
    bad('Invalid condition "state<x" (only = and != can be used'
        ' with "state")',lambda : v.where('state<x'))
    bad('Invalid condition "crystal=maybe" (value "maybe" is not yes or'
        ' no)',lambda : v.where('crystal=maybe'))
    bad('Invalid condition "crystalsystem=bla" ("crystalsystem" must be'
        ' one of: triclinic, monoclinic, orthorhombic, tetragonal, trigonal,'
        ' hexagonal, cubic)',lambda : v.where('crystalsystem=bla'))
    print('Errors OK')

def test_lazy_loading():
    #Physics is loaded on demand, once per factory, shared by derived
    #browsers (the progress hook is internal):
    nstdlib = evaluate_query(['util','browsedb'])['stdlib']
    progress = []
    b = nb.DataBrowser( _progress = lambda a,b : progress.append(b) )
    assert not progress
    s1 = b.match('virtual::').where('crystal=yes')
    assert set(progress) == {4}#only the virtual factory loaded
    s2 = b.match('virtual::').where('density>1')
    assert set(progress) == {4} and s1 and s2
    b.match('stdlib::Al_sg225').where('crystal=yes')
    assert set(progress) == {4,nstdlib}
    #Physics is also loaded on demand when accessed via an entry:
    progress.clear()
    b = nb.DataBrowser( _progress = lambda a,b : progress.append(b) )
    e = first(b.match('virtual::crystal'))
    assert not progress
    assert e.props.sg == 225 and e.error is None and progress
    progress.clear()
    assert first(b.match('virtual::broken')).error and not progress#cached
    print('Lazy loading OK')

def test_lazlau():
    #Non-NCMAT data, both in-memory and on-disk (in current directory):
    import pathlib

    from NCTestUtils.common import work_in_tmpdir
    NC.registerInMemoryFileData('mem.laz',lazlau('laz'))
    NC.registerInMemoryFileData('mem.lau',lazlau('lau'))
    NC.enableRelativePaths(True)
    with work_in_tmpdir():
        pathlib.Path('disk.laz').write_text(lazlau('laz'))
        pathlib.Path('disk.lau').write_text(lazlau('lau'))
        #relpath entries only with factory="relpath":
        assert not nb.DataBrowser().match('relpath::')
        for e in ( list( nb.DataBrowser( '.la', factory = 'virtual' ) )
                   + list( nb.DataBrowser( factory = 'relpath' ) ) ):
            print(f'  {e.fullkey}: comments={e.comments} error={e.error}')
            p = e.props
            print(f'     sg={p.sg} crystalsystem={p.crystalsystem}'
                  f' formula={p.formula} natoms={p.natoms}'
                  f' braggthreshold={p.braggthreshold:.6g}'
                  f' vdos={p.vdos} scatknl={p.scatknl}')
            assert e.load().info.hasStructureInfo()
            print('     load OK',flush=True)
        #Staged flushed markers: this area was observed to hang on the
        #Windows Debug CI legs, and these localise the next such hang:
        print('lazlau: leaving tmpdir',flush=True)
    print('lazlau: tmpdir cleaned up',flush=True)
    NC.enableRelativePaths(False)
    print('lazlau: relative paths disabled again',flush=True)
    print('LAZ/LAU data OK')

def test_elements():
    #Natural elements and isotopes never match each other:
    NC.registerInMemoryFileData('b10.ncmat',
                                'NCMAT v7\n@DENSITY\n  1 g_per_cm3\n'
                                '@DYNINFO\n  element B10\n  fraction 1\n'
                                '  type freegas\n')
    b = nb.DataBrowser()
    for lbls in (('B',),('B10',),('D',),('H',),('He','B10'),('Al','Xx')):
        try:
            print(f'with_elements{lbls}:',b.with_elements(*lbls).names()[:8])
        except NC.NCBadInput as e:
            print(f'with_elements{lbls}: NCBadInput({e})')
    assert 'stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat' in b.with_elements('D').names()
    assert 'stdlib::LiquidHeavyWaterD2O_T293.6K.ncmat' not in b.with_elements('H').names()
    assert 'virtual::b10.ncmat' in b.with_elements('B10').names()
    assert 'virtual::b10.ncmat' not in b.with_elements('B').names()
    print('Elements OK')

def test_yesno():
    v = nb.DataBrowser(factory='virtual')
    for val in ('yes','YES','true','On','1','no','False','off','0'):
        print(f'crystal={val}:',v.where(f'crystal={val}').names())
    print('Yes/no values OK')

def test_threads():
    #Factory threads are only changed temporarily (NB: must be last, since
    #it configures factory threads):
    assert nthreads() == 1
    nb.DataBrowser(factory='stdlib').where('crystal=yes')
    assert nthreads() == 1
    NC.enableFactoryThreads(2)
    nb.DataBrowser(factory='stdlib').where('crystal=yes')
    assert nthreads() == ( 2 if evaluate_query(['util','factorythreads'])
                           ['threads_available'] else 1 )
    NC.enableFactoryThreads(1)
    print('Factory thread handling OK')

def main():
    NC.removeAllDataSources()
    NC.enableStandardDataLibrary()
    NC.registerInMemoryFileData('crystal.ncmat',_crystal)
    NC.registerInMemoryFileData('gas.ncmat',_gas)
    NC.registerInMemoryFileData('broken.ncmat',
                                'NCMAT v7\n@DENSITY\n -1 g_per_cm3\n')
    #Hides stdlib::Al_sg225.ncmat:
    NC.registerInMemoryFileData('Al_sg225.ncmat',_gas)

    print('Properties:',[ n for n,d in nb.physics_props_doc() ])
    test_entries()
    test_selection()
    test_output()
    test_errors()
    test_lazy_loading()
    test_lazlau()
    test_elements()
    test_yesno()
    test_threads()

if __name__ == '__main__':
    main()
