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

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.stabilise_ncpprint # noqa F401
import NCrystalDev.cli as nc_cli
from NCrystalDev.exceptions import NCBadInput
from NCTestUtils.env import ncsetenv
from NCTestUtils.common import ensure_error
from argparse import ArgumentError
import shlex
import json
import pathlib
from NCrystalDev.misc import evaluate_query as ncquery

def rawlist( *items ):
    #List argument passing the items verbatim (cf. the sgsym and cif queries):
    return [ 'RAWDATA', '@' + '@'.join( items ) ]

def safefmt(s):
    return str(s).encode('unicode_escape').decode('utf-8')

def test(*query):
    print()

    print(f">>> Sending query: {[safefmt(e) for e in query]!r}")
    r = ncquery(query)
    import pprint
    pprint.pprint(r)

def test_cli( *args ):
    print(f"============= CLI >>{shlex.join(safefmt(a) for a in args)}"
          "<< ====================")
    nc_cli.run('query',*args)
    print("===========================================")

def test_matcfg():
    #Raw data can follow a "RAWDATA" argument, exempting it from the
    #character restrictions on query strings:
    from NCrystalDev.misc import evaluate_query
    data = ( 'NCMAT v7\n# Comment with forbidden chars "#<>&;" and a NCRYSTALMATCFG[temp=20K]\n'
             '@CELL\n  cubic 4.04958\n@SPACEGROUP\n  225\n'
             '@ATOMPOSITIONS\n  Al 0 0 0\n  Al 0 0.5 0.5\n'
             '  Al 0.5 0 0.5\n  Al 0.5 0.5 0\n'
             '@DEBYETEMPERATURE\n  Al   410.4\n' )
    for name in ('','my #data.ncmat'):
        d = evaluate_query(['matcfg','decodedata','ncmat','vdoslux=1',
                            'RAWDATA',name,'RAWDATA',data])
        print('decodedata:',d['data_name'],d['textdata_name'],d['pars'])
    with ensure_error(NCBadInput,
                      'Forbidden character "#" in query! Problem found in:'
                      ' "#x"'):
        evaluate_query(['matcfg','decodedata','ncmat','','RAWDATA','',
                        'NOTRAWDATA','#x'])
    with ensure_error(NCBadInput,
                      'Forbidden character "#" in query! Problem found in:'
                      ' "#"'):
        evaluate_query(['matcfg','decode','RAWDATA','x','#'])
    with ensure_error(NCBadInput,
                      'Invalid matcfg JSON query: ["matcfg","decode"]'
                      ' (correct usage: ["matcfg","decode",CFGSTR])'):
        evaluate_query(['matcfg','decode'])
    with ensure_error(NCBadInput,
                      'Invalid matcfg JSON query: ["matcfg","doc","x"]'
                      ' (correct usage: ["matcfg","doc",MODE] with MODE'
                      ' "full", "short" or "json")'):
        evaluate_query(['matcfg','doc','x'])

def main():
    from NCrystalDev._common import print
    ncsetenv('FIX_QUERY_VERSION_FOR_TESTS','1')

    with ensure_error(NCBadInput,
                      'Missing or empty query'):
        test()
    with ensure_error(NCBadInput,
                      'Invalid JSON query key: "foo"'):
        test('foo')

    with ensure_error(NCBadInput,
                      'Invalid JSON query key: ""'):
        test('')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc",""]'):
        test('mmc','')

    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: "-bla"'):
        test('mmc','-bla')
    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: "-bla"'):
        test('mmc','-bla')
    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: " -bla"'):
        test('mmc',' -bla')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","bla-bla"]'):
        test('mmc','bla-bla')

    with ensure_error(NCBadInput,
                      "Invalid query (not all entries are strings):"
                      " ('foo', 117)"):
        test('foo',117)

    #An argument following "RAWDATA" is exempt from the character checks (so
    #here the query fails later, on the argument count):
    with ensure_error(NCBadInput,
                      'Forbidden character "#" in query! Problem found in:'
                      ' "#"'):
        test('test','wl2ekin','#')
    with ensure_error(NCBadInput,'Invalid test query: ["test","wl2ekin",'
                      '"RAWDATA","#"] (correct usage: ["test","wl2ekin",'
                      'VALUE_ANGSTROM])'):
        test('test','wl2ekin','RAWDATA','#')

    for a in ('bar'+'\x07','\x07','\x07\x07','\x07'+'bar'):
        with ensure_error(NCBadInput,
                          'Can not use character 0x07 (ASCII BEL) in JSON'
                          ' query strings'):
            test('foo',a)

    test_cli('--help')
    test_cli('version')
    test_cli('version','-j')
    fn = 'tmp_test_file.json'
    if pathlib.Path(fn).exists():
        pathlib.Path(fn).unlink()
    test_cli('version','-o',fn)
    print(f"Contents of {fn}:",'='*40)
    print(pathlib.Path(fn).read_text())
    with ensure_error(RuntimeError,
                      'File already exists (use --force to overwrite):'
                      ' tmp_test_file.json'):
        test_cli('list','-o',fn)

    test_cli('list','-o',fn,'--force')
    print(f"Contents of {fn}:",'='*40)
    print(pathlib.Path(fn).read_text())

    with ensure_error(NCBadInput,
                      'Invalid JSON query key: "foo"'):
        test_cli('foo')

    test_cli('list')
    test_cli('mmc','list')
    test_cli('util','list')
    test_cli('test','list')
    test_cli('sab','list')
    test_cli('vdos','list')
    test_cli('matcfg','list')
    test_cli('matcfg','normalise','stdlib::Al_sg225.ncmat ; temp = 20C')
    test_matcfg()
    for k in ('sab','vdos'):
        with ensure_error(NCBadInput,
                          f'Invalid {k.upper()} JSON query: ["{k}","list","x"]'
                          f' (no arguments should come after: ["{k}","list"])'):
            test(k,'list','x')
    with ensure_error(NCBadInput,
                      'Invalid SAB JSON query: ["sab","surveyor","@1@2@x","@1@2"]'
                      ' (correct usage: ["sab","surveyor",ALPHAGRID,BETAGRID],'
                      " with grids encoded as '@val1@..@valn')"):
        test('sab','surveyor','@1@2@x','@1@2')
    test_atomsites()
    test_cif()
    test_uniquestructprov()
    test_cli('test','wl2ekin','1.8')
    test_cli('test','ekin2wl','0.025')
    test_cli('test','eqrefl','relaxed','149','1,2,3','1,1,0')
    test_cli('test','eqrefl','relaxed','150','1,2,3','1,1,0')
    test_cli('test','eqrefl','strict','166:r','1,2,3','1,1,0')
    test_cli('sgsym','info','strict','166:r')
    test_cli('sgsym','info','strict','62:cab')
    test_cli('sgsym','info','relaxed','1')
    test_cli('sgsym','info','relaxed','227')#(default setting)
    with ensure_error(NCBadInput,'Invalid sgsym query: ["sgsym","info"]'
                      ' (correct usage: ["sgsym","info","strict|relaxed",SG])'):
        test_cli('sgsym','info')
    with ensure_error(NCBadInput,'Invalid space group specification: "227" (space'
                      ' group 227 has several settings, so one of them must be'
                      ' specified: 227:1 (default), 227:2)'):
        test_cli('sgsym','info','strict','227')
    test_cli('sgsym','list')
    #Settings from operations and symbols:
    test_cli('sgsym','fromops','1',*rawlist('x,y,z','-x,-y,-z'))
    test_cli('sgsym','fromops','1',*rawlist('x,y,z','-x+1/2,-y+1/2,-z+1/2'))#(origin shift)
    test_cli('sgsym','fromops','0',*rawlist('x,y,z','x+1/5,y,z'))#(not a multiple of 1/24)
    #An operation which can not be parsed, with DOTHROW=0 and 1:
    test_cli('sgsym','fromops','0',*rawlist('x,y,z','x,y,q'))
    test_cli('sgsym','fromops','0',*rawlist('x,y,z','x+y,y,z'))
    with ensure_error(NCBadInput,'Invalid symmetry operation "x,y,q":'
                      ' unexpected character'):
        test_cli('sgsym','fromops','1',*rawlist('x,y,z','x,y,q'))
    test_cli('sgsym','fromhall','RAWDATA','-P 2ybc')
    test_cli('sgsym','fromhm','RAWDATA','P 21/c')
    test_cli('sgsym','fromhm','RAWDATA','R -3 m')
    test_cli('sgsym','fromhm','RAWDATA','C m c e')
    test_cli('sgsym','fromhm','RAWDATA','C -1')#(no tabulated setting)
    test_cli('sgsym','fromnumber','227')
    test_cli('sgsym','fromnumber','14')
    #A setting, or with DOTHROW=0 the reason it is none:
    test_cli('sgsym','setting','0','RAWDATA','70:1')
    test_cli('sgsym','setting','1','RAWDATA','225')
    test_cli('sgsym','setting','0','RAWDATA','70:9')
    test_cli('sgsym','setting','0','RAWDATA','70')
    test_cli('sgsym','setting','0','RAWDATA','F d d d')
    with ensure_error(NCBadInput,'Invalid space group specification: "70:9"'
                      ' (not a setting of space group 70, the settings are:'
                      ' 70:1 (default), 70:2)'):
        test_cli('sgsym','setting','1','RAWDATA','70:9')
    with ensure_error(NCBadInput,'Invalid sgsym query: ["sgsym","setting","2",'
                      '"RAWDATA","70:1"] (correct usage: ["sgsym","setting",'
                      '0|1,"RAWDATA",SG])'):
        test_cli('sgsym','setting','2','RAWDATA','70:1')
    test_cli('sgsym','fromnumber','225')
    for bad in ( '0', '231', 'x', '1.0' ):
        with ensure_error( NCBadInput, 'Invalid sgsym query: ["sgsym",'
                           f'"fromnumber","{bad}"] (correct usage: ["sgsym",'
                           '"fromnumber",N], with N in 1..230)' ):
            ncquery( [ 'sgsym', 'fromnumber', bad ] )
    test( 'sgsym','fromops','1','RAWDATA','|x,y,z|-x,-y,-z|x+1/2,y,z' )
    for bad in ( [ 'x,y,z' ], [ '1', 'x,y,z' ], [ '1', 'RAWDATA' ],
                 [ '1', 'RAWDATA', '' ], [ '1', 'RAWDATA', '@x,y,z', '@-x,-y,-z' ] ):
        with ensure_error(NCBadInput,'Invalid sgsym query: '
                          + json.dumps( [ 'sgsym', 'fromops' ] + bad,
                                        separators = ( ',', ':' ) )
                          + ' (correct usage: ["sgsym","fromops",0|1,"RAWDATA",'
                          'OPLIST])'):
            ncquery( [ 'sgsym', 'fromops', *bad ] )
    #Cell parameters snapped to the constraints of a setting (printed values
    #with rounding errors with lenient tolerances):
    def cellargs( *v ):
        return [ e for x in v for e in ( 'RAWDATA', x ) ]
    test_cli('sgsym','cell','lenient','0','225',*cellargs('4.123(2)','4.124','4.1234','90','90.0','90'))
    test_cli('sgsym','cell','strict','1','194',*cellargs('3','3','5','90','90','120'))
    test_cli('sgsym','cell','strict','0','225',*cellargs('4.123','4.124','4.123','90','90','90'))
    with ensure_error(NCBadInput,'Lattice parameters a and b must be equal for'
                      ' this space group (got a=4.123 and b=4.124)'):
        test_cli('sgsym','cell','strict','1','225',*cellargs('4.123','4.124','4.123','90','90','90'))
    with ensure_error(NCBadInput,'Invalid sgsym query: ["sgsym","cell","lenient",'
                      '"0","225","4.1","4.1","4.1","90","90","90"] (correct usage:'
                      ' ["sgsym","cell","strict|lenient",0|1,SG,"RAWDATA",a,'
                      '"RAWDATA",b,"RAWDATA",c,"RAWDATA",alpha,"RAWDATA",beta,'
                      '"RAWDATA",gamma])'):
        test_cli('sgsym','cell','lenient','0','225','4.1','4.1','4.1','90','90','90')
    with ensure_error(NCBadInput,'Invalid sgsym query: ["sgsym","cell","lenient",'
                      '"0","225","RAWDATA","4.1","RAWDATA","4.1","RAWDATA","4.1",'
                      '"RAWDATA","90","RAWDATA","90","RAWDATA","x"] (correct usage:'
                      ' ["sgsym","cell","strict|lenient",0|1,SG,"RAWDATA",a,'
                      '"RAWDATA",b,"RAWDATA",c,"RAWDATA",alpha,"RAWDATA",beta,'
                      '"RAWDATA",gamma])'):
        test_cli('sgsym','cell','lenient','0','225',*cellargs('4.1','4.1','4.1','90','90','x'))
    with ensure_error(NCBadInput,'Invalid sgsym query: ["sgsym","cell","lenient",'
                      '"2","225","RAWDATA","4.1","RAWDATA","4.1","RAWDATA","4.1",'
                      '"RAWDATA","90","RAWDATA","90","RAWDATA","90"] (correct usage:'
                      ' ["sgsym","cell","strict|lenient",0|1,SG,"RAWDATA",a,'
                      '"RAWDATA",b,"RAWDATA",c,"RAWDATA",alpha,"RAWDATA",beta,'
                      '"RAWDATA",gamma])'):
        test_cli('sgsym','cell','lenient','2','225',*cellargs('4.1','4.1','4.1','90','90','90'))
    for bad in ( [], [ 'P 21/c' ] ):
        with ensure_error(NCBadInput,'Invalid sgsym query: '
                          + json.dumps( [ 'sgsym', 'fromhm' ] + bad,
                                        separators = ( ',', ':' ) )
                          + ' (correct usage: ["sgsym","fromhm","RAWDATA",'
                          'SYMBOL])'):
            ncquery( [ 'sgsym', 'fromhm', *bad ] )
    with ensure_error(NCBadInput,'Invalid sgsym query: ["sgsym","foo"]'
                      ' (unknown key)'):
        test_cli('sgsym','foo')

    with ensure_error(ArgumentError,
                      'the following arguments are required: STR'):
        test_cli()
    with ensure_error(NCBadInput,
                      'Invalid JSON query key: "foo"'):
        test_cli('foo')
    with ensure_error(NCBadInput,
                      'Invalid JSON query key: ""'):
        test_cli('')

    with ensure_error(NCBadInput,
                      'Invalid JSON query key: ""'):
        test_cli('')
    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc",""]'):
        test_cli('mmc','')

    with ensure_error(NCBadInput,
                      'Query items can not start with a "-" character.'
                      ' Problem found in: " -bla"'):
        test_cli('mmc',' -bla')

    with ensure_error(NCBadInput,
                      'Invalid MiniMC JSON query: ["mmc","bla-bla"]'):
        test_cli('mmc','bla-bla')

    for a in ('bar'+'\x07','\x07','\x07\x07','\x07'+'bar'):
        with ensure_error(NCBadInput,
                          'Can not use character 0x07 (ASCII BEL) in JSON'
                          ' query strings'):
            test_cli('foo',a)

    test_huge_arrays()
    test_info_queries()

def test_info_queries():
    #Object queries on an Info object (cached, shared results):
    import pprint

    import NCrystalDev as NC
    info = NC.createInfo('Al_sg225.ncmat;dcutoff=0.5')
    print(">>> Info query ['list']:")
    pprint.pprint( info._query('list') )
    print(">>> Info query ['structure']:")
    pprint.pprint( info._query('structure') )
    assert info._query('structure') is info.structure_info
    assert info._query('structure') is info.getStructureInfo()
    diamond2 = [ '1/8 1/8 1/8', '1/8 5/8 5/8', '5/8 1/8 5/8', '5/8 5/8 1/8',
                 '7/8 7/8 7/8', '7/8 3/8 3/8', '3/8 7/8 3/8', '3/8 3/8 7/8' ]
    for name, sg, pos, angles in [ ( 'origin2', 227, diamond2, '90 90 90' ),
                                   ( 'rhombo', 166, ['0 0 0'], '60 60 60' ) ]:
        atoms = '\n'.join( f'  Si {p}' for p in pos )
        NC.registerInMemoryFileData( f'{name}.ncmat', f"""NCMAT v7
@CELL
  lengths 5 5 5
  angles {angles}
@SPACEGROUP
  {sg}
@ATOMPOSITIONS
{atoms}
@DYNINFO
  element Si
  fraction 1
  type vdosdebye
  debye_temp 400
""")
        si = NC.createInfo(f'{name}.ncmat;inelas=0').structure_info
        print( f">>> {name}: spacegroup {si['spacegroup']}"
               f" setting {si['spacegroup_setting']}" )
    gas = NC.createInfo('freegas::Ar/1e-3perAa3')
    assert not gas.hasStructureInfo() and gas._query('structure') is None
    with ensure_error(NCBadInput,
                      'Invalid Info JSON query key: "foo"'):
        info._query('foo')
    with ensure_error(NCBadInput,
                      'Invalid Info JSON query (no arguments should'
                      ' follow key "structure")'):
        info._query('structure','bar')

def test_huge_arrays():
    #In the mode where huge arrays are transferred separately from the JSON
    #data, arrays which are never collected must be discarded when caches are
    #cleared, and when the mode is enabled the next time.
    import json

    from NCrystalDev._chooks import _get_raw_cfcts
    from NCrystalDev.core import clearCaches
    rawfct = _get_raw_cfcts()
    query = [ 'sab', 'sglcell', '@40.0', '@158.28456094789',
              '@160.60010716438572', '@-2.2701714000978403e-14',
              '@-2.120250108329995e-17', '@0.02833179339754119',
              '@0.027532855700519738', '@0.027532855700519426',
              '@0.026748257342144074', '3', '12462' ]
    def find_keys( data ):
        if isinstance( data, dict ):
            data = list( data.values() )
        if isinstance( data, list ):
            return [ k for e in data for k in find_keys(e) ]
        is_key = ( isinstance( data, str )
                   and data.startswith('__ncrystal__dblarray::') )
        return [ data ] if is_key else []
    def query_without_collecting():
        rawfct['enablejsonarr'](True)
        res = rawfct['jsonquery'](query)
        rawfct['enablejsonarr'](False)
        keys = find_keys( json.loads( res ) )
        assert len(keys) >= 2
        return keys
    print()
    keys = query_without_collecting()
    assert len( rawfct['getjsonarr']( keys[0] ) ) > 0
    print("Collected huge array")
    clearCaches()
    def errmsg( key ):
        return 'Invalid key for JSON Huge Vector: "{}"'.format(
            key.split('::')[1] )
    with ensure_error(NCBadInput,errmsg(keys[1])):
        rawfct['getjsonarr']( keys[1] )
    keys = query_without_collecting()
    rawfct['enablejsonarr'](True)
    rawfct['enablejsonarr'](False)
    with ensure_error(NCBadInput,errmsg(keys[0])):
        rawfct['getjsonarr']( keys[0] )
    res = ncquery( query, huge_arrays = True )
    assert not find_keys( res )
    print("Normal query with huge arrays worked afterwards")


def test_atomsites():
    #["sgsym","positions",...] expands sites, ["sgsym","sites",...]
    #groups positions into sites (diamond in origin choice 2):
    site = 'Si,0.125,0.125,0.125'
    cell = [ 'RAWDATA', '5.43,5.43,5.43,90,90,90' ]
    r = ncquery( [ 'sgsym', 'positions','strict','227:2','strict',*cell,
                   *rawlist(site)] )
    s, = r['sites']
    print( f'227:2 site {site} -> {len(s["positions"])} positions' )
    assert s['multiplicity'] == len( s['positions'] ) == 8
    assert s['split'] is None and s['uij'] is None and s['uijs'] is None
    pos = [ f'Si,{x!r},{y!r},{z!r}' for x,y,z in s['positions'] ]
    test( 'sgsym','sites','strict','227:2','strict',*cell,*rawlist(*pos) )
    #Relaxed: the setting is deduced from the positions:
    r = ncquery( [ 'sgsym', 'sites','relaxed','227','strict',*cell,
                   *rawlist(*pos)] )
    print( 'relaxed 227 -> setting', r['spacegroup_setting'] )
    assert r['spacegroup_setting'] == '227:2'
    with ensure_error( NCBadInput, 'Atom positions are not consistent with'
                       ' space group 227:1 (largest deviation of a symmetry'
                       ' image from an atom of the same type: 0.25 in'
                       ' fractional coordinates)' ):
        ncquery( [ 'sgsym', 'sites','strict','227:1','strict',*cell,
                   *rawlist(*pos)] )
    with ensure_error( NCBadInput, 'Invalid sgsym query: ["sgsym","sites",'
                       '"227:2","Si,0.125,0.125,0.125"] (correct usage:'
                       ' ["sgsym","sites","strict|relaxed",SG,"strict|lenient",'
                       '"RAWDATA","a,b,c,alpha,beta,gamma","RAWDATA",'
                       '"@lbl,x,y,z[,occupancy][,u11,u22,u33,u12,u13,u23]@..."])' ):
        ncquery( [ 'sgsym', 'sites','227:2',site] )
    with ensure_error( NCBadInput, 'Invalid sgsym query: ["sgsym","positions",'
                       '"strict","227:2","strict","RAWDATA",'
                       '"5.43,5.43,5.43,90,90,90","RAWDATA",'
                       '"@Si,0.125,0.125,0.125,1,2"] (correct usage:'
                       ' ["sgsym","positions","strict|relaxed",SG,'
                       '"strict|lenient","RAWDATA","a,b,c,alpha,beta,gamma",'
                       '"RAWDATA",'
                       '"@lbl,x,y,z[,occupancy][,u11,u22,u33,u12,u13,u23]@..."])' ):#(6 fields)
        ncquery( [ 'sgsym', 'positions','strict','227:2','strict',*cell,
                   *rawlist('Si,0.125,0.125,0.125,1,2')] )
    #A split site (two positions 0.054 Aa apart, about a mirror plane of Pm):
    pm = [ 'strict','6:b','lenient','RAWDATA','5,5.4,6,90,90,90' ]
    r = ncquery( [ 'sgsym', 'positions', *pm, *rawlist('O,0.1,0.005,0.2,0.5') ] )
    print( 'split site:', rnd( r['sites'] ) )
    with ensure_error( NCBadInput, 'At most one of the 2 symmetry-equivalent'
                       ' positions of the site (0.1, 0.005, 0.2) within 0.5 Aa'
                       ' (the closest 0.054 Aa apart) can be occupied, so its'
                       ' occupancy must be at most 1/2 (got 1). For sites on'
                       ' special positions, specify coordinates with'
                       ' sufficient precision, preferably as exact fractions'
                       ' like 1/3.' ):
        ncquery( [ 'sgsym', 'positions', *pm, *rawlist('O,0.1,0.005,0.2') ] )
    test_uij()

def test_cif():
    test( 'cif', 'list' )
    test( 'cif', 'numbers', *rawlist( '0.0283', '2344', '1.234e-2(6)', '3400e-6',
                                  '-0.000', '0', '0.000(2)', '?', '.', '1,5',
                                  '1.2(3)e-4', '1e999' ) )
    test( 'cif', 'elements', 'typesymbol', *rawlist( 'O2-', 'Fe3+', 'CO', 'Co',
                                                 'D', 'T+', 'Mn3.614+', 'Fe0',
                                                 'O-2', 'D+1', 'Ag+1', '?',
                                                 'Xx', 'O2+-', 'O+-2',
                                                 'O-2-' ) )
    test( 'cif', 'elements', 'label', 'RAWDATA', 'Ca Co0.5 O3 C2 H Cl',
          *rawlist( 'CaM1', 'OA1', 'Fe12', 'D1', 'Ow1', 'X1',
                                            'C', 'Co1', 'CO1', 'Cl(1)', 'Cs1',
                                            'H01a', 'CL1', 'CA1' ) )
    test( 'cif', 'elements', 'label', 'RAWDATA', 'Al2O3' )
    test( 'cif', 'elements', 'label', 'RAWDATA', 'not a formula',
          *rawlist( 'Al1' ) )
    test( 'cif', 'numbers' )
    for bad in ( [], [ 'label' ], [ 'label', 'Al2O3' ], [ 'bla', 'X' ],
                 [ 'typesymbol', 'Al' ], [ 'typesymbol', 'RAWDATA', '' ],
                 [ 'typesymbol', *rawlist( 'Al' ), 'Fe' ] ):
        with ensure_error( NCBadInput, 'Invalid cif query: '
                           + json.dumps( [ 'cif', 'elements' ] + bad,
                                         separators = ( ',', ':' ) )
                           + ' (correct usage: ["cif","elements","typesymbol",'
                           '"RAWDATA",LIST] or ["cif","elements","label",'
                           '"RAWDATA",FORMULA,"RAWDATA",LIST])' ):
            ncquery( [ 'cif', 'elements', *bad ] )
    for bad in ( [ '1.0' ], [ 'RAWDATA', '' ] ):
        with ensure_error( NCBadInput, 'Invalid cif query: '
                           + json.dumps( [ 'cif', 'numbers' ] + bad,
                                         separators = ( ',', ':' ) )
                           + ' (correct usage: ["cif","numbers","RAWDATA",'
                           'LIST])' ):
            ncquery( [ 'cif', 'numbers', *bad ] )
    with ensure_error( NCBadInput,
                       'Invalid cif query: ["cif","bla"] (unknown key)' ):
        ncquery( [ 'cif', 'bla' ] )

def test_uniquestructprov():
    #Items of the CIF block with a crystal structure, as [text,line] (cf. the
    #unit test app_cifitems):
    cif = """data_info
_journal_name_full 'From data_info'
data_s
_cell_length_a 4.05(1)
_cell_length_b 4.05(1)
_cell_length_c 4.05(1)
_cell_angle_beta 90.0
_space_group_name_H-M_alt 'F m -3 m'
_symmetry_space_group_name_H-M 'F m 3 m'
_space_group_IT_number 225
_chemical_formula_sum 'Al'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
Al1 Al 0 0 0 ?
Al2 Al 0.5 0.5 0.5 0.01
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
Al1 ? ? ? ? ? ?
Al2 0.01 0.01 0.01 0 0 0
_publ_section_title 'A title'
loop_ _publ_author_name 'Doe, J.' 'Roe, R.'
_cod_database_code 1234567
"""
    test( 'cif', 'uniquestructprov', 'permissive', 'RAWDATA', 'q.cif',
          'RAWDATA', cif )
    for name, shown in ( ( 'q.cif', 'q.cif' ), ( '', 'raw CIF data' ) ):
        with ensure_error( NCBadInput, f'Invalid CIF data in "{shown}"'
                           ' (data_s, line 29): anisotropic displacements of'
                           ' "Al1" are all "?" or ".". Use permissive mode to'
                           ' ignore such entries.' ):
            ncquery( [ 'cif', 'uniquestructprov', 'strict', 'RAWDATA', name,
                       'RAWDATA', cif ] )
    import NCrystalDev as NC
    NC.registerInMemoryFileData( 'query_test.cif', """data_x
_cell_length_a 3
_cell_length_b 3
_cell_length_c 3
_atom_site_label Fe1
_atom_site_fract_x 0
_atom_site_fract_y 0
_atom_site_fract_z 0
""" )
    test( 'cif', 'uniquestructprov', 'strict', 'query_test.cif' )
    #Synthetic file with the publication only in a data_global block:
    from NCTestUtils.dirs import test_data_dir
    r = ncquery( [ 'cif', 'uniquestructprov', 'strict',
                   str( test_data_dir / 'cif_globalblock.cif' ) ] )
    print( 'cif_globalblock.cif: block', r['block'], 'provenance_block',
           r['provenance_block'] )
    import pprint
    pprint.pprint( r['provenance'] )
    for bad in ( [], [ 'strict' ], [ '1', 'x.cif' ],
                 [ 'strict', 'RAWDATA' ], [ 'strict', 'RAWDATA', 'x' ],
                 [ 'strict', 'RAWDATA', 'n', 'x', 'y' ],
                 [ 'strict', 'x', 'y' ] ):
        with ensure_error( NCBadInput, 'Invalid cif query: '
                           + json.dumps( [ 'cif', 'uniquestructprov' ] + bad,
                                         separators = ( ',', ':' ) )
                           + ' (correct usage: ["cif","uniquestructprov",'
                           'PERMISSIVE,NAME] or ["cif","uniquestructprov",'
                           'PERMISSIVE,"RAWDATA",NAME,"RAWDATA",TEXT], with'
                           ' PERMISSIVE "strict" or "permissive")' ):
            ncquery( [ 'cif', 'uniquestructprov', *bad ] )

def rnd( o ):
    #Round floats (also nested) to 10 digits, for platform-stable printouts:
    if isinstance( o, float ):
        return float( f'{o:.10g}' )
    if isinstance( o, dict ):
        return { k : rnd(v) for k,v in o.items() }
    if isinstance( o, (list,tuple) ):
        return [ rnd(v) for v in o ]
    return o

def test_uij():
    #Anisotropic displacements (uij, as Ustd): on the mirror plane of Pm
    #(y=0), U12 and U23 must vanish:
    pm = [ 'strict','6:b','strict','RAWDATA','5,5.4,6,90,95,90' ]
    sym = 'O,0.1,0,0.2,0.012,0.009,0.015,0,0.003,0'
    r = ncquery( [ 'sgsym', 'positions', *pm, *rawlist( sym ) ] )
    s, = r['sites']
    print( 'uij on the mirror plane:', rnd( s['uij'] ), 'multiplicity', s['multiplicity'] )
    assert s['multiplicity'] == 1 and s['uijs'] == [ s['uij'] ]
    #Slightly off (U12 2% of the largest component) is snapped, more is an error:
    r = ncquery( [ 'sgsym', 'positions', *pm,
                   *rawlist( 'O,0.1,0,0.2,0.012,0.009,0.015,0.0003,0.003,0' ) ] )
    print( 'U12=0.0003 snapped to:', rnd( r['sites'][0]['uij'] ) )
    assert r['sites'][0]['uij'][3] == 0.0
    with ensure_error( NCBadInput, 'The uij values of the site "O" (0.1, 0, 0.2)'
                       ' do not obey the site symmetry' ):
        ncquery( [ 'sgsym', 'positions', *pm,
                   *rawlist( 'O,0.1,0,0.2,0.012,0.009,0.015,0.0012,0.003,0' ) ] )
    #With lenient tolerances, the rounding errors of printed numbers are also
    #accepted (here U12=0.0012 printed with rounding error 0.0005 each):
    pml = [ 'strict','6:b','lenient','RAWDATA','5,5.4,6,90,95,90' ]
    r = ncquery( [ 'sgsym', 'positions', *pml,
                   *rawlist( 'O,0.1,0,0.2,0.012,0.009,0.015,0.001,0.003,0' ) ] )
    print( 'lenient, printed U12=0.001 snapped to:', rnd( r['sites'][0]['uij'] ) )
    with ensure_error( NCBadInput, 'The uij values of the site "O" (0.1, 0, 0.2)'
                       ' do not obey the site symmetry' ):
        ncquery( [ 'sgsym', 'positions', *pml,
                   *rawlist( 'O,0.1,0,0.2,0.0120,0.0090,0.0150,0.0012,0.0030,0.0000' ) ] )
    with ensure_error( NCBadInput, 'The uij values of the site "O" (0.1, 0, 0.2)'
                       ' do not form a positive definite matrix' ):
        ncquery( [ 'sgsym', 'positions', *pm, *rawlist( 'O,0.1,0,0.2,0.012,-0.009,0.015,0,0.003,0' ) ] )
    #A general position: the uij at the mirror image (y -> -y), and back:
    r = ncquery( [ 'sgsym', 'positions', *pm,
                   *rawlist( 'O,0.1,0.3,0.2,0.5,0.012,0.009,0.015,0.002,0.003,0.001' ) ] )
    s, = r['sites']
    print( 'uij at the two positions of a general site:', rnd( s['uijs'] ) )
    assert s['uijs'][1][3] < 0.0 and s['uijs'][1][5] < 0.0
    pos = [ 'O,' + ','.join( repr(v) for v in p ) + ',0.5,'
            + ','.join( repr(v) for v in u ) for p, u in zip( s['positions'], s['uijs'] ) ]
    r = ncquery( [ 'sgsym', 'sites', *pm, *rawlist( *pos ) ] )
    print( 'grouped back into sites:', rnd( r['sites'] ) )
    assert r['sites'][0]['multiplicity'] == 2
    #Tensors at the two positions which are not symmetry images of each other:
    bad = [ pos[0], pos[1].replace( ',-0.002', ',0.002' ) ]
    assert bad[1] != pos[1]
    with ensure_error( NCBadInput, 'The uij values of the symmetry-equivalent'
                       ' positions of "O" (0.1, 0.3, 0.2) are not symmetry'
                       ' images of each other' ):
        ncquery( [ 'sgsym', 'sites', *pm, *rawlist( *bad ) ] )
    #Uij must be given for all or none of the entries of a label:
    with ensure_error( NCBadInput, 'Anisotropic displacements (uij) must be'
                       ' given for all or none of the entries of O' ):
        ncquery( [ 'sgsym', 'sites', *pm,
                   *rawlist( pos[0], 'O,0.1,0.7,0.2,0.5' ) ] )
    #Printed CIF numbers with uncertainties are accepted:
    r = ncquery( [ 'sgsym', 'positions', *pml,
                   *rawlist( 'O,0.1000(3),0,0.2000(2),0.0120(4),0.0090(3),'
                             '0.0150(5),0,0.0030(2),0' ) ] )
    print( 'printed CIF numbers:', rnd( r['sites'][0]['uij'] ) )
    #Ueq and its rounding error from the printed tensor and cell (here the
    #angles 90 printed without uncertainty are exact):
    for c in ( '5,5.4,6,90,95,90', '5.0(1),5.4(2),6.0(1),90,95.0(3),90',
               '5.0(1),5.4(2),6.0(1),90.0(1),95.0(3),90' ):
        r = ncquery( [ 'sgsym', 'positions', 'strict', '6:b', 'lenient',
                       'RAWDATA', c,
                       *rawlist( 'O,0.1,0.3,0.2,0.0120(4),0.0090(3),0.0150(5),'
                                 '0.0020(2),0.0030(2),0.0010(1)' ) ] )
        print( 'ueq in cell', c, ':', [ float( f'{v:.6g}' )
                                         for v in r['sites'][0]['ueq'] ] )
    r = ncquery( [ 'sgsym', 'positions', 'strict', '6:b', 'strict',
                   'RAWDATA', '5,5.4,6,90,95,90', *rawlist( 'O,0.1,0.3,0.2' ) ] )
    assert r['sites'][0]['ueq'] is None

if __name__ == '__main__':
    main()
