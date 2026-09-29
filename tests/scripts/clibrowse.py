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

# Test the "ncrystal browse" command-line tool.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
import NCrystalDev.cli as nc_cli
from NCrystalDev._common import capture_print_ctxmgr
from NCTestUtils.common import ensure_error
import os
import shlex

_ncmat_withcomments = """NCMAT v7
#
#   My test material (for testing only).
#   Second line of description which is long enough that it will have to be
#   truncated in the listing.
#
#   Some more comments mentioning Togo and VDOS.
#
@DENSITY
  1 g_per_cm3
@DYNINFO
  element H
  fraction 1
  type vdosdebye
  debye_temp 300
"""

_ncmat_nocomments = """NCMAT v7
@DENSITY
  2 g_per_cm3
@DYNINFO
  element O
  fraction 1
  type vdosdebye
  debye_temp 400
"""

_ncmat_crystal = """NCMAT v7
# A small crystal.
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

def run( *args, show = True ):
    print(f"============= CLI >>browse {shlex.join(args)}<< =============")
    with capture_print_ctxmgr() as cap:
        nc_cli.run('browse',*args)
    out = ''.join(cap.data)
    out = out.replace('\x1b','<ESC>')#make color codes visible in log
    if show:
        print(out,end='')
    return out

_color_envvars = ('NO_COLOR','FORCE_COLOR','CLICOLOR_FORCE',
                  'GREP_COLORS','TERM')

def setenv( **kw ):
    for k,v in kw.items():
        if v is None:
            os.environ.pop(k,None)
        else:
            os.environ[k] = v

def test_colors():
    #NB: Output is captured, so "auto" mode means no colors unless forced.
    run('-s','togo','-s','my','-f','virtual')
    run('-s','togo','-s','my','-f','virtual','--color=always')
    run('-s','togo','-f','virtual','--colour=yes')
    run('-s','togo','-f','virtual','--color=never')
    run('-s','togo','--color','mytestmat')#plain --color means auto
    setenv( FORCE_COLOR = '1' )
    run('-s','togo','-f','virtual')
    setenv( NO_COLOR = '1' )
    run('-s','togo','-f','virtual')
    run('-s','togo','-f','virtual','--color=always')
    setenv( NO_COLOR = None, FORCE_COLOR = None,
            GREP_COLORS = 'sl=1:ms=01;32:ln=35' )
    run('-s','togo','-f','virtual','--color=always')
    setenv( GREP_COLORS = None )
    #Highlighting, with overlapping matches merged:
    run('-s','tog','-s','ogo','-f','virtual','--color=always')
    import argparse
    with ensure_error(argparse.ArgumentError,
                      "argument --color/--colour: invalid choice: 'blue'"
                      " (choose from 'always', 'auto', 'force', 'if-tty',"
                      " 'never', 'no', 'none', 'tty', 'yes')"):
        run('--color=blue')

def test_lazlau():
    #Non-NCMAT data, both in-memory and on-disk (in current directory):
    import pathlib

    from NCrystalDev.mcstasutils import cfgstr_2_hkl
    from NCTestUtils.common import work_in_tmpdir
    def lazlau( fmt ):
        return '\n'.join( cfgstr_2_hkl( cfgstr = ( 'stdlib::Al_sg225.ncmat;'
                                                   'dcutoff=1.0' ),
                                        tgtformat = fmt, verbose = False,
                                        fp_format = '%.8g' ) ) + '\n'
    NC.registerInMemoryFileData('mem.laz',lazlau('laz'))
    NC.registerInMemoryFileData('mem.lau',lazlau('lau'))
    NC.enableRelativePaths(True)
    with work_in_tmpdir():
        pathlib.Path('disk.laz').write_text(lazlau('laz'))
        pathlib.Path('disk.lau').write_text(lazlau('lau'))
        run('.la')#relpath entries only with -f relpath
        run('.la','-f','relpath')
        run('-s','disk','-f','relpath')
        run('-f','relpath','sg=225')
        run('-f','virtual','sg==225')
        #Staged flushed markers: this area was observed to hang on the
        #Windows Debug CI legs, and these localise the next such hang:
        print('lazlau: leaving tmpdir',flush=True)
    print('lazlau: tmpdir cleaned up',flush=True)
    NC.enableRelativePaths(False)
    print('lazlau: relative paths disabled again',flush=True)

def test_physics():
    run('formula=o','-f','virtual')#case-insensitive
    run('crystal=yes','sg=225','natoms=4','-f','virtual')
    run('sg > 200','-f','virtual')#missing values never match
    run('crystal=no','density<1.5','-f','virtual')
    run('absxs>0.2','nelements=1','state=solid','-f','virtual')
    run('formula=D2O','-f','stdlib','Liquid')
    run('state=liquid','nelements=2','-f','stdlib','Liquid')
    run('formula!=O','nphases=1','-f','virtual')
    run('crystalsystem=CUBIC','braggthreshold>4','-f','virtual')
    run('vdosdebye=YES','crystal=off','-f','virtual')#yes/no spellings
    #Element names also select data containing them (in addition to names):
    run('Al','-f','virtual')
    run('O','nelements=1','-f','virtual')
    run('al','-f','virtual')#not an element name (lower case)
    #Suggestions:
    run('-f','virtual','mycrystl')
    run('-f','virtual','mycrystal','absxs > 100')#no suggestions
    run('-f','virtual','qwertyzzz')
    #JSON:
    import json
    out = run('-f','virtual','mycrystal','--json',show=False)
    import NCTestUtils.stabilise_ncpprint # noqa F401
    import NCrystalDev._common as nc_common
    nc_common.ncpprint( json.loads(out) )#FP precision clipped
    #No truncation:
    run('-f','virtual','mytestmat','-s','truncated')#never shortened
    import argparse
    def bad_cond( cond, errmsg ):
        with ensure_error(argparse.ArgumentError,errmsg):
            run(cond)
    from NCrystalDev.browse import physics_props_doc
    propnames = ', '.join( n for n,d in physics_props_doc( conditions = True ) )
    bad_cond('foo>1', ('Invalid property "foo" in condition "foo>1" (must be'
                       ' one of: '+propnames+')'))
    bad_cond('COHSL<0', ('Invalid property "COHSL" in condition "COHSL<0"'
                         ' (names are case-sensitive, did you mean "cohsl"?)'))
    bad_cond('absxs>x', 'Invalid condition "absxs>x" (value "x" is not a'
             ' number)')
    bad_cond('elements=B', ('Invalid condition "elements=B" ("elements" can'
                            ' not be used in conditions)'))
    bad_cond('crystalsystem=bla', ('Invalid condition "crystalsystem=bla"'
                                   ' ("crystalsystem" must be one of:'
                                   ' triclinic, monoclinic, orthorhombic,'
                                   ' tetragonal, trigonal, hexagonal,'
                                   ' cubic)'))
    bad_cond('crystal=maybe', ('Invalid condition "crystal=maybe"'
                               ' (value "maybe" is not yes or no)'))

def test_examples():
    #All examples in the help must work:
    import shlex

    from NCrystalDev._cli_browse import _cond_examples, _examples
    for e,_ in _examples + _cond_examples:
        out = run( *shlex.split(e), show = False )
        print(f'  -> {len(out.splitlines())} lines of output')

def main():
    #Colors in output must not depend on the environment of the test:
    orig_env = { k: os.environ.get(k) for k in _color_envvars }
    setenv( **{ k: None for k in _color_envvars } )
    try:
        main_impl()
        test_colors()
    finally:
        setenv( **orig_env )

def main_impl():
    NC.removeAllDataSources()
    NC.enableStandardDataLibrary()
    NC.registerInMemoryFileData('mytestmat.ncmat',_ncmat_withcomments)
    NC.registerInMemoryFileData('othermat.ncmat',_ncmat_nocomments)
    NC.registerInMemoryFileData('notncmat.laz','whatever')
    NC.registerInMemoryFileData('mycrystal.ncmat',_ncmat_crystal)
    #Hide stdlib file by in-memory file with same name:
    NC.registerInMemoryFileData('Al_sg225.ncmat',_ncmat_nocomments)

    run('--help')
    run('--full-help')
    run('--full')#abbreviation
    run('--cond-help')
    run('Al','help')
    run('-f','help')
    run('-f','virtual')
    run('Al_sg225')
    run('stdlib::Be')
    run('-f','stdlib','sg229')
    run('-s','togo','-f','virtual')
    run('-s','TOGO','-s','vdos','-f','virtual')
    run('-s','togo','-s','nonexistentword')
    run('-s','togo','-f','stdlib','Al_sg225')
    #Search words and patterns are matched literally:
    run('-s','togo|nomatch','-f','virtual')
    run('my.*mat')
    run('mytest','othermat')#several patterns: any of them
    test_physics()
    test_lazlau()
    test_examples()
    run('-x','mytestmat.ncmat')
    with ensure_error(NC.NCFileNotFound,
                      'Could not find data: "nonexistent.ncmat"'):
        run('-x','nonexistent.ncmat')
    import argparse
    with ensure_error(argparse.ArgumentError,
                      '--extract can not be combined with other options.'):
        run('-x','mytestmat.ncmat','-f','virtual')

if __name__ == '__main__':
    main()
