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

# Test the "ncrystal atomdb" command-line tool.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev.cli as nc_cli
from NCrystalDev._common import capture_print_ctxmgr
from NCTestUtils.common import ensure_error
import argparse
import shlex

def run( *args, show = True ):
    print(f"============= CLI >>atomdb {shlex.join(args)}<< =============")
    with capture_print_ctxmgr() as cap:
        nc_cli.run('atomdb',*args)
    out = ''.join(cap.data)
    if show:
        print(out,end='')
    return out

def bad( errmsg, *args ):
    with ensure_error(argparse.ArgumentError,errmsg):
        run(*args)

def main():
    run('--help')
    run('--cond-help')
    run('He','help')
    run('He','b10')
    run('absxs > 1000','A=0')
    run('Li')
    run('H','--json')
    run('He','N15','--ncmat')
    run('D','T')#aliases of H2 and H3
    run('absxs > 1e9')
    #All examples in the help must work:
    from NCrystalDev._cli_atomdb import _cond_examples, _examples
    for e,_ in _examples + _cond_examples:
        out = run( *shlex.split(e), show = False )
        print(f'  -> {len(out.splitlines())} lines of output')
    bad(('Invalid property "sg" in condition "sg>1" (must be one of: Z, A,'
         ' mass, cohsl, cohxs, incohxs, scatxs, absxs)'),'sg>1')
    bad('Invalid condition "scatxs>foo" (value "foo" is not a number)',
        'scatxs>foo')
    bad(('Invalid property "COHSL" in condition "COHSL<0" (names are'
         ' case-sensitive, did you mean "cohsl"?)'),'COHSL<0')
    bad(('Invalid condition "label=He" ("label" can not be used in'
         ' conditions)'),'label=He')
    bad(('Unknown element or isotope "foo" (specify elements or isotopes,'
         ' e.g. He B10)'),'He','foo')
    bad('Do not specify both --json and --ncmat.','--json','--ncmat')
    bad(('Invalid property "a" in condition "a=0" (names are case-sensitive,'
         ' did you mean "A"?)'),'a=0')

if __name__ == '__main__':
    main()
