
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

#Checks that the [tool.pyodide.build] flags in ncrystal_core/pyproject.toml are
#Pyodide's current default flags, except for the optimisation level, include
#and library paths, and our own additions. Run from the repository root, with
#pyodide-build and its cross-build environment installed.

import pathlib
import subprocess
import tempfile

ours_extra = { '-fvisibility=hidden', '-fvisibility-inlines-hidden' }

def get_flags( key, cwd ):
    #pyodide-build reads [tool.pyodide.build] from pyproject.toml in cwd:
    rv = subprocess.run( ['pyodide','config','get',key], cwd = cwd,
                         capture_output = True, text = True, check = True )
    return rv.stdout.strip().splitlines()[-1].replace('-s ','-s')

def normalise( flags, extra = () ):
    return { f for f in flags.split()
             if not f.startswith(('-O','-I','-L')) and f not in extra }

def main():
    ourdir = pathlib.Path('ncrystal_core').absolute()
    with tempfile.TemporaryDirectory() as tmpdir:
        defaults = { k : get_flags( k, tmpdir )
                     for k in ('cflags','cxxflags','ldflags') }
    ours = { k : get_flags( k, ourdir )
             for k in ('cflags','cxxflags','ldflags') }
    #C++ compilation uses cflags+cxxflags:
    checks = [ ( 'cflags', defaults['cflags'], ours['cflags'] ),
               ( 'cflags+cxxflags',
                 defaults['cflags']+' '+defaults['cxxflags'],
                 ours['cflags']+' '+ours['cxxflags'] ),
               ( 'ldflags', defaults['ldflags'], ours['ldflags'] ) ]
    ok = True
    for name, d, o in checks:
        print(f'{name}:\n  Pyodide default: {d}\n  ours:            {o}')
        dn, on = normalise(d), normalise(o,ours_extra)
        if dn != on:
            print(f'  ERROR: differ in {sorted(dn^on)}')
            ok = False
        if o.split()[-1] != '-O3':
            print('  ERROR: ours does not end with -O3')
            ok = False
    if not ok:
        raise SystemExit('Pyodide flags in pyproject.toml are out of sync'
                         ' with the Pyodide defaults')
    print('OK: Pyodide flags in sync')

if __name__ == '__main__':
    main()
