
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

# Windows-only CI test of how the ncrystal-core Python module locates
# NCrystal.dll (_ncrystal_core/info.py: RECORD lookup, directory heuristic,
# sysconfig scheme search). Each scenario installs the ncrystal-core wheel in
# a different way into a temporary directory and checks which DLL is found.
#
# Usage: python windll_lookup.py <dir-with-ncrystal-core-wheel> <srcdir>

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

wheeldir = pathlib.Path(sys.argv[1]).absolute()
srcdir = pathlib.Path(sys.argv[2]).absolute()
wheels = list(wheeldir.glob('ncrystal_core-*.whl'))
assert len(wheels) == 1, wheels
wheel = wheels[0]
tmp = pathlib.Path(tempfile.mkdtemp(prefix='ncwindll_'))
failures = []

_query = ( 'import _ncrystal_core.info as m;'
           'print(m.libpath())' )

def norm( p ):
    return os.path.normcase(os.path.realpath(str(p)))

def run( cmd, env_extra = None, check = True ):
    env = dict(os.environ)
    env.pop('PYTHONPATH',None)
    env.pop('PYTHONUSERBASE',None)
    env.update( env_extra or {} )
    print('  $', ' '.join(str(c) for c in cmd), flush = True)
    r = subprocess.run( [str(c) for c in cmd], env = env, check = False,
                        capture_output = True, text = True )
    if check and r.returncode != 0:
        print(r.stdout)
        print(r.stderr)
        raise RuntimeError(f'Command failed: {cmd}')
    return r

def mkvenv( name ):
    d = tmp / name
    run( [ sys.executable, '-m', 'venv', d ] )
    return d, d / 'Scripts' / 'python.exe'

def pipinstall( py, *args, env_extra = None ):
    run( [ py, '-m', 'pip', 'install', '--no-deps', '--no-index',
           *args, wheel ], env_extra = env_extra )

def check( label, py, expected, env_extra = None ):
    r = run( [ py, '-c', _query ], env_extra = env_extra, check = False )
    got = r.stdout.strip()
    ok = ( r.returncode == 0 and norm(got) == norm(expected) )
    print(f'{"OK  " if ok else "FAIL"} {label}: got {got or r.stderr.strip()[-300:]}')
    if not ok:
        print(f'       expected {expected}')
        failures.append(label)

def check_cfg( label, cfgexe, expected, env_extra = None ):
    r = run( [ cfgexe, '--show', 'shlibpath' ], env_extra = env_extra,
             check = False )
    got = r.stdout.strip()
    ok = ( r.returncode == 0 and norm(got) == norm(expected) )
    print(f'{"OK  " if ok else "FAIL"} {label}: got {got or r.stderr.strip()[-300:]}')
    if not ok:
        failures.append(label)

def sitepkgs_of( py, env_extra = None ):
    return pathlib.Path( run( [ py, '-c', ( 'import sysconfig;'
                                  'print(sysconfig.get_path("purelib"))' ) ],
                              env_extra = env_extra ).stdout.strip() )

#1) Standard venv installation (RECORD lists ..\..\Scripts\NCrystal.dll):
venv, vpy = mkvenv('venv')
pipinstall( vpy )
check( 'venv install', vpy, venv / 'Scripts' / 'NCrystal.dll' )
check_cfg( 'venv install (ncrystal-config)', venv / 'Scripts' / 'ncrystal-config.exe',
           venv / 'Scripts' / 'NCrystal.dll' )

#2) User installation (scripts in USERBASE\PythonXY\Scripts), with a decoy
#   DLL in USERBASE\Scripts which the directory heuristic would pick first:
ub = tmp / 'userbase'
ubenv = { 'PYTHONUSERBASE' : str(ub) }
pipinstall( sys.executable, '--user', env_extra = ubenv )
ub_scripts = pathlib.Path( run( [ sys.executable, '-c', ( 'import sysconfig;'
                                  'print(sysconfig.get_path("scripts","nt_user"))' ) ],
                                env_extra = ubenv ).stdout.strip() )
decoy = ub / 'Scripts' / 'NCrystal.dll'
decoy.parent.mkdir( parents = True, exist_ok = True )
decoy.write_bytes(b'decoy')
check( 'user install, decoy in USERBASE\\Scripts', sys.executable,
       ub_scripts / 'NCrystal.dll', env_extra = ubenv )

#3) --target installation, imported from a venv which also has its own
#   installation (whose DLL must not be picked up). pip installs into a
#   temporary home scheme and moves the result into DIR, so the DLL ends up in
#   DIR\Scripts while RECORD still says ../../Scripts/NCrystal.dll, which we
#   make point to a decoy:
tgt = tmp / 't1' / 't2' / 'target'
pipinstall( vpy, '--target', tgt )
stale = tgt / '..' / '..' / 'Scripts' / 'NCrystal.dll'
stale.parent.mkdir( parents = True, exist_ok = True )
stale.write_bytes(b'decoy')
check( '--target install, decoy at stale RECORD path', vpy,
       tgt / 'Scripts' / 'NCrystal.dll',
       env_extra = { 'PYTHONPATH' : str(tgt) } )

#4) Without RECORD, the directory heuristic must find the DLL:
venv2, vpy2 = mkvenv('venv_norecord')
pipinstall( vpy2 )
for rec in sitepkgs_of( vpy2 ).glob('ncrystal_core-*.dist-info/RECORD'):
    rec.unlink()
check( 'venv install without RECORD', vpy2, venv2 / 'Scripts' / 'NCrystal.dll' )

#5) Module copied elsewhere (no RECORD, heuristic misses): the sysconfig
#   scheme search must find the DLL of the running venv:
elsewhere = tmp / 'a' / 'b' / 'elsewhere'
shutil.copytree( sitepkgs_of( vpy ) / '_ncrystal_core',
                 elsewhere / '_ncrystal_core' )
check( 'copied module, scheme fallback', vpy, venv / 'Scripts' / 'NCrystal.dll',
       env_extra = { 'PYTHONPATH' : str(elsewhere) } )

#6) Copied module in an environment without any DLL: must fail clearly.
venv3, vpy3 = mkvenv('venv_empty')
r = run( [ vpy3, '-c', _query ], env_extra = { 'PYTHONPATH' : str(elsewhere) },
         check = False )
ok = r.returncode != 0 and 'Unable to locate shared library' in r.stderr
print(f'{"OK  " if ok else "FAIL"} no DLL anywhere gives a clear error')
if not ok:
    print(r.stdout, r.stderr)
    failures.append('no DLL anywhere')

#7) Build isolation: building the example plugin in a venv without
#   ncrystal-core installed, so ncrystal-config only exists in pip's temporary
#   build environment (this failed before the fix):
venv4, vpy4 = mkvenv('venv_buildiso')
#Pin ncrystal-core to the tested wheel (PIP_CONSTRAINT also applies inside the
#build environment, unlike -c):
constraints = tmp / 'constraints.txt'
constraints.write_text(f'ncrystal-core=={wheel.name.split("-")[1]}\n')
r = run( [ vpy4, '-m', 'pip', 'wheel', '--no-deps', '--wheel-dir', tmp / 'pw',
           '--find-links', wheeldir, srcdir / 'examples' / 'plugin' ],
         env_extra = { 'PIP_CONSTRAINT' : str(constraints) }, check = False )
ok = r.returncode == 0
print(f'{"OK  " if ok else "FAIL"} plugin build with build isolation')
if not ok:
    print(r.stdout[-3000:], r.stderr[-3000:])
    failures.append('build isolation')

if failures:
    for f in failures:
        #GitHub annotation, visible without access to the full logs:
        print(f'::error title=windll_lookup::Scenario failed: {f}')
    raise SystemExit('FAILURES: ' + ', '.join(failures))
print('ALL OK')
