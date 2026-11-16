
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

from contextlib import contextmanager as _ctxmgr


@_ctxmgr
def work_in_tmpdir():
    """Context manager for working in a temporary directory (automatically
    created+cleaned) and then switching back"""
    import os
    import tempfile
    the_cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            os.chdir(tmpdir)
            yield
        finally:
            os.chdir(the_cwd)#Important to leave tmpdir *before* deletion, to
                             #avoid PermissionError on Windows.

def prepend_to_path_var( env, varname, val ):
    import platform
    val = str(val)
    varsep = ':' if platform.system()!='Windows' else ';'
    if varname in env:
        env[varname] = f'{val}{varsep}{env[varname]}'
    else:
        env[varname] = val


def run_test( script, verbose, pycmd ):
    name = script.stem
    print(f" .. running test {name}")
    pypath = script.parent.parent.joinpath('pypath')
    assert pypath.is_dir()
    import os
    import subprocess
    env = os.environ.copy()
    prepend_to_path_var( env, 'PYTHONPATH', pypath )
    #Needed for windows, leaving them on all the time for now:
    env['PYTHONIOENCODING']='UTF-8'
    env['PYTHONLEGACYWINDOWSSTDIO']='UTF-8'
    #Python 3.13+'s argparse/traceback output colours can break tests:
    env['PYTHON_COLORS']='0'
    #Make sure we ignore any plugins in the environment, or we might get
    #spurious test failures:
    env['NCRYSTAL_DISABLE_DYNLOAD']='1'
    with work_in_tmpdir():
        rv = subprocess.run( [pycmd, str(script)],
                             env = env,
                             capture_output = True,
                             check = False )
        def print_output():
            snip = not verbose
            print_lines_with_snipping(rv.stdout,prefix='stdout: ', snip = snip)
            print_lines_with_snipping(rv.stderr,prefix='stderr: ', snip = snip)

        if rv.returncode != 0:
            print(" .. failed !")
            print_output()
            return False
        else:
            logfile = script.parent.joinpath(f'{name}.log')
            if logfile.exists():
                diff = calc_diff_output( logfile.read_bytes(),
                                         rv.stdout,
                                         'REFERENCE',
                                         'ACTUAL' )
                if diff:
                    print('Log file differs!')
                    print(diff)
                    return False
            if verbose:
                print_output()
            print(" .. success")
    return True #success

def normalise_eolmarkers( s ):
    return ( s[:-2]+'\n' ) if s.endswith('\r\n') else s

def as_lines( s ):
    #Safe decode with EOL normalisation to avoid spurious failures on windows:
    lines = ( s.decode('utf8','backslashreplace')
              if isinstance(s,bytes)
              else s ).splitlines(keepends=True)
    if not any( e.endswith('\r\n') for e in lines ):
        return lines
    return [ normalise_eolmarkers(e) for e in lines ]

def print_lines_with_snipping( b, prefix, snip = True ):
    if not b:
        return
    lines = as_lines(b)
    n1,n2 = (20,80) if snip else (999999999, 999999999)
    if len(lines)>n1+n2+20:
        lines = ( lines[0:n1]
                  + [f'<SNIPPED {len(lines)-n1-n2} lines>']
                  +lines[-n2:] )
    print(prefix + prefix.join(lines))

def calc_diff_output( a, b, beforetxt, aftertxt ):
    if a == b:
        return ''
    import difflib
    ud = difflib.unified_diff( as_lines(a), as_lines(b),
                               fromfile=beforetxt, tofile=aftertxt )
    return ''.join(ud)

def extract_needs_statement( pyfile ):
    with pyfile.open('rt') as fh:
        for line in fh:
            if line.startswith('# NEEDS: '):
                return set(line[len('# NEEDS: '):].split())
    return None

def try_cfg_with_cmake_project(verbose):
    #Return true if able to use CMake with C and CXX compilers.
    import os
    import pathlib
    import shlex
    import shutil
    import subprocess
    cmake_cmd = shutil.which('cmake')
    if not cmake_cmd:
        return False
    with work_in_tmpdir():
        td = pathlib.Path('.')
        src=td.joinpath('src')
        src.mkdir()
        bld=td.joinpath('bld')
        bld.mkdir()
        src.joinpath('main.cpp').write_text('int main() { return 0; }\n')
        src.joinpath('CMakeLists.txt').write_text(
            'cmake_minimum_required(VERSION 3.16...3.31)\n'
            'project( NCrystalTestCompilers LANGUAGES C CXX )\n'
            'add_executable( myapp "main.cpp")\n'
        )
        print(end='',flush=True)
        rv = subprocess.run(([cmake_cmd,'-S','src','-B','bld']
                             +shlex.split(os.environ.get('CMAKE_ARGS',''))),
                            capture_output=not verbose,check=False)
        print(end='',flush=True)
    return rv.returncode == 0

#Code run with the Python command of the tests (via "-c"), to check which of
#the dependencies are present there (with the same version checks as used by
#NCrystal). Prints a JSON dict: dep -> [present, note], plus sys.platform:
_depcheck_code = r"""
import importlib, importlib.metadata, json, re, sys
def vtuple( s ):
    return tuple( int(e) for e in re.findall(r'\d+',s)[:3] )
res = {}
for dep, ( minv, modname, pypiname ) in json.loads( sys.argv[1] ).items():
    if dep == 'subprocess':
        res[dep] = [ sys.platform != 'emscripten', '' ]
        continue
    try:
        mod = importlib.import_module( modname )
    except ModuleNotFoundError:
        res[dep] = [ False, '' ]
        continue
    v = getattr( importlib.import_module( modname.split('.')[0] ),
                 '__version__', None )
    if v is None:
        try:
            v = importlib.metadata.version( pypiname )
        except importlib.metadata.PackageNotFoundError:
            pass
    if minv and v and vtuple(v) < vtuple(minv):
        res[dep] = [ False, f' (version {v} is too old, {minv} or later needed)' ]
    else:
        res[dep] = [ True, '' ]
print( json.dumps( { 'deps' : res, 'platform' : sys.platform } ) )
"""

def check_deps_with_pycmd( deps, pycmd ):
    import json
    import subprocess

    from ._depdb import deps as depdb
    #"subprocess" is not a package, just the ability to start processes:
    table = { d : ( [None,None,None] if d == 'subprocess' else depdb[d] )
              for d in deps }
    rv = subprocess.run( [ pycmd, '-c', _depcheck_code, json.dumps(table) ],
                         capture_output = True, check = False )
    if rv.returncode != 0:
        raise SystemExit(f'Failed to check dependencies with {pycmd}:\n'
                         + rv.stderr.decode('utf8','backslashreplace'))
    return json.loads( rv.stdout.decode().strip().splitlines()[-1] )

def prepare_needs( scripts, verbose, pycmd ):
    script2dep = {}
    alldeps = set()
    for script in scripts:
        name = script.stem
        assert name not in script2dep
        needs = extract_needs_statement( script )
        if needs:
            script2dep[name] = needs
            alldeps.update( needs )
    res = check_deps_with_pycmd( alldeps - {'cmake'}, pycmd )
    absent = set()
    for dep in sorted(alldeps):
        deppretty = 'cmake (with compilers)' if dep=='cmake' else dep
        print(f"Checking for presence of {deppretty}...",end='',flush=True)
        if dep == 'cmake':
            #special, must try to tiny cmake project with C and C++ (not for
            #Emscripten, where the tests do not use the native compilers):
            ok = ( res['platform'] != 'emscripten'
                   and try_cfg_with_cmake_project(verbose) )
            note = ''
        else:
            ok, note = res['deps'][dep]
        if ok:
            print(' yes')
        else:
            absent.add( dep )
            print(' no'+note)
    return { 'absent': absent,
             'script2dep': script2dep }


def parse_args():
    import textwrap
    from argparse import ArgumentParser, RawTextHelpFormatter
    def wrap(t,w=59):
        return textwrap.fill( ' '.join(t.split()), width=w )

    descr = """Validate NCrystal installation by launching several tests."""
    parser = ArgumentParser( description=wrap(descr,79),
                             formatter_class = RawTextHelpFormatter )
    parser.add_argument('-m','--allow-missing', action='append',
                        metavar='DEP',
                        help=wrap(
                            """Skip tests which do not have the listed
                            optional dependencies installed (the special name
                            "all" means all of them). Any test missing a
                            dependency not listed with this flag will be marked
                            as a failure. Specify this flag multiple times or
                            use comma-separation to list more than one""" )
                        )

    parser.add_argument('-f','--filter', action='append',
                        metavar='PATTERN', dest='filters',
                        help=wrap(
                            """Select tests based on fnmatch'ing againt their
                            names, and only run these tests (default is to run
                            all). Specify this flag multiple times or use
                            comma-separation to list more than one pattern""" )
                        )

    parser.add_argument('-v','--verbose', action='store_true',
                        help=wrap(
                            """Always print all output of tests.""" )
                        )

    parser.add_argument('-l','--long', action='store_true',dest='longtests',
                        help=wrap(
                            """Enable longer tests as well.""" )
                        )

    parser.add_argument('--python', metavar='CMD', default=None,
                        help=wrap(
                            """Python command for running the tests (default:
                            the current Python). For example the python of a
                            Pyodide virtual environment, to test NCrystal in
                            WebAssembly.""" )
                        )

    args = parser.parse_args()

    def flatten_to_set( list_of_strs ):
        m = set()
        for a in list_of_strs or []:
            m.update({e.strip() for e in a.split(',')})
        return m

    args.allow_missing = flatten_to_set( args.allow_missing )
    args.filters = flatten_to_set( args.filters )

    return args

def main():
    args = parse_args()
    import os
    import sys
    pycmd = args.python or sys.executable
    if os.sep in pycmd or '/' in pycmd:
        #Absolute, since the tests run in temporary directories:
        pycmd = os.path.abspath( pycmd )
    filterfct = None
    if args.filters:
        import fnmatch
        def filterfct( name ):
            return any( fnmatch.fnmatch( name, f )
                        for f in args.filters )


    import pathlib
    datadir = pathlib.Path(__file__).parent.joinpath('data')
    scripts = sorted( datadir.joinpath('scripts').glob('*.py'),
                      key = lambda p : (p.name, p) )
    nhidden = 0
    if not args.longtests:
        _ = []
        for s in scripts:
            if not s.name.startswith('long_'):
                _.append(s)
        nhidden = len(scripts)-len(_)
        assert nhidden > 0, "no long tests defined"
        scripts = _

    depinfo = prepare_needs( scripts, args.verbose, pycmd )
    jobs = []

    for script in scripts:
        name = script.stem
        needs = depinfo['script2dep'].get(name,set())
        missing = needs.intersection( depinfo['absent'] )
        jobs.append( (name, script, missing ) )

    failures = 0
    successes = 0
    dep_skipped_ok = 0
    dep_skipped_bad = 0
    hidden = 0
    ntot = 0
    for name, script, missing in jobs:
        ntot += 1
        if filterfct and not filterfct(name):
            hidden += 1
            continue
        if missing:
            ok = ( 'all' in args.allow_missing
                   or not (missing-args.allow_missing) )
            if ok:
                dep_skipped_ok += 1
                print(" .. skipping test {} (missing: {})".format(name,' '.join(missing)))
            else:
                dep_skipped_bad += 1
                print(" .. can not run test {} (missing: {})".format(name,' '.join(missing)))
        else:
            if not run_test(script,args.verbose,pycmd):
                failures += 1
            else:
                successes += 1
    assert len(jobs) == ntot

    print()
    print('Test summary:')
    print()
    stats = [
        ('successful tests', successes ),
        ('failed tests', failures ),
        ('tests w/o deps (ok)', dep_skipped_ok),
        ('tests w/o deps (not ok)', dep_skipped_bad),
        ('tests hidden by filters', hidden ),
    ]
    stats.append( ( 'TOTAL', ntot) )

    w1 = max( len(s) for s,n in stats )
    w2 = max( len(str(n)) for s,n in stats )
    for s,n in stats:
        print(f'   {s.rjust(w1)}:  {str(n).rjust(w2)}')
    print()

    if nhidden:
        print("Note: Additional tests can be enabled with the --long flag.\n")

    assert ntot == ( successes
                     + failures
                     + dep_skipped_ok
                     + dep_skipped_bad
                     + hidden )
    assert sum( n for s,n in stats ) == 2*ntot

    nbad = failures + dep_skipped_bad
    if nbad:
        print(f'ERROR: Problems with {nbad} tests.')
        raise SystemExit(1)
    if not successes:
        print('ERROR: Requires at least a '
              'single test to run succesfully.')
        raise SystemExit(1)
    print('All selected tests ran successfully.')
