
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

#Local version of the plugindb.yml workflow: build the NCrystal wheels once,
#then install and test each plugin in its own fresh venv.

def short_description():
    return ('Install and test plugins from the plugin database against'
            ' wheels built from this repo')

def main( parser ):
    parser.init( short_description() )
    parser.add_argument(
        'plugins', nargs='*', metavar='PLUGIN',
        help="""Plugins to test (default: all in the plugin database)."""
    )
    parser.add_argument(
        '--workdir', metavar='DIR',
        help="""Directory for all the work (wheels, plugin sources, venvs and
        logs). Default: a new temporary directory. Wheels already present in
        DIR/wheels are reused, so a rerun for a subset of plugins is quick."""
    )
    parser.add_argument(
        '--rebuild', action='store_true',
        help="""Rebuild the NCrystal wheels even if DIR/wheels has them."""
    )
    parser.add_argument(
        '--python', metavar='PY', default=None,
        help="""Python interpreter for the venvs (default: the one running
        this script)."""
    )
    parser.add_argument(
        '--no-pypi', action='store_true',
        help="""Skip installing and testing the PyPI versions of plugins."""
    )
    args = parser.parse_args()

    import pathlib
    import sys
    import tempfile

    from ._climode_plugindb import load_and_check_data
    data = load_and_check_data()
    plugins = args.plugins or sorted(data.keys())
    for p in plugins:
        if p not in data:
            parser.error(f'Unknown plugin "{p}" (choices:'
                         f' {", ".join(sorted(data.keys()))})')
    workdir = pathlib.Path( args.workdir or
                            tempfile.mkdtemp(prefix='ncplugindbtest_') )
    workdir = workdir.absolute()
    workdir.mkdir( parents = True, exist_ok = True )
    python = args.python or sys.executable
    print(f'Work directory: {workdir}')

    wheeldir = _build_wheels( workdir, python, args.rebuild )
    results = []
    for p in plugins:
        results.append( _test_plugin( p, data[p], workdir, wheeldir, python,
                                      skip_pypi = args.no_pypi ) )
    _print_summary( results, workdir )
    if any( not r['ok'] for r in results ):
        raise SystemExit(1)

def _run( cmd, logfh, cwd = None, env = None ):
    import subprocess
    cmd = [str(c) for c in cmd]
    logfh.write('\n>>> ' + ' '.join(cmd) + '\n')
    logfh.flush()
    rv = subprocess.run( cmd, cwd = cwd, env = env, stdout = logfh,
                         stderr = subprocess.STDOUT, check = False )
    return rv.returncode == 0

def _venv_env( venvdir, extra = None ):
    #Environment as if venvdir was activated, without anything from the
    #environment of this script which could affect NCrystal or Python:
    import os
    bindir = _venv_bin( venvdir, 'python' ).parent
    old_venv = os.environ.get('VIRTUAL_ENV')
    path = [ p for p in os.environ.get('PATH','').split(os.pathsep)
             if p and not ( old_venv and p.startswith(old_venv) ) ]
    env = { k : v for k, v in os.environ.items()
            if not ( k.startswith( ('NCRYSTAL','PYTHON') )
                     or k in ('VIRTUAL_ENV','VIRTUAL_ENV_PROMPT') ) }
    env['PATH'] = os.pathsep.join( [ str(bindir) ] + path )
    env['VIRTUAL_ENV'] = str(venvdir)
    env.update( extra or {} )
    return env

def _venv_bin( venvdir, name ):
    import os
    return ( venvdir / 'Scripts' / f'{name}.exe' if os.name == 'nt'
             else venvdir / 'bin' / name )

def _new_venv( path, python, logfh ):
    import shutil
    if path.exists():
        shutil.rmtree( path )
    if not _run( [ python, '-m', 'venv', path ], logfh ):
        return None
    py = _venv_bin( path, 'python' )
    if not _run( [ py, '-m', 'pip', 'install', '--upgrade', 'pip' ], logfh,
                 env = _venv_env( path ) ):
        return None
    return py

def _version():
    from .dirs import reporoot
    return reporoot.joinpath('VERSION').read_text().strip()

def _build_wheels( workdir, python, rebuild ):
    from .dirs import reporoot
    wheeldir = workdir / 'wheels'
    pkgs = [ 'ncrystal_core', 'ncrystal_python', 'ncrystal_metapkg',
             'ncrystal_pypluginmgr' ]
    v = _version()
    if ( not rebuild and wheeldir.is_dir()
         and sum( 1 for f in wheeldir.glob('*.whl') if f'-{v}-' in f.name )
         >= len(pkgs) ):
        print(f'Reusing NCrystal {v} wheels in {wheeldir}')
        return wheeldir
    import shutil
    if wheeldir.exists():
        shutil.rmtree( wheeldir )
    wheeldir.mkdir()
    logfile = workdir / 'build_wheels.log'
    print(f'Building NCrystal {v} wheels (log: {logfile})')
    with logfile.open('wt') as fh:
        py = _new_venv( workdir / 'venv_build', python, fh )
        env = _venv_env( workdir / 'venv_build' )
        ok = py is not None and _run( [ py, '-m', 'pip', 'install', 'build' ],
                                      fh, env = env )
        for pkg in pkgs:
            ok = ok and _run( [ py, '-m', 'build', '--wheel', '--outdir',
                                wheeldir, reporoot / pkg ], fh, env = env )
    if not ok:
        raise SystemExit(f'ERROR: Failed to build wheels (see {logfile})')
    return wheeldir

def _get_plugin_src( name, info, workdir, logfh ):
    import shutil

    from .dirs import reporoot
    srcdir = workdir / 'src' / name
    if srcdir.exists():
        shutil.rmtree( srcdir )
    srcdir.parent.mkdir( parents = True, exist_ok = True )
    subdir = info['repo_subdir']
    if info['github_repo_key'] == 'mctools/ncrystal':
        shutil.copytree( reporoot / subdir, srcdir )
        return srcdir
    clonedir = workdir / 'src' / f'{name}_reporoot'
    if clonedir.exists():
        shutil.rmtree( clonedir )
    cmd = [ 'git', 'clone', '--quiet', info['repourl'], clonedir ]
    if not _run( cmd, logfh ):
        return None
    if info['gitref'] and not _run( [ 'git', 'checkout', '--quiet',
                                      info['gitref'] ], logfh,
                                    cwd = clonedir ):
        return None
    if not subdir:
        clonedir.rename( srcdir )
        return srcdir
    shutil.copytree( clonedir / subdir, srcdir )
    return srcdir

def _test_plugin( name, info, workdir, wheeldir, python, skip_pypi ):
    logdir = workdir / 'logs'
    logdir.mkdir( exist_ok = True )
    logfile = logdir / f'{name}.log'
    res = { 'name' : name, 'ok' : False, 'failed_step' : None,
            'log' : logfile, 'notes' : [] }
    print(f'Testing plugin {name} (log: {logfile}) ...', flush = True)
    constraints = workdir / 'constraints.txt'
    v = _version()
    constraints.write_text( ''.join( f'{p}=={v}\n' for p in
                                     ( 'ncrystal', 'ncrystal-core',
                                       'ncrystal-python',
                                       'ncrystal-verify' ) ) )
    pkgname = f'ncrystal-plugin-{name}'
    do_tests = not info['disable_tests']
    if not do_tests:
        res['notes'].append('tests disabled in plugin database')
    with logfile.open('wt') as fh:
        def step( label, ok ):
            if not ok and res['failed_step'] is None:
                res['failed_step'] = label
            return ok
        srcdir = _get_plugin_src( name, info, workdir, fh )
        if not step( 'get plugin source', srcdir is not None ):
            return _done( res )
        venvdir = workdir / 'venvs' / name
        py = _new_venv( venvdir, python, fh )
        if not step( 'create venv', py is not None ):
            return _done( res )
        env = _venv_env( venvdir )
        pip = [ py, '-m', 'pip' ]
        c = [ '-c', constraints ]
        steps = [
            ( 'install build backends',
              pip + [ 'install' ] + c + [ 'numpy', 'setuptools>=77.0',
                                          'scikit-build-core>=0.11.1' ] ),
            ( 'install ncrystal',
              pip + [ 'install' ] + c + [ '--no-index', '--find-links',
                                          wheeldir, 'ncrystal',
                                          'ncrystal-pypluginmgr' ] ),
            ( 'nctool --test', [ _venv_bin( venvdir, 'nctool' ), '--test' ] ),
            ( 'install plugin',
              pip + [ 'install' ] + c + [ srcdir, '--no-build-isolation' ] ),
            ( 'pip check', pip + [ 'check' ] ),
            ( 'plugin package name', pip + [ 'show', pkgname ] ),
        ]
        pm = _venv_bin( venvdir, 'ncrystal-pluginmanager' )
        if do_tests:
            steps.append( ( 'test plugin', [ pm, '--test', name ] ) )
        steps.append( ( 'uninstall plugin',
                        pip + [ 'uninstall', '-y', pkgname ] ) )
        for label, cmd in steps:
            if not step( label, _run( cmd, fh, env = env ) ):
                return _done( res )
        #Users typically install with build isolation (PIP_CONSTRAINT also
        #applies to the isolated build environment, unlike -c):
        env_iso = _venv_env( venvdir, { 'PIP_CONSTRAINT' : str(constraints) } )
        if not step( 'install plugin with build isolation',
                     _run( pip + [ 'install', '--find-links', wheeldir,
                                   srcdir ], fh, env = env_iso ) ):
            return _done( res )
        if do_tests and not step( 'test plugin (build isolation)',
                                  _run( [ pm, '--test', name ], fh,
                                        env = env ) ):
            return _done( res )
        if not step( 'uninstall plugin (build isolation)',
                     _run( pip + [ 'uninstall', '-y', pkgname ], fh,
                           env = env ) ):
            return _done( res )
        if info['is_on_pypi'] and not skip_pypi:
            if not step( 'install plugin from PyPI',
                         _run( pip + [ 'install' ] + c
                               + [ pkgname, '--no-build-isolation' ], fh,
                               env = env ) ):
                return _done( res )
            if do_tests and not step( 'test PyPI plugin',
                                      _run( [ pm, '--test', name ], fh,
                                            env = env ) ):
                return _done( res )
    res['ok'] = True
    return _done( res )

def _done( res ):
    print( '  -> ' + ( 'OK' if res['ok'] else
                       f'FAILED in step "{res["failed_step"]}"' ),
           flush = True )
    return res

def _print_summary( results, workdir ):
    print()
    print('Summary:')
    w = max( len(r['name']) for r in results )
    for r in results:
        status = 'OK' if r['ok'] else f'FAILED ({r["failed_step"]})'
        notes = f' [{"; ".join(r["notes"])}]' if r['notes'] else ''
        print(f'  {r["name"].ljust(w)} : {status}{notes}')
    nfail = sum( 1 for r in results if not r['ok'] )
    print()
    print(f'Tested {len(results)} plugin(s), {nfail} failed.')
    for r in results:
        if not r['ok']:
            print(f'  Log for {r["name"]}: {r["log"]}')
    print(f'Work directory: {workdir}')
