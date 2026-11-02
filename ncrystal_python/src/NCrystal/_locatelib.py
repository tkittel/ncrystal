
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

"""Internal module for locating the compiled NCrystal shared library and
determining its namespace"""

import pathlib

_cache = [None]
def get_libpath_and_namespace():
    if _cache[0] is None:
        _cache[0] = _search()
    return _cache[0]

def _detect_monolithic_installation():
    #Determine if these python modules were installed via "pip install
    #<ncrystalrepo>" (monolithic package with bore core and python parts), or if
    #using a standard ncrystal-python package. Also detect conflicts in case a
    #user did both, and thus messed up their environment.

    d = pathlib.Path(__file__).parent
    is_std = d.joinpath('_is_std.py').is_file()
    is_mono = d.joinpath('_is_monolithic.py').is_file()
    assert is_mono or is_std
    if not ( is_mono and is_std ):
        return is_mono

    import textwrap
    conflict_msg = textwrap.dedent("""

    ERROR: Broken environment detected.

    The current environment has traces of both a normal NCrystal installation
    (where ncrystal-core and ncrystal-python modules are installed separately),
    and a monolithic installation (most likely from running "pip install"
    directly on the root of the ncrystal source tree). This latter method is
    primarily intended for developers, and should never be performed in an
    environment where NCrystal is already installed through other conventional
    methods.

    Most likely your environment is now broken. You can attempt to fix it by
    uninstalling all NCrystal packages (with the same package manager that did
    the original installation of that package), and then installing NCrystal
    again from one source only. Or you might simply wish to remove and recreate
    the environment from scratch.
    """)
    #We can not simply emit ImportError, since it might be taken by clients as
    #"ncrystal not installed". We can choose RuntimeError or SystemExit. To be
    #100% this triggers a warning, we pick SystemExit.
    raise SystemExit(conflict_msg)

def _search():
    #Try to determine the location of the NCrystal shared library. In order of
    #preference, we try:
    #
    #  1) The NCrystal_LIB environment variable (intended for specialised usage,
    #     like CTests). Normally this variable is also enough for us to be able
    #     to decode the NCrystal namespace protection used, if any. But if it
    #     fails, one can use the NCRYSTAL_LIB_NAMESPACE_PROTECTION variable to
    #     specify it explicitly. If set, it is always used.
    #
    #  2) Attempt to get the location via a python module providing it (in case
    #     ncrystal-core was installed via PyPI for instance).
    #
    #  3) Invoke "ncrystal-config --show ..." for the information, first the one
    #     in the environment of the current Python interpreter (which might not
    #     be activated, i.e. not in the PATH), then the one in the PATH.
    #
    #  4) Last resort: Let the system dynamic loader find the library by its
    #     SONAME (libNCrystal.so.<major>). This is for installations where the
    #     library is in the default library search path but ncrystal-config is
    #     not available (e.g. a system package with ncrystal-config in a
    #     separate development package). Only non-namespaced libraries can be
    #     found like this.
    #
    #  Candidates from 2-4 with a version different from that of the Python
    #  modules are skipped. If simplebuild devel mode we allow env var
    #  overrides, but apart from that always go straight to the sb_nccmd_config
    #  command in the PATH.

    # Always invoke _detect_monolithic_installation() since it also detects
    # broken installations.
    import os
    verbose = 'NCRYSTAL_DEBUG_LIBSEARCH' in os.environ
    if verbose:
        from ._common import print
        print('NCrystal._locatelib: Starting search for'
              ' NCrystal shared library')
    is_monolithic = _detect_monolithic_installation()
    if verbose:
        print('NCrystal._locatelib: monolithic installation'
              f' = {"yes" if is_monolithic else "no"}')

    v = _search_env_overrides()
    if v:
        if verbose:
            print('NCrystal._locatelib: Succesfully searched via method:'
                  ' env vars')
        return _finalise( v[0], v[1], verbose )

    is_simplebuild_devel = ( pathlib.Path(__file__).parent
                             .joinpath('_is_sblddevel.py').is_file() )

    from . import __version__ as _nc_version
    rejected = []#candidates with wrong version
    for method, fct, by_name in _candidate_searches( is_simplebuild_devel,
                                                     is_monolithic ):
        v = fct()
        if not v:
            continue
        lib, namespace, version = v
        if verbose:
            print(f'NCrystal._locatelib: Found via method {method}:'
                  f' "{lib}" (version {version})')
        if version is not None and version.strip() != _nc_version.strip():
            rejected.append( ( method, lib, version.strip() ) )
            continue
        if verbose:
            print('NCrystal._locatelib: Succesfully searched via method:'
                  f' {method}')
        return _finalise( lib, namespace, verbose, by_name )

    if not rejected:
        from .exceptions import NCFileNotFound
        raise NCFileNotFound('Could not locate the NCrystal shared library'
                             ' file. Have you installed the ncrystal-core'
                             ' package?')

    #Found something, but only with version mismatches!
    import textwrap
    found = '\n'.join( f'  {ver} ({meth}: {lib})'
                       for meth, lib, ver in rejected )
    msg = textwrap.dedent(f"""
    ERROR: Inconsistent environment detected.

    The ncrystal Python API has version {_nc_version}, but only core ncrystal
    libraries with other versions were found. This is not supported. Found:
    """) + found + '\n'
    #The environment is not necessarily completely broken, but the NCrystal
    #packages can not be used. Hence, we emit an ImportError here. However, we
    #also emit a warning, in case the client code is hiding ImportError's, on
    #the assumption it simply means NCrystal is not available.
    from ._common import warn
    warn(msg)
    raise ImportError(msg)

def _finalise( lib, namespace, verbose, by_name = False ):
    #by_name: lib is a name for the system dynamic loader, not a path.
    if verbose:
        from ._common import print
        print(f'NCrystal._locatelib: namespace = "{namespace or ""}"')
        print(f'NCrystal._locatelib: lib = "{lib}"')
    if by_name:
        return lib, ( namespace or '' )
    lib = pathlib.Path(lib)
    if not lib.is_file():
        from .exceptions import NCFileNotFound
        raise NCFileNotFound('Problems locating the NCrystal shared library'
                             f' file. Tried "{lib}" unsuccesfully.')
    return ( lib.absolute().resolve(), ( namespace or '' ) )

def _candidate_searches( is_simplebuild_devel, is_monolithic ):
    #Yields (methodname, fct, by_name) in order of preference, where fct()
    #returns (lib,namespace,version) or None, and by_name indicates that lib is
    #a name for the system dynamic loader rather than a path.
    if is_simplebuild_devel:
        yield ( 'sb_nccmd_config',
                lambda : _search_nccfgapp('sb_nccmd_config'), False )
        return
    yield 'pymod', lambda : _search_core_info_mod( is_monolithic ), False
    for cmd in _nccfgapp_cmds():
        yield ( f'ncrystal-config ({cmd})',
                lambda c=cmd: _search_nccfgapp(c), False )
    yield 'system loader (SONAME)', _search_system_soname, True

def _nccfgapp_cmds():
    #ncrystal-config commands to try: first the one in the environment of the
    #current Python interpreter (which might not be activated), then the one in
    #the PATH (if different).
    import os
    import shutil
    import sys
    import sysconfig
    dirs = [ sysconfig.get_path('scripts') ]
    if sys.platform.startswith('win'):
        #conda on Windows:
        dirs.append( str(pathlib.Path(sys.prefix).joinpath('Library','bin')) )
    cmds = [ shutil.which( 'ncrystal-config', path = d ) for d in dirs if d ]
    cmds.append( shutil.which( 'ncrystal-config' ) )
    seen = set()
    for c in cmds:
        key = os.path.normcase( os.path.abspath( c ) ) if c else None
        if key and key not in seen:
            seen.add( key )
            yield c

def _search_env_overrides():
    import os
    #NB: The next two variables can NOT be namespaced. E.g. it will always be
    #NCRYSTAL_LIB, and never e.g. NCRYSTAL<NAMESPACEHERE>_LIB:
    lib = os.environ.get('NCRYSTAL_LIB')

    #This second variable is hopefully not needed, since the namespace can
    #hopefully be inferred from the shared library name, but we allow it as an
    #ultimate fall-back option in case one is for some reason running with
    #non-standard naming of the shared libraries:
    ns = os.environ.get('NCRYSTAL_LIB_NAMESPACE_PROTECTION','')
    if lib:
        if not ns:
            #Try to infer the namespace from the library file name (so it is
            #enough to set NCRYSTAL_LIB).
            ll = pathlib.Path(lib).name.split('.')[0]
            if 'NCrystal-' in ll:
                ll = ll.split('NCrystal-')[-1]
                if ll and 'NCrystal-' not in ll and '.' not in ll:
                    ns = ll
        lib = pathlib.Path(lib)
        if not lib.exists() or lib.is_dir():
            from .exceptions import NCFileNotFound
            raise NCFileNotFound('NCRYSTAL_LIB environment variable is set'
                                 f' ("{lib}") but does not point'
                                 ' to an actual file.')
        return lib, ns, None

def _search_core_info_mod( is_monolithic ):
    #Look for the standard module, installed by the ncrystal-core package.
    try:
        if is_monolithic:
            import _ncrystal_core_monolithic.info as mod
        else:
            import _ncrystal_core.info as mod
    except ModuleNotFoundError:
        mod = None
    if mod:
        return mod.libpath(), mod.namespace(), mod.version()

def _search_system_soname():
    #Try to load libNCrystal.so.<major> (libNCrystal.<major>.dylib on macOS)
    #via the system dynamic loader, to read its version. Not attempted on
    #Windows. A library with another version stays loaded (C++ libraries can in
    #general not be unloaded), which is harmless: this is the last candidate, so
    #the import then fails, and with RTLD_LOCAL its symbols can not interfere
    #with other libraries:
    import platform
    sysname = platform.system()
    if sysname == 'Windows':
        return None
    from . import __version__ as _nc_version
    major = _nc_version.split('.')[0]
    name = ( f'libNCrystal.{major}.dylib' if sysname == 'Darwin'
             else f'libNCrystal.so.{major}' )
    import ctypes
    try:
        lib = ctypes.CDLL( name, mode = ctypes.RTLD_LOCAL )
        f = lib.ncrystal_version_str
    except (OSError, AttributeError):
        return None
    f.restype = ctypes.c_char_p
    f.argtypes = []
    return name, '', f().decode('utf8')

def _search_nccfgapp( cmdname ):
    #Try to query ncrystal-config script:
    import subprocess
    try:
        res = subprocess.run( [cmdname,'--show','shlibpath',
                               'namespace','version'],
                              capture_output=True, check=False )
    except OSError:
        #Not found, not executable, or e.g. built for another architecture:
        return None

    if res.returncode == 0:
        lines = res.stdout.decode('utf8').splitlines()
        if len(lines) == 3:
            p = pathlib.Path(lines[0])
            if p.is_file():
                return p, ( lines[1].strip() or None ), lines[2].strip()
