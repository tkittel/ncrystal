
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


"""Internal module for importing optional third-party modules, with a check
of their version against the minimum in _depdb, and consistent errors telling
users what to install."""

_cache = {}

def _version_tuple( s ):
    import re
    return tuple( int(e) for e in re.findall(r'\d+',s or '')[:3] )

def find_optional( name ):
    """Returns (module,None) with the top-level module of the optional
    dependency "name" (a key in _depdb.deps), or (None,problem) if it is
    missing, broken, or older than the minimum version. Results are cached."""
    if name in _cache:
        return _cache[name]
    import importlib

    from ._depdb import deps
    minv, modname, pypiname, _ = deps[name]
    problem, mod = None, None
    try:
        importlib.import_module( modname )
        mod = importlib.import_module( modname.split('.')[0] )
    except ImportError as e:
        import importlib.util
        if importlib.util.find_spec( modname.split('.')[0] ) is None:
            problem = 'is not installed'
        else:
            problem = f'could not be imported ({e})'
    if mod is not None and minv:
        v = getattr( mod, '__version__', None )
        if v is None:
            import importlib.metadata
            try:
                v = importlib.metadata.version( pypiname )
            except importlib.metadata.PackageNotFoundError:
                v = None
        if v and _version_tuple(v) < _version_tuple(minv):
            problem = f'is too old (found version {v})'
    _cache[name] = ( None, problem ) if problem else ( mod, None )
    return _cache[name]

def optional_error_msg( name, problem, purpose = None ):
    """Error message for a problem reported by find_optional, saying what to
    install."""
    from ._depdb import deps
    minv, _, pypiname, condaname = deps[name]
    vreq = f'>={minv}' if minv else ''
    return ( f'The Python package "{pypiname}"'
             + ( f', needed {purpose},' if purpose else '' )
             + f' {problem}. NCrystal requires {pypiname}{vreq}, which is'
             ' available on PyPI ("python3 -m pip install'
             f' \'{pypiname}{vreq}\'") and conda-forge ("conda install -c'
             f' conda-forge \'{condaname}{vreq}\'").' )

def import_optional( name, *, purpose = None, sysexit = False ):
    """Imports and returns the top-level module of the optional dependency
    "name" (a key in _depdb.deps). If it is missing or older than the minimum
    version, raises ImportError (SystemExit if sysexit=True) with installation
    instructions. The optional purpose (e.g. "to process CIF files") is
    included in the error message."""
    mod, problem = find_optional( name )
    if problem:
        msg = optional_error_msg( name, problem, purpose )
        if sysexit:
            raise SystemExit( f'ERROR: {msg}' )
        raise ImportError( msg )
    return mod
