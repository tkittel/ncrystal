
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

"""Module with internal utilities used by several NCrystal modules"""

#The print function used for all output of NCrystal (allowing one to capture
#and redirect it, cf. also print(..) and NCrystal.msg below):

def _builtin_print():
    import builtins
    return builtins.print

_stdprint = [_builtin_print()]

def set_ncrystal_print_fct( fct ):
    oldfct = _stdprint[0]
    _stdprint[0] = fct or _builtin_print()
    return oldfct

def get_ncrystal_print_fct():
    return _stdprint[0]

_cache_ns = [None]
def get_namespace():
    if _cache_ns[0] is None:
        from ._locatelib import get_libpath_and_namespace as f
        _p,_ns = f()
        _cache_ns[0] = _ns
    return _cache_ns[0]

_cache_nsev = [None]
def has_namespaced_envvars():
    if _cache_nsev[0] is None:
        import pathlib
        _cache_nsev[0] = ( pathlib.Path(__file__).parent
                           .joinpath('_do_namespace_envvars.py').is_file() )
    return _cache_nsev[0]

def expand_envname( name ):
    assert not name.startswith('NCRYSTAL_')#common mistake
    if has_namespaced_envvars():
        return 'NCRYSTAL' + get_namespace().upper()+'_' + name
    else:
        return 'NCRYSTAL_' + name

def ncgetenv( name, defval = None ):
    import os
    return os.environ.get( expand_envname(name), defval )

def ncgetenv_bool( name ):
    v = ncgetenv( name )
    return not (v is None or v=='0')

def ncgetenv_int( name, defval ):
    assert isinstance( defval, int)
    v = ncgetenv( name )
    if v is None:
        return defval
    try:
        v=int(v)
    except ValueError:
        n = expand_envname( name )
        raise ValueError(f"ERROR: {n} environment variable must"
                         " be integral value if set")
    return v

def ncgetenv_int_nonneg( name, defval ):
    assert isinstance( defval, int) and defval >= 0
    v = ncgetenv_int( name, defval )
    if v < 0:
        n = expand_envname( name )
        raise ValueError(f"ERROR: {n} environment variable must"
                         " be non-negative value if set")
    return v

def ncsetenv( name, val ):
    import os
    varname = expand_envname( name )
    if val is None:
        if varname in os.environ:
            del os.environ[varname]
    else:
        os.environ[varname] = str(val)

#Messages: notices via infomsg(..), warnings via warn(..) and raw output via
#print(..) (the counterparts of NCRYSTAL_MSG, NCRYSTAL_WARN and NCRYSTAL_RAWOUT
#in C++, whose messages also end up here via _msg.py) are first passed to all
#active NCrystal.msg.MsgCtrl contexts (which may collect them), and then
#printed via the print function (with a prefix for notices and warnings), or
#sent to a logger (cf. NCrystal.msg.use_logging and MsgCtrl).

_msg_parts = ( 'info', 'warning', 'raw' )
_msg_global = { 'logger' : None, 'include_raw' : False }
_msg_default_prefix = { 'info' : 'NCrystal: ', 'warning' : 'NCrystal WARNING: ' }
_active_msgctrls = []#active MsgCtrl objects, outermost first

def infomsg(msg):
    #Emit a notice (printed as "NCrystal: <msg>" by default).
    _emit( 'info', str(msg) )

def warn(msg):
    #Emit a warning (printed as "NCrystal WARNING: <msg>" by default).
    _emit( 'warning', str(msg) )

def print(*args,**kwargs):
    #Raw output (like the builtin print, but via the print function of
    #NCrystal and active MsgCtrl contexts).
    if not _active_msgctrls and not ( _msg_global['logger']
                                      and _msg_global['include_raw'] ):
        _stdprint[0](*args,**kwargs)
        return
    sep, end = kwargs.get('sep',' '), kwargs.get('end','\n')
    text = ( ( ' ' if sep is None else sep ).join( str(a) for a in args )
             + ( '\n' if end is None else end ) )
    _emit( 'raw', text, ( args, kwargs ) )

def _emit( part, text, rawcall = None ):
    ctrls = _active_msgctrls
    show = True
    for c in ctrls:
        if not c._take( part, text, 'NCrystalWarning' ):
            show = False
    if not show or any( part in c._block for c in ctrls ):
        return
    #Settings of the innermost context which has them:
    logger, include_raw = _msg_global['logger'], _msg_global['include_raw']
    prefix = _msg_default_prefix.get( part )
    for c in ctrls:
        if c._logger is not None:
            logger = c._logger or None
        if c._include_raw is not None:
            include_raw = c._include_raw
        if c._prefix.get( part ) is not None:
            prefix = c._prefix[part]
    if logger is not None and any( part in c._ignore_logging for c in ctrls ):
        logger = None
    if part == 'raw':
        if logger is not None and include_raw:
            logger.info( text.removesuffix('\n') )
        else:
            args, kwargs = rawcall
            _stdprint[0]( *args, **kwargs )
    elif logger is not None:
        ( logger.info if part == 'info' else logger.warning )( text )
    else:
        _stdprint[0]( f'{prefix}{text}' )

def _fmt_pywarning( message, category ):
    #Python warning as (message,category) strings:
    catname = None
    if category is not None:
        catname = str(getattr(category,'__name__',category)).strip()
    elif isinstance(message,Warning):
        catname = str(message.__class__.__name__).strip()
    if not catname:
        catname='Warning'
    return ( ( str(message).strip() or '<no message>'), catname )

_orig_showwarning = [None]

def _msgctrl_showwarning( message, category, *args, **kwargs ):
    #warnings.showwarning while a MsgCtrl handling warnings is active:
    msg_str, cat_str = _fmt_pywarning( message, category )
    ctrls = _active_msgctrls
    show = True
    for c in ctrls:
        if not c._take( 'warning', msg_str, cat_str ):
            show = False
    if show and not any( 'warning' in c._block for c in ctrls ):
        return _orig_showwarning[0]( message, category, *args, **kwargs )

def warn_deprecated( msg, stacklevel = 3 ):
    #Emit an NCrystalFutureWarning about deprecated functionality via the
    #warnings module, attributed to the caller of the function calling this
    #(or higher up, with a larger stacklevel).
    import warnings

    from .exceptions import NCrystalFutureWarning
    warnings.warn( NCrystalFutureWarning( str(msg) ), stacklevel = stacklevel )

def _split_trailing_digit( s ):
    if not s or not s[-1].isdigit():
        return s,None
    n=0
    while (n+1)<len(s) and s[-(n+1)].isdigit():
        n += 1
    if n:
        return s[0:-n], int( s[-n:] )
    else:
        return s, None

def check_elem_or_isotope_marker( s ):
    """If input is of form "Al", "O16", ..., return it, otherwise return
       None. Ignores excess whitespace in input. Will also map "H2"->"D" and
       "H3"->"T".
    """
    s=' '.join(str(s).split())
    if not s:
        return None
    elem_name, isotope_val = _split_trailing_digit( s )

    #special support for D,T:
    if not isotope_val and elem_name in ('D','T'):
        return elem_name
    if elem_name == 'H' and isotope_val in (2,3):
        return 'D' if isotope_val == 2 else 'T'

    from .atomdata import isElementName

    if isElementName( elem_name) and ( isotope_val is None or ( 0 <= isotope_val < 999 ) ):
        return elem_name if not isotope_val else f'{elem_name}{isotope_val}'


def _classifySG(sgno):
    assert 1<=sgno<=230
    ll=[(195,'cubic'),(168,'hexagonal'),(143,'trigonal'),
        (75,'tetragonal'),(16,'orthorombic'),(3,'monoclinic'),(1,'triclinic')]
    for thr,nme in ll:
        if sgno>=thr:
            return nme
    assert False

#colors inspired by http://www.mulinblog.com/a-color-palette-optimized-for-data-visualization/
#supposedly from Stephen Few's book, "Show Me the Numbers":
_palette_Few = {'red': "#F15854",
                'blue': "#5DA5DA",
                'orange': "#FAA43A",
                'green': "#60BD68",
                'brown': "#B2912F",
                'purple': "#B276B2",
                'yellow': "#DECF3F",
                'pink': "#F17CB0",
                'gray': "#4D4D4D"}

def _grid_is_linspace( grid, tol = 1e-6 ):
    if len(grid)<=2:
        return len(grid)==2
    from ._numpy import _ensure_numpy, _np
    _ensure_numpy()
    g = _np.asarray(grid,dtype=float)
    bws = g[1:] - g[:-1]
    bmin, bmax = bws.min(), bws.max()
    if ( not ( bmax >= bmin > 0.0 ) ) or _np.isinf(bmin) or _np.isinf(bmax):
        return False
    return ( bmax - bmin ) < tol * ( bmax + bmin )

def extract_path( s ):
    """Returns pathlib.Path path object from argument if it is a path, otherwise
    None. Strings with no newlines will be assumed to be paths."""
    import pathlib
    if hasattr( s, '__fspath__' ):
        return pathlib.Path(s)
    if isinstance( s, str ) and '\n' not in s:
        res_try = _lookup_existing_file(s)
        if res_try:
            return res_try
        return pathlib.Path(s)

def download_url( url,
                  decode_as_utf8_str = True,
                  wrap_exception = True,
                  timeout = None,
                  quiet_network_fail = False ):
    """Download the provided url. The timeout value is in seconds if given.
    Returns None on network failure if quiet_network_fail is True.
    """
    import urllib.error
    import urllib.request
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout = timeout ) as response:
            data = response.read()
    except urllib.error.URLError as e:
        if quiet_network_fail:
            return None
        if wrap_exception:
            from .exceptions import NCException
            raise NCException(f'Error downloading url "{url}": {e}') from e
        else:
            raise
    if decode_as_utf8_str:
        try:
            data = data.decode('utf8')
        except UnicodeDecodeError as e:
            if wrap_exception:
                from .exceptions import NCException
                raise NCException('Error decoding url to utf8 data'
                                  f' "{url}": {e}') from e
            else:
                raise
    return data

def _decodeflt(s):
    import numbers
    if s is None:
        return None
    if isinstance(s,numbers.Real):
        return float(s)
    try:
        x = float( s )
    except (TypeError,ValueError):
        return None
    return x

def _calc_md5hexdigest( str_or_bytes, / ):
    import hashlib
    if hasattr(str_or_bytes,'encode'):
        data = str_or_bytes.encode('utf8',errors='backslashreplace')
    else:
        data = str_or_bytes
    return hashlib.md5( data ).hexdigest()

def write_text( path, content ):
    """Like path.write_text(content) but forcing some global settings
    to ensure consistent NCrystal behaviour on any platform. Specifically
    encoding='utf8' and newline='\n'"""
    #The newline parameter was only added for pathlib.Path.write_text in python
    #3.10, so we do instead:
    import pathlib
    with pathlib.Path(path).open( 'wt',
                                  encoding = 'utf8',
                                  newline='\n' ) as fh:
        fh.write(content)

def _lookup_existing_file( path ):
    """Will wrap in pathlib.Path, except in simplebuild development mode where
    it will try $SBLD_DATA_DIR/+path when the path is a string with a single
    path separator. Returns None if unable to find an existing file at the
    location, otherwise an absolute path to the file.

    """
    import pathlib
    p = pathlib.Path( path ).expanduser()
    if p.exists():
        return p.absolute()
    if ( isinstance( path, str )
         and path.count('/')==1
         and ( pathlib.Path(__file__).parent
             .joinpath('_is_sblddevel.py').is_file() ) ):
        import os
        dd = os.environ.get('SBLD_DATA_DIR')
        if dd:
            p = pathlib.Path(dd).joinpath(path)
            if p.exists():
                return p.absolute()

_override_datetime_now = [ None ]
def _datetime_now():
    """Returns datetime.datetime.now() but possibly intercepted for unit tests"""
    import datetime
    n = datetime.datetime.now(datetime.timezone.utc)
    return _override_datetime_now[0] or n

class FixedFakeDatetimeNow:
    # Context manager to be used in unit tests to modify the returned value from
    # the datetime_now function below, to a fixed value.
    def __enter__(self):
        import datetime
        self.__orig = _override_datetime_now[0]
        _override_datetime_now[0] = datetime.datetime( 2017, 8, 29, 13, 40,
                                                       tzinfo
                                                       = datetime.timezone.utc )
    def __exit__(self,*args,**kwargs):
        _override_datetime_now[0] = self.__orig

def fixed_fake_datetime_now(f):
    #To be used as a decorator
    def fw(*a, **kw):
        with FixedFakeDatetimeNow():
            return f(*a, **kw)
    return fw

def flex_load_json( data, copy = False ):
    """If data is not a string or buffer, simply return it (possibly after a
    deep copy if copy=True).  Otherwise, it is assumed to be JSON data encoded
    in either a string, a utf8-encoded buffer, or a utf8-encoded and
    gzip-compressed buffer. In that case, the JSON data is loaded and the
    resulting object returned.

    The function also accepts file paths from which to read data (like
    pathlib.Path objects), but str or buffer objects will be assumed to contain
    JSON data.

    Thus, this function can be used where we need to accept flexible input data
    and return a decoded object (normally a dictionary or a list).
    """
    if isinstance( data, bytes ):
        if data.startswith(b'\x1f\x8b'):
            import gzip
            data = gzip.decompress(data)
        data = data.decode('utf8')
    if isinstance( data, str ):
        import json
        data = json.loads( data )
    if hasattr( data, '__fspath__' ):
        import pathlib
        return flex_load_json(pathlib.Path(data).read_bytes(),copy=copy)
    if copy:
        import copy
        data = copy.deepcopy(data)
    return data

def copy_and_deobjectify_data( data ):
    """Traverse lists, tuples and dicts of data and replace any objects inside
    it that has a ._to_json_compat_object() method with the result of calling
    that method.
    """
    import copy
    def o( d ):
        if hasattr( d, '_to_json_compat_object'):
            return d._to_json_compat_object()
        if isinstance(d,list):
            return [o(e) for e in d]
        if isinstance(d,tuple):
            return tuple([o(e) for e in d])
        if isinstance(d,dict):
            return { copy.deepcopy(k): o(v) for k,v in d.items() }
        #something else, just pass through:
        return copy.deepcopy( d )
    return o( data )

class ReadOnlyDict(dict):
    """Read-only dict: behaves like a dict for all lookups, printing, JSON
    encoding and isinstance(..,dict) checks, but any attempt to modify it
    raises a TypeError. Nested dicts are also read-only and nested lists are
    tuples. Use .copy() (or dict(x)) for an editable shallow copy, and
    copy.deepcopy(x) for an editable copy of the whole structure."""

    def __readonly(self, *a, **kw):
        raise TypeError('read-only dict (use .copy() for an editable copy)')
    __setitem__ = __delitem__ = clear = pop = popitem = __readonly
    setdefault = update = __ior__ = __readonly

    def copy(self):
        """Editable shallow copy as a plain dict."""
        return dict(self)
    __copy__ = copy

    def __deepcopy__(self, memo):
        return _create_mutable_copy(self)

    def __reduce__(self):
        return (dict, (dict(self),))

def create_read_only_view(data):
    #Create read-only view of data, assuming only dict/list/tuple + basic
    #immutable types are used.
    if isinstance(data, dict):
        return ReadOnlyDict( (k, create_read_only_view(v))
                             for k,v in data.items() )
    elif isinstance(data, (list, tuple)):
        return tuple(create_read_only_view(v) for v in data)
    else:
        return data

def _create_mutable_copy(data):
    #Inverse of create_read_only_view: dicts become plain dicts, tuples lists.
    if isinstance(data, dict):
        return { k: _create_mutable_copy(v) for k,v in data.items() }
    elif isinstance(data, (list, tuple)):
        return [ _create_mutable_copy(v) for v in data ]
    else:
        return data

def ncpprint(obj, *, do_sort=False):
    #pprint.pprint-like, but via our print (dicts sorted only if do_sort):
    import pprint
    print(pprint.pformat(obj, sort_dicts=bool(do_sort)))
