
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


"""Internal module providing optional access to numpy (as _np, which evaluates
to False if numpy is unavailable) along with a few numpy-related helper
functions"""

__all__ = [
    '_ensure_numpy',
    '_np',
    '_np_geomspace',
    '_np_linspace',
    '_np_logspace',
    '_np_trapezoid',
]

_numpy_mod = []#loaded module (or None if unavailable), once tried

class _LazyNumpy:
    """Stand-in for the numpy module, which is only imported upon first use
    (saving ~30ms when importing NCrystal for tasks not needing numpy). It
    evaluates to False if numpy is unavailable. Attributes are cached on the
    instance, so later lookups bypass __getattr__ (faster than the module)."""
    @staticmethod
    def _load():
        if not _numpy_mod:
            try:
                import numpy
            except ImportError:
                numpy = None
            _numpy_mod.append( numpy )
        return _numpy_mod[0]
    def __bool__( self ):
        return _LazyNumpy._load() is not None
    def __getattr__( self, name ):
        m = _LazyNumpy._load()
        if m is None:
            if name.startswith('__'):
                #Python's own probing (hasattr, copy, pickle) expects this:
                raise AttributeError( name )
            _ensure_numpy()
        v = getattr( m, name )
        setattr( self, name, v )
        return v

_np = _LazyNumpy()

def _ensure_numpy():
    if not _np:
        from .exceptions import NCException
        raise NCException("Numpy not available - array "
                          "based functionality is unavailable")
    return _np

def _np_linspace(start,stop,num=50):
    """linspace with reproducible endpoint value"""
    _ensure_numpy()
    assert num >= 2
    ll = _np.linspace(start,stop,num)
    ll[0] = start
    ll[-1] = stop
    return ll

def _np_geomspace(start,stop,num=50):
    """geomspace with reproducible endpoint value"""
    _ensure_numpy()
    assert num >= 2
    ll = _np.geomspace(start,stop,num)
    ll[0] = start
    ll[-1] = stop
    return ll

def _np_logspace(start,stop,num=50):
    """logspace with reproducible endpoint value"""
    _ensure_numpy()
    assert num >= 2
    ll = _np.logspace(start,stop,num)
    ll[0] = 10.0**start
    ll[-1] = 10.0**stop
    return ll

def _np_trapezoid( *args, **kwargs ):
    _ensure_numpy()
    if hasattr(_np,'trapezoid'):
        return _np.trapezoid( *args, **kwargs )
    return _np.trapz( *args, **kwargs )
