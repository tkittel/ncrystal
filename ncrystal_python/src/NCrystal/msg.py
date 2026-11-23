
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

"""Control of the messages from NCrystal: notices, warnings and raw output.

By default, notices and warnings (from both the Python and the C++ code) are
printed with the prefixes "NCrystal: " and "NCrystal WARNING: ", and raw output
(e.g. dumps) is printed as it is. With use_logging(), notices and warnings are
instead sent to the logger "NCrystal". The MsgCtrl context manager can
temporarily change this, and collect or suppress messages.

"""

__all__ = [ 'MsgCtrl', 'use_logging' ]

def use_logging( enable = True, include_raw = False ):
    """Send notices and warnings from NCrystal (from both the Python and C++
    code) to the logger named "NCrystal" (as logging.INFO and logging.WARNING
    records), instead of printing them, which is the default. If include_raw
    is True, raw output (e.g. dumps) is also sent there (as logging.INFO
    records). After that, they are controlled by the logging configuration of
    the application, e.g.:

       import logging
       import NCrystal
       NCrystal.use_logging()
       logging.basicConfig(level=logging.INFO)

    Calling use_logging(False) switches back to printing. MsgCtrl contexts can
    override this temporarily."""
    from ._common import _msg_global
    if enable:
        import logging
        _msg_global['logger'] = logging.getLogger('NCrystal')
        _msg_global['include_raw'] = bool(include_raw)
    else:
        _msg_global['logger'] = None
        _msg_global['include_raw'] = False

class MsgCtrl:
    """Context manager controlling the messages from NCrystal while active.
    Messages come in three parts: 'info' (notices), 'warning' (warnings from
    NCrystal, and Python warnings via the warnings module) and 'raw' (raw
    output, e.g. dumps). Parameters left as None do not change the settings of
    enclosing contexts:

      info_prefix, warning_prefix : Prefixes of printed notices and warnings
                                    (default "NCrystal: " and
                                    "NCrystal WARNING: ").
      logger : False (print), True (the logger "NCrystal") or a
               logging.Logger, for notices and warnings (cf. use_logging).
      include_raw : Whether a logger also gets the raw output.
      spy : Parts to collect, available as the lists .info, .warnings (of
            (message,category) tuples, with category "NCrystalWarning" for
            those from NCrystal) and .raw (the printed texts).
      block : Parts not shown (they are still collected).
      ignore_logging : Parts printed even if a logger is in use.
      dedup : Parts on which repeated messages are neither shown nor
              collected again while the context is active.

    Parts are given as a list (or a single string). Prefixes, logger and
    include_raw are those of the innermost context setting them, while block,
    ignore_logging and dedup apply if any active context has them (so inner
    contexts can not undo them). Each context with spy collects its parts,
    regardless of the settings of other contexts. While a context handles
    warnings (in spy, block or dedup), the filters of the warnings module are
    disabled, so all Python warnings are seen.

    Example:

       with NCrystal.msg.MsgCtrl( spy = ['warning'], block = ['warning'] ) as mc:
           mat = NCrystal.load('Al_sg225.ncmat')
       print( mc.warnings )
    """

    def __init__( self, *, info_prefix = None, warning_prefix = None,
                  logger = None, include_raw = None, spy = None, block = None,
                  ignore_logging = None, dedup = None ):
        from ._common import _msg_parts
        from .exceptions import NCBadInput
        def parts( x, name ):
            if x is None:
                return frozenset()
            if isinstance( x, str ):
                x = [ x ]
            res = frozenset( x )
            bad = sorted( res - set( _msg_parts ) )
            if bad:
                raise NCBadInput( f'Invalid MsgCtrl {name} part(s): {bad}'
                                  f' (must be among {list(_msg_parts)})' )
            return res
        import logging
        if logger is True:
            logger = logging.getLogger('NCrystal')
        elif not ( logger is None or logger is False
                   or isinstance( logger, logging.Logger ) ):
            raise NCBadInput( 'MsgCtrl logger must be None, False, True or a'
                              ' logging.Logger' )
        self._prefix = { 'info' : info_prefix, 'warning' : warning_prefix }
        self._logger = logger
        self._include_raw = None if include_raw is None else bool(include_raw)
        self._spy = parts( spy, 'spy' )
        self._block = parts( block, 'block' )
        self._ignore_logging = parts( ignore_logging, 'ignore_logging' )
        self._dedup = parts( dedup, 'dedup' )
        self._seen = set()
        self._data = { 'info' : [], 'warning' : [], 'raw' : [] }
        self._warnctx = None

    @property
    def info( self ):
        """Collected notices (if 'info' in spy)."""
        return self._data['info']

    @property
    def warnings( self ):
        """Collected warnings as (message,category) tuples (if 'warning' in
        spy)."""
        return self._data['warning']

    @property
    def raw( self ):
        """Collected raw output texts (if 'raw' in spy)."""
        return self._data['raw']

    def _take( self, part, text, category ):
        #Collect the message (unless a repeat with dedup) and return whether
        #it may be shown as far as this context is concerned:
        if part in self._dedup:
            key = ( part, text, category if part == 'warning' else None )
            if key in self._seen:
                return False
            self._seen.add( key )
        if part in self._spy:
            self._data[part].append( ( text, category ) if part == 'warning'
                                     else text )
        return True

    def __enter__( self ):
        from . import _common
        if 'warning' in ( self._spy | self._block | self._dedup ):
            import warnings
            self._warnctx = warnings.catch_warnings()
            self._warnctx.__enter__()
            if warnings.showwarning is not _common._msgctrl_showwarning:
                _common._orig_showwarning[0] = warnings.showwarning
                warnings.showwarning = _common._msgctrl_showwarning
            warnings.simplefilter('always')
        _common._active_msgctrls.append( self )
        return self

    def __exit__( self, *args, **kwargs ):
        from . import _common
        assert _common._active_msgctrls and _common._active_msgctrls[-1] is self
        _common._active_msgctrls.pop()
        if self._warnctx is not None:
            self._warnctx.__exit__( *args, **kwargs )
            self._warnctx = None
        return False
