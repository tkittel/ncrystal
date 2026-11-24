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

# Test of NCrystal.msg: notices, warnings and raw output (from both Python and
# C++), printed by default, logged after use_logging(), and controlled and
# collected by (nested) MsgCtrl contexts. Also msg_style of cli.run.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev as NC
import NCrystalDev.cli as nc_cli
from NCrystalDev._common import infomsg, warn, set_ncrystal_print_fct
from NCrystalDev._common import print as ncprint
from NCrystalDev.msg import MsgCtrl, use_logging
from NCTestUtils.common import ensure_error
import logging
import sys
import warnings

#Debye temperature + VDOS triggers a deterministic NCRYSTAL_WARN on load:
data = """NCMAT v7
@CELL
 cubic 4.04958
@SPACEGROUP
 225
@ATOMPOSITIONS
 Al 0 1/2 1/2
 Al 0 0 0
 Al 1/2 1/2 0
 Al 1/2 0 1/2
@DEBYETEMPERATURE
 Al 400
@DYNINFO
 element Al
 fraction 1
 type vdos
 vdos_egrid 0.01 0.04
 vdos_density 1 4 9 16 25 20 5
"""

nload = [0]
def cxxwarn():
    #Unique comment to avoid any caching:
    nload[0] += 1
    NC.directLoad( data + f'#{nload[0]}\n',
                   doScatter = False, doAbsorption = False )

def emit_all():
    infomsg( 'a notice' )
    warn( 'a warning' )
    warn( 'a warning' )#(repeated, not deduplicated without dedup)
    ncprint( 'raw output' )
    cxxwarn()

def show( title, mc ):
    print( f'   {title}: info={mc.info} warnings={mc.warnings} raw={mc.raw}' )

class LogToStdout:
    #Logger "NCrystal" (or other) printing records to stdout:
    def __init__( self, name = 'NCrystal' ):
        self.logger = logging.getLogger( name )
        self.h = logging.StreamHandler( sys.stdout )
        self.h.setFormatter( logging.Formatter(
            'LOG %(levelname)s %(name)s: %(message)s' ) )
    def __enter__( self ):
        self.logger.addHandler( self.h )
        self.logger.setLevel( logging.INFO )
        self.logger.propagate = False
        return self.logger
    def __exit__( self, *a ):
        self.logger.removeHandler( self.h )

def test_default():
    print('==> Default: printed')
    emit_all()
    print('==> Captured by the print function:')
    captured = []
    orig = set_ncrystal_print_fct( lambda *a, **kw : captured.append(
        ' '.join( str(e) for e in a ) ) )
    emit_all()
    set_ncrystal_print_fct( orig )
    for c in captured:
        print( f'   captured: {c}' )

def test_logging():
    with LogToStdout():
        print('==> use_logging(): notices and warnings logged')
        use_logging()
        emit_all()
        print('==> use_logging(include_raw=True): raw output logged too')
        use_logging( include_raw = True )
        emit_all()
        print('==> MsgCtrl collects regardless of logging:')
        with MsgCtrl( spy = ['info','warning','raw'],
                      block = ['info','warning','raw'] ) as mc:
            emit_all()
        show( 'collected', mc )
        print('==> MsgCtrl(logger=False) prints despite use_logging():')
        with MsgCtrl( logger = False ):
            emit_all()
        print('==> ignore_logging=["warning"] (and inner contexts can not undo it):')
        use_logging()
        with MsgCtrl( ignore_logging = ['warning'] ), MsgCtrl( ignore_logging = ['info'] ):
            emit_all()
        use_logging( False )
    with LogToStdout( 'myapp' ) as mylogger:
        print('==> MsgCtrl(logger=<custom logger>, include_raw=True):')
        with MsgCtrl( logger = mylogger, include_raw = True ):
            emit_all()
        print('==> Inner MsgCtrl(logger=False) wins:')
        with MsgCtrl( logger = mylogger ), MsgCtrl( logger = False ):
            emit_all()
    print('==> After use_logging(False), printed again:')
    emit_all()

def test_prefixes():
    print('==> Prefixes (each context its own, innermost wins):')
    with MsgCtrl( info_prefix = '', warning_prefix = 'WARNING: ' ):
        emit_all()
        with MsgCtrl( info_prefix = 'INNER: ' ):
            emit_all()
        warn( 'outer again' )

def test_spy_block_dedup():
    print('==> Nested: outer blocks and collects, inner dedups and collects:')
    with MsgCtrl( spy = ['info','warning','raw'],
                  block = ['info','warning','raw'] ) as outer:
        with MsgCtrl( spy = ['info','warning','raw'],
                      dedup = ['info','warning','raw'] ) as inner:
            emit_all()
            emit_all()
        warn( 'after inner' )
    show( 'inner', inner )
    show( 'outer', outer )
    print('==> block is cumulative (inner can not unblock):')
    with MsgCtrl( block = ['info'] ), MsgCtrl( block = ['raw'] ):
        emit_all()
    print('==> Inner blocks, outer still collects:')
    with MsgCtrl( spy = ['warning'] ) as outer, \
         MsgCtrl( spy = ['warning'], block = ['warning'] ) as inner:
        warn( 'blocked by inner' )
        warnings.warn( 'a python warning' )
    show( 'inner', inner )
    show( 'outer', outer )
    print('==> Python warnings are shown unless blocked:')
    with MsgCtrl( spy = ['warning'] ) as mc:
        warnings.warn( 'a python warning, shown' )
    show( 'mc', mc )
    print('==> Each context starts with a clean slate:')
    for i in range(2):
        with MsgCtrl( spy = ['warning'], dedup = ['warning'],
                      block = ['warning'] ) as mc:
            warn( 'same' )
            warn( 'same' )
        show( f'context {i}', mc )
    print('==> Single string accepted as parts:')
    with MsgCtrl( spy = 'raw', block = 'raw' ) as mc:
        ncprint( 'x' )
    show( 'mc', mc )

def pywarn( msg ):
    warnings.warn( msg )#always from this same location

def test_pywarn_filters():
    print('==> Python warning filters and deduplication are respected:')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with MsgCtrl( spy = ['warning'] ) as mc:
            pywarn( 'ignored by filter' )
        show( 'ignore filter', mc )
    with warnings.catch_warnings():
        warnings.simplefilter('default')
        filters_before = list( warnings.filters )
        with MsgCtrl( spy = ['warning'], block = ['warning'] ) as mc:
            pywarn( 'once per location' )
            pywarn( 'once per location' )
        show( 'default filter', mc )
        #Not shown again after the context ended:
        with MsgCtrl( spy = ['warning'], block = ['warning'] ) as mc:
            pywarn( 'once per location' )
        show( 'next context', mc )
        print( '   filters unchanged:', warnings.filters == filters_before )

def test_errors():
    print('==> Errors:')
    with ensure_error( NC.NCBadInput, "Invalid MsgCtrl spy part(s): ['bla']"
                       " (must be among ['info', 'warning', 'raw'])" ):
        MsgCtrl( spy = ['info','bla'] )
    with ensure_error( NC.NCBadInput, 'MsgCtrl logger must be None, False,'
                       ' True or a logging.Logger' ):
        MsgCtrl( logger = 'NCrystal' )

def test_cli_msgstyle():
    for style in ( 'api', 'cli' ):
        print( f'==> cli.run with msg_style="{style}" (dump output blocked):' )
        NC.registerInMemoryFileData( f'msgtest_{style}.ncmat', data )
        with MsgCtrl( block = ['raw'] ):
            nc_cli.run( 'nctool', '--dump', f'msgtest_{style}.ncmat;inelas=0',
                        msg_style = style )
    with ensure_error( NC.NCBadInput, 'Invalid msg_style: "bla" (must be "api"'
                       ' or "cli")' ):
        nc_cli.run( 'nctool', '--help', msg_style = 'bla' )

def main():
    test_default()
    test_logging()
    test_prefixes()
    test_spy_block_dedup()
    test_pywarn_filters()
    test_errors()
    test_cli_msgstyle()

if __name__ == '__main__':
    main()
