
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

"""

Internal implementation details for NCrystal output handling (see also
functions in _common.py).

"""

def _default_pymsghandler( msg, msgtype ):
    #Messages from the C++ library: notices (NCRYSTAL_MSG), warnings
    #(NCRYSTAL_WARN) and raw output (NCRYSTAL_RAWOUT):
    from . import _common
    if msgtype == 0:
        _common.infomsg(msg)
    elif msgtype == 1:
        _common.warn(msg)
    else:
        assert msgtype == 2
        _common.print(msg,end='')

#TODO: Overlaps somewhat with the (py-only) set_ncrystal_print_fct from _common.py:

#NB: This next function could become part of a public API, allowing e.g. a GUI
#to redirect all NCrystal output to appropriate text boxes, etc.:
_was_set = [False]
_current = [None]#None means the C++ default handler
def _setMsgHandler( handler ):
    from ._chooks import _get_raw_cfcts
    _rawfct = _get_raw_cfcts()
    _rawfct['setmsghandler'](handler)
    _was_set[0] = True
    _current[0] = handler

def _setDefaultPyMsgHandlerIfNotSet():
    if not _was_set[0]:
        _setMsgHandler(_default_pymsghandler)
