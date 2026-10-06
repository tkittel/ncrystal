
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


"""Utilities for tests which must work both with and without thread support
in NCrystal (cf. the NCRYSTAL_ENABLE_THREADS CMake option)."""

def has_threads():
    """Whether NCrystal was built with thread support."""
    from NCrystalDev.misc import evaluate_query
    res = evaluate_query(['util','has_threads'],unpack=False)
    assert res in ('[true]','[false]')
    return res == '[true]'

_nothreads_warning = 'NCrystal installation does not support threads.'

def hide_nothreads_warnings():
    """For tests requesting MiniMC simulations with several threads: in
    builds without thread support, hide the warnings about running in a
    single thread instead, so the output is the same as with threads."""
    if has_threads():
        return
    from NCrystalDev._msg import _default_pymsghandler, _setMsgHandler
    def handler( msg, msgtype ):
        if not ( msgtype == 1 and msg.startswith(_nothreads_warning) ):
            _default_pymsghandler( msg, msgtype )
    _setMsgHandler( handler )
