
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

Utilities for NCrystal cfg-strings, like normalising or decoding them, and
generating documentation about the available cfg-string parameters.

"""

def normaliseCfg(cfgstr):
    """Returns normalised version of cfgstr. This is done behind the scenes by
       loading the specified cfg-string into a C++ MatCfg object and then
       re-encoding it as a string.
    """
    from .misc import evaluate_query
    return evaluate_query( ['matcfg','normalise',cfgstr] )

def decodeCfg( material, cfg_params = None, *, asJSONStr = False ):
    """Decodes the cfg of a material and returns it as a Python data structure
       (a dictionary). The material can be a cfg-string, raw NCMAT data, a
       TextData or NCMATComposer object, or anything else accepted by
       NCrystal.misc.MaterialSource, and cfg_params are extra cfg parameters
       to apply. The "allpars" entry lists all parameters, including those at
       default values. The format of this data structure should be mostly
       self-evident by inspection, and is not guaranteed to stay the same
       across NCrystal versions. If asJSONStr=true, the data structure will
       be returned as a JSON-encoded string, instead of a Python dictionary.
    """
    from .misc import MaterialSource
    _js = MaterialSource( material, cfg_params = cfg_params )._decodecfg_json()
    if asJSONStr:
        return _js
    import json
    return json.loads(_js)

def generateCfgStrDoc( mode = "print" ):
    """Generates documentation about the available cfg-str variables. Mode can
    either be 'print' (print detailed explanation to stdout), 'txt_full' (return
    detailed explanations as string), 'txt_short' (return short explanations as
    string), 'json' (return json-encoded string), or 'python' (return python
    data structures).
    """
    modemap={'print':'full','txt_full':'full','txt_short':'short',
             'json':'json','python':'json'}
    qmode = modemap.get(mode,None)
    if qmode is None:
        from .exceptions import NCBadInput
        raise NCBadInput(f'mode must be one of {sorted(modemap.keys())}')
    from .misc import evaluate_query
    _=evaluate_query( ['matcfg','doc',qmode], unpack = ( mode != 'json' ) )
    if mode == 'print':
        from ._common import print
        print(_)
    else:
        return _

def decodecfg_vdoslux( material, cfg_params = None ):
    """Extract vdoslux value from the cfg of a material (accepting the same
    arguments as decodeCfg)."""
    def get( d ):
        if not d['ismultiphase']:
            return { v for n,v,_ in d['allpars'] if n == 'vdoslux' }
        return set().union( *( get(ph) for _,ph in d['phases'] ) )
    vals = get( decodeCfg( material, cfg_params ) )
    if len(vals) != 1:
        from .exceptions import NCCalcError
        raise NCCalcError( 'Could not determine unique value of parameter'
                           ' "vdoslux" on multiphase material (different'
                           ' values found in different phases).' )
    return int( vals.pop() )
