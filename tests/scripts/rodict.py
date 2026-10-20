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

# Read-only views of JSON-like data (dicts/lists/tuples of basic immutable
# values), as returned by cached queries: dicts become _ReadOnlyDict (a dict
# subclass with the mutating methods disabled, so repr/pprint/json/isinstance
# behave as for a dict), lists become tuples. Copies are plain mutable objects.

import copy
import json
import operator
import pickle
import pprint

################################################################################
# Code for NCrystal/_common.py:
from NCrystalDev._common import _create_mutable_copy, create_read_only_view

################################################################################
# Tests:

def ensure_readonly( f ):
    try:
        f()
    except TypeError as e:
        print('  blocked:', e)
        return
    raise RuntimeError('mutation not blocked')

data = { 'a': 4.04958, 'spacegroup': 225, 'spacegroup_setting': '225',
         'sites': [ { 'position': [0.0, 0.0, 0.0], 'multiplicity': 4 } ],
         'nothing': None, 'flag': True }
ro = create_read_only_view( data )

print('type:', type(ro).__name__, '| isinstance dict:', isinstance(ro,dict))
print('repr:', repr(ro))
pprint.pprint(ro)
print('json:', json.dumps(ro, sort_keys=True))
print('access:', ro['a'], ro['sites'][0]['position'], ro.get('missing','dflt'),
      len(ro), sorted(ro.keys()), 'sites' in ro)
print('nested types:', type(ro['sites']).__name__,
      type(ro['sites'][0]).__name__,
      type(ro['sites'][0]['position']).__name__)
#NB: ro == data is False, since nested lists became tuples (dict equality
#compares values); the mutable copy compares equal:
print('equality with original dict:', ro == data,
      '| mutable copy equal:', _create_mutable_copy(ro) == data)

print('mutations:')
ensure_readonly( lambda : operator.setitem(ro, 'a', 1.0) )
ensure_readonly( lambda : operator.delitem(ro, 'a') )
ensure_readonly( lambda : ro.clear() )
ensure_readonly( lambda : ro.pop('a') )
ensure_readonly( lambda : ro.popitem() )
ensure_readonly( lambda : ro.setdefault('new', 1) )
ensure_readonly( lambda : ro.update( {'new':1} ) )
ensure_readonly( lambda : ro.__ior__( {'new':1} ) )
ensure_readonly( lambda : operator.setitem(ro['sites'][0], 'multiplicity', 8) )
ensure_readonly( lambda : operator.setitem(ro['sites'], 0, None) )
assert _create_mutable_copy(ro) == data

print('copies:')
for name, c in [ ('.copy()', ro.copy()),
                 ('dict()', dict(ro)),
                 ('copy.copy', copy.copy(ro)),
                 ('pickle', pickle.loads(pickle.dumps(ro))) ]:
    assert type(c) is dict and c == ro
    c['a'] = 1.0 # editable at top level
    print(f'  {name}: plain dict, nested still read-only:',
          type(c['sites']).__name__, type(c['sites'][0]).__name__)
dc = copy.deepcopy(ro)
mc = _create_mutable_copy(ro)
assert type(dc) is dict and dc == data and mc == data
print('  _create_mutable_copy: same as deepcopy:', mc == dc)
dc['sites'][0]['position'][0] = 0.5 # editable all the way down
print('  copy.deepcopy: fully mutable:', type(dc['sites']).__name__,
      type(dc['sites'][0]).__name__, type(dc['sites'][0]['position']).__name__)
assert ro['sites'][0]['position'][0] == 0.0

print('view of a view and of None:', type(create_read_only_view(ro)).__name__,
      create_read_only_view(None))
print('all ok')
