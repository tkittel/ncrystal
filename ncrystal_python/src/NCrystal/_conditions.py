
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

"""Internal conditions like "absxs>100", for selecting entries by property
values."""

def _docs_for( docs, kinds, conditions ):
    return [ (n,d) for n,d in docs
             if not ( conditions and kinds.get(n) == 'nocond' ) ]

def _bad_name_msg( what, name, names, context = '' ):
    #Error message for an unknown name, with a hint for wrong casing (e.g.
    #"SG" or "a"):
    m = [ n for n in names if n.lower() == name.lower() ]
    hint = ( f'names are case-sensitive, did you mean "{m[0]}"?' if m
             else 'must be one of: ' + ', '.join( names ) )
    return f'Invalid {what} "{name}"{context} ({hint})'

_bool_values = { 'yes': True, 'true': True, 'on': True, '1': True,
                 'no': False, 'false': False, 'off': False, '0': False }

def _condition_re():
    #Conditions like "absxs>100":
    import re
    return re.compile( r'\s*([A-Za-z_][A-Za-z0-9_]*)\s*'
                       r'(<=|>=|!=|==|=|<|>)\s*(.*?)\s*' )

def _is_condition( s ):
    #True if s has the form of a condition (like "absxs>100"):
    return _condition_re().fullmatch( s ) is not None

_condition_ops = { '=': lambda a,b: a == b, '!=': lambda a,b: a != b,
                   '<': lambda a,b: a < b, '<=': lambda a,b: a <= b,
                   '>': lambda a,b: a > b, '>=': lambda a,b: a >= b }

class _Condition:
    #Where condition like "sg=225" or "absxs>100", called with a dict of
    #values:

    def __init__( self, cond, names, kinds ):
        import re

        from .exceptions import NCBadInput
        def bad( msg ):
            return NCBadInput(f'Invalid condition "{cond}" ({msg})')
        m = _condition_re().fullmatch( str(cond) )
        if not m:
            raise bad('must be like "absxs>100" or sg==225')
        key, op, val = m.groups()
        op = '=' if op == '==' else op
        if key not in names:
            usable = [ n for n in names if kinds.get( n ) != 'nocond' ]
            raise NCBadInput( _bad_name_msg( 'property', key, usable,
                                             f' in condition "{cond}"' ) )
        if len(val) >= 2 and val[0] == val[-1] and val[0] in '"\'':
            val = val[1:-1]
        if not val:
            raise bad('missing value')
        kind = kinds.get( key, 'num' )
        allowed = None
        if isinstance( kind, tuple ):
            allowed, kind = kind, 'str'
        if kind == 'nocond':
            raise bad(f'"{key}" can not be used in conditions')
        if kind != 'num' and op not in ('=','!='):
            raise bad(f'only = and != can be used with "{key}"')
        if kind == 'num':
            if not re.fullmatch( r'[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)'
                                 r'([eE][+-]?[0-9]+)?', val ):
                raise bad(f'value "{val}" is not a number')
            val = float( val )
        elif kind == 'bool':
            if val.lower() not in _bool_values:
                raise bad(f'value "{val}" is not yes or no')
            val = _bool_values[val.lower()]
        else:
            val = val.lower()
            if allowed and val not in allowed:
                raise bad(f'"{key}" must be one of: '+', '.join(allowed))
        self.__key, self.__op, self.__kind, self.__val = key, op, kind, val

    def __call__( self, d ):
        v = d.get( self.__key )
        if v is None:
            return False
        if self.__kind == 'str':
            v = v.lower()
        elif self.__kind == 'bool':
            v = bool( v )
        return _condition_ops[self.__op]( v, self.__val )
