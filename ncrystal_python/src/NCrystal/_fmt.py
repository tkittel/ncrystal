
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

"""Internal utilities for formatting values."""

def fmt( value ):
    #like our C++ NC::fmt()
    s = repr( float( value ) or 0.0 )
    return s.removesuffix('.0')

def _hexdoubles( values ):
    #The values (flattened) as a text with 16 hex digits for each, the 8
    #bytes of the IEEE 754 double in little-endian order, for passing many
    #numbers to C++ queries exactly and much faster than as decimal texts:
    if hasattr( values, 'ravel' ):
        from ._numpy import _np
        a = _np.ascontiguousarray( values, dtype = '<f8' )
        return a.ravel().tobytes().hex()
    import struct
    v = _flat_values( values )
    return struct.pack( f'<{len(v)}d', *v ).hex()

def _strquery( query ):
    #Query returning a JSON string without escapes, so no need for a JSON
    #decoder (e.g. texts separated by spaces):
    from ._miscimpl import evalquery
    raw = evalquery( query, unpack = False, readonly = False,
                     huge_arrays = False )
    assert raw[0] == '"' and raw[-1] == '"' and '\\' not in raw
    return raw[1:-1]

def _parse_text( t ):
    if '/' in t:
        p, q = t.split( '/' )
        return float( p ) / float( q )
    return float( t )

def fmt_vector_rtol( values, rtol, *, is_grid = False,
                     strip_leading_zero = True, max_denominator = 0 ):
    """Texts for values, rounded to the significant digits for the relative
    tolerance rtol (the fewest n with 10^(1-n) <= rtol, or exact if more than
    15 would be needed, cf. the C++ numFmtDigits), so they read back within
    rtol/2 of the values, and formatting the read back values gives the same
    texts. Zeros are written as "0". With max_denominator, values in (0,1)
    which are fractions p/q with q up to that within the precision are written
    as such (then within rtol). For a grid, more digits are used where needed
    so neighbouring values also keep their differences within rtol. All values
    are formatted in a single call to the C++ fmtRelTolVector."""
    hx = _hexdoubles( values )
    if not hx:
        return []
    return _strquery( [ 'util', 'numfmt', repr( float( rtol ) ),
                        '1' if is_grid else '0',
                        '1' if strip_leading_zero else '0',
                        str( int( max_denominator ) ), hx ] ).split( ' ' )

def _flat_values( values ):
    out = []
    for v in values:
        if hasattr( v, '__len__' ) and not isinstance( v, str ):
            out += _flat_values( v )
        else:
            out.append( v )
    return out

def fmt_value_rtol( x, rtol, *, strip_leading_zero = False,
                    max_denominator = 0 ):
    """Text for x as from fmt_vector_rtol."""
    return fmt_vector_rtol( [ x ], rtol,
                            strip_leading_zero = strip_leading_zero,
                            max_denominator = max_denominator )[0]

def find_fraction( x, tol = 1e-15, max_denom = 1000000 ):
    import fractions
    f = fractions.Fraction.from_float( x )
    f = f.limit_denominator(max_denom)
    if hasattr(f,'as_integer_ratio'):
        #py 3.8+
        a,b = f.as_integer_ratio()
    else:
        a,b = f.numerator,f.denominator
    if abs( a / b - x ) < tol:
        return a,b
    else:
        return None

def _hill_sort( chemform ):
    #takes chemform like [('Al',2'),('O',3)]  and sorts in order of Hill system notation

    #Remap H2/H3 and remove duplicates:
    remap = {'H2':'D','H3':'T'}
    if any(k in remap for k,v in chemform) or len({k for k,v in chemform})!=len(chemform):
        d={}
        for k,v in chemform:
            k = remap.get(k,k)
            d[k] = d.setdefault(k,0) + v
        return _hill_sort( list(d.items()) )

    has_carbon = any( en=='C' for en,c in chemform )
    if not has_carbon:
        return sorted( chemform )
    def hillsortkey( e ):
        #if not has carbon, then all in alphabetical order
        #first carbon, then H/D/T, then in alphabetical order
        if e[0]=='C':
            return ( -999999, e )
        is_hydrogen = ( e[0] in ('H','D','T') or ( e[0].startswith('H') and e[0][1:].isdigit() ) )
        return ( -1, e ) if is_hydrogen else (0, e )
    return sorted(chemform,key=hillsortkey)

def _gcd( *vals ):
    import math
    if len(vals)<2:
        return vals[0] if vals else None
    gcd = math.gcd(vals[0],vals[1])
    for e in vals[2:]:
        gcd = math.gcd( gcd, e )
    return gcd

def format_chemform( chemform, *, allow_rescaling = True ):
    #takes chemform like [('Al',2),('O',3)] and returns nicely formatted string
    #"Al2O3", with no duplicated element names and sorted according to the Hill
    #system of notation. Additionally, integral counts are divided by greatest
    #common divisor, and an attempt is made to scale up non-integral counts
    #where it makes sense (e.g. Al0.5O1.5 becomes AlO3).
    if len(chemform)==1:
        return str(chemform[0][0])
    cf = _hill_sort(chemform)
    def is_near_int(v):
        return abs(int(v)-v)<1e-15
    def get_non_ints(_l):
        return [ v for k,v in _l if not is_near_int(v) ]

    non_ints = get_non_ints( cf )
    def _find_frac(x):
        return find_fraction( x, max_denom = 100 )
    non_int_ffractions = [ _find_frac(x) for x in non_ints ]
    denoms = [ ff[1] for ff in non_int_ffractions if ff is not None ]
    cf_alt = None
    if allow_rescaling and denoms:
        import math
        factor = 1
        while denoms:
            _d = denoms.pop(0)
            _g = math.gcd( factor, _d )
            factor *= ( _d // _g )
        _cf = []
        for k,v in cf:
            if is_near_int(v):
                _cf.append( (k,v*factor ) )
            else:
                ff = _find_frac( v )
                if ff is None:
                    _cf.append( (k,v*factor ) )
                else:
                    assert factor % ff[1] == 0
                    _cf.append( (k,ff[0] * ( factor // ff[1] ) ) )
        if max( v for k,v in _cf ) < 100:
            cf_alt = _cf

    #find and apply gcd of all integers:
    def final_format( the_cf ):
        ll=[ int(v) for k,v in the_cf if is_near_int(v) ]
        gcd = _gcd( *ll ) if ( ll and allow_rescaling ) else 1
        def wrapiso( x ):
            #nb: these curly braces are not great for filenames...:
            return x if not x[-1].isdigit() else f'{{{x}}}'
        the_cf = [ (wrapiso(en),(count//gcd if is_near_int(count) else count/gcd)) for en,count in the_cf ]
        def fmtcount( c ):
            return f'{int(c) if c==int(c) else c:g}'
        return ''.join( ( en if count==1 else en+fmtcount(count) )
                        for en,count in the_cf )

    f1 = final_format(cf)
    f2 = final_format(cf_alt) if cf_alt else f1
    return f1 if len(f1) < len(f2) else f2


def _frexp10(x):
    import math
    exp = math.floor(math.log10(abs(x)))
    return x / 10**exp, exp

def _latex_format(x):
    if len(f'{x:f}')<6:
        return f'{x:f}'
    b,e = _frexp10(x)
    e=f'10^{e}'
    if b==1:
        return e
    return f'{b:g}'+r'\cdot'+e
