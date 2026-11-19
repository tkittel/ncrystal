
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

def _fmt_ndigits( x, nd, strip_leading_zero ):
    #x (nonzero) with nd significant digits (17: exactly), e.g. "0.25"
    #(".25" with strip_leading_zero):
    if nd >= 17:
        s = repr( x )
    else:
        s = f'{x:.{nd}g}'
        if 'e' in s:
            alt = repr( float( s ) )#same value, sometimes shorter ("300")
            if len( alt.removesuffix('.0') ) < len( s ):
                s = alt
    s = s.removesuffix('.0')
    if strip_leading_zero:
        if s.startswith('0.'):
            s = s[1:]
        elif s.startswith('-0.'):
            s = '-' + s[2:]
    return s

def fmt_vector_rtol( values, rtol, *, is_grid = False,
                     strip_leading_zero = True ):
    """Texts for values, each with the fewest significant digits which read
    back within the relative tolerance rtol of the value (exactly if rtol is
    0). Zeros are written as "0". For a grid, more digits are used where needed
    so neighbouring values also keep their differences within rtol."""
    from ._numpy import _ensure_numpy, _np
    _ensure_numpy()
    x = _np.asarray( values, dtype = float ).ravel()
    nonzero = x != 0.0
    #Estimate the needed number of digits (17: exact), and flag values where
    #one digit less is within float noise of the tolerance:
    nd = _np.full( x.shape, 17, dtype = int )
    try_fewer = _np.zeros( x.shape, dtype = bool )
    if rtol > 0.0 and nonzero.any():
        ax = _np.abs( x[nonzero] )
        e = _np.floor( _np.log10( ax ) )
        e -= ( 10.0**e > ax )
        e += ( 10.0**( e + 1 ) <= ax )
        def rounding_error( d ):
            unit = 10.0**( e - d + 1 )
            return _np.abs( _np.round( ax / unit ) * unit - ax )
        ndnz = _np.full( ax.shape, 17, dtype = int )
        for d in range( 16, 0, -1 ):
            ndnz[ rounding_error( d ) <= rtol * ax ] = d
        nd[nonzero] = ndnz
        try_fewer[nonzero] = ( ( ndnz > 1 ) & ( rounding_error( ndnz - 1 )
                                                <= ( rtol + 1e-15 ) * ax ) )
    def fmt( i, d = None ):
        return ( _fmt_ndigits( float( x[i] ), int( nd[i] if d is None else d ),
                               strip_leading_zero )
                 if x[i] != 0.0 else '0' )
    res = [ fmt( i ) for i in range( len( x ) ) ]
    def within( idx, texts ):
        y = _np.array( [ float( t ) for t in texts ], dtype = float )
        return _np.abs( y - x[idx] ) <= rtol * _np.abs( x[idx] )
    def fix_values( bad_indices ):
        for i in bad_indices:
            nd[i] = min( 17, nd[i] + 1 )
            res[i] = fmt( i )
    #The estimate can be one off near rounding boundaries. Add digits where
    #needed, then remove where possible:
    idx = _np.nonzero( nonzero )[0]
    while len( idx ):
        idx = idx[ ~within( idx, [ res[i] for i in idx ] ) ]
        fix_values( idx )
    idx = _np.nonzero( try_fewer & ( nd > 1 ) & ( nd < 17 ) )[0]
    while len( idx ):
        texts = [ fmt( i, nd[i] - 1 ) for i in idx ]
        ok = within( idx, texts )
        for i, t in zip( idx[ok], _np.asarray( texts, dtype = object )[ok] ):
            nd[i] -= 1
            res[i] = t
        idx = idx[ ok & ( nd[idx] > 1 ) ]
    if is_grid and rtol > 0.0:
        while True:
            y = _np.array( [ float( r ) for r in res ] )
            dx = _np.diff( x )
            bad = _np.nonzero( _np.abs( _np.diff( y ) - dx )
                               > rtol * _np.abs( dx ) )[0]
            if not len( bad ):
                break
            fix_values( set( bad ) | set( bad + 1 ) )
    return res

def fmt_value_rtol( x, rtol, *, strip_leading_zero = False ):
    """Text for x as from fmt_vector_rtol."""
    return fmt_vector_rtol( [ x ], rtol,
                            strip_leading_zero = strip_leading_zero )[0]

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

def prettyFmtValue(x):
    """Recognises common fractions in truncated values like '0.33333333' or
    '0.4166667' and formats them as '1/3' and '5/12'. This not only saves space,
    but can also reinject accuracy in atom positions (which seems to be a clear
    improvement in all cases tested so far). Additionally also performs
    optimisations for storage efficiency, in particular omitting leading 0s
    with no information value (e.g. "0.2341" becomes ".2341").

    """
    assert 0.0 <= x <= 1.0
    if x==0.0:
        return '0'
    if x==1.0:
        return '1'
    if x==0.5:
        return '1/2'
    def stripleading0(s):
        return ( s[1:] if (len(s) > 2 and s.startswith('0.')) else s )
    xfmt = stripleading0(f'{x:.13g}')#.14g leads to irreproducibility issues in
    #our tests but there are sooo many numbers ending with 3333.... or
    #666667.... that we can safely "snap" these to their correct values:
    if xfmt[0]=='.' and xfmt.endswith('33333') and len(xfmt)<18:
        xfmt += '3'*(18-len(xfmt))
    if xfmt[0]=='.' and xfmt.endswith('66667') and len(xfmt)<18:
        xfmt = xfmt[:-1]+'6'*(18-len(xfmt))+'7'
    ff = find_fraction( x, tol = 1e-14, max_denom = 40 )
    if ff is not None and ff[1]!=1:
        assert ff[1] > 0
        v = ff[0] / ff[1]
        v_dbfmt = f'{ff[0]}/{ff[1]}'
        #check if credible that the numbers are identical (must be close AND the str
        #representations must be consistent with truncation):
        expected_at_prec = stripleading0((f'%.{len(xfmt)-1}g')%v)
        if abs(v-x)<1e-7 and expected_at_prec == xfmt:
            return v_dbfmt
        if abs(v-x)<1e-14:
            return v_dbfmt
    #No fraction found:
    if xfmt.isdigit() and float( xfmt ) != x:
        #abort, we don't want 0.99999999999 to print as '1' (which would then be
        #mapped to 0 in a unit cell):
        return stripleading0(f'{x:.19g}')
    return xfmt

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
