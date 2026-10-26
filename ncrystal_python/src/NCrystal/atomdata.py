
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

Utilities for accessing NCrystal's database of atomic data, with information
about atomic masses, scattering lengths, etc. Also contains a few other related
utilities, like a list of all element names, and the AtomDBBrowser class for
browsing and searching the database.

"""

__atomdb={}
def atomDB(Z,A=None,throwOnErrors=True):
    """Access internal database with data for isotopes and natural elements.

    If A is provided, both A and Z must be integers, thus defining a specific isotope.

    If Z is an integer and A is 0 or None, the corresponding natural element is provided.

    Finally, the function can be called with a string identifying either natural
    elements or isotopes: atomDB("Al"), atomDB("He3"), ...

    In all cases, in case of errors or missing entries in the database, either
    an NCBadInput exception is thrown (throwOnErrors==True) or None is
    returned (when throwOnErrors==False).
    """
    import numbers
    if isinstance(Z,numbers.Integral):
        Z=int(Z)
        key=(Z,int(A or 0))
        strkey=False
    else:
        assert A is None,"Do not supply two arguments unless the first argument is an integer"
        assert isinstance(Z,str),"The first argument to the function must either be of int or str type"
        key=Z
        strkey=True
    obj=__atomdb.get(key,None)
    if obj:
        return obj
    from ._chooks import _get_raw_cfcts, _str2cstr
    _rawfct = _get_raw_cfcts()
    if strkey:
        rawatomdata=_rawfct['ncrystal_create_atomdata_fromdbstr'](_str2cstr(key))
    else:
        rawatomdata=_rawfct['ncrystal_create_atomdata_fromdb'](*key)
    if not _rawfct['ncrystal_valid'](rawatomdata):
        if not throwOnErrors:
            return None
        if strkey:
            s=f'key="{key}"'
        else:
            if key[1]==0:
                s=f'Z={int(key[0])}'
            else:
                s=f'Z={int(key[0])},A={int(key[1])}'
        from .exceptions import NCBadInput
        raise NCBadInput(f'atomDB: Could not find entry for key ({s})')
    from .core import AtomData
    ad = AtomData(rawatomdata)
    assert ad.isElement()
    Z,A = ad.Z(), (ad.A() if ad.isSingleIsotope() else 0)
    keys=[ (Z,A)]
    if Z==1 and A==2:
        keys+=['H2','D']
    elif Z==1 and A==3:
        keys+=['H3','T']
    else:
        assert ad.isNaturalElement() or ad.isSingleIsotope()
        keys += [ ad.description(False) ]#guaranteed to give just symbol for natelem/singleisotope!
    assert key in keys#Should always be true unless we forgot some keys above
    assert ad.description(False) in keys#Should also be true, given guarantees for AtomData::description(false)
    for k in keys:
        __atomdb[k] = ad
    return ad

def iterateAtomDB(objects=True):
    """Iterate over all entries in the internal database with data for isotopes and
       natural elements. If objects=True, AtomData objects are returned. If
       objects=False, (Z,A) values are returned (A=0 indicates a natural
       element)."""
    from ._chooks import _get_raw_cfcts
    _rawfct = _get_raw_cfcts()
    for z,a in _rawfct['atomdb_getall_za']():
        yield atomDB(z,a) if objects else (int(z),int(a))

__all_element_names=("H", "He", "Li", "Be", "B", "C",
                     "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
                     "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
                     "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
                     "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
                     "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
                     "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
                     "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                     "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
                     "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
                     "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh",
                     "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts",
                     "Og")

__elem2z = [None]
def elementNameToZValue( element_name, allow_isotopes = False ):
    """Return Z value for element name (e.g. element_name="C" returns
    6). Returns None if argument is not an element name. If
    allow_isotopes==True, isotope names like "B10", "D", etc. will be decoded as
    well and return the Z value of the corresponding element.
    """
    if allow_isotopes:
        if element_name in ('D','T'):
            return 1
        while element_name and element_name[-1].isdigit():
            element_name = element_name[:-1]
    if __elem2z[0] is None:
        __elem2z[0] = { e: i+1 for i,e in enumerate(__all_element_names) }
    return __elem2z[0].get(element_name, None )

def elementZToName( Z ):
    """
    Return element name corresponding to given Z value, or None in case of invalid or out of range argument.
    """
    import numbers
    i = int(Z)-1 if isinstance(Z,numbers.Integral) else -1
    return __all_element_names[i] if 0 <= i < len(__all_element_names) else None

def allElementNames():
    """Like knownElementNames, but also returns names of elements for which
    NCrystal has no data. See also knownElementNames()."""
    return __all_element_names

def isElementName( label ):
    """Determine whether label is a an element name like "H" or "He". See also
    isKnownElement(..).

    """
    return ( ( elementNameToZValue(label) is not None )
             if __elem2z[0] is None else ( label in __elem2z[0] ) )

__all_known_element_names=[None,None]#list,set
def knownElementNames():
    """Returns tuple of all element names for which NCrystal has data
    values. The elements are ordered by increasing Z value. See also
    allElementNames()."""
    if __all_known_element_names[0] is None:
        ll=[]
        for a in iterateAtomDB():
            if a.isNaturalElement():
                ll.append( (a.Z(),a.elementName()) )
        ll = tuple( name for z,name in sorted(ll) )
        __all_known_element_names[0] = ll
        __all_known_element_names[1] = frozenset(ll)
    return __all_known_element_names[0]

def isKnownElement( label ):
    """Determine whether or not label is an element name for which NCrystal has
    data values. See also isElementName(..).
    """
    if __all_known_element_names[1] is None:
        knownElementNames()#trigger init
    return label in __all_known_element_names[1]

def atomdb_fields_doc( conditions = False ):
    """List of (name,description) of the AtomDBBrowser fields (or only those
    usable in conditions)."""
    from ._conditions import _docs_for
    return _docs_for( _atomdb_docs, _atomdb_kinds, conditions )

class AtomDBBrowser:
    """Browse NCrystal's database of isotopes and natural elements. The
    selection methods return new AtomDBBrowser objects, and iteration gives
    the entries as dicts (see atomdb_fields_doc())."""

    def __init__( self, *patterns, where = () ):
        """Browse all entries, applying any patterns and where conditions."""
        rows = [ { k: d[k] for k,_ in _atomdb_docs }
                 for d in _atomdb_data() ]
        self._rows = tuple( rows )
        b = self.match( *patterns ).where(
            *( [where] if isinstance( where, str ) else where ) )
        self._rows = b._rows

    def _derived( self, rows ):
        b = object.__new__( AtomDBBrowser )
        b._rows = tuple( rows )
        return b

    def match( self, *patterns ):
        """Select entries matching any of the patterns: an element (e.g. "He",
        also selecting its isotopes) or an isotope (e.g. "He3", with D and T
        meaning H2 and H3). Matching is case-insensitive."""
        if not patterns:
            return self
        #Element names by Z (natural elements have A=0):
        elem = { d['Z']: d['label'].lower()
                 for d in _atomdb_data() if d['A'] == 0 }
        def ok( d, p ):
            p = { 'd': 'h2', 't': 'h3' }.get( p.lower(), p.lower() )
            return p in ( d['label'].lower(), elem.get( d['Z'], '' ) )
        return self._derived( [ d for d in self._rows
                                if any( ok( d, p ) for p in patterns ) ] )

    def where( self, *conditions ):
        """Select entries fulfilling all conditions like "absxs>100" or "A=0"
        (see atomdb_fields_doc(conditions=True))."""
        from ._conditions import _Condition
        names = [ n for n,_ in _atomdb_docs ]
        conds = [ _Condition( c, names, _atomdb_kinds ) for c in conditions ]
        return self._derived( [ d for d in self._rows
                                if all( c( d ) for c in conds ) ] )

    def names( self ):
        """List of the labels of the entries."""
        return [ d['label'] for d in self._rows ]

    def __len__( self ):
        return len( self._rows )

    def __bool__( self ):
        return bool( self._rows )

    def __iter__( self ):
        return ( dict(d) for d in self._rows )

    def __repr__( self ):
        n = len(self)
        return f'AtomDBBrowser({n} entr{"y" if n==1 else "ies"})'

    def table( self ):
        """Text table of the entries, with stored values as stored and derived
        ones (cohxs and scatxs) rounded at their uncertainty."""
        cols = [ n for n,_ in _atomdb_docs if n != 'label' ]
        header = [ 'label' ] + cols
        rows = [ _atomdb_row_strs( d, header ) for d in self._rows ]
        return _fmt_text_table( [ h.upper() for h in header ], rows,
                                leftcols = {0} )

    def atomdb_lines( self, align = True ):
        """Lines in NCMAT @ATOMDB format (e.g. "He3 3.01602932265u 5.74fm 1.6b
        5333b"), with aligned columns if align=True."""
        lines = { (d['Z'],d['A']): d['atomdbline'] for d in _atomdb_data() }
        res = [ lines[(d['Z'],d['A'])] for d in self._rows ]
        if not align or not res:
            return res
        rows = [ ll.split() for ll in res ]
        w = [ max( len( r[k] ) for r in rows ) for k in range(5) ]
        return [ ' '.join( ( v.ljust(w[k]) if k == 0 else v.rjust(w[k]) )
                           for k,v in enumerate(r) ) for r in rows ]

    def to_json( self ):
        """JSON string with all fields of the entries."""
        import json
        return json.dumps( [ dict(d) for d in self._rows ], indent = 1 ) + '\n'

_atomdb_docs = [
    ('Z', 'atomic number'),
    ('A', 'mass number (0 for natural elements)'),
    ('label', 'e.g. Al or He3'),
    ('mass', 'atomic mass [amu]'),
    ('cohsl', 'coherent scattering length [fm]'),
    ('cohxs', 'coherent scattering cross section [barn]'),
    ('incohxs', 'incoherent scattering cross section [barn]'),
    ('scatxs', 'scattering cross section [barn]'),
    ('absxs', 'absorption cross section at 2200m/s [barn]'),
]

_atomdb_kinds = { 'label': 'nocond' }

_atomdb_cache = [None]
def _atomdb_data():
    #Entries of the "atomdb" query, which is only needed once since the atom
    #database is static (NB: callers must not modify the returned dicts):
    if _atomdb_cache[0] is None:
        from .misc import evaluate_query
        _atomdb_cache[0] = tuple( evaluate_query( ['util','atomdb'] ) )
    return _atomdb_cache[0]

def _atomdb_row_strs( d, header ):
    #Text values of an atomdb entry. Values stored in the database (mass,
    #cohsl, incohxs, absxs) are formatted with at most 15 significant digits,
    #which reproduces their stored decimal form. Their uncertainty is taken
    #as half a unit in their last digit, and propagated to the derived
    #cohxs=4*pi*cohsl^2 and scatxs=cohxs+incohxs, which are rounded at the
    #decimal place of the leading digit of their uncertainty (so rounding adds
    #at most half of it).
    stored = { k: _dbval_str( d[k] ) for k in ('mass','cohsl','incohxs','absxs') }
    dcoh = ( d['cohxs'] * 2.0 * _half_ulp( stored['cohsl'] ) / abs( d['cohsl'] )
             if d['cohsl'] else 0.0 )
    dinc = _half_ulp( stored['incohxs'] ) if d['incohxs'] else 0.0
    vals = dict( stored, cohxs = _fmt_with_uncertainty( d['cohxs'], dcoh ),
                 scatxs = _fmt_with_uncertainty( d['scatxs'], dcoh + dinc ) )
    return [ vals[c] if c in vals else str( d[c] ) for c in header ]

def _dbval_str( v ):
    #Stored decimal form in plain notation (e.g. "0.00005", not "5e-05"):
    s = format( v, '.15g' )
    if 'e' not in s:
        return s
    import decimal
    return format( decimal.Decimal( s ), 'f' )

def _half_ulp( s ):
    #Half a unit in the last digit of the plain number string s (e.g. 0.005
    #for "5.74"):
    return 0.5 * 10.0**( -len( s.partition('.')[2] ) )

def _fmt_with_uncertainty( v, dv ):
    #Round v at the decimal place of the leading digit of its uncertainty dv:
    if v == 0.0 or dv <= 0.0:
        return _dbval_str( v )
    import math
    p = math.floor( math.log10( dv ) )
    return f'{v:.{-p}f}' if p < 0 else f'{round( v, -p ):.0f}'

def _fmt_text_table( header, rows, *, leftcols ):
    #NB: Rows are never truncated (no data should be lost).
    widths = [ max( len(r[k]) for r in rows + [header] )
               for k in range(len(header)) ]
    def fmt( r ):
        return '  '.join( ( v.ljust(w) if k in leftcols else v.rjust(w) )
                          for k,(v,w) in enumerate(zip(r,widths)) ).rstrip()
    rule = '  '.join( '-'*w for w in widths )
    return ''.join( ll + '\n' for ll in [ fmt(header), rule ]
                    + [ fmt(r) for r in rows ] )
