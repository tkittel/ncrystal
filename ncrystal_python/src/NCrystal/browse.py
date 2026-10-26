
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

"""Utilities for browsing and searching the data available to NCrystal (e.g.
files in the standard data library or in the current directory, in-memory
data, or data created on-demand like "solid::B4C/2.52gcm3"), optionally based
on the physics content of the materials they describe.

The main tool is the DataBrowser class. For example:

  import NCrystal.browse as nb
  b = nb.DataBrowser( factory = 'stdlib' )
  b.match('Al').dump()                    #names and short descriptions
  sel = b.search('phonopy').where('crystal=yes','absxs<0.1')
  for e in sel:
      print( e.fullkey, e.props.density )
  js = sel.to_json()

Physics properties (see physics_props_doc()) are loaded when first needed,
in parallel by the C++ layer.
"""

__all__ = [ 'DataBrowser', 'DataEntry', 'PhysicsProps', 'physics_props_doc' ]

from ._conditions import _Condition, _docs_for


def physics_props_doc( conditions = False ):
    """List of (name,description) of the PhysicsProps properties (or only
    those usable in conditions)."""
    return _docs_for( _propdocs, _propkinds, conditions )

class PhysicsProps:
    """Physics properties of a loaded material as attributes (e.g.
    props.absxs, see physics_props_doc())."""

    def __init__( self, info ):
        #info: the "info" dict from the browsedb query.
        from ._fmt import format_chemform
        from .atomdata import elementZToName
        compos = tuple( ( Z, tuple( (A,fr) for A,fr in isotopes ) )
                        for Z, isotopes in info['composition'] )
        atomfracs = {}
        for Z, isotopes in compos:
            en = elementZToName(Z)
            for A, fr in isotopes:
                lbl = ( en if A == 0 else
                        { (1,2):'D', (1,3):'T' }.get( (Z,A), f'{en}{A}' ) )
                atomfracs[lbl] = atomfracs.get(lbl,0.0) + fr
        elements = frozenset( elementZToName(Z) for Z,_ in compos )
        self.__d = {
            'elements': elements,
            'nelements': len( elements ),
            'formula': format_chemform( sorted( atomfracs.items() ) ),
            'mass': info['mass'],
            'cohsl': info['cohsl'],
            'cohxs': info['xsect_coh'],
            'incohxs': info['xsect_incoh'],
            'scatxs': info['xsect_free'],
            'absxs': info['xsect_absorption'],
            'density': info['density'],
            'numdens': info['numberdensity'],
            'state': info['stateofmatter'].lower(),
            'crystal': info['crystalline'],
            'crystalsystem': _crystal_system( info['spacegroup'] ),
            'sg': info['spacegroup'],
            'natoms': info['natoms_unitcell'],
            'braggthreshold': info['braggthreshold'],
            **{ t: ( t in info['dyninfo_types'] ) for t in _dyninfo_types },
            'nphases': info['nphases'],
        }
        self.__compos = compos

    def _as_dict( self, json_compatible = False ):
        #The properties as a (new) dictionary (with sets as sorted lists if
        #json_compatible=True):
        if not json_compatible:
            return dict( self.__d )
        return { k: ( sorted(v) if isinstance(v,frozenset) else v )
                 for k,v in self.__d.items() }

    @property
    def composition( self ):
        """Composition as ((Z,((A,fraction),...)),...), with A=0 for natural
        elements."""
        return self.__compos

    def __getattr__( self, name ):
        d = self.__dict__.get('_PhysicsProps__d')
        if d is not None and name in d:
            return d[name]
        raise AttributeError(f'PhysicsProps has no attribute "{name}"')

    def __str__( self ):
        s = ', '.join( f'{k}={_fmt_prop(v)}' for k,v in self.__d.items() )
        return f'PhysicsProps({s})'

    def __repr__( self ):
        return str(self)

class DataEntry:
    """A data entry available to NCrystal (e.g. a file), as provided by
    DataBrowser objects."""

    def __init__( self, data, _source = None ):
        #data: a dict from the browsedb query, _source: for loading physics.
        self._d = data
        c = data.get('comments')
        self.__comments = tuple(c) if c is not None else None
        self.__props = ( PhysicsProps( data['info'] ) if 'info' in data
                         else None )
        loaded = 'info' in data or 'error' in data
        self.__source = None if loaded else _source

    def _ensure_loaded( self ):
        src, self.__source = self.__source, None
        if src is not None:
            e = src.with_physics( [ self ] )[0]
            self._d, self.__props = e._d, e.__props

    @property
    def name( self ):
        """Name of the entry (e.g. a file name)."""
        return self._d['name']

    @property
    def fullkey( self ):
        """The name with its factory, as "<factory>::<name>"."""
        return self._d['fullkey']

    @property
    def factory( self ):
        """Name of the factory providing the entry."""
        return self._d['factory']

    @property
    def comments( self ):
        """Header comment lines of NCMAT data (without "#"), or None if not
        NCMAT data."""
        return self.__comments

    @property
    def props( self ):
        """PhysicsProps of the material (loaded when needed), or None if
        loading failed (see .error)."""
        self._ensure_loaded()
        return self.__props

    @property
    def error( self ):
        """Error message if loading the material failed, otherwise None."""
        self._ensure_loaded()
        return self._d.get('error')

    def load( self, cfg_params = '' ):
        """Load the material, optionally with extra cfg parameters (e.g.
        ";temp=20K")."""
        from .core import load
        return load( self.fullkey + cfg_params )

    def _as_dict( self ):
        #JSON-compatible dictionary with all information about the entry:
        p, d = self.props, self._d
        return { 'name': self.name, 'fullkey': self.fullkey,
                 'factory': self.factory, 'source': d['source'],
                 'priority': d['priority'], 'hidden': d['hidden'],
                 'datatype': d.get('datatype'), 'path': d.get('path'),
                 'description': _short_descr( self.__comments ),
                 'comments': ( list(self.__comments)
                               if self.__comments is not None else None ),
                 'props': ( p._as_dict( json_compatible = True )
                            if p is not None else None ),
                 'error': self.error }

    def __str__( self ):
        return f'DataEntry({self.fullkey})'

    def __repr__( self ):
        return str(self)

class DataBrowser:
    """Browse and search the data available to NCrystal. The selection
    methods return new DataBrowser objects, and iteration gives the selected
    entries (DataEntry objects). Physics properties are loaded when first
    needed."""

    def __init__( self, *patterns, factory = None, search = (), where = (),
                  _progress = None ):
        """Browse all data, or that of a factory (files relative to the
        current directory only with factory="relpath"), applying any
        patterns, search words, and where conditions."""
        src = _Source( factory, progress = _progress )
        self._init( src, [ e for e in src.entries()
                           if factory == 'relpath' or e.factory != 'relpath' ] )
        b = self.match( *patterns ).search( *_aslist(search) )
        self._init( src, b.where( *_aslist(where) )._entries )

    def _init( self, src, entries ):
        self._src = src
        self._entries = tuple( entries )

    def _derived( self, entries ):
        b = object.__new__( DataBrowser )
        b._init( self._src, entries )
        return b

    def _filter( self, fct ):
        return self._derived( [ e for e in self._entries if fct(e) ] )

    def _with_physics( self ):
        return self._derived( self._src.with_physics( self._entries ) )

    def match( self, *patterns ):
        """Select entries whose names contain any of the patterns
        (case-insensitive). Patterns with "::" are matched against the full
        key (e.g. "stdlib::Al")."""
        if not patterns:
            return self
        return self._filter( lambda e : any( _name_matches( e, p )
                                             for p in patterns ) )

    def search( self, *words ):
        """Select entries with all the words in their name or NCMAT header
        comments (case-insensitive)."""
        lw = [ w.lower() for w in words ]
        def ok( e ):
            txt = '\n'.join( [ e.name ] + list( e.comments or [] ) ).lower()
            return all( w in txt for w in lw )
        return self._filter( ok ) if lw else self

    def with_elements( self, *labels ):
        """Select materials containing any of the elements or isotopes (e.g.
        "B", "B10", or "D"). Natural elements and isotopes never match each
        other (e.g. "B" does not select materials with only B10)."""
        if not labels:
            return self
        za = { _atom_za( lbl ) for lbl in labels }
        def ok( e ):
            return e.props is not None and any(
                (Z,A) in za for Z,isotopes in e.props.composition
                for A,_ in isotopes )
        return self._with_physics()._filter( ok )

    def where( self, *conditions ):
        """Select materials fulfilling all conditions like "sg=225" or
        "absxs>100" (see physics_props_doc(conditions=True)). Text comparisons
        ignore case, and missing values (e.g. sg of a liquid) never match."""
        conds = [ _Condition( c, _propnames(), _propkinds ) for c in conditions ]
        if not conds:
            return self
        return self._with_physics()._filter(
            lambda e : e.props is not None
            and all( c( e.props._as_dict() ) for c in conds ) )

    def names( self ):
        """List of the full keys of the entries."""
        return [ e.fullkey for e in self._entries ]

    def __len__( self ):
        return len( self._entries )

    def __bool__( self ):
        return bool( self._entries )

    def __iter__( self ):
        return iter( self._entries )

    def __repr__( self ):
        n = len(self)
        return f'DataBrowser({n} entr{"y" if n==1 else "ies"})'

    def dump( self ):
        """Print the full key and short description of each entry."""
        from ._common import print
        w = min( 40, max( ( len(e.fullkey) for e in self ), default = 0 ) )
        for e in self:
            print( f'{e.fullkey.ljust(w)}  {_short_descr(e.comments)}'.rstrip() )

    def to_json( self ):
        """JSON string with all information about the entries."""
        import json
        return json.dumps( [ e._as_dict() for e in self._with_physics() ],
                           indent = 1 ) + '\n'

###############################################################################
# Implementation details:

_progress_chunk = 16

_crystal_systems = ( 'triclinic', 'monoclinic', 'orthorhombic', 'tetragonal',
                     'trigonal', 'hexagonal', 'cubic' )

_dyninfo_types = ( 'vdos', 'vdosdebye', 'scatknl', 'freegas', 'sterile' )

_propdocs = [
    ('elements', 'elements, e.g. Al,O'),
    ('nelements', 'number of elements'),
    ('formula', 'chemical formula, e.g. Al2O3'),
    ('mass', 'average atomic mass [amu]'),
    ('cohsl', 'average coherent scattering length [fm]'),
    ('cohxs', 'coherent scattering cross section [barn/atom]'),
    ('incohxs', 'incoherent scattering cross section [barn/atom]'),
    ('scatxs', 'free scattering cross section [barn/atom]'),
    ('absxs', 'absorption cross section at 2200m/s [barn/atom]'),
    ('density', 'density [g/cm3]'),
    ('numdens', 'number density [atoms/Aa3]'),
    ('state', 'solid, liquid, gas, or unknown'),
    ('crystal', 'yes or no'),
    ('crystalsystem', ', '.join(_crystal_systems[:-1]) + ', or cubic'),
    ('sg', 'space group number'),
    ('natoms', 'number of atoms in the unit cell'),
    ('braggthreshold', 'Bragg threshold [Aa]'),
    *( ( t, 'whether some atoms have this type of dynamics (yes or no)' )
       for t in _dyninfo_types ),
    ('nphases', 'number of phases'),
]

#Kinds of properties, which are numbers unless listed here (a tuple lists the
#allowed values of a text property, and "nocond" marks properties which can
#not be used in conditions):
_propkinds = { 'elements': 'nocond', 'formula': 'str',
               'state': ('solid','liquid','gas','unknown'),
               'crystalsystem': _crystal_systems, 'crystal': 'bool',
               **{ t: 'bool' for t in _dyninfo_types } }

def _propnames():
    return [ n for n,d in _propdocs ]

def _crystal_system( sg ):
    if not sg:
        return None
    #First space group number of each crystal system:
    for sgmin, name in zip( (195,168,143,75,16,3,1),
                            reversed(_crystal_systems) ):
        if sg >= sgmin:
            return name

def _q( *args ):
    from .misc import evaluate_query
    return evaluate_query( ['util'] + list(args) )

def _aslist( x ):
    return [ x ] if isinstance( x, str ) or callable( x ) else list( x )

class _nullctx:
    def __enter__( self ):
        pass
    def __exit__( self, *a ):
        pass

def _sortkey( d ):
    #Sort key for dicts from browsedb query:
    p = d['priority']
    rank = p if isinstance( p, int ) else ( -1 if p == 'OnlyOnExplicitRequest'
                                            else -2 )
    return ( -rank, d['factory'], d['source'], d['name'] )

def _query( facts, counts, *, physics, nthreads, quiet, progress ):
    #Query given factories (with or without physics), returning sorted list
    #of dicts:
    ntotal = sum( counts[f] for f in facts )
    extra = ( ( f'nthreads={nthreads}', ) if physics else ('cheap',) )
    data, ndone = [], 0
    from ._msg import _suppress_msgs_ctx
    msgctx = _suppress_msgs_ctx() if ( physics and quiet ) else _nullctx()
    with msgctx:
        for f in facts:
            #Smaller chunks if progress reporting is needed:
            nchunks = ( 1 if progress is None
                        else max( 1, -(-counts[f]//_progress_chunk) ) )
            for i in range( nchunks ):
                if progress is not None:
                    progress( ndone, ntotal )
                d = _q( 'browsedb', f, str(i), str(nchunks), *extra )
                data += d
                ndone += len(d)
    if progress is not None:
        progress( ndone, ntotal )
    data.sort( key = _sortkey )
    return data

class _Source:
    #Data shared by a DataBrowser and all browsers derived from it, in
    #particular caching physics properties loaded on demand (per factory).
    def __init__( self, factory, *, progress ):
        self.counts = _q('browsedb')
        if factory is not None and factory not in self.counts:
            from .exceptions import NCBadInput
            raise NCBadInput(f'Unknown TextData factory: "{factory}"')
        self.facts = [ f for f,n in self.counts.items()
                       if n and ( factory is None or f == factory ) ]
        self.progress = progress
        self.cheap = None
        self.loaded = {}
        self.loaded_facts = set()

    def entries( self ):
        if self.cheap is None:
            self.cheap = [ DataEntry( d, self ) for d in
                           _query( self.facts, self.counts, physics = False,
                                   nthreads = 1, quiet = True,
                                   progress = None ) ]
        return self.cheap

    def with_physics( self, entries ):
        needed = []
        for e in entries:
            if e.factory not in self.loaded_facts and e.factory not in needed:
                needed.append( e.factory )
        if needed:
            for d in _query( needed, self.counts, physics = True,
                             nthreads = 'auto', quiet = True,
                             progress = self.progress ):
                self.loaded[d['fullkey']] = DataEntry(d)
            self.loaded_facts.update( needed )
        return [ self.loaded.get( e.fullkey, e ) for e in entries ]

def _name_matches( entry, pattern ):
    s = entry.fullkey if '::' in pattern else entry.name
    return pattern.lower() in s.lower()

def _short_descr( comments ):
    #First paragraph of the header comments, as a single line.
    words = []
    for line in ( comments or [] ):
        if not line.strip():
            if words:
                break
            continue
        #Skip ascii-art rulers like "----" or "=====":
        if not any( c.isalnum() for c in line ):
            if words:
                break
            continue
        words += line.split()
    return ' '.join(words)

def _fmt_prop( v ):
    if v is None:
        return 'None'
    if isinstance( v, bool ):
        return 'yes' if v else 'no'
    if isinstance( v, frozenset ):
        return '{' + ','.join( repr(e) for e in sorted(v) ) + '}'
    if isinstance( v, float ):
        return f'{v:g}'
    return repr(v)

def _atom_labels():
    #Labels of all atom database entries, and the D and T aliases:
    from .atomdata import _atomdb_data
    return { d['label'] for d in _atomdb_data() } | { 'D', 'T' }

def _atom_za( label ):
    #(Z,A) of an element or isotope label (A=0 for natural elements):
    if label not in _atom_labels():
        from .exceptions import NCBadInput
        raise NCBadInput(f'Unknown element or isotope: "{label}"')
    from .atomdata import atomDB
    a = atomDB( label )
    return ( a.Z(), a.A() )
