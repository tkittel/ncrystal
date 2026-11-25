
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

Internal implementation details for the sites of CIF data in cifutils.py:
finding and grouping the final sites of the atoms (from C++, cf. read_cif),
refined with spglib if requested. Symmetry and tensor operations are done by
the sgsym queries.

"""

__all__ = []

from . import _common as _nc_common
from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core

_rawlist = _nc_sgsymimpl._rawlist


def _query( *args ):
    from .misc import evaluate_query
    return evaluate_query( list( args ) )

def _query_sites( c, labelled_positions ):
    #Sites (one per orbit) of labelled positions (label, position, occupancy)
    #in the cell c (cellsg dict), from the sgsym sites query:
    q = [ 'sgsym', 'sites', 'relaxed',
          _nc_sgsymimpl._sgspec( c['spacegroup']['number'], c['alpha'] ),
          'lenient', 'RAWDATA', _nc_sgsymimpl._cellarg( c ) ]
    q += _rawlist( f'{lbl},{float(x)!r},{float(y)!r},{float(z)!r},'
                   f'{float(o)!r}' for lbl,(x,y,z),o in labelled_positions )
    return _query( *q )

def _orbit_positions( cellsg, setting, pos, occupancy ):
    #All positions of the orbit of pos, in the space group setting:
    x, y, z = ( float(v) for v in pos )
    site, = _query( 'sgsym', 'positions', 'strict', setting, 'lenient',
                    'RAWDATA', _nc_sgsymimpl._cellarg( cellsg ),
                    *_rawlist( [ ( f'X,{x!r},{y!r},{z!r},'
                                   f'{float(occupancy)!r}' ) ] ) )['sites']
    return [ tuple( _nc_sgsymimpl._remap_fract_pos(v) for v in p )
             for p in site['positions'] ]

def _mksite( a, cif_labels, pos, positions, aniso ):
    #Site dict (cf. CIFLoader.sites) of the CIF atom a at position pos, with
    #all positions of its orbit:
    return { 'cif_labels': tuple( sorted( set( cif_labels ) ) ),
             'composition': tuple( a['composition'] ),
             'occupancy': a['occupancy'],
             'position': tuple( _nc_sgsymimpl._remap_fract_pos(v) for v in pos ),
             'multiplicity': len( positions ),
             'uiso': a['uiso'],
             'aniso': aniso,
             '_positions': positions }

def _impl_find_sites( cellsg, cif_atoms, setting, *, tols, refine, permissive,
                      quiet, warnfct ):
    #Returns (cellsg, sites, setting) for the cif_atoms (cf. merge_atoms) in
    #the cell cellsg and space group setting: one site per orbit (cf.
    #_mksite), after refining the structure with spglib if refine is True
    #(which can change the cell and setting).
    from collections import Counter
    sgs = _nc_sgsymimpl
    if not refine:
        return cellsg, [ _mksite( a, a['cif_labels'], a['positions'][0],
                                  a['positions'], a['aniso'] )
                         for a in cif_atoms ], setting

    #All positions, with their tensors:
    entries = []
    for i, a in enumerate( cif_atoms ):
        for j, x in enumerate( a['positions'] ):
            entries.append( sgs.SymSite( key = tuple( a['composition'] ),
                                         pos = x, occu = a['occupancy'],
                                         uiso = a['uiso'],
                                         ustd = ( a['uijs'][j] if a['uijs']
                                                  else None ),
                                         tag = i ) )
    r = sgs._refine_sites( cellsg, entries, tols )
    for w in r['warnings']:
        warnfct(w)
    if not quiet:
        for m in r['msgs']:
            _nc_common.infomsg(m)
    nc = r['cellsg']
    cellsg = { k: nc[k] for k in ('a','b','c','alpha','beta','gamma') }
    cellsg['spacegroup'] = { 'number': nc['spacegroup'], 'hm': r['sgsymb_hm'] }
    q = _query_sites( cellsg, [ ( f't{i}', p, entries[i].occu )
                                for i, p in r['allpos'] ] )
    setting = q['spacegroup_setting']
    def labels_of( indices ):
        return [ lbl for i in indices
                 for lbl in cif_atoms[ entries[i].tag ]['cif_labels'] ]
    types = r['types']
    rsites = r['sites']
    if rsites is not None:
        #Check that the sites are one per orbit of the final setting:
        type2rep = { types[i]: i for i, _ in r['allpos'] }
        if ( Counter( s['label'] for s in q['sites'] )
             != Counter( f't{type2rep[ types[e["index"]] ]}' for e in rsites ) ):
            rsites = None
    if rsites is None:
        if any( e.ustd for e in entries ):
            msg = ( 'The CIF sites could not be followed through the'
                    ' refinement of the structure by spglib, so their'
                    ' anisotropic displacements can not be used' )
            if not permissive:
                raise _nc_core.NCBadInput( msg + ' (try permissive mode to'
                                           ' discard them)' )
            warnfct( msg + ' (discarded due to permissive mode)' )
        #Sites from the query (without tensors), all CIF labels of each type:
        res = []
        for s in q['sites']:
            pos = s['position']
            irep = int( s['label'][1:] )
            members = [ i for i in range( len( entries ) )
                        if types[i] == types[irep] ]
            res.append( _mksite( cif_atoms[ entries[irep].tag ],
                                 labels_of( members ), pos,
                                 _orbit_positions( cellsg, setting, pos,
                                                   cif_atoms[ entries[irep].tag ]['occupancy'] ),
                                 None ) )
        return cellsg, res, setting
    res, dropped = [], []
    for e in rsites:
        a = cif_atoms[ entries[ e['index'] ].tag ]
        if e['aniso_dropped']:
            dropped += labels_of( e['members'] )
        res.append( _mksite( a, labels_of( e['members'] ), e['pos'],
                             _orbit_positions( cellsg, setting, e['pos'],
                                               a['occupancy'] ),
                             e['aniso'] ) )
    if dropped:
        lbls = ', '.join( f'"{x}"' for x in sorted( set( dropped ) ) )
        msg = ( 'The anisotropic displacements of the symmetry equivalent CIF'
                f' sites {lbls} are not symmetry images of each other' )
        if not permissive:
            raise _nc_core.NCBadInput( msg + ' (try permissive mode to'
                                       ' discard them)' )
        warnfct( msg + ' (discarded due to permissive mode)' )
    return cellsg, res, setting

def _impl_group_sites( sites, *, tols, merge_equiv ):
    #Group the sites into atoms (if merge_equiv is True, sites with equal
    #composition and occupancies and uiso values within tols are merged).
    #Returns (atoms, site2atom): the atoms (sorted) and the atom index of each
    #site. NB: removes the '_positions' entries from the sites (as added by
    #_impl_find_sites).
    sgs = _nc_sgsymimpl
    def opt_close( a, b, tol ):
        return ( ( a is None ) == ( b is None )
                 and ( a is None or sgs._relclose( a, b, tol ) ) )
    def same( s1, s2 ):
        return ( s1['composition'] == s2['composition']
                 and ( s1['aniso'] is None ) == ( s2['aniso'] is None )
                 and sgs._relclose( s1['occupancy'], s2['occupancy'],
                                    tols['occu'] )
                 and opt_close( s1['uiso'], s2['uiso'], tols['uij'] ) )
    groups = ( sgs._group_indices( sites, same ) if merge_equiv
               else list( range( len( sites ) ) ) )
    atoms = {}
    for g, s in zip( groups, sites ):
        if g not in atoms:
            atoms[g] = { 'cif_labels': [], 'composition': s['composition'],
                         'equivalent_positions': [],
                         'occupancy': s['occupancy'], 'uiso': s['uiso'] }
        atoms[g]['cif_labels'] += list( s['cif_labels'] )
        atoms[g]['equivalent_positions'] += list( s['_positions'] )
    def finalise( a ):
        return { 'cif_labels': tuple( sorted( set( a['cif_labels'] ) ) ),
                 'composition': a['composition'],
                 'equivalent_positions':
                 tuple( sorted( a['equivalent_positions'] ) ),
                 'occupancy': a['occupancy'], 'uiso': a['uiso'] }
    def _atom_sort( a ):
        #well-defined sort order for reproducible output
        return ( tuple(sorted(a['composition'])),
                 len(a['equivalent_positions']),
                 a['uiso'] if a['uiso'] is not None else -1.0,
                 a['cif_labels'], a['equivalent_positions'] )
    order = sorted( atoms, key = lambda g : _atom_sort( finalise( atoms[g] ) ) )
    newidx = { g: i for i, g in enumerate( order ) }
    for s in sites:
        del s['_positions']
    return ( [ finalise( atoms[g] ) for g in order ],
             [ newidx[g] for g in groups ] )
