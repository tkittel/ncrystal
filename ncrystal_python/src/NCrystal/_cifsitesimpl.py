
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

Internal implementation of finding and grouping the sites of CIF data, for
the cifutils module.

"""

__all__ = []

from . import _common as _nc_common
from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core

_uij_keys = ( 'u11', 'u22', 'u33', 'u12', 'u13', 'u23' )

def _impl_find_sites( cellsg, cif_atoms, *, tols, refine, permissive, quiet ):
    """Sites (one per orbit) of the CIF atoms (dicts with 'positions' and the
    'rotations' of the symmetry operations producing them from the first),
    after refining the structure with spglib if requested. Returns the
    (possibly refined) cellsg, the sites, and the space group setting."""
    from .misc import evaluate_query
    sgs = _nc_sgsymimpl
    def sgspec_of( c ):
        sg = c['spacegroup']['number']
        rhomb = sgs._uses_rhombohedral_axes( sg, c['alpha'] )
        return f'{sg}:r' if rhomb else str(sg)
    def query_sites( c, labelled_positions ):
        q = [ 'util', 'atomsites', 'relaxed', sgspec_of(c) ]
        q += [ f'{lbl},{float(x)!r},{float(y)!r},{float(z)!r}'
               for lbl,(x,y,z) in labelled_positions ]
        try:
            return evaluate_query( q )
        except _nc_core.NCBadInput as e:
            raise _nc_core.NCBadInput( 'The atom positions do not have the'
                                       ' symmetry of space group'
                                       f' {c["spacegroup"]["number"]}: {e}' )
    def orbit_positions( setting, pos ):
        x, y, z = ( float(v) for v in pos )
        q = [ 'util', 'atompositions', 'strict', setting, f'X,{x!r},{y!r},{z!r}' ]
        return [ tuple( sgs._remap_fract_pos(v) for v in p )
                 for _, p in evaluate_query( q ) ]
    def mksite( a, cif_labels, pos, setting, aniso ):
        allpos = orbit_positions( setting, pos )
        return { 'cif_labels': tuple( sorted( set( cif_labels ) ) ),
                 'composition': tuple( a['composition'] ),
                 'occupancy': a['occupancy'],
                 'position': tuple( sgs._remap_fract_pos(v) for v in pos ),
                 'multiplicity': len( allpos ),
                 'uiso': a['uiso'],
                 'aniso': aniso,
                 '_positions': allpos }

    if not refine:
        q = query_sites( cellsg, [ ( f'a{i}', p )
                                   for i, a in enumerate( cif_atoms )
                                   for p in a['positions'] ] )
        setting = q['spacegroup_setting']
        return cellsg, [ mksite( a, a['cif_labels'], a['positions'][0], setting,
                                 a['aniso'] ) for a in cif_atoms ], setting

    #All positions, with tensors rotated accordingly:
    entries = []
    for i, a in enumerate( cif_atoms ):
        uf = sgs._cif_to_frac( a['aniso'], cellsg ) if a['aniso'] else None
        for x, R in zip( a['positions'], a['rotations'] ):
            u = a['aniso']
            if uf is not None and not sgs._np_is_identity( R ):
                u = sgs._frac_to_cif( sgs._transform_frac_tensor( R, uf ),
                                      cellsg )
            entries.append( { 'key': tuple( a['composition'] ), 'pos': x,
                              'occu': a['occupancy'], 'uiso': a['uiso'],
                              'aniso': u, 'src': i } )
    r = sgs._refine_sites( cellsg, entries, tols )
    for w in r['warnings']:
        _nc_common.warn(w)
    if not quiet:
        for m in r['msgs']:
            _nc_common.print(m)
    nc = r['cellsg']
    cellsg = { k: nc[k] for k in ('a','b','c','alpha','beta','gamma') }
    cellsg['spacegroup'] = { 'number': nc['spacegroup'], 'hm': r['sgsymb_hm'] }
    q = query_sites( cellsg, [ ( f't{i}', p ) for i, p in r['allpos'] ] )
    setting = q['spacegroup_setting']
    def labels_of( indices ):
        return [ lbl for i in indices
                 for lbl in cif_atoms[ entries[i]['src'] ]['cif_labels'] ]
    types = r['types']
    def count( ll ):
        res = {}
        for e in ll:
            res[e] = res.get( e, 0 ) + 1
        return res
    rsites = r['sites']
    if rsites is not None:
        #Check that the sites are one per orbit of the final setting:
        type2rep = { types[i]: i for i, _ in r['allpos'] }
        if ( count( lbl for lbl, _, _ in q['sites'] )
             != count( f't{type2rep[ types[e["index"]] ]}' for e in rsites ) ):
            rsites = None
    has_aniso = any( e['aniso'] for e in entries )
    if rsites is None:
        if has_aniso:
            msg = ( 'The CIF sites could not be followed through the'
                    ' refinement of the structure by spglib, so their'
                    ' anisotropic displacements can not be used' )
            if not permissive:
                raise _nc_core.NCBadInput( msg + ' (try permissive mode to'
                                           ' discard them)' )
            _nc_common.warn( msg + ' (discarded due to permissive mode)' )
        res = []
        for lbl, pos, _ in q['sites']:
            irep = int( lbl[1:] )
            members = [ i for i in range( len( entries ) )
                        if types[i] == types[irep] ]
            res.append( mksite( cif_atoms[ entries[irep]['src'] ],
                                labels_of( members ), pos, setting, None ) )
        return cellsg, res, setting
    res, dropped = [], []
    for e in rsites:
        a = cif_atoms[ entries[ e['index'] ]['src'] ]
        if e['aniso_dropped']:
            dropped += labels_of( e['members'] )
        res.append( mksite( a, labels_of( e['members'] ), e['pos'], setting,
                            e['aniso'] ) )
    if dropped:
        msg = ( 'The anisotropic displacements of the symmetry equivalent CIF'
                f' sites {", ".join( sorted( set( dropped ) ) )} are not'
                ' symmetry images of each other' )
        if not permissive:
            raise _nc_core.NCBadInput( msg + ' (try permissive mode to'
                                       ' discard them)' )
        _nc_common.warn( msg + ' (discarded due to permissive mode)' )
    return cellsg, res, setting

def _impl_group_sites( sites, *, tols, merge_equiv ):
    """Group sites into atoms (with merge_equiv, sites with equal composition
    and occupancies and uiso values within the tolerances are merged). Returns
    the atoms (sorted) and the atom index of each site."""
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
    groups = ( sgs._group_by_tol( sites, same ) if merge_equiv
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
