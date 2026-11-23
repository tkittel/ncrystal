
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
collecting the atoms of the sites (elements, positions, a single displacement
magnitude from Uiso and Ueq), merging those at the same positions, and finding
and grouping the final sites (refined with spglib if requested). Symmetry and
tensor operations are done by the sgsym queries.

"""

__all__ = []

from . import _common as _nc_common
from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core
from ._cifreadimpl import CIFContent

_rawlist = _nc_sgsymimpl._rawlist


def _query( *args ):
    from .misc import evaluate_query
    return evaluate_query( list( args ) )

def _cellarg( cell ):
    #Cell (dict of printed values or numbers) as argument of the sgsym queries:
    return ','.join( str( cell[k] ) for k in ( 'a', 'b', 'c', 'alpha',
                                               'beta', 'gamma' ) )

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

def _uij_dict( u ):
    return dict( zip( _nc_sgsymimpl._uij_keys, u ) ) if u else None

def _expand( setting, cell, entries ):
    #The sgsym positions query for entries (label, x, y, z, occupancy, uij),
    #with lenient tolerances (values as printed strings or numbers):
    return _query( 'sgsym', 'positions', 'strict', setting, 'lenient',
                   'RAWDATA', _cellarg( cell ),
                   *_rawlist( ','.join( [ str( e ) for e in en[:5] ]
                                        + ( list( en[5] ) if en[5] else [] ) )
                              for en in entries ) )['sites']

def _site_numbers( sites ):
    #The printed numbers of the sites (string -> (value, rounding error) or
    #None), from one cif query:
    strs = sorted( { s[k] for s in sites
                     for k in ( 'x', 'y', 'z', 'occ', 'uiso', 'biso' )
                     if s[k] is not None }
                   | { v for s in sites for v in ( s['uij'] or [] ) } )
    return dict( zip( strs, ( tuple( v ) if v else None for v in
                              _query( 'cif', 'numbers', *_rawlist( strs ) ) ) ) )

def _element_markers( sites, formula ):
    #Element (or isotope, e.g. "H2" for deuterium) of each CIF site (dicts,
    #cf. read_cif), from its type symbol, or else from its label checked
    #against the chemical formula, or without a formula only from labels like
    #"Al1" (with one warning; cf. elementFromCIFLabel):
    from .atomdata import elementZToName
    def marker( e ):
        return elementZToName( e[0] ) + ( str( e[1] ) if e[1] else '' )
    elems = _query( 'cif', 'elements', 'typesymbol',
                    *_rawlist( s['type_symbol'] or '?' for s in sites ) )
    missing = [ i for i, e in enumerate( elems ) if e is None ]
    if missing:
        fromlbl = _query( 'cif', 'elements', 'label', 'RAWDATA', formula or '',
                          *_rawlist( sites[i]['label'] for i in missing ) )
        examples = {}
        for i, e in zip( missing, fromlbl ):
            lbl = sites[i]['label']
            if e is None:
                raise _nc_core.NCBadInput(
                    f'CIF site "{lbl}" has no usable _atom_site_type_symbol,'
                    + ( ' and its label does not unambiguously start with an'
                        f' element of the chemical formula ("{formula}")'
                        if formula else
                        ' and without a _chemical_formula_sum its label must'
                        ' be an element name followed by nothing but digits'
                        ' and parentheses (like "Al", "Fe2" or "O(1)")' ) )
            elems[i] = e
            examples.setdefault( marker( e ), lbl )
        ex = [ f'{m} from {lbl}' for m, lbl in examples.items() ]
        _nc_common.warn( 'Some CIF sites have no _atom_site_type_symbol and'
                         ' their element content was guessed from their label'
                         f' ({", ".join( ex[:4] )}{", ..." if len(ex) > 4 else ""})' )
    return [ marker( e ) for e in elems ]

def _resolve_uiso_ueq( lbl, uiso, ueq, *, tols, permissive ):
    #Uiso of a CIF site with a tensor, from its Uiso and the Ueq of its tensor
    #(both as (value, rounding error), uiso None if not given). Differences
    #explained by the rounding errors of the printed values plus tol_uij are
    #accepted:
    if not uiso:
        return ueq[0]
    allow = uiso[1] + ueq[1] + tols['uij'] * max( ueq[0], uiso[0] )
    if abs( ueq[0] - uiso[0] ) > allow:
        msg = ( f'CIF site "{lbl}" has Uiso={uiso[0]:g} which'
                f' differs from the Ueq={ueq[0]:g} of its anisotropic'
                ' displacements by more than the rounding errors of'
                ' the printed values explain' )
        if not permissive:
            raise _nc_core.NCBadInput( msg + ' (try permissive mode'
                                       ' to use the Uiso value)' )
        _nc_common.warn( msg + ' (using the Uiso value due to'
                         ' permissive mode)' )
    return uiso[0]

def collect_atoms( cif_content, setting, cell, *, tols, permissive ):
    #Returns the sites with positive occupancy of cif_content, expanded in the
    #space group setting and cell (dict of printed values), as dicts with all
    #their 'positions' (the first the symmetrised site), the 'element_marker'
    #(e.g. "Al" or "H2"), 'occupancy', 'cif_label', 'split' (cf. the sgsym
    #positions query), and 'uiso' (with the tensors 'aniso' and those at each
    #position 'uijs' rescaled to it).
    assert isinstance( cif_content, CIFContent )
    import math
    nums = _site_numbers( cif_content.sites )
    biso_to_uiso = 1.0 / ( 8.0 * math.pi * math.pi )
    def num( s ):
        return None if s is None else nums[s]
    sites, entries, uisos = [], [], []
    for s in cif_content.sites:
        lbl = s['label']
        occ = ( 1.0, 0.0 ) if s['occ'] is None else num( s['occ'] )
        if occ is None:
            raise _nc_core.NCBadInput( f'CIF site "{lbl}" has an invalid'
                                       f' occupancy ("{s["occ"]}")' )
        if not occ[0] > 0.0:
            continue
        if occ[0] > 1.0 + 1e-6:
            raise _nc_core.NCBadInput( f'CIF site "{lbl}" has an occupancy'
                                       f' above 1 ({occ[0]:g})' )
        xyz = [ num( s[k] ) for k in ( 'x', 'y', 'z' ) ]
        if any( v is None for v in xyz ):
            raise _nc_core.NCBadInput( f'CIF site "{lbl}" has missing or'
                                       ' invalid fractional coordinates' )
        uij = s['uij']
        if uij and any( num( v ) is None for v in uij ):
            raise _nc_core.NCBadInput( f'CIF site "{lbl}" has invalid'
                                       ' anisotropic displacement values' )
        if uij and all( num( v )[0] == 0.0 for v in uij ):
            uij = None#(all zero: a placeholder for no tensor)
        u = num( s['uiso'] )
        if not ( u and u[0] > 0.0 ):
            b = num( s['biso'] )
            u = ( ( b[0] * biso_to_uiso, b[1] * biso_to_uiso )
                  if b and b[0] > 0.0 else None )
        sites.append( s )
        uisos.append( u )
        #(CIF labels in the query, so they appear in its error messages, made
        #unique and without commas):
        qlbl = ( lbl or '?' ).replace( ',', ';' )
        while any( e[0] == qlbl for e in entries ):
            qlbl += "'"
        entries.append( ( qlbl, s['x'], s['y'], s['z'],
                          repr( min( occ[0], 1.0 ) ), uij ) )
    if not sites:
        raise _nc_core.NCBadInput('No atom sites found in CIF data.')
    markers = _element_markers( sites, cif_content.formula )
    res = []
    for s, en, u, marker, q in zip( sites, entries, uisos, markers,
                                    _expand( setting, cell, entries ) ):
        uiso, aniso, uijs = ( u[0] if u else None ), None, None
        if q['uij']:
            ueq = tuple( q['ueq'] )
            uiso = _resolve_uiso_ueq( s['label'], u, ueq, tols = tols,
                                      permissive = permissive )
            aniso = _nc_sgsymimpl._rescale_ustd( _uij_dict( q['uij'] ), ueq[0],
                                                uiso )
            uijs = [ _nc_sgsymimpl._rescale_ustd( _uij_dict( v ), ueq[0], uiso )
                     for v in q['uijs'] ]
        res.append( { 'positions': [ tuple( p ) for p in q['positions'] ],
                      'element_marker': marker,
                      'occupancy': float( en[4] ),
                      'cif_label': s['label'],
                      'split': q['split'],
                      'uiso': uiso,
                      'aniso': aniso,
                      'uijs': uijs,
                      '_entry': en } )
    return res

def merge_atoms( collected_atoms, expected_tot_atom_in_orig_cell, setting,
                 cell ):
    #Merge the collected_atoms (cf. collect_atoms) which occupy the same
    #positions into atoms with a composition (fractions of the total occupancy)
    #and occupancy, checking the total number of positions against
    #expected_tot_atom_in_orig_cell (if not None). Merged atoms are expanded
    #again in the setting and cell with their total occupancy (cf. the split
    #site rule of the sgsym positions query). Returns the merged atoms.
    from math import fsum as _math_fsum
    list_ro = tuple
    pos_tolerance = 0.0001#NB: Best if value matches the one in NCInfoBuilder.cc
    def has_same_sites( atom1, atom2 ):
        l1,l2 = atom1['positions'],atom2['positions']
        return len(l1) == len(l2) and any( ( _nc_sgsymimpl._unit_cell_point_dist(c2,l1[0])<pos_tolerance) for c2 in l2 )

    atoms_to_process, final_atoms, reexpand = collected_atoms, [], []
    while atoms_to_process:
        atom1 = atoms_to_process[0]
        atoms_with_same_pos = [ atom1 ]
        atoms_with_different_pos = []
        for atom2 in atoms_to_process[1:]:
            if has_same_sites(atom1,atom2):
                if atom1['uiso'] != atom2['uiso'] or atom1['aniso'] != atom2['aniso']:
                    raise _nc_core.NCBadInput('atoms labelled "{}" and "{}" have same cell '.format(atom1['cif_label'],atom2['cif_label'])
                                              +'positions but different uiso information. This is not supported.')
                atoms_with_same_pos.append( atom2 )
            else:
                atoms_with_different_pos.append( atom2 )

        cif_labels = [ a['cif_label'] for a in atoms_with_same_pos ]
        cif_labels_str = '"{}"'.format('", "'.join(cif_labels))

        composition = [ ( a['occupancy'], a['element_marker'] ) for a in atoms_with_same_pos ]
        occupancy = _math_fsum( fr for fr, elem in composition )
        assert occupancy > 0.0
        if occupancy != 1.0:
            composition = [ (fr/occupancy, elem) for fr, elem in composition ]
            if occupancy > (1.0 + 1e-6) :
                m=f'Error: Too high total occupancy ({occupancy} which is >1) for atoms with CIF labels: {cif_labels_str}'
                if abs(int(occupancy+0.5)-occupancy)<1e-10:
                    m+='. A possible cause could be multiple entries in the CIF file using atomic positions that are actually symmetrically equivalent!'
                raise _nc_core.NCBadInput(m)
            occupancy = min( 1.0, occupancy )

        a0 = atoms_with_same_pos[0]
        final_atoms.append( { 'cif_labels': list_ro( cif_labels ),
                              'positions': list( a0['positions'] ),
                              'uijs': a0['uijs'],
                              'composition': list_ro( sorted( composition ) ),
                              'occupancy': occupancy,
                              'split': a0['split'],
                              'uiso': a0['uiso'],
                              'aniso': a0['aniso'] } )
        if len( atoms_with_same_pos ) > 1:
            reexpand.append( ( final_atoms[-1], a0['_entry'][:4]
                               + ( repr( occupancy ), None ) ) )
        atoms_to_process = atoms_with_different_pos

    if reexpand:
        for ( a, _ ), q in zip( reexpand, _expand( setting, cell,
                                                   [ e for _, e in reexpand ] ) ):
            a['split'] = q['split']

    if expected_tot_atom_in_orig_cell:
        n_actual = sum( len(a['positions']) for a in final_atoms )
        def rd( x,y ):
            return abs(x-y)/(max(1e-300,abs(x)+abs(y)))
        if rd( n_actual, expected_tot_atom_in_orig_cell ) > 0.01:
            raise _nc_core.NCBadInput(f'Expanded number of atoms per cell ({n_actual}) is different '
                                      'from what is stated explicitly in the CIF data'
                                      f' ({expected_tot_atom_in_orig_cell}). Wrong spacegroup (or'
                                      ' choice/setting improperly specified)?')
    return final_atoms

def check_split_sites( atoms, *, permissive ):
    #Split sites (cf. the sgsym positions query) are an error, unless
    #permissive (then a warning).
    for a in atoms:
        sp = a['split']
        if not sp:
            continue
        labels = '", "'.join( a['cif_labels'] )
        msg = ( f'Split site (CIF label "{labels}"): {sp["positions"]}'
                ' symmetry-equivalent positions'
                f' {sp["distance"]:.2g} Aa apart (occupancy'
                f' {a["occupancy"]:g}). NCrystal currently treats all'
                ' vacancies as uncorrelated, so incoherent elastic scattering'
                ' is overestimated for this material' )
        if not permissive:
            raise _nc_core.NCBadInput( msg + ' (try permissive mode to load'
                                       ' it anyway)' )
        _nc_common.warn( msg + ' (accepted due to permissive mode)' )

def _impl_find_sites( cellsg, cif_atoms, setting, *, tols, refine, permissive,
                      quiet ):
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
        _nc_common.warn(w)
    if not quiet:
        for m in r['msgs']:
            _nc_common.print(m)
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
            _nc_common.warn( msg + ' (discarded due to permissive mode)' )
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
        _nc_common.warn( msg + ' (discarded due to permissive mode)' )
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
