
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

def _guess_elem_name( site_label ):
    #Element name from a CIF label like "CaM1" or "OA1" (None if unclear):
    from .atomdata import allElementNames
    sl, traildigit = _nc_common._split_trailing_digit( site_label )
    if traildigit is None:
        return
    #first check element names with 2 chars:
    candidates = [ e for e in allElementNames() if ( sl.startswith(e) and len(e)==2 ) ]
    if len(candidates) > 1:
        return
    #Then those with 1 char if no hit already (also check for 'D' and 'T'):
    if not candidates:
        candidates = [ e for e in allElementNames() if ( sl.startswith(e) and len(e)==1 ) ]
        candidates += [ e for e in ('D','T') if ( sl.startswith(e) and len(e)==1 ) ]
    if len(candidates)!=1:
        return
    c = candidates[0]
    assert site_label.startswith(c)
    leftover = site_label[len(c):].strip()
    while leftover and leftover[0].isupper()  and leftover[0].isalpha():
        leftover = leftover[1:].strip()
    if not leftover or leftover.isdigit():
        return c

def collect_atoms( content, choice, cellsg, *, tols, permissive ):
    """The CIF sites with positive occupancy (those of the symmetry choice if
    already expanded, else those of the content, cf. decide_symmetry and
    read_cif), as dicts with all their positions in the cell ('expanded_coords',
    'rotations'), element, occupancy, label, and a single magnitude of their
    displacements ('uiso', with 'aniso' rescaled to it)."""
    from . import _cifsymimpl as sym
    sgs = _nc_sgsymimpl
    cell = content['cell']
    fractcoord_approx_1angstrom = 1.0 / ( (cell['a']+cell['b']+cell['c'])/3.0 )
    uiso_printed, ueq_printed = content['uiso_printed'], content['ueq_printed']
    p1_sites = choice['p1_sites']
    collected_atoms = []
    for site in ( p1_sites if p1_sites is not None else content['sites'] ):
        if not site['occ'] > 0.0:
            continue
        pos0 = sym.wrap_position( choice, site['fract'] )
        expanded_coords = sym.expand_position( choice, pos0,
                                               fractcoord_approx_1angstrom )

        aniso = site['aniso']
        u_iso = site['u_iso']
        if not u_iso > 0.0:
            u_iso = None

        if aniso:
            #Resolve Uiso versus Ueq here, so each site carries one magnitude
            #(the tensor is rescaled to it). Accepted differences are those
            #explained by the rounding errors of the printed values, plus
            #tol_uij:
            ueq = sgs._cif_uij_ueq( aniso, cellsg )
            lbl = site['label']
            allow = ( uiso_printed.get( lbl, ( None, 0.0 ) )[1]
                      + ueq_printed.get( lbl, ( None, 0.0 ) )[1]
                      + tols['uij'] * max( ueq, u_iso or 0.0 ) )
            if not u_iso:
                u_iso = ueq
            elif abs( ueq - u_iso ) > allow:
                msg = ( f'CIF site "{site["label"]}" has Uiso={u_iso:g} which'
                        f' differs from the Ueq={ueq:g} of its anisotropic'
                        ' displacements by more than the rounding errors of'
                        ' the printed values explain' )
                if not permissive:
                    raise _nc_core.NCBadInput( msg + ' (try permissive mode'
                                               ' to use the Uiso value)' )
                _nc_common.warn( msg + ' (using the Uiso value due to'
                                 ' permissive mode)' )
            aniso = sgs._rescale_uij( aniso, ueq, u_iso )

        site_en = site['element_name']
        if site_en == 'X' or not site_en:
            #Issue seen in codid::9005777, which had no atomic symbols but
            #labels like "CaM1", "CaM2", "OA1".."OA3", "OB1", "OB2", "OC1",
            #"OC3",... Correct if beginning of string starts with exactly one
            #element name, followed by either a digit or an upper case string
            #followed by digits.
            elem_name = _guess_elem_name( site['label'] )
            if elem_name:
                _nc_common.warn(f'Assuming atomic symbol "{elem_name}" based on CIF label "{site["label"]}"')
            else:
                raise _nc_core.NCBadInput(f'Neither Gemmi nor NCrystal\'s custom'
                                          ' code could deduce an atomic symbol'
                                          f' for the CIF label "{site["label"]}"')
            elem_Z = None
        else:
            elem_name, elem_Z = site_en, site['element_Z']

        elem_A = None
        if elem_name == 'D':
            elem_name, elem_Z, elem_A = 'H',1,2
        if elem_name == 'T':
            elem_name, elem_Z, elem_A = 'T',1,3
        elem_marker = elem_name + ( str(elem_A) if elem_A else '' )

        from .atomdata import elementNameToZValue
        z2name = elementNameToZValue( elem_name, allow_isotopes = False )
        if z2name is None:
            raise _nc_core.NCBadInput(f'Not an element name: "{elem_name}"')

        if elem_Z is None:
            elem_Z = z2name
        if elem_Z != z2name:
            raise _nc_core.NCBadInput(f'Wrong atomic number ({elem_Z}) provided'
                                      f' for element "{elem_name}"')

        collected_atoms.append( { 'expanded_coords': expanded_coords,
                                  'rotations':
                                  sym.image_rotations( choice,
                                                       expanded_coords ),
                                  'element_marker': elem_marker,
                                  'occupancy': site['occ'],
                                  'cif_label': site['label'],
                                  'uiso': u_iso,
                                  'aniso': aniso } )

    if not collected_atoms:
        raise _nc_core.NCBadInput('No atom sites found in CIF data.')

    if p1_sites is None:
        #(already checked when expanding p1_sites)
        mult = content['multiplicities']
        sym._check_multiplicities( [ ( a['cif_label'],
                                       mult.get( a['cif_label'] ),
                                       len( a['expanded_coords'] ) )
                                     for a in collected_atoms ],
                                   choice['nops'],
                                   strict = choice['require_mult'] )
    return collected_atoms

def merge_atoms( collected_atoms, expected_tot_atom_in_orig_cell ):
    """Merge the collected atoms (cf. collect_atoms) which occupy the same
    positions into atoms with a composition (fractions of the total occupancy)
    and occupancy, checking the total number of positions against the
    expected one (if not None)."""
    from math import fsum as _math_fsum
    list_ro = tuple
    pos_tolerance = 0.0001#NB: Best if value matches the one in NCInfoBuilder.cc
    def has_same_sites( atom1, atom2 ):
        l1,l2 = atom1['expanded_coords'],atom2['expanded_coords']
        assert len(l1)>0 and len(l2)>0
        return len(l1) == len(l2) and any( ( _nc_sgsymimpl._unit_cell_point_dist(c2,l1[0])<pos_tolerance) for c2 in l2 )

    atoms_to_process, final_atoms = collected_atoms, []
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

        final_atoms.append( { 'cif_labels': list_ro( cif_labels ),
                              'positions': list( atoms_with_same_pos[0]['expanded_coords'] ),
                              'rotations': atoms_with_same_pos[0]['rotations'],
                              'composition': list_ro( sorted( composition ) ),
                              'occupancy': occupancy,
                              'uiso': atoms_with_same_pos[0]['uiso'],
                              'aniso': atoms_with_same_pos[0]['aniso'] } )
        atoms_to_process = atoms_with_different_pos

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
