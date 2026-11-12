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

# NEEDS: numpy spglib gemmi

# CIFLoader with symmetry equivalent sites given separately (P1, also in a
# primitive cell): symmetry detection with site data, .sites with anisotropic
# displacements, the warning and dropping of tensors which are not symmetry
# images of each other, and the tolerances. Results are compared (via Bragg
# reflections, which depend on the tensors) with the same structures given in
# their space group.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.reprint_escaped_warnings # noqa F401
import numpy as np
import NCrystalDev as NC
import NCrystalDev._sgsymimpl as sgs
from NCrystalDev.cifutils import CIFLoader
from NCrystalDev.misc import evaluate_query

keys = ( 'u11', 'u22', 'u33', 'u12', 'u13', 'u23' )

def cif_text( cell, hm, sites, ops = None, mults = None, header = '' ):
    #sites: (label, element, pos, occu, uiso, aniso dict or None)
    ll = ( [ header ] if header else [] ) + [ 'data_test' ]
    ll += [ f'_cell_length_{k} {cell[k]!r}' for k in 'abc' ]
    ll += [ f'_cell_angle_{k} {cell[k]!r}' for k in ('alpha','beta','gamma') ]
    ll += [ f"_symmetry_space_group_name_H-M '{hm}'" ]
    if ops:
        ll += [ 'loop_', '_symmetry_equiv_pos_as_xyz' ]
        ll += [ f"'{o}'" for o in ops ]
    ll += [ 'loop_', '_atom_site_label', '_atom_site_type_symbol',
            '_atom_site_fract_x', '_atom_site_fract_y', '_atom_site_fract_z',
            '_atom_site_occupancy', '_atom_site_U_iso_or_equiv' ]
    ll += [ '_atom_site_symmetry_multiplicity' ] if mults else []
    for i, ( lbl, el, p, occ, uiso, _ ) in enumerate( sites ):
        ll.append( f'{lbl} {el} {p[0]!r} {p[1]!r} {p[2]!r} {occ!r} {uiso!r}'
                   + ( f' {mults[i]}' if mults else '' ) )
    if any( s[5] for s in sites ):
        ll += [ 'loop_', '_atom_site_aniso_label' ]
        ll += [ f'_atom_site_aniso_U_{k[1:]}' for k in keys ]
        for lbl, _, _, _, _, u in sites:
            if u:
                ll.append( lbl + ' ' + ' '.join( repr( u[k] ) for k in keys ) )
    return '\n'.join( ll ) + '\n'

def expand( setting, cell, asym ):
    #All positions of the asymmetric unit, with rotated tensors:
    ops = evaluate_query( [ 'util', 'sgsym', 'strict', setting ] )['ops']
    res = []
    for lbl, el, p, occ, uiso, u in asym:
        e = sgs._expand_sites( ops, cell, [ { 'pos': p, 'aniso': u } ] )
        for i, x in enumerate( e ):
            res.append( ( f'{lbl}_{i}', el, x['pos'], occ, uiso, x['aniso'] ) )
    return res

def load( txt, **kw ):
    return CIFLoader( txt, quiet = True, **kw )

def expect_error( descr, txt, **kw ):
    try:
        load( txt, **kw )
    except NC.NCBadInput as e:
        print( f'   {descr}: rejected: {e}' )
    else:
        raise SystemExit( f'{descr}: not rejected' )

def bragg( lc ):
    info = lc.create_ncmat_composer( uiso_temperature = 300,
                                     quiet = True ).load( 'dcutoff=0.5' ).info
    return sorted( ( round( h.d, 9 ), round( h.f2 * h.mult, 6 ) )
                   for h in info.hklObjects() )

def same_bragg( b1, b2 ):
    return len( b1 ) == len( b2 ) and all(
        abs( x[0] - y[0] ) < 1e-6
        and abs( x[1] - y[1] ) < 1e-5 * max( 1.0, x[1] )
        for x, y in zip( b1, b2 ) )

def show( name, lc ):
    print( f'== {name}: sg={lc.cellsg["spacegroup"]["number"]}'
           f' natoms={len(lc.atoms)}' )
    for s in lc.sites:
        u = s['aniso']
        us = ( ' aniso=' + ' '.join( f'{u[k]:.6g}' for k in keys )
               if u else '' )
        print( f'   {",".join(s["cif_labels"]):24} mult={s["multiplicity"]}'
               f' occu={s["occupancy"]:g} uiso={s["uiso"]:.6g}{us}' )
    for a in lc.atoms:
        assert set( a.keys() ) == { 'cif_labels', 'composition',
                                    'equivalent_positions', 'occupancy',
                                    'uiso' }

rutile = { 'a': 4.594, 'b': 4.594, 'c': 2.959,
           'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0 }
u_ti = { 'u11': 0.006, 'u22': 0.006, 'u33': 0.004,
         'u12': 0.0005, 'u13': 0.0, 'u23': 0.0 }
u_o = { 'u11': 0.009, 'u22': 0.009, 'u33': 0.006,
        'u12': -0.003, 'u13': 0.0, 'u23': 0.0 }
asym_rutile = [ ( 'Ti1', 'Ti', ( 0.0, 0.0, 0.0 ), 1.0,
                  sgs._cif_uij_ueq( u_ti, rutile ), u_ti ),
                ( 'O1', 'O', ( 0.305, 0.305, 0.0 ), 1.0,
                  sgs._cif_uij_ueq( u_o, rutile ), u_o ) ]

def test_rutile():
    ref = load( cif_text( rutile, 'P 42/m n m', asym_rutile ) )
    show( 'rutile, P42/mnm', ref )
    bref = bragg( ref )
    p1 = expand( '136', rutile, asym_rutile )
    lc = load( cif_text( rutile, 'P 1', p1 ) )
    show( 'rutile, P1', lc )
    assert lc.cellsg['spacegroup']['number'] == 136
    assert all( s['aniso'] for s in lc.sites )
    assert same_bragg( bref, bragg( lc ) )
    print( '   Bragg reflections identical to P42/mnm' )
    #Unrotated copies (as in some exported files):
    p1b = [ ( lbl, el, p, occ, uiso, ( u_ti if el == 'Ti' else u_o ) )
            for lbl, el, p, occ, uiso, _ in p1 ]
    txt = cif_text( rutile, 'P 1', p1b )
    expect_error( 'rutile, P1, unrotated tensor copies', txt )
    lc = load( txt, permissive = True )
    show( 'rutile, P1, unrotated tensor copies, permissive', lc )
    for w, _ in lc.warnings:
        print( '   warning:', w )
    assert not any( s['aniso'] for s in lc.sites )
    #Occupancies of the two Ti positions differ, which lowers the symmetry:
    p1c = [ ( lbl, el, p, ( 0.46 if lbl == 'Ti1_0' else
                            0.54 if lbl == 'Ti1_1' else occ ), uiso, u )
            for lbl, el, p, occ, uiso, u in p1 ]
    lc = load( cif_text( rutile, 'P 1', p1c ) )
    show( 'rutile, P1, Ti occupancies 0.46/0.54', lc )
    assert lc.cellsg['spacegroup']['number'] != 136
    #Tiny differences are within the tolerances, larger ones not:
    p1d = [ ( lbl, el, p, ( 0.9997 if lbl == 'Ti1_1' else occ ), uiso, u )
            for lbl, el, p, occ, uiso, u in p1 ]
    for tol in ( None, 1e-4 ):
        lc = load( cif_text( rutile, 'P 1', p1d ), tol = tol )
        show( f'rutile, P1, Ti occupancies 1/0.9997, tol={tol}', lc )

def test_symmetry_information():
    #Listed operations and multiplicities, consistent or not:
    ops = evaluate_query( [ 'util', 'sgsym', 'strict', '136' ] )['ops']
    good = cif_text( rutile, 'P 42/m n m', asym_rutile, ops = ops,
                     mults = [ 2, 4 ] )
    lc = load( good )
    show( 'rutile, with listed operations and multiplicities', lc )
    expect_error( 'wrong multiplicity', cif_text( rutile, 'P 42/m n m',
                                                  asym_rutile, ops = ops,
                                                  mults = [ 2, 2 ] ) )
    #An unparsable operation: rejected, but in permissive mode the symbol is
    #used:
    bad = good.replace( f"'{ops[1]}'", "'x+500/999,y,z'" )
    expect_error( 'unparsable operation', bad )
    lc = load( bad, permissive = True )
    show( 'rutile, unparsable operation, permissive', lc )
    for w, _ in lc.warnings:
        print( '   warning:', w )
    #Cell slightly inconsistent with the space group (accepted since spglib
    #finds the same space group), or more than that (rejected):
    #(Uiso values are the Ueq of the tensors in the given cell. When clearly
    #distorted, the tensors of the symmetry images no longer have that Ueq,
    #which is rejected first, so it is also tested without tensors):
    for gamma, with_uij in ( ( 90.002, True ), ( 90.5, True ),
                             ( 90.5, False ) ):
        c = dict( rutile, gamma = gamma )
        asym = [ ( lbl, el, pos, occ, sgs._cif_uij_ueq( u, c ),
                   u if with_uij else None )
                 for lbl, el, pos, occ, _, u in asym_rutile ]
        txt = cif_text( c, 'P 42/m n m', asym )
        descr = f'rutile, gamma={gamma}' + ( '' if with_uij else ', no uij' )
        if gamma < 90.1:
            lc = load( txt )
            show( descr, lc )
            for w, _ in lc.warnings:
                print( '   warning:', w )
            assert lc.cellsg['spacegroup']['number'] == 136
        else:
            expect_error( descr, txt )

def test_primitive():
    #C2/c with tensors on both sites, given in P1 in a primitive cell:
    conv = { 'a': 9.7, 'b': 8.9, 'c': 5.3, 'alpha': 90.0, 'beta': 105.6,
             'gamma': 90.0 }
    u_si = { 'u11': 0.012, 'u22': 0.016, 'u33': 0.010, 'u12': 0.002,
             'u13': 0.003, 'u23': -0.0015 }
    u_o = { 'u11': 0.020, 'u22': 0.014, 'u33': 0.017, 'u12': 0.0,
            'u13': 0.004, 'u23': 0.0 }
    asym = [ ( 'Si1', 'Si', ( 0.11, 0.23, 0.37 ), 1.0,
               sgs._cif_uij_ueq( u_si, conv ), u_si ),
             ( 'O1', 'O', ( 0.0, 0.41, 0.25 ), 1.0,
               sgs._cif_uij_ueq( u_o, conv ), u_o ) ]
    ref = load( cif_text( conv, 'C 1 2/c 1', asym ) )
    show( 'C2/c', ref )
    #Site multiplicities, either the number of positions (as in the CIF
    #dictionary) or the order of the site symmetry (as written by SHELXL), but
    #consistently so:
    for mults, ok in ( ( [ 8, 4 ], True ), ( [ 1, 2 ], True ),
                       ( [ 1, 4 ], False ), ( [ 8, 2 ], False ),
                       ( [ 4, 4 ], False ) ):
        txt = cif_text( conv, 'C 1 2/c 1', asym, mults = mults )
        if ok:
            load( txt )
            print( f'   multiplicities {mults}: accepted' )
        else:
            expect_error( f'multiplicities {mults}', txt )
    allc = expand( '15:b1', conv, asym )
    L = np.array( sgs._cellparams_to_spglib_lattice( conv ) ).T
    M = np.array( [ [ 0.5, 0.5, 0 ], [ -0.5, 0.5, 0 ], [ 0, 0, 1 ] ] )
    Lp = L @ M
    v = [ Lp[:,i] for i in range(3) ]
    n = [ float( np.linalg.norm( e ) ) for e in v ]
    def ang( i, j ):
        return float( np.degrees( np.arccos( np.dot( v[i], v[j] )
                                             / ( n[i] * n[j] ) ) ) )
    prim = { 'a': n[0], 'b': n[1], 'c': n[2], 'alpha': ang(1,2),
             'beta': ang(0,2), 'gamma': ang(0,1) }
    Minv = np.linalg.inv( M )
    sites, seen = [], []
    for lbl, el, p, occ, uiso, u in allc:
        x = tuple( float( q ) % 1.0 for q in Minv @ np.array( p ) )
        if any( np.max( np.abs( ( np.array( x ) - y + 0.5 ) % 1.0 - 0.5 ) )
                < 1e-9 for y in seen ):
            continue
        seen.append( np.array( x ) )
        uf = Minv @ sgs._cif_to_frac( u, conv ) @ Minv.T
        sites.append( ( lbl, el, x, occ, uiso, sgs._frac_to_cif( uf, prim ) ) )
    lc = load( cif_text( prim, 'P 1', sites ) )
    show( 'C2/c given in P1 in a primitive cell', lc )
    assert lc.cellsg['spacegroup']['number'] == 15
    assert all( s['aniso'] for s in lc.sites )
    bref = bragg( ref )
    assert same_bragg( bref, bragg( lc ) )
    print( '   Bragg reflections identical to C2/c' )

    #As written by pymatgen for primitive cells: the asymmetric unit with
    #symmetry operations in the primitive basis, but the standard symbol:
    from fractions import Fraction
    def fmtop( R, t ):
        comps = []
        for i in range(3):
            c = ''
            for j, v in enumerate( 'xyz' ):
                k = round( R[i][j] )
                if k:
                    c += ( '+' if k > 0 and c else '' ) + ( '-' if k < 0 else '' )
                    c += ( f'{abs(k)}*' if abs(k) != 1 else '' ) + v
            f = Fraction( float( t[i] ) % 1.0 ).limit_denominator( 12 )
            c += f'+{f}' if f else ''
            comps.append( c )
        return ','.join( comps )
    ops = []
    for o in evaluate_query( [ 'util', 'sgsym', 'strict', '15:b1' ] )['ops']:
        R, t = sgs._parse_xyz_op( o )
        Rp, tp = Minv @ R @ M, ( Minv @ t ) % 1.0
        s = fmtop( Rp, tp )
        if s not in ops:
            ops.append( s )
    asym_p = []
    for lbl, el, p, occ, uiso, u in asym:
        x = tuple( float( q ) % 1.0 for q in Minv @ np.array( p ) )
        uf = Minv @ sgs._cif_to_frac( u, conv ) @ Minv.T
        asym_p.append( ( lbl, el, x, occ, uiso, sgs._frac_to_cif( uf, prim ) ) )
    txt = cif_text( prim, 'C 1 2/c 1', asym_p, ops = ops, mults = [ 4, 2 ],
                    header = '# generated using pymatgen' )
    lc = load( txt )
    show( f'C2/c, primitive cell with {len(ops)} listed operations', lc )
    assert lc.cellsg['spacegroup']['number'] == 15
    assert same_bragg( bref, bragg( lc ) )
    print( '   Bragg reflections identical to C2/c' )
    #Operations missing (caught by the multiplicities):
    for bad, descr in ( ( txt.replace( f"'{ops[1]}'\n", '' ),
                          'missing operation' ),
                        ( txt.replace( f"'{ops[1]}'\n",
                                       "'x+500/999,y,z'\n" ),
                          'unparsable operation' ),
                        ( txt.replace( f"'{ops[1]}'\n", '' ).replace(
                              f"'{ops[3]}'\n", '' ),
                          'operations of a smaller group' ),
                        ( txt.replace( f"'{ops[2]}'\n", '' ).replace(
                              f"'{ops[3]}'\n", '' ),
                          'operations of a tabulated smaller group' ) ):
        for permissive in ( False, True ):
            expect_error( f'{descr}, permissive={permissive}', bad,
                          permissive = permissive )

test_rutile()
test_symmetry_information()
test_primitive()
