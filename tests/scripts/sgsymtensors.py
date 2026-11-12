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

# NEEDS: numpy spglib

# Unit tests of the basis change and displacement tensor helpers in _sgsymimpl,
# against Cartesian calculations (mostly non-cubic cells, which unlike cubic
# ones do not hide mix-ups of axes or transpositions).

import NCTestUtils.enable_fpe # noqa F401
import numpy as np
import NCrystalDev._sgsymimpl as sgs

def latt( c ):
    #Columns are the lattice vectors (a along x, b in the xy plane):
    a, b, cc = c['a'], c['b'], c['c']
    ca, cb, cg = ( np.cos( np.radians( c[k] ) )
                   for k in ( 'alpha', 'beta', 'gamma' ) )
    sg = np.sin( np.radians( c['gamma'] ) )
    cx = cc * cb
    cy = cc * ( ca - cb * cg ) / sg
    return np.array( [ [ a, b*cg, cx ], [ 0.0, b*sg, cy ],
                       [ 0.0, 0.0, np.sqrt( cc*cc - cx*cx - cy*cy ) ] ] )

def params( L ):
    #Cell parameters of lattice with vectors as columns:
    v = [ L[:,i] for i in range(3) ]
    n = [ np.linalg.norm( e ) for e in v ]
    def ang( i, j ):
        return float( np.degrees( np.arccos( np.dot( v[i], v[j] )
                                             / ( n[i]*n[j] ) ) ) )
    return { 'a': float(n[0]), 'b': float(n[1]), 'c': float(n[2]),
             'alpha': ang(1,2), 'beta': ang(0,2), 'gamma': ang(0,1) }

def ucart_from_cif( aniso, c ):
    #Directly from the CIF definition: U_cart = A.N.U.N.A^T
    L = latt( c )
    N = np.diag( np.linalg.norm( np.linalg.inv( L ), axis = 1 ) )
    u = np.array( [ [ aniso['u11'], aniso['u12'], aniso['u13'] ],
                    [ aniso['u12'], aniso['u22'], aniso['u23'] ],
                    [ aniso['u13'], aniso['u23'], aniso['u33'] ] ] )
    return L @ N @ u @ N @ L.T

def cif_from_ucart( ucart, c ):
    L = latt( c )
    N = np.linalg.norm( np.linalg.inv( L ), axis = 1 )
    uf = np.linalg.inv( L ) @ ucart @ np.linalg.inv( L ).T
    u = uf / np.outer( N, N )
    return { 'u11': u[0,0], 'u22': u[1,1], 'u33': u[2,2],
             'u12': u[0,1], 'u13': u[0,2], 'u23': u[1,2] }

def fmt( aniso ):
    def f( v ):
        return 0.0 if abs( v ) < 1e-12 else v
    return ' '.join( f'{k}={f(aniso[k]):+.6f}'
                     for k in ( 'u11','u22','u33','u12','u13','u23' ) )

def maxdiff( a1, a2 ):
    return max( abs( a1[k] - a2[k] ) for k in a1 )

def change_basis( name, cell, aniso, Pinv, *, expect = None ):
    #New basis vectors are the columns of L.Pinv, so x' = P.x. Checks the
    #mapped tensor against the Cartesian calculation (which is
    #rotation-free here, since both cells share the same Cartesian frame),
    #and the round trip, eigenvalues and Ueq.
    print( f'== {name}' )
    Pinv = np.array( Pinv, dtype = float )
    P = np.linalg.inv( Pinv )
    newcell = params( latt( cell ) @ Pinv )
    uf = sgs._cif_to_frac( aniso, cell )
    uf2 = sgs._transform_frac_tensor( P, uf )
    res = sgs._frac_to_cif( uf2, newcell )
    #The new cell parameters in Cartesian form might be rotated relative to
    #L.Pinv (a along x convention), so compare via the rotation-free route:
    ucart = ucart_from_cif( aniso, cell )
    L2 = latt( cell ) @ Pinv
    uf_ref = np.linalg.inv( L2 ) @ ucart @ np.linalg.inv( L2 ).T
    N2 = np.linalg.norm( np.linalg.inv( L2 ), axis = 1 )
    u = uf_ref / np.outer( N2, N2 )
    ref = { 'u11': u[0,0], 'u22': u[1,1], 'u33': u[2,2],
            'u12': u[0,1], 'u13': u[0,2], 'u23': u[1,2] }
    print( '   old:', fmt( aniso ) )
    print( '   new:', fmt( res ) )
    assert maxdiff( res, ref ) < 1e-12
    back = sgs._frac_to_cif( sgs._transform_frac_tensor( Pinv, uf2 ), cell )
    assert maxdiff( back, aniso ) < 1e-12
    e1 = sgs._cart_eigenvalues( uf, cell )
    e2 = sgs._cart_eigenvalues( uf2, newcell )
    e_ref = sorted( np.linalg.eigvalsh( ucart ) )
    assert max( abs( x - y ) for x, y in zip( e1, e_ref ) ) < 1e-12
    assert max( abs( x - y ) for x, y in zip( e2, e_ref ) ) < 1e-12
    assert abs( sgs._cif_uij_ueq( res, newcell )
                - sgs._cif_uij_ueq( aniso, cell ) ) < 1e-12
    if expect:
        expect( res )
    return res

def test_basis_changes():
    #Monoclinic (beta=103, nonzero U13), unique axis b, cell choice 1 -> 3
    #(a' = c, c' = -a-c):
    mono = { 'a': 5.1, 'b': 6.2, 'c': 7.3,
             'alpha': 90, 'beta': 103, 'gamma': 90 }
    u_mono = { 'u11': 0.011, 'u22': 0.017, 'u33': 0.023,
               'u12': 0.0, 'u13': 0.004, 'u23': 0.0 }
    def mono_ok( r ):
        #Still the unique axis b: U12 and U23 remain zero.
        assert abs( r['u12'] ) < 1e-14 and abs( r['u23'] ) < 1e-14
        assert abs( r['u22'] - u_mono['u22'] ) < 1e-14
    change_basis( 'monoclinic, cell choice change', mono, u_mono,
                  [ [ 0, 0, -1 ], [ 0, 1, 0 ], [ 1, 0, -1 ] ],
                  expect = mono_ok )

    #Hexagonal site tensor (U11=U22=2*U12) to the C-centred orthohexagonal
    #cell (a,a+2b,c), where it must be diagonal with U11'=U22':
    hexa = { 'a': 3.2, 'b': 3.2, 'c': 5.2,
             'alpha': 90, 'beta': 90, 'gamma': 120 }
    u_hex = { 'u11': 0.012, 'u22': 0.012, 'u33': 0.021,
              'u12': 0.006, 'u13': 0.0, 'u23': 0.0 }
    def ortho_ok( r ):
        for k in ( 'u12', 'u13', 'u23' ):
            assert abs( r[k] ) < 1e-14
        assert abs( r['u11'] - r['u22'] ) < 1e-14
    change_basis( 'hexagonal to orthohexagonal', hexa, u_hex,
                  [ [ 1, 1, 0 ], [ 0, 2, 0 ], [ 0, 0, 1 ] ],
                  expect = ortho_ok )

    #Rhombohedral axes (site on the 3-fold axis: U11=U22=U33, U12=U13=U23) to
    #hexagonal axes (obverse), where U11=U22=2*U12, U13=U23=0:
    rho = { 'a': 5.0, 'b': 5.0, 'c': 5.0,
            'alpha': 55.0, 'beta': 55.0, 'gamma': 55.0 }
    u_rho = { 'u11': 0.015, 'u22': 0.015, 'u33': 0.015,
              'u12': 0.004, 'u13': 0.004, 'u23': 0.004 }
    def hex_ok( r ):
        assert abs( r['u13'] ) < 1e-14 and abs( r['u23'] ) < 1e-14
        assert abs( r['u11'] - r['u22'] ) < 1e-14
        assert abs( r['u11'] - 2 * r['u12'] ) < 1e-14
    change_basis( 'rhombohedral to hexagonal axes', rho, u_rho,
                  [ [ 1, 0, 1 ], [ -1, 1, 1 ], [ 0, -1, 1 ] ],
                  expect = hex_ok )

    #Triclinic, general tensor, unimodular basis change:
    tri = { 'a': 4.3, 'b': 5.7, 'c': 6.1,
            'alpha': 81.0, 'beta': 97.0, 'gamma': 107.0 }
    u_tri = { 'u11': 0.014, 'u22': 0.019, 'u33': 0.011,
              'u12': 0.003, 'u13': -0.002, 'u23': 0.0015 }
    change_basis( 'triclinic, unimodular basis change', tri, u_tri,
                  [ [ 1, 1, 0 ], [ 0, 1, 1 ], [ 1, 0, 2 ] ] )

def test_transforms():
    t1 = ( np.array( [ [ 0, 1, 0 ], [ -1, 1, 0 ], [ 0, 0, 1 ] ] ),
           np.array( [ 0.25, 0.0, 0.5 ] ) )
    t2 = ( np.array( [ [ 1, 0, 1 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ] ),
           np.array( [ 0.0, 0.5, 0.125 ] ) )
    x = np.array( [ 0.1, 0.2, 0.3 ] )
    t = sgs._compose_transforms( t1, t2 )
    a = sgs._transform_position( t2, sgs._transform_position( t1, x ) )
    b = sgs._transform_position( t, x )
    assert np.max( np.abs( a - b ) ) < 1e-15
    t = sgs._compose_transforms( sgs._identity_transform(), t1 )
    assert np.array_equal( t[0], t1[0] ) and np.array_equal( t[1], t1[1] )
    print( 'transforms OK' )

def refine_primitive( name, conv, frac, ucart_site, prim_pinv, expect_sg ):
    #A conventional structure (cell, fractional positions, Cartesian site
    #tensor applying to all atoms) given in a primitive cell with tensors,
    #refined by spglib. The tracked transform must give the tensors of the
    #conventional cell.
    print( f'== {name}' )
    Lc = latt( conv )
    Lp = Lc @ np.array( prim_pinv, dtype = float )
    pcell = params( Lp )
    Lp2 = latt( pcell )#same cell in the a-along-x frame
    R = Lp2 @ np.linalg.inv( Lp )#rotation between the frames
    pos = []
    for f in frac:
        xp = np.linalg.solve( Lp, Lc @ np.array( f ) ) % 1.0
        if not any( np.allclose( xp, q ) for q in pos ):
            pos.append( xp )
    aniso_p = cif_from_ucart( R @ ucart_site @ R.T, pcell )
    d = sgs._spglib_refine_cell( ( Lp2.T, pos, [ 0 ] * len( pos ) ) )
    print( '   space group:', d['sgno'], d['sgsymb_hm'],
           '  can_keep_aniso:', d['can_keep_anisotropic_properties'] )
    assert d['sgno'] == expect_sg
    t = d['transform']
    assert t is not None
    newcell = d['cellparams_snapped']
    res = sgs._frac_to_cif( sgs._transform_frac_tensor(
        t[0], sgs._cif_to_frac( aniso_p, pcell ) ), newcell )
    print( '   primitive:', fmt( aniso_p ) )
    print( '   refined:  ', fmt( res ) )
    print( '   refined cell:', ', '.join( f'{k}={newcell[k]:g}' for k in
                                         ( 'a','b','c','alpha','beta',
                                           'gamma' ) ) )
    return res, newcell

def test_primitive_to_conventional():
    #Body-centred tetragonal (I4/mmm, atom at origin, site symmetry 4/mmm
    #so the tensor is diag(u1,u1,u3)):
    u1, u3 = 0.012, 0.027
    conv = { 'a': 3.1, 'b': 3.1, 'c': 5.3,
             'alpha': 90, 'beta': 90, 'gamma': 90 }
    res, c = refine_primitive( 'body-centred tetragonal', conv,
                               [ ( 0, 0, 0 ) ], np.diag( [ u1, u1, u3 ] ),
                               [ [ -0.5, 0.5, 0.5 ], [ 0.5, -0.5, 0.5 ],
                                 [ 0.5, 0.5, -0.5 ] ], 139 )
    assert abs( c['c'] - 5.3 ) < 1e-6 and abs( c['a'] - 3.1 ) < 1e-6
    for k, v in ( ( 'u11', u1 ), ( 'u22', u1 ), ( 'u33', u3 ),
                  ( 'u12', 0 ), ( 'u13', 0 ), ( 'u23', 0 ) ):
        assert abs( res[k] - v ) < 1e-6, ( k, res[k], v )

    #C-centred orthorhombic (Cmmm, atom at origin, tensor diag(u1,u2,u3)
    #along a,b,c with all lengths different):
    u = ( 0.011, 0.019, 0.025 )
    conv = { 'a': 3.3, 'b': 5.9, 'c': 4.4,
             'alpha': 90, 'beta': 90, 'gamma': 90 }
    res, c = refine_primitive( 'C-centred orthorhombic', conv,
                               [ ( 0, 0, 0 ) ], np.diag( u ),
                               [ [ 0.5, 0.5, 0 ], [ -0.5, 0.5, 0 ],
                                 [ 0, 0, 1 ] ], 65 )
    #spglib might permute the axes, so check by axis length:
    for i, k in enumerate( ( 'u11', 'u22', 'u33' ) ):
        L = c['abc'[i]]
        j = [ abs( L - conv[q] ) < 1e-6 for q in 'abc' ].index( True )
        assert abs( res[k] - u[j] ) < 1e-6, ( k, res[k], u[j] )
    for k in ( 'u12', 'u13', 'u23' ):
        assert abs( res[k] ) < 1e-6

def test_periodic_distance():
    #Squared Cartesian distance via the nearest periodic image:
    for cell in ( { 'a': 3.0, 'b': 4.0, 'c': 5.0, 'alpha': 90.0, 'beta': 90.0,
                    'gamma': 90.0 },
                  { 'a': 5.1, 'b': 6.3, 'c': 7.2, 'alpha': 84.0, 'beta': 105.0,
                    'gamma': 97.0 } ):
        L = latt( cell )#(columns are the cell vectors)
        d2 = sgs._periodic_distsq_fn( L.T )
        for x, y in ( ( ( 0.01, 0.5, 0.5 ), ( 0.99, 0.5, 0.5 ) ),
                      ( ( 0.1, 0.2, 0.3 ), ( 0.3, 0.1, 0.9 ) ),
                      ( ( 0.0, 0.0, 0.0 ), ( 1.0, -1.0, 2.0 ) ) ):
            d = ( np.asarray( x ) - np.asarray( y ) )
            d = d - np.round( d )
            ref = float( np.dot( L @ d, L @ d ) )
            assert abs( d2( x, y ) - ref ) < 1e-12 and d2( x, y ) == d2( y, x )
            print( f'periodic distance {x} {y}: {d2( x, y ) ** 0.5:.6f} Aa' )
        assert d2( ( 0.2, 0.2, 0.2 ), ( 0.2, 0.2, 0.2 ) ) == 0.0
    #Across the cell boundary along a (3 Aa): 0.02*3 Aa, not 0.98*3 Aa:
    d2 = sgs._periodic_distsq_fn( latt( { 'a': 3.0, 'b': 4.0, 'c': 5.0,
                                          'alpha': 90.0, 'beta': 90.0,
                                          'gamma': 90.0 } ).T )
    assert abs( d2( ( 0.01, 0, 0 ), ( 0.99, 0, 0 ) ) ** 0.5 - 0.06 ) < 1e-12

test_basis_changes()
test_transforms()
test_primitive_to_conventional()
test_periodic_distance()
