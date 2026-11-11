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

# Unit tests of _sgsymimpl._find_site_orbits: symmetry detection with site data
# (occupancies, uiso and tensor eigenvalues) in the spglib types, and the check
# that tensors of equivalent sites are each other's symmetry images. Test
# structures are expanded to P1 with the operations of a space group.

import NCTestUtils.enable_fpe # noqa F401
import numpy as np
import spglib
import NCrystalDev._sgsymimpl as sgs

def sg_ops( sgno ):
    #Operations of the first (standard) setting of the space group:
    for hall in range( 1, 531 ):
        if spglib.get_spacegroup_type( hall ).number == sgno:
            d = spglib.get_symmetry_from_database( hall )
            return list( zip( d['rotations'], d['translations'] ) )

def expand( sgno, cellsg, asym ):
    #P1 sites from the asymmetric unit, with tensors rotated as they should:
    res = []
    for a in asym:
        pos = []
        uf = sgs._cif_to_frac( a['aniso'], cellsg ) if a.get('aniso') else None
        for R, t in sg_ops( sgno ):
            x = ( np.dot( R, a['pos'] ) + t ) % 1.0
            if any( np.max( np.abs( ( x - p + 0.5 ) % 1.0 - 0.5 ) ) < 1e-9
                    for p in pos ):
                continue
            pos.append( x )
            s = dict( a, pos = tuple( float( e ) for e in x ) )
            if uf is not None:
                uf2 = sgs._transform_frac_tensor( R, uf )
                s['aniso'] = sgs._frac_to_cif( uf2, cellsg )
            res.append( s )
    return res

def show( name, cellsg, sites, **tolkw ):
    r = sgs._find_site_orbits( cellsg, sites,
                               sgs._resolve_tolerances( **tolkw ) )
    def fmt( e ):
        bad = ' BAD-ANISO' if e['aniso_ok'] is False else ''
        return f"{sites[e['members'][0]]['key']}{e['members']}{bad}"
    o = '  '.join( fmt( e ) for e in r['orbits'] )
    print( f'{name:46} sg={r["sgno"]:<3} {o}' )
    return r

rutile = { 'a': 4.594, 'b': 4.594, 'c': 2.959,
           'alpha': 90, 'beta': 90, 'gamma': 90 }
u_ti = { 'u11': 0.006, 'u22': 0.006, 'u33': 0.004,
         'u12': 0.0005, 'u13': 0.0, 'u23': 0.0 }
u_o = { 'u11': 0.009, 'u22': 0.009, 'u33': 0.006,
        'u12': -0.003, 'u13': 0.0, 'u23': 0.0 }

def rutile_sites( **kw ):
    return expand( 136, rutile,
                   [ dict( key = 'Ti', pos = ( 0, 0, 0 ), occu = 1.0,
                           uiso = 0.0053, aniso = u_ti, **kw ),
                     dict( key = 'O', pos = ( 0.305, 0.305, 0 ), occu = 1.0,
                           uiso = 0.008, aniso = u_o, **kw ) ] )

def test_rutile():
    s = rutile_sites()
    assert len( s ) == 6
    r = show( 'rutile', rutile, s )
    assert r['sgno'] == 136 and all( e['aniso_ok'] for e in r['orbits'] )

    #Tensors copied without rotating them:
    s2 = [ dict( e, aniso = ( u_ti if e['key'] == 'Ti' else u_o ) ) for e in s ]
    r = show( 'rutile, unrotated tensor copies', rutile, s2 )
    assert r['sgno'] == 136
    assert not any( e['aniso_ok'] for e in r['orbits'] )

    #Element-only data (no tensors etc.):
    s3 = [ { 'key': e['key'], 'pos': e['pos'] } for e in s ]
    r = show( 'rutile, positions only', rutile, s3 )
    assert r['sgno'] == 136 and all( e['aniso_ok'] is None
                                     for e in r['orbits'] )

    #Small noise on tensors of non-representative sites:
    def noisy( rel ):
        res = []
        for i, e in enumerate( s ):
            a = dict( e['aniso'] )
            if i not in ( 0, 2 ):
                a['u33'] *= 1.0 + rel
            res.append( dict( e, aniso = a ) )
        return res
    r = show( 'rutile, 5e-4 tensor noise, tol=1e-3', rutile, noisy( 5e-4 ),
              tol = 1e-3 )
    assert r['sgno'] == 136 and all( e['aniso_ok'] for e in r['orbits'] )
    #With a tighter tolerance the eigenvalues differ, giving lower symmetry:
    r = show( 'rutile, 5e-4 tensor noise, tol=1e-4', rutile, noisy( 5e-4 ),
              tol = 1e-4 )
    assert r['sgno'] != 136

def test_occupancies():
    #Different occupancies on the two Ti positions (like COD 4117533):
    s = rutile_sites()
    def ti_occu( o1, o2 ):
        res, n = [], 0
        for e in s:
            if e['key'] == 'Ti':
                e = dict( e, occu = ( o1 if n == 0 else o2 ) )
                n += 1
            res.append( e )
        return res
    r = show( 'rutile, Ti occupancies 0.46/0.54', rutile,
              ti_occu( 0.46, 0.54 ) )
    assert r['sgno'] != 136
    r = show( 'rutile, Ti occu 0.5/0.5003, tol=1e-3', rutile,
              ti_occu( 0.5, 0.5003 ), tol = 1e-3 )
    assert r['sgno'] == 136
    r = show( 'rutile, Ti occu 0.5/0.5003, tol=1e-4', rutile,
              ti_occu( 0.5, 0.5003 ), tol = 1e-4 )
    assert r['sgno'] != 136
    #Same via uiso, with tol_uij overriding tol:
    res = [ dict( e, uiso = 0.0054 ) if i == 1 else e
            for i, e in enumerate( s ) ]
    r = show( 'rutile, one Ti uiso differs 1.9%', rutile, res )
    assert r['sgno'] != 136
    r = show( 'rutile, one Ti uiso differs 1.9%, tol_uij=0.05', rutile, res,
              tol_uij = 0.05 )
    assert r['sgno'] == 136

def test_eigenvalues():
    #Two of the O sites with differently shaped tensors (same uiso):
    s = rutile_sites()
    u2 = dict( u_o, u33 = 0.0066 )
    res, n = [], 0
    for e in s:
        if e['key'] == 'O':
            if n < 2:
                e = dict( e, aniso = u2 )
            n += 1
        res.append( e )
    r = show( 'rutile, O eigenvalues differ', rutile, res )
    assert r['sgno'] != 136

def test_monoclinic():
    #C2/c, general position (8 sites), beta != 90, tensor with all components:
    cell = { 'a': 9.7, 'b': 8.9, 'c': 5.3, 'alpha': 90, 'beta': 105.6,
             'gamma': 90 }
    u = { 'u11': 0.012, 'u22': 0.016, 'u33': 0.010, 'u12': 0.002,
          'u13': 0.003, 'u23': -0.0015 }
    s = expand( 15, cell, [ { 'key': 'Si', 'pos': ( 0.11, 0.23, 0.37 ),
                              'uiso': 0.0125, 'aniso': u },
                            { 'key': 'O', 'pos': ( 0.0, 0.41, 0.25 ),
                              'uiso': 0.02 } ] )
    assert len( s ) == 12
    r = show( 'C2/c', cell, s )
    assert r['sgno'] == 15 and r['orbits'][0]['aniso_ok']
    #Swap two tensors (keeps eigenvalues, breaks orientation):
    s2 = list( s )
    s2[1] = dict( s[1], aniso = s[2]['aniso'] )
    s2[2] = dict( s[2], aniso = s[1]['aniso'] )
    r = show( 'C2/c, two tensors swapped', cell, s2 )
    assert r['sgno'] == 15 and not r['orbits'][0]['aniso_ok']

test_rutile()
test_occupancies()
test_eigenvalues()
test_monoclinic()
