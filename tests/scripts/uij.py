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

# NEEDS: numpy

# Anisotropic displacements (uij) on @ATOMSITES lines: transformation to the
# orbit representative, snapping to the site symmetry, normalisation to the
# msd, and the rejection of invalid input.

import NCTestUtils.enable_fpe # noqa F401
import numpy as np
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query

_cnt = [0]
def load( cell, sg, sitelines, debye_temps = None, fractions = None ):
    elems = sorted( { ll.split()[0] for ll in sitelines } )
    debye_temps = debye_temps or {}
    fractions = fractions or {}
    dyninfo = ''.join( f"""@DYNINFO
  element {e}
  fraction {fractions.get(e,f'1/{len(elems)}')}
  type vdosdebye
  debye_temp {debye_temps.get(e,300)}
""" for e in elems )
    sites = ''.join( f'  {ll}\n' for ll in sitelines )
    _cnt[0] += 1
    fn = f'uijtest{_cnt[0]}.ncmat'
    NC.registerInMemoryFileData( fn, f"""NCMAT v8
@CELL
  {cell}
@SPACEGROUP
  {sg}
@ATOMSITES
{sites}{dyninfo}""" )
    return NC.createInfo( f'{fn};inelas=0;dcutoff=1' )

def require_error( descr, *args, **kwargs ):
    try:
        load( *args, **kwargs )
    except NC.NCBadInput as e:
        print( f'{descr}: {e}' )
    else:
        raise SystemExit(f'{descr}: did not raise')

def metric( si ):
    a, b, c = ( si[k] for k in ('a','b','c') )
    ca, cb, cg = ( np.cos( np.radians( si[k] ) )
                   for k in ('alpha','beta','gamma') )
    g = np.array( [ [ a*a, a*b*cg, a*c*cb ],
                    [ a*b*cg, b*b, b*c*ca ],
                    [ a*c*cb, b*c*ca, c*c ] ] )
    return g, np.sqrt( np.diag( np.linalg.inv( g ) ) )

def to_mat( u ):
    return np.array( [ [ u[0], u[3], u[4] ],
                       [ u[3], u[1], u[5] ],
                       [ u[4], u[5], u[2] ] ] )

def from_mat( m ):
    return [ m[0,0], m[1,1], m[2,2], m[0,1], m[0,2], m[1,2] ]

def ueq( si, u ):
    g, r = metric( si )
    ufrac = np.outer( r, r ) * to_mat( u )
    return np.sum( ufrac * g ) / 3.0

def ops_of( setting ):
    #Parse the CIF-like operation strings ("-y,x-y,z+1/3") of the sgsym query:
    ops = []
    for s in evaluate_query( ['util','sgsym','strict',setting] )['ops']:
        rot = np.zeros((3,3))
        tr = np.zeros(3)
        for i, comp in enumerate( s.replace(' ','').split(',') ):
            for term in comp.replace('-','+-').split('+'):
                if not term:
                    continue
                sign = -1.0 if term.startswith('-') else 1.0
                t = term.lstrip('-')
                if t in 'xyz':
                    rot[i,'xyz'.index(t)] += sign
                else:
                    num, _, den = t.partition('/')
                    tr[i] += sign * float(num) / ( float(den) if den else 1.0 )
        ops.append( (rot,tr) )
    return ops

def fmtu( u ):
    return '(' + ', '.join( f'{x:.6g}' for x in u ) + ')'

def check_normalisation( info ):
    si = info.structure_info
    for ai in info.atominfos:
        tot = sum( s['multiplicity'] * ueq( si, s['uij'] ) for s in ai.sites )
        avg = tot / sum( s['multiplicity'] for s in ai.sites )
        assert abs( avg - ai.msd ) < 1e-14 * ai.msd, ( avg, ai.msd )

def test_hexagonal():
    #Site 2c of P6_3/mmc (site symmetry -6m2): U11=U22=2*U12, U13=U23=0.
    cell = 'lengths 3.21 3.21 5.21\n  angles 90 90 120'
    keep = []
    def mg( uij ):
        keep.append( load( cell, 194, [ f'Mg 1/3 2/3 1/4 uij {uij}' ] ) )
        return keep[-1]
    info = mg( '0.0100 0.0100 0.0120 0.0050 0 0' )
    check_normalisation( info )
    u = info.atominfos[0].sites[0]['uij']
    print( 'Mg 2c, uij:', fmtu( u ) )
    assert u[0] == u[1] and u[0] == 2*u[3] and u[4] == 0.0 and u[5] == 0.0
    #Arbitrary normalisation, same result:
    u2 = mg( '1 1 1.2 0.5 0 0' ).atominfos[0].sites[0]['uij']
    assert np.allclose( u, u2, rtol = 1e-14, atol = 0 )
    #Rounded values are snapped:
    u3 = mg( '0.0100 0.0101 0.0120 0.0049 0 0' ).atominfos[0].sites[0]['uij']
    print( 'Mg 2c (rounded input), uij:', fmtu( u3 ) )
    assert abs( u3[0] - u3[1] ) <= 1e-17 and abs( u3[0] - 2*u3[3] ) <= 1e-17
    assert u3[4] == 0.0 and u3[5] == 0.0
    require_error( 'Broken site symmetry', cell, 194,
                   [ 'Mg 1/3 2/3 1/4 uij 0.0100 0.0100 0.0120 0.0050 0.002 0' ] )
    require_error( 'Broken site symmetry', cell, 194,
                   [ 'Mg 1/3 2/3 1/4 uij 0.0100 0.0100 0.0120 0 0 0' ] )

def test_transformation():
    #General position of P2_1/c with a non-orthogonal cell: giving the site
    #at any position of the orbit, with the tensor transformed accordingly,
    #gives the same result.
    setting = '14:b1'
    cell = 'lengths 5.1 6.3 7.2\n  angles 90 105 90'
    pos = np.array( [ 0.21, 0.13, 0.37 ] )
    u = [ 0.012, 0.009, 0.015, 0.002, -0.003, 0.001 ]
    ref = None
    si = None
    for rot, tr in ops_of( setting ):
        p = ( rot @ pos + tr ) % 1.0
        if si is None:
            si = load( cell, setting, [ 'O 0.21 0.13 0.37 uij 1 1 1 0 0 0' ] ).structure_info
        _, r = metric( si )
        ufrac = np.outer( r, r ) * to_mat( u )
        ut = from_mat( ( rot @ ufrac @ rot.T ) / np.outer( r, r ) )
        info = load( cell, setting, [ 'O {} {} {} uij {}'.format(
            *( repr(float(x)) for x in p ), ' '.join( repr(float(x)) for x in ut ) ) ] )
        check_normalisation( info )
        s = info.atominfos[0].sites[0]
        if ref is None:
            ref = s
            print( 'O in 14:b1 at', fmtu( s['position'] ), 'uij:', fmtu( s['uij'] ) )
        assert np.allclose( s['position'], ref['position'], rtol = 0, atol = 1e-12 )
        assert np.allclose( s['uij'], ref['uij'], rtol = 1e-12, atol = 0 )

def test_normalisation():
    #Two sites of the same element: the multiplicity-weighted average Ueq
    #equals the msd, and the ratio between the sites is kept:
    cell = 'lengths 3.21 3.21 5.21\n  angles 90 90 120'
    info = load( cell, 194, [ 'Mg 0 0 0 uij 1 1 1 0.5 0 0',
                              'Mg 1/3 2/3 1/4 uij 2 2 2 1 0 0',
                              'Zn 0 0 1/4 uij 1 1 3 0.5 0 0' ],
                 debye_temps = { 'Zn' : 200 },
                 fractions = { 'Mg' : '2/3', 'Zn' : '1/3' } )
    check_normalisation( info )
    for ai in info.atominfos:
        for s in ai.sites:
            print( f'{ai.atomData.displayLabel()} x{s["multiplicity"]}'
                   f' uij: {fmtu(s["uij"])}' )
    mg = info.atominfos[0] if info.atominfos[0].atomData.displayLabel() == 'Mg' else info.atominfos[1]
    u1, u2 = ( s['uij'] for s in mg.sites )
    assert np.allclose( np.asarray(u2), 2 * np.asarray(u1), rtol = 1e-14, atol = 0 )

def test_errors():
    cubic = 'lengths 4.05 4.05 4.05\n  angles 90 90 90'
    require_error( 'Not positive definite', cubic, 1,
                   [ 'Al 0 0 0 uij 0.01 0.01 0.01 0.011 0 0' ] )
    require_error( 'Not positive definite', cubic, 1,
                   [ 'Al 0 0 0 uij -0.01 -0.01 -0.01 0 0 0' ] )
    require_error( 'Only some sites', cubic, 1,
                   [ 'Al 0 0 0 uij 0.01 0.01 0.01 0 0 0', 'Al 1/2 1/2 1/2' ] )
    require_error( 'Too few numbers', cubic, 1,
                   [ 'Al 0 0 0 uij 0.01 0.01 0.01 0 0' ] )
    require_error( 'Bad keyword', cubic, 1,
                   [ 'Al 0 0 0 uji 0.01 0.01 0.01 0 0 0' ] )
    require_error( 'Bad number', cubic, 1,
                   [ 'Al 0 0 0 uij 0.01 0.01 0.01 0 0 abc' ] )
    require_error( 'Not finite', cubic, 1,
                   [ 'Al 0 0 0 uij 0.01 0.01 0.01 0 0 inf' ] )

def test_dump():
    cell = 'lengths 3.21 3.21 5.21\n  angles 90 90 120'
    info = load( cell, 194, [ 'Mg 1/3 2/3 1/4 uij 0.0100 0.0100 0.0120 0.0050 0 0' ] )
    lines = info.dump_str( verbose = 0 ).splitlines()
    i = lines.index('Atomic sites:')
    print( '\n'.join( lines[i:i+2] ) )

def fsq_python( info ):
    #|F|^2 of the hkl families of the info, calculated independently: the
    #tensor at each position is that of its site, transformed with the
    #operation mapping the site to the position.
    si = info.structure_info
    ops = ops_of( si['spacegroup_setting'] )
    _, r = metric( si )
    rr = np.outer( r, r )
    contribs = []
    for ai in info.atominfos:
        b = ai.atomData.coherentScatLen() * ai.occupancy
        for p in ai.positions:
            u = None
            for s in ai.sites:
                for rot, tr in ops:
                    q = rot @ np.asarray( s['position'] ) + tr - np.asarray( p )
                    if np.max( np.abs( q - np.round( q ) ) ) < 1e-6:
                        u = rot @ ( rr * to_mat( s['uij'] ) ) @ rot.T
                        break
                if u is not None:
                    break
            assert u is not None
            contribs.append( ( b, np.asarray( p ), u ) )
    res = []
    for h, k, ll, mult, d, f2 in info.hklList():
        hkl = np.array( [ h, k, ll ], dtype = float )
        f = sum( b * np.exp( -2.0 * np.pi**2 * ( hkl @ u @ hkl ) )
                 * np.exp( 2j * np.pi * ( hkl @ p ) ) for b, p, u in contribs )
        res.append( ( f2, abs( f )**2 ) )
    return res

def test_bragg():
    setting = '14:b1'
    cell = 'lengths 5.1 6.3 7.2\n  angles 90 105 90'
    sites = [ 'O 0.21 0.13 0.37 uij 0.012 0.009 0.015 0.002 -0.003 0.001',
              'Mg 0 0 0 uij 0.006 0.010 0.008 0.001 0.002 -0.001' ]
    fr = { 'O' : '2/3', 'Mg' : '1/3' }
    info = load( cell, setting, sites, fractions = fr )
    vals = fsq_python( info )
    maxrd = max( abs( a - b ) / b for a, b in vals )
    print( f'Bragg |F|^2 with anisotropic displacements: {len(vals)}'
           ' families agree with independent calculation' )
    assert maxrd < 1e-10, maxrd
    #Isotropic tensors (U_cart proportional to the identity) give the same
    #results as no tensors:
    si = info.structure_info
    _, r = metric( si )
    ginv = np.linalg.inv( metric( si )[0] )
    uiso = from_mat( ginv / np.outer( r, r ) )
    isosites = [ ' '.join( ll.split()[:4] + [ 'uij' ]
                           + [ repr(float(x)) for x in uiso ] ) for ll in sites ]
    plain = load( cell, setting, [ ' '.join( ll.split()[:4] ) for ll in sites ],
                  fractions = fr )
    iso = load( cell, setting, isosites, fractions = fr )
    a = [ e[5] for e in plain.hklList() ]
    b = [ e[5] for e in iso.hklList() ]
    assert len(a) == len(b)
    maxrd = max( abs( x - y ) / x for x, y in zip( a, b ) )
    print( 'Bragg |F|^2 with isotropic tensors: same as without tensors' )
    assert maxrd < 1e-10, maxrd

def main():
    test_hexagonal()
    test_transformation()
    test_normalisation()
    test_errors()
    test_dump()
    test_bragg()

if __name__ == '__main__':
    main()
