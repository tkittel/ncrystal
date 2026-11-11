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
import NCTestUtils.reprint_escaped_warnings # noqa F401
import NCrystalDev as NC
from NCTestUtils.common import ensure_error


def main():
    print('\n\n  ================> Al example 1\n\n')

    c_Al = NC.NCMATComposer()
    c_Al.set_cellsg_cubic( 4.05 )
    c_Al.set_atompos( [ ('Al',0,0,0),
                        ('Al',0,1/2,1/2),
                        ('Al',1/2,0,1/2),
                       ('Al',1/2,1/2,0)])
    c_Al.allow_fallback_dyninfo()
    with ensure_error(NC.NCBadInput,
                      "Must provide a space group number (or invoke"
                      " .refine_crystal_structure()) before it is"
                      " possible to verify a crystal structure"):
        c_Al.create_ncmat()

    c_Al.refine_crystal_structure()

    print( c_Al() )
    c_Al.load().dump()

    print('\n\n  ================> Al example 2\n\n')
    c_Al2 = NC.NCMATComposer()
    c_Al2.set_cellsg_cubic( 4.05, spacegroup=225 )
    c_Al2.set_atompos( [('Al',0,0,0),('Al',0,1/2,1/2),('Al',1/2,0,1/2),('Al',1/2,1/2,0)])
    c_Al2.allow_fallback_dyninfo()
    c_Al2.set_composition('Al','0.99 Al 0.01 Cr')
    print(c_Al2())

    print('\n\n  ================> Al example 3\n\n')

    c_Al3 = NC.NCMATComposer()
    c_Al3.set_cellsg_cubic( 4.05 )
    c_Al3.set_atompos( [('tight_atom',0,0,0),('loose_atom',0,1/2,1/2),('loose_atom',1/2,0,1/2),('loose_atom',1/2,1/2,0)])
    c_Al3.set_dyninfo_msd('tight_atom',msd=0.005, temperature=200)
    c_Al3.set_dyninfo_msd('loose_atom',msd=0.02, temperature=200)
    c_Al3.set_composition('tight_atom','Al')
    c_Al3.set_composition('loose_atom','Al')
    c_Al3.refine_crystal_structure()#Detect spacegroup
    print(c_Al3())

    print('\n\n  ================> Al SANS example\n\n')

    c=NC.NCMATComposer('Al_sg225.ncmat') #<--- Can init from cfg-string.
    c.add_secondary_phase(0.01,'void.ncmat')
    c.add_hard_sphere_sans_model(50)
    print(c())
    c.load().dump()


def test_verify_origin_and_triclinic():
    #verify_crystal_structure must not reject valid structures just because
    #spglib standardises with a different origin (arbitrary for P1 and polar
    #space groups) or to a different (reduced) triclinic basis:
    print('\n\n  ================> verify_crystal_structure\n\n')
    def verify( cell, atoms, sg ):
        c = NC.NCMATComposer()
        a, b, cc, al, be, ga = cell
        c.set_cellsg( a=a, b=b, c=cc, alpha=al, beta=be, gamma=ga,
                      spacegroup = sg )
        c.set_atompos( atoms )
        c.allow_fallback_dyninfo()
        c.verify_crystal_structure( quiet = True )
        #(creating the data converts the positions to sites, which also
        #requires the symmetry to hold):
        txt = c.create_ncmat()
        sgline = txt.splitlines()[ txt.splitlines().index('@SPACEGROUP') + 1 ].strip()
        nsites = txt.count('\n', txt.index('@ATOMSITES'), txt.index('@DYNINFO') ) - 1
        return sgline, nsites
    tricl = (5.1,6.3,7.2,81,97,103)#mixed acute/obtuse angles
    p1 = [('Al',0.13,0.71,0.29),('O',0.62,0.05,0.93),('Fe',0.41,0.38,0.57)]
    verify( tricl, p1, 1 )
    pm1 = [('Al',0.1,0.2,0.3),('Al',0.9,0.8,0.7),
           ('O',0.35,0.05,0.6),('O',0.65,0.95,0.4)]
    verify( tricl, pm1, 2 )
    p = 4.05/2**0.5 #cell where spglib keeps the basis but shifts the origin:
    verify( (p,p,p,60,60,60), [('Ni',0.0424,0.8626,0.7802),
                               ('Al',0.0126,0.2354,0.7916)], 1 )
    fcc = [('Al',0,0,0),('Al',0,.5,.5),('Al',.5,0,.5),('Al',.5,.5,0)]
    #An fcc structure shifted off the standard origin: spglib verifies it
    #(it finds the symmetry for any origin), but NCrystal requires the
    #standard origin of the setting, so the data can not be created:
    c = NC.NCMATComposer()
    c.set_cellsg( a=4.05, b=4.05, c=4.05, alpha=90, beta=90, gamma=90,
                  spacegroup = 225 )
    c.set_atompos( [ (e,(x+.123)%1,(y+.456)%1,(z+.789)%1) for e,x,y,z in fcc ] )
    c.allow_fallback_dyninfo()
    c.verify_crystal_structure( quiet = True )
    with ensure_error( NC.NCBadInput ):
        c.create_ncmat()
    print('Origin-shifted fcc verified by spglib but rejected at creation')
    #Non-standard settings are accepted when they are settings of ITVB (the
    #data is then written with the setting declared, here: permuted
    #orthorhombic axes, and monoclinic with another unique axis), and rejected
    #otherwise (tetragonal with the 4-fold axis along a):
    import random

    import NCTestUtils.cifgen as cg
    rng = random.Random( 1 )
    cyclic = [[0,0,1],[1,0,0],[0,1,0]]
    def verify_struct( s, sg ):
        return verify( s.cellpars, [ (e,*p) for e,p in s.atoms ], sg )
    print( 'Permuted SG-19 written as:', *verify_struct( cg.transformed_basis(
        cg.random_structure( rng, 19, nsites = 2 ), cyclic ), 19 ) )
    print('Valid structures verified OK')
    #Wrong space groups or settings must be rejected:
    for cell, atoms, sg in ( ( tricl, p1, 2 ),
                             ( (4.05,)*3+(90,)*3, fcc, 229 ) ):
        with ensure_error(NC.NCBadInput):
            verify( cell, atoms, sg )
    s = cg.transformed_basis( cg.random_structure( rng, 4, nsites = 3 ),
                              [[0,1,0],[0,0,1],[1,0,0]] )
    sgline, nsites = verify_struct( s, 4 )
    print( f'Monoclinic with a non-b unique axis written as @SPACEGROUP'
           f' {sgline}, {nsites} sites' )
    assert sgline in ( '4:a', '4:c' )
    s = cg.transformed_basis( cg.random_structure( rng, 76, nsites = 3 ),
                              cyclic )
    try:
        verify_struct( s, 76 )
    except NC.NCBadInput as e:
        assert 'space group 76' in str(e).lower(), str(e)
        print('Tetragonal with the 4-fold axis along a rejected')
    else:
        raise RuntimeError('Incompatible setting was accepted')

def test_refine_corrections():
    #refine_crystal_structure must only report corrections when spglib
    #actually moves atoms, not for mere basis changes (e.g. triclinic
    #reduction, primitive to centred cells) or origin shifts. Anisotropic
    #properties can only be kept if the basis is unchanged.
    print('\n\n  ================> refine corrections\n\n')
    from NCrystalDev._sgsymimpl import _spglib_refine_cell
    def refine( cell, atoms ):
        c = NC.NCMATComposer()
        a, b, cc, al, be, ga = cell
        c.set_cellsg( a=a, b=b, c=cc, alpha=al, beta=be, gamma=ga )
        c.set_atompos( atoms )
        d = _spglib_refine_cell( c.as_spglib_cell()[0] )
        return d['sgno'], d['warnings'], d['can_keep_anisotropic_properties']
    p = 4.05/2**0.5
    fcc = [('Al',0,0,0),('Al',0,.5,.5),('Al',.5,0,.5),('Al',.5,.5,0)]
    for descr, cell, atoms, expect in (
            ( 'P1, triclinic basis change', (5.1,6.3,7.2,81,97,103),
              [('Al',0.13,0.71,0.29),('O',0.62,0.05,0.93),
               ('Fe',0.41,0.38,0.57)], (1,[],False) ),
            ( 'P1, repeated reductions', (p,p,p,60,60,60),
              [('Al',0.24,0.54,0.38),('O',0.62,0.06,0.94),
               ('Fe',0.40,0.83,0.11)], (1,[],False) ),
            ( 'Imm2, from primitive cell', (p,p,p,60,60,60),
              [('Al',.3,.1,.7),('O',.9,.5,.2)], (44,[],False) ),
            ( 'fcc, shifted origin', (4.05,)*3+(90,)*3,
              [ (e,(x+.1)%1,(y+.2)%1,(z+.3)%1) for e,x,y,z in fcc ],
              (225,[],True) ),
            ( 'fcc, displaced atom', (4.05,)*3+(90,)*3,
              [('Al',0.001,0,0)]+fcc[1:],
              (225,[('Structure received minor corrections by spglib'
                     ' (at the 0.1% level)')],True) ) ):
        res = refine( cell, atoms )
        print(f'{descr}: SG-{res[0]} warnings={res[1]} keep_aniso={res[2]}')
        assert res == expect, (descr,res)

def test_v8_scatknl_keywords():
    #The effective_temperature/msd/debye_temp parameters of
    #set_dyninfo_scatknl emit NCMAT v8 keywords, flowing to the
    #extension/elastic machinery on load. A non-crystalline material
    #without a declared state of matter is refused already at
    #create_ncmat time (NCMAT v8 requires it):
    import NCrystalDev as NC
    import numpy as np
    from NCrystalDev.misc import evaluate_query as q
    def mkc( **kwargs ):
        c = NC.NCMATComposer()
        c.set_composition('Al','Al')
        c.set_density( 2.7, 'g/cm3' )
        a = np.linspace(1e-6,100.0,10)
        b = np.linspace(0.0,100.0,10)#symmetric-scaled form
        c.set_dyninfo_scatknl( 'Al', alphagrid=a, betagrid=b,
                               temperature=300.0,
                               sab_scaled=np.ones(100), **kwargs )
        return c
    try:
        mkc( effective_temperature = 450.0 ).create_ncmat()
    except NC.NCBadInput as e:
        print(f'=> correctly refused stateless v8: {e}')
    c = mkc( effective_temperature = 450.0, msd = 0.05 )
    c.set_state_of_matter('solid')
    txt = c.create_ncmat()
    assert txt.startswith('NCMAT v8')
    assert '\n  effective_temperature 450\n' in txt
    assert '\n  msd .05\n' in txt
    r = q( ['sab','analyse',c.register_as('comp_v8.ncmat',''),''] )
    r = r['sabanalyse']
    assert r['teff_extender'] == 450.0
    print('=> v8 keywords roundtrip (teff_extender=450): OK')
    #no keywords -> also v8 (the composer only emits NCMAT v8):
    c2 = mkc()
    c2.set_state_of_matter('liquid')
    assert c2.create_ncmat().startswith('NCMAT v8')
    print('=> plain scatknl composer output is NCMAT v8: OK')

def test_atomsites():
    #set_atomsites lists one position per site, expanded with the given space
    #group setting (which must be unambiguous), and the NCMAT data lists
    #exactly the given sites:
    import NCrystalDev as NC
    print('\n\n  ================> set_atomsites\n\n')
    def mk( sg, sites, *, use_sites = True ):
        c = NC.NCMATComposer()
        c.set_cellsg_cubic( 5.43, spacegroup = sg )
        if use_sites:
            c.set_atomsites( sites )
        else:
            c.set_atompos( sites )
        c.set_dyninfo_debyetemp( 'Si', 400 )
        return c
    def natoms( c ):
        return c.load().info.structure_info['n_atoms_average']
    c = mk( '227:2', [ ('Si',1/8,1/8,1/8) ] )
    txt = c.create_ncmat()
    assert '@SPACEGROUP\n  227:2\n@ATOMSITES\n  Si 1/8 1/8 1/8\n' in txt
    print('diamond Si from a single site (227:2):',natoms(c),'atoms')
    print('same site with origin choice 1:',
          natoms( mk( '227:1', [ ('Si',0,0,0) ] ) ),'atoms')
    try:
        mk( 227, [ ('Si',0,0,0) ] ).create_ncmat()
    except NC.NCBadInput as e:
        print(f'=> correctly refused ambiguous setting: {e}')
    try:
        mk( '227:3', [ ('Si',0,0,0) ] )
    except NC.NCBadInput as e:
        print(f'=> correctly refused invalid setting: {e}')
    #A site occupancy (written in @ATOMOCCUPANCIES):
    c = mk( '227:2', [ ('Si',1/8,1/8,1/8,0.9) ] )
    assert '@ATOMOCCUPANCIES\n  Si .9\n' in c.create_ncmat()
    print('with occupancy 0.9:',natoms(c),'atoms')
    #set_atompos with a given setting checks exactly that setting:
    allpos = [ (lbl,p[0],p[1],p[2]) for lbl,p in
               NC.misc.evaluate_query(['util','atompositions','strict',
                                       '227:2','Si,0.125,0.125,0.125']) ]
    print('all positions with 227:2:',
          natoms( mk( '227:2', allpos, use_sites = False ) ),'atoms')
    try:
        mk( '227:1', allpos, use_sites = False ).create_ncmat()
    except NC.NCBadInput as e:
        print(f'=> correctly refused positions in wrong setting: {e}')

def test_atomsites_uij():
    #Anisotropic displacements given with set_atomsites:
    import NCrystalDev as NC
    print('\n\n  ================> set_atomsites with uij\n\n')
    def mk( sites ):
        c = NC.NCMATComposer()
        c.set_cellsg( a = 3.21, b = 3.21, c = 5.21, alpha = 90, beta = 90,
                      gamma = 120, spacegroup = 194 )
        c.set_atomsites( sites )
        c.set_dyninfo_debyetemp( 'Mg', 400 )
        return c
    u = ( 0.0100, 0.0100, 0.0120, 0.0050, 0, 0 )
    for sites in ( [ ('Mg',1/3,2/3,1/4,u) ], [ ('Mg',1/3,2/3,1/4,0.9,u) ] ):
        c = mk( sites )
        txt = c.create_ncmat()
        i = txt.index('@ATOMSITES')
        print( txt[i:txt.index('@DYNINFO')], end = '' )
        info = c.load().info
        ai = info.atominfos[0]
        print( 'loaded uij:', ', '.join( f'{x:.6g}' for x in ai.sites[0]['uij'] ),
               f'(msd={ai.msd:.6g})' )
    for descr, sites in [
            ( 'uij for only some sites',
              [ ('Mg',1/3,2/3,1/4,u), ('Mg',0,0,0) ] ),
            ( 'too few uij values', [ ('Mg',1/3,2/3,1/4,u[:5]) ] ),
            ( 'non-finite uij', [ ('Mg',1/3,2/3,1/4,u[:5]+(float('nan'),)) ] ) ]:
        try:
            mk( sites )
        except NC.NCBadInput as e:
            print(f'=> correctly refused {descr}: {e}')
        else:
            raise SystemExit(f'{descr}: not refused')
    #Refinement keeps sites and tensors (here unchanged):
    c = mk( [ ('Mg',1/3,2/3,1/4,u) ] )
    c.refine_crystal_structure( quiet = True )
    txt = c.create_ncmat()
    print( txt[txt.index('@CELL'):txt.index('@DYNINFO')], end = '' )

def test_symprec_alias():
    #symprec is an obsolete alias of tol_symprec (warning on first use):
    c = NC.NCMATComposer()
    c.set_cellsg_cubic( 4.05, spacegroup = 225 )
    c.set_atompos( [ ('Al',0,0,0), ('Al',0,1/2,1/2), ('Al',1/2,0,1/2),
                     ('Al',1/2,1/2,0) ] )
    c.refine_crystal_structure( tol_symprec = 0.02, quiet = True )
    c.verify_crystal_structure( tol_symprec = 0.02, quiet = True )
    c.verify_crystal_structure( tol = 2e-3, quiet = True )
    for i in range(2):
        c.refine_crystal_structure( symprec = 0.02, quiet = True )
        c.verify_crystal_structure( symprec = 0.02, quiet = True )
    with ensure_error( NC.NCBadInput,
                       'Do not specify both symprec and tol_symprec' ):
        c.verify_crystal_structure( symprec = 0.02, tol_symprec = 0.02 )
    #All parameters are keyword-only:
    for f in ( c.refine_crystal_structure, c.verify_crystal_structure ):
        try:
            f( 0.02 )
        except TypeError:
            pass
        else:
            raise SystemExit('positional parameter accepted')
    with ensure_error( NC.NCBadInput, 'Invalid tol_symprec value: 0.6 (must'
                       ' not exceed 0.5Aa)' ):
        c.refine_crystal_structure( tol_symprec = 0.6 )
    print( 'symprec alias OK' )

def bragg( info ):
    return sorted( ( round( h.d, 9 ), round( h.f2 * h.mult, 6 ) )
                   for h in info.hklObjects() )

def test_refine_uij():
    #A C2/c structure with tensors given as P1 in a primitive cell (so spglib
    #must change the basis) is refined to C2/c sites with tensors, giving the
    #same Bragg reflections (d, mult*F^2, which depend on the tensors) as the
    #structure given directly in C2/c:
    import NCrystalDev._sgsymimpl as sgs
    import numpy as np
    from NCrystalDev.misc import evaluate_query
    conv = { 'a': 9.7, 'b': 8.9, 'c': 5.3, 'alpha': 90, 'beta': 105.6,
             'gamma': 90 }
    u_si = ( 0.012, 0.016, 0.010, 0.002, 0.003, -0.0015 )
    u_o = ( 0.020, 0.014, 0.017, 0.0, 0.004, 0.0 )
    asym = [ ( 'Si', ( 0.11, 0.23, 0.37 ), u_si ),
             ( 'O', ( 0.0, 0.41, 0.25 ), u_o ) ]
    keys = ( 'u11', 'u22', 'u33', 'u12', 'u13', 'u23' )
    def mk( cellsg, sites ):
        c = NC.NCMATComposer()
        c.set_cellsg( **cellsg )
        c.set_atomsites( sites )
        c.set_dyninfo_debyetemp( 'Si', 500 )
        c.set_dyninfo_debyetemp( 'O', 400 )
        return c
    ref = mk( dict( conv, spacegroup = '15:b1' ),
              [ ( lbl, *p, u ) for lbl, p, u in asym ] )
    ref = ref.load( 'dcutoff=0.5' ).info
    #All positions in the conventional cell:
    ops = evaluate_query( [ 'util', 'sgsym', 'strict', '15:b1' ] )['ops']
    allc = sgs._expand_sites( ops, conv,
                              [ { 'key': lbl, 'pos': p,
                                  'aniso': dict( zip( keys, u ) ) }
                                for lbl, p, u in asym ] )
    #Primitive cell (a-b)/2, (a+b)/2, c:
    def latt( c ):
        return np.array( sgs._cellparams_to_spglib_lattice( c ) ).T
    M = np.array( [ [ 0.5, 0.5, 0 ], [ -0.5, 0.5, 0 ], [ 0, 0, 1 ] ] )
    Lp = latt( conv ) @ M
    v = [ Lp[:,i] for i in range(3) ]
    n = [ float( np.linalg.norm( e ) ) for e in v ]
    def ang( i, j ):
        return float( np.degrees( np.arccos( np.dot( v[i], v[j] )
                                             / ( n[i] * n[j] ) ) ) )
    prim = { 'a': n[0], 'b': n[1], 'c': n[2], 'alpha': ang(1,2),
             'beta': ang(0,2), 'gamma': ang(0,1) }
    Minv = np.linalg.inv( M )
    def to_prim( e, unrotated = None ):
        x = tuple( float( q ) % 1.0 for q in Minv @ np.array( e['pos'] ) )
        a = unrotated or e['aniso']
        uf = Minv @ sgs._cif_to_frac( a, conv ) @ Minv.T
        u = sgs._frac_to_cif( uf, prim )
        return ( e['key'], *x, tuple( u[k] for k in keys ) )
    def same( s, q ):
        return s[0] == q[0] and max( min( abs( a - b ), 1 - abs( a - b ) )
                                     for a, b in zip( s[1:4], q[1:4] ) ) < 1e-9
    def prim_sites( unrotated = False ):
        res = []
        for e in allc:
            s = to_prim( e, dict( zip( keys, u_si ) ) if ( unrotated and
                                                          e['key'] == 'Si' )
                         else None )
            if not any( same( s, q ) for q in res ):
                res.append( s )
        return res
    sites = prim_sites()
    print( f'P1 primitive cell: {len(sites)} sites' )
    c = mk( dict( prim, spacegroup = 1 ), sites )
    c.refine_crystal_structure( quiet = True )
    txt = c.create_ncmat()
    print( txt[txt.index('@CELL'):txt.index('@DYNINFO')], end = '' )
    b1, b2 = bragg( ref ), bragg( c.load( 'dcutoff=0.5' ).info )
    assert len( b1 ) == len( b2 ) and len( b1 ) > 100
    assert all( abs( x[0] - y[0] ) < 1e-6
                and abs( x[1] - y[1] ) < 1e-6 * max( 1.0, x[1] )
                for x, y in zip( b1, b2 ) )
    print( f'Bragg reflections identical to direct C2/c: {len(b1)} families' )

    #Unrotated copies of the Si tensor (general position, so images under
    #the 2-fold axis have U12 and U23 of opposite sign): those are dropped
    #(with a warning), the O tensors are kept:
    c = mk( dict( prim, spacegroup = 1 ), prim_sites( unrotated = True ) )
    c.refine_crystal_structure( quiet = False )
    txt = c.create_ncmat()
    print( txt[txt.index('@CELL'):txt.index('@DYNINFO')], end = '' )

if __name__ == '__main__':
    main()
    test_verify_origin_and_triclinic()
    test_refine_corrections()
    test_v8_scatknl_keywords()
    test_atomsites()
    test_atomsites_uij()
    test_refine_uij()
    test_symprec_alias()
