
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

Internal implementation of space group and symmetry utilities shared by
NCMATComposer and CIFLoader (spglib refinement, cell geometry, displacement
tensors).

"""

__all__ = []

import copy
import math

from . import _common as _nc_common
from . import core as _nc_core


def _lattice_params_to_vectors( a, b, c, alpha, beta, gamma ):
    if max(max(alpha,beta),gamma)<=math.pi:
        raise _nc_core.NCBadInput("angles are expected in degrees not radians")
    def _sincos( x ):
        if x == 0:
            return 0.0, 1.0
        if x == 90.0:
            return 1.0, 0.0
        if x == 120.0:
            return 0.86602540378443864676372317075293618347140262690519, -0.5
        return math.sin( x*(math.pi/180) ), math.cos( x*(math.pi/180) )
    _sa,ca = _sincos( alpha )
    sb,cb = _sincos( beta )
    sg,cg = _sincos( gamma )
    assert sg > 0.0
    vect_a = tuple( a*e for e in ( 1, 0.0, 0.0 ) )
    vect_b = tuple( b*e for e in ( cg, sg, 0.0 ) )
    ux = cb
    assert sg > 0.0
    assert sg**2 > 0.0
    uy = (ca-cb*cg)/sg
    uzsq = sb**2 - uy**2
    if uzsq < -1e-99:
        raise _nc_core.NCBadInput('Inconsistent cell parameters detected')
    vect_c = tuple( c*e for e in ( ux, uy, math.sqrt(max(0.0,uzsq)) ) )#sign of third component positive since we want the three vectors to form a right handed basis.
    return ( vect_a, vect_b, vect_c )

def _lattice_vectors_to_params( vect_a, vect_b, vect_c, as_dict = True ):
    def mag( v ):
        return math.sqrt( math.fsum(e**2 for e in v) )
    def dot( v1,v2 ):
        return math.fsum(e1*e2 for e1,e2 in zip(v1,v2))
    def ang( v1,v2 ):
        return math.acos( dot(v1,v2)/(mag(v1)*mag(v2)) )*180/math.pi
    a,b,c = mag(vect_a), mag(vect_b), mag(vect_c)
    alpha,beta,gamma = ang(vect_b,vect_c), ang(vect_a,vect_c), ang(vect_a,vect_b)
    if as_dict:
        return { 'a': a, 'b': b, 'c': c,
                 'alpha': alpha, 'beta': beta, 'gamma': gamma }
    else:
        return a, b, c, alpha, beta, gamma

def _snap_lattice_params( params_dict, relative_tolerance = 1e-6 ):
    log = []
    p = copy.deepcopy(params_dict)
    def rd( x,y ):
        return abs(x-y)/(max(1e-300,abs(x)+abs(y)))
    super_tolerance = 1e-14#almost 64 bit FP prec
    assert relative_tolerance > super_tolerance
    def is_near( x,y ):
        return rd(x,y) < relative_tolerance
    def is_super_near( x,y ):
        return rd(x,y) < super_tolerance

    #"snap" b,c to a,b and angles to 60/90/120:
    for l1,l2 in (('a','b'),('a','c'),('b','c')):
        if is_near(p[l1],p[l2]):
            if not is_super_near(p[l1],p[l2]):
                log.append(f'Snapping lattice parameter {l2} to {l1} since it was very near.')
            p[l2] = p[l1]
    for commonangle in (45,60,90,120):
        for angname in ['alpha','beta','gamma']:
            if is_near(p[angname],commonangle):
                if not is_super_near(p[angname],commonangle):
                    log.append(f'Snapping lattice angle {angname} to {commonangle} since it was very near.')
                p[angname] = commonangle
    return log, p

def _cellparams_to_spglib_lattice( cellsg ):
    return _lattice_params_to_vectors( cellsg['a'], cellsg['b'], cellsg['c'],
                                       cellsg['alpha'], cellsg['beta'], cellsg['gamma'] )

def _reldiff_cellparams( c1, c2 ):
    #compares a,b,c,alpha,beta,gamma in the two cellsg arrays.
    def rd( x,y ):
        return abs(x-y)/(max(1e-300,abs(x)+abs(y)))
    def rdf( s ):
        return rd(c1[s],c2[s])
    return max( rdf(s) for s in ('a','b','c','alpha','beta','gamma') )

def _reldiff_atompos( spglib_cell1, spglib_cell2 ):
    #input in spglib cell format, returns largest dist between positions, or
    #None if not the same type+number of atoms. A common translation of all
    #positions (i.e. a different choice of origin, which is arbitrary for
    #e.g. P1 or polar space groups) is not considered a difference.

    for c in ( spglib_cell1, spglib_cell2 ):
        assert c is not None
        assert len(c) >= 3
        assert c[0] is not None
    _,allpos1,idxlist1 = spglib_cell1
    _,allpos2,idxlist2 = spglib_cell2
    assert len(idxlist1) == len(allpos1)
    assert len(idxlist2) == len(allpos2)
    if sorted(idxlist1) != sorted(idxlist2) or len(idxlist1)!=len(idxlist2):
        return None

    def calc_maxdistsq( l1, l2, shift, abort_above ):
        #Max over l1 (shifted) of the squared distance to the nearest point in
        #l2. Returns None as soon as it is known to be above abort_above.
        distsqmax = 0.0
        for e1 in l1:
            e1 = ( e1[0]+shift[0], e1[1]+shift[1], e1[2]+shift[2] )
            dsqmin_e1 = min( _unit_cell_point_distsq(e1,e2) for e2 in l2 )
            if abort_above is not None and dsqmin_e1 > abort_above:
                return None
            distsqmax = max( dsqmin_e1, distsqmax )
        return distsqmax

    bytype = {}
    for idx in sorted(set(idxlist1)):
        bytype[idx] = tuple( [ tuple(float(x) for x in p)
                               for i,p in zip(idxlist,allpos) if i == idx ]
                             for idxlist, allpos in ( (idxlist1,allpos1),
                                                      (idxlist2,allpos2) ) )
    #Candidate shifts: those mapping one atom of the least abundant type onto
    #an atom of the same type:
    l1r, l2r = min( bytype.values(), key = lambda ll : len(ll[0]) )
    best = None
    for q in l2r:
        shift = tuple( e2-e1 for e1,e2 in zip(l1r[0],q) )
        worst = 0.0
        for l1, l2 in bytype.values():
            d = calc_maxdistsq( l1, l2, shift, best )
            if d is None:
                break
            worst = max( worst, d )
        else:
            if best is None or worst < best:
                best = worst
    return math.sqrt(best) if best is not None else None

def _spglib_std_transform( spglib_cell, symprec ):
    #Get the transformation (P,p,centring translations) used by spglib to
    #standardise a cell (x_std = P.x + p), or None if it can not be handled
    #(e.g. a supercell as input).
    from itertools import product as itertools_product

    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    ds = _import_spglib().get_symmetry_dataset( spglib_cell, symprec = symprec )
    if not ds:
        return None
    def get( k ):
        return getattr( ds, k ) if hasattr( ds, k ) else ds[k]
    P = np.asarray( get('transformation_matrix'), dtype = float )
    p = np.asarray( get('origin_shift'), dtype = float )
    detinv = abs( np.linalg.det( np.linalg.inv( P ) ) )
    ncentr = round( detinv )
    if ncentr < 1 or abs( detinv - ncentr ) > 1e-6:
        return None
    #Lattice points of the input lattice inside the standard cell:
    centr = []
    for n in itertools_product( range(ncentr), repeat = 3 ):
        t = np.dot( P, np.asarray(n,dtype=float) ) % 1.0
        if not any( _unit_cell_point_distsq(t,c) < 1e-12 for c in centr ):
            centr.append( t )
    if len(centr) != ncentr:
        return None
    return P, p, centr

def _spglib_apply_std_transform( spglib_cell, transform ):
    #Express cell in the basis and origin given by the transformation.
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    P, p, centr = transform
    lattice, positions, types = spglib_cell
    pos, typ = [], []
    for x, i in zip( positions, types ):
        xs = np.dot( P, np.asarray( x, dtype = float ) ) + p
        for t in centr:
            pos.append( tuple( float(e) for e in (xs+t) % 1.0 ) )
            typ.append( i )
    lattice_std = np.dot( np.linalg.inv( P ).T,
                          np.asarray( lattice, dtype = float ) )
    return lattice_std, pos, typ

def _import_spglib( *, sysexit = False ):
    try:
        import spglib  #both available on pypi and conda-forge
    except ImportError:
        m = ( 'Could not import spglib module needed to standardise and verify crystal structures.'
              +' The spglib package is available on both PyPI ("python3 -mpip install'
              +' spglib") and conda ("conda install -c conda-forge spglib")' )
        if sysexit:
            raise SystemExit(m)
        else:
            raise ImportError(m)
    if hasattr(spglib,'spglib'):
        #spglib version < 2.7
        return spglib.spglib
    if hasattr(spglib,'OLD_ERROR_HANDLING'):
        #For version 2.7 (and likely not later), we can avoid a deprecation
        #warning by this trick:
        try:
            import spglib.error
            spglib.error.OLD_ERROR_HANDLING = False
        except (ImportError,AttributeError):
            pass
    return spglib

def _spglib_extractsg( spglib_symdata ):
    sd = spglib_symdata
    sgno = ( sd.number
             if hasattr(sd,'number')
             else sd['number'] )
    sgsymb_hermann_mauguin = ( sd.international
                               if hasattr(sd,'international')
                               else sd['international'] )
    c = ( sd.choice
          if hasattr(sd,'choice')
          else sd.get('choice',None))
    if c:
        sgsymb_hermann_mauguin += f':{c}'
    return sgno, sgsymb_hermann_mauguin

def _uses_rhombohedral_axes( sgnumber, alpha ):
    #Rhombohedral (R) space groups in rhombohedral rather than hexagonal axes
    #(like usesRhombohedralAxes in NCLatticeUtils.hh):
    return ( int(sgnumber) in ( 146, 148, 155, 160, 161, 166, 167 )
             and alpha != 90 )

def _check_cell_sg_consistency( sgnumber, a, b, c, alpha, beta, gamma ):
    sgclass = _nc_common._classifySG(int(sgnumber))
    if sgclass in ( 'orthorhombic', 'tetragonal', 'cubic'):
        if sgclass=='cubic' and ( a!=b or a!=c ):
            raise _nc_core.NCBadInput(f'Cubic space group {sgnumber} requires cell parameters a==b==c (found a={a:g}, b={b:g}, c={c:g})')
        if not all( (_ == 90) for _ in (alpha, beta, gamma) ):
            raise _nc_core.NCBadInput(f'Space group {sgnumber} requires cell parameters alpha=beta=gamma=90 '
                                +f'(found alpha={alpha}, beta={beta}, gamma={gamma})')
    elif _uses_rhombohedral_axes( sgnumber, alpha ):
        if not ( a == b == c and alpha == beta == gamma and alpha < 120 ):
            raise _nc_core.NCBadInput(f'Space group {sgnumber} in rhombohedral'
                                      ' axes requires cell parameters a=b=c'
                                      ' and alpha=beta=gamma<120 (found'
                                      f' a={a:g}, b={b:g}, c={c:g},'
                                      f' alpha={alpha}, beta={beta},'
                                      f' gamma={gamma})')
        return
    elif sgclass in ('trigonal','hexagonal'):
        if not ( alpha == 90 and beta == 90 and gamma == 120 ):
            raise _nc_core.NCBadInput(f'Space group {sgnumber} requires cell parameters alpha=beta=90 and gamma=120 '
                                +f'(found alpha={alpha}, beta={beta}, gamma={gamma})')
    elif sgclass == 'monoclinic':
        #(which two angles are 90 depends on the unique axis of the setting,
        #which is verified precisely when the data is loaded)
        if not ( all( (0<_<180) for _ in ( alpha, beta, gamma ) )
                 and sum( (_ == 90) for _ in ( alpha, beta, gamma ) ) >= 2 ):
            raise _nc_core.NCBadInput(f'Space group {sgnumber} requires two'
                                      ' cell angles of 90 degrees (found'
                                      f' alpha={alpha}, beta={beta},'
                                      f' gamma={gamma})')
    else:
        if not all( (0<a<180) for a in ( alpha, beta, gamma ) ):
            raise _nc_core.NCBadInput('Unit cell must have all angles between 0 and 180'
                                +f' (found alpha={alpha}, beta={beta}, gamma={gamma})')
    if int(sgnumber) >= 75 and a!=b:
        raise _nc_core.NCBadInput(f'Space group {sgnumber} requires cell parameters a==b (found a={a:g}, b={b:g}, c={c:g})')

def _spglib_refine_cell( spglib_cell, symprec = 0.01 ):
    assert spglib_cell
    p_orig = _lattice_vectors_to_params( *spglib_cell[0], as_dict = True )

    warnings, msgs = [], []
    trackinfo = { 'cell' : spglib_cell, 'same_basis' : True }
    refined_cell, symdata = _impl_spglib_refine_cell( spglib_cell,
                                                      symprec=symprec,
                                                      warnings=warnings,
                                                      msgs=msgs,
                                                      trackinfo = trackinfo )

    #check difference from first to final, to warn/msg about corrections
    #(comparing in the same basis, since spglib might transform the cell):
    p_final = _lattice_vectors_to_params( *refined_cell[0], as_dict = True )
    _,p_final = _snap_lattice_params( p_final )
    sgno, sgsymb_hm = _spglib_extractsg( symdata )
    cell_cmp, same_basis = trackinfo['cell'], trackinfo['same_basis']
    if cell_cmp is None:
        cell_cmp = spglib_cell#could not track, compare as-is
    else:
        p_orig = _lattice_vectors_to_params( *cell_cmp[0], as_dict = True )
    rdl = _reldiff_cellparams( p_orig, p_final )
    rda = _reldiff_atompos( cell_cmp, refined_cell )
    rd = None if (rdl is None or rda is None) else max(rdl,rda)
    #Anisotropic properties refer to the original axes:
    discard_aniso = not same_basis
    if not same_basis:
        msgs.append('Unit cell was transformed to the'
                    ' standard setting by spglib')
    if rd is None or rd > 1e-2:
        discard_aniso = True
        _ = ( 'possibly due to conversion from primitive cell'
              if rd is None else f'at the {100.0*rd:g}% level' )
        warnings.append(f'Structure received major corrections by spglib ({_})')
    elif rd > 1e-6:
        warnings.append('Structure received minor corrections by '
                        f'spglib (at the {100.0*rd:g}% level)')
    else:
        msgs.append('Self-consistency of structure was verified by spglib')

    #remove all but the last of any kind of error message (otherwise we might
    #get many repetitions of 'snapping b=a' kind of messages, due to the
    #recursion.

    warnrev = []
    for w in warnings[::-1]:
        if w not in warnrev:
            warnrev.append(w)
    warnings = warnrev[::-1]

    return { 'refined_cell': refined_cell,
             'cellparams_snapped': p_final,
             'sgno': sgno,
             'sgsymb_hm': sgsymb_hm,
             'warnings': warnings,
             'msgs': msgs,
             'can_keep_anisotropic_properties': not discard_aniso,
             'reldiff' : rd }

def _impl_spglib_refine_cell( spglib_cell, *, symprec, warnings, msgs,
                              trackinfo, nrepeat = 0 ):
    #trackinfo['cell'] is the original cell, which is transformed along with
    #spglib_cell (so it can be compared with the final cell in the same basis).
    spglib = _import_spglib()

    orig_cell = spglib_cell
    if trackinfo['cell'] is not None:
        tr = _spglib_std_transform( orig_cell, symprec )
        if tr is None:
            trackinfo['cell'] = None
        else:
            trackinfo['cell'] = _spglib_apply_std_transform( trackinfo['cell'],
                                                             tr )
            trackinfo['same_basis'] = ( trackinfo['same_basis']
                                        and _np_is_identity( tr[0] ) )

    orig_cell_copy = copy.deepcopy( orig_cell )#copy is just a safeguard
    refined_cell = spglib.standardize_cell( orig_cell_copy, symprec = symprec )
    symdata = spglib.get_symmetry_dataset( orig_cell_copy )
    if not refined_cell or not symdata or refined_cell[0] is None:
        raise _nc_core.NCBadInput('Could not standardize provided'
                                  ' crystal structure with spglib')

    p = _lattice_vectors_to_params( *refined_cell[0], as_dict = True )

    snaplog, p = _snap_lattice_params( p )
    if snaplog or nrepeat == 0:
        #This is either the first try, of we snapped some parameter
        #(e.g. alpha=89.99999 -> alpha=90.0). After this a second refinement
        #might find a different (higher) symmetry:
        warnings += snaplog[:]
        refined_cell = tuple( [ _cellparams_to_spglib_lattice( p ) ]
                              + list( refined_cell[1:] ) )
        if nrepeat <= 3:
            return _impl_spglib_refine_cell( refined_cell,
                                             symprec = symprec,
                                             warnings=warnings,
                                             msgs=msgs,
                                             trackinfo = trackinfo,
                                             nrepeat = nrepeat + 1 )

    return refined_cell, symdata


def _np_is_identity( m ):
    #nb: spglib transformation matrices can have tiny numerical noise
    return all( abs( m[i][j] - (1.0 if i==j else 0.0) ) < 1e-9
                for i in range(3) for j in range(3) )

def _remap_fract_pos_pt(xyz):
    return ( _remap_fract_pos(xyz[0]),
             _remap_fract_pos(xyz[1]),
             _remap_fract_pos(xyz[2]) )

def _remap_fract_pos(x):
    x = float(x)
    if -1.0<=x<0.0:
        x += 1.0
    if 1.0<=x<2.0:
        x -= 1.0
    return x

def _unit_cell_point_dist( p1, p2 ):
    return math.sqrt( _unit_cell_point_distsq( p1, p2 ) )

def _unit_cell_point_distsq( p1, p2 ):
    #take into account unit cell periodicity, so e.g. (a,b,epsilon) and
    #(a,b,1.0-epsilon) are actually only 2*epsilon from each other, not
    #1.0-epsilon.
    def remap( x ):
        return x-math.floor(x)
    d = [ abs(remap(e1)-remap(e2)) for e1,e2 in zip(p1,p2)]
    return math.fsum(  min(e,1.0-e)**2 for e in d )

def _cif_uij_ueq( aniso, cellsg ):
    """Equivalent isotropic displacement, Ueq = trace(U_cart)/3, of CIF uij
    values (dict with keys u11,...,u23, relative to the normalised reciprocal
    axes), in the cell given by a dict with keys a,b,c,alpha,beta,gamma (in
    angstrom and degrees)."""
    import math
    a, b, c = ( cellsg[k] for k in 'abc' )
    ca, cb, cg = ( 0.0 if cellsg[k] == 90 else math.cos(math.radians(cellsg[k]))
                   for k in ('alpha','beta','gamma') )
    g = [ [ a*a, a*b*cg, a*c*cb ],
          [ a*b*cg, b*b, b*c*ca ],
          [ a*c*cb, b*c*ca, c*c ] ]
    #Reciprocal axis lengths from the diagonal of the inverse metric:
    cof = [ g[1][1]*g[2][2] - g[1][2]**2,
            g[0][0]*g[2][2] - g[0][2]**2,
            g[0][0]*g[1][1] - g[0][1]**2 ]
    det = ( g[0][0]*cof[0] - g[0][1]*( g[0][1]*g[2][2] - g[1][2]*g[0][2] )
            + g[0][2]*( g[0][1]*g[1][2] - g[1][1]*g[0][2] ) )
    r = [ math.sqrt( x / det ) for x in cof ]
    u = [ [ aniso['u11'], aniso['u12'], aniso['u13'] ],
          [ aniso['u12'], aniso['u22'], aniso['u23'] ],
          [ aniso['u13'], aniso['u23'], aniso['u33'] ] ]
    return math.fsum( u[i][j] * r[i] * r[j] * g[i][j]
                      for i in range(3) for j in range(3) ) / 3.0

################################################################################
#Tolerances for symmetry detection and for merging site data (positions,
#displacements and occupancies). One relative tolerance, tol, applies to the
#site data and sets the position tolerance (spglib symprec) proportionally:

_TOL_DEFAULT = 1e-3
_TOL_SYMPREC_AT_DEFAULT = 0.01 #Aa (symprec for tol=_TOL_DEFAULT)

def _resolve_tolerances( tol = None, *, tol_symprec = None, tol_uij = None,
                         tol_occu = None ):
    """Tolerances as dict with keys symprec (positions, in Aa), uij
    (displacements, uiso values included) and occu (occupancies). The latter
    two are relative. Values not given follow from tol (default _TOL_DEFAULT):
    uij=occu=tol, and symprec=0.01Aa*tol/_TOL_DEFAULT."""
    def check( name, v, relative ):
        try:
            v = float( v )
        except (TypeError, ValueError):
            pass
        if not ( isinstance( v, float ) and math.isfinite( v ) and v > 0.0 ):
            raise _nc_core.NCBadInput( f'Invalid {name} value: {v!r} (must be'
                                       ' a positive number)' )
        if relative and not v < 1.0:
            raise _nc_core.NCBadInput( f'Invalid {name} value: {v!r} (a'
                                       ' relative tolerance must be less'
                                       ' than 1)' )
        return v
    tol = _TOL_DEFAULT if tol is None else check( 'tol', tol, True )
    def get( name, v, relative, default ):
        return default if v is None else check( name, v, relative )
    return { 'symprec': get( 'tol_symprec', tol_symprec, False,
                             _TOL_SYMPREC_AT_DEFAULT * ( tol / _TOL_DEFAULT ) ),
             'uij': get( 'tol_uij', tol_uij, True, tol ),
             'occu': get( 'tol_occu', tol_occu, True, tol ) }

def _relclose( a, b, tol ):
    """Whether a and b agree within the relative tolerance tol."""
    return abs( a - b ) <= tol * max( abs( a ), abs( b ) )

def _eigclose( e1, e2, tol ):
    """Whether two sets of (sorted) tensor eigenvalues agree within the
    tolerance tol, relative to the largest of them."""
    scale = max( max( abs( x ) for x in e1 ), max( abs( x ) for x in e2 ) )
    return all( abs( x1 - x2 ) <= tol * scale for x1, x2 in zip( e1, e2 ) )

def _group_by_tol( items, same ):
    """Group index of each item. Each item joins the first group whose first
    item (its representative) it is the same as, according to same(rep,item),
    and otherwise starts a new group. Items are only compared with
    representatives, so groups never drift through chains of near-equal
    items."""
    reps, res = [], []
    for it in items:
        for i, r in enumerate( reps ):
            if same( r, it ):
                res.append( i )
                break
        else:
            res.append( len( reps ) )
            reps.append( it )
    return res
