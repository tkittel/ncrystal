
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
from . import _stablemath as _nc_stablemath
from . import core as _nc_core


def _lattice_params_to_vectors( a, b, c, alpha, beta, gamma ):
    if max(max(alpha,beta),gamma)<=math.pi:
        raise _nc_core.NCBadInput("angles are expected in degrees not radians")
    def _sincos( x ):
        return _nc_stablemath.sin_deg( x ), _nc_stablemath.cos_deg( x )
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
    trackinfo = { 'cell' : spglib_cell, 'same_basis' : True,
                  'transform' : None }
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
             'reldiff' : rd,
             'transform' : trackinfo['transform'] }

def _impl_spglib_refine_cell( spglib_cell, *, symprec, warnings, msgs,
                              trackinfo, nrepeat = 0 ):
    #trackinfo['cell'] is the original cell, which is transformed along with
    #spglib_cell (so it can be compared with the final cell in the same basis),
    #and trackinfo['transform'] the total (P,p) applied to it so far.
    spglib = _import_spglib()

    orig_cell = spglib_cell
    if trackinfo['cell'] is not None:
        tr = _spglib_std_transform( orig_cell, symprec )
        if tr is None:
            trackinfo['cell'] = None
            trackinfo['transform'] = None
        else:
            trackinfo['cell'] = _spglib_apply_std_transform( trackinfo['cell'],
                                                             tr )
            t0 = trackinfo['transform'] or _identity_transform()
            trackinfo['transform'] = _compose_transforms( t0, tr[:2] )
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
    #Snap numerical noise near the cell edges (e.g. 7e-17 or 1-1e-16) to 0:
    return 0.0 if ( abs( x ) < 1e-12 or abs( x - 1.0 ) < 1e-12 ) else x

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

def _cell_metric( cellsg ):
    """Metric tensor (3x3 nested lists) and reciprocal axis lengths of the cell
    given by a dict with keys a,b,c,alpha,beta,gamma (angstrom and degrees)."""
    a, b, c = ( cellsg[k] for k in 'abc' )
    ca, cb, cg = ( _nc_stablemath.cos_deg( cellsg[k] )
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
    return g, [ math.sqrt( x / det ) for x in cof ]

def _cif_uij_ueq( aniso, cellsg ):
    """Equivalent isotropic displacement, Ueq = trace(U_cart)/3, of CIF uij
    values (dict with keys u11,...,u23, relative to the normalised reciprocal
    axes), in the cell given by a dict with keys a,b,c,alpha,beta,gamma (in
    angstrom and degrees)."""
    g, r = _cell_metric( cellsg )
    u = [ [ aniso['u11'], aniso['u12'], aniso['u13'] ],
          [ aniso['u12'], aniso['u22'], aniso['u23'] ],
          [ aniso['u13'], aniso['u23'], aniso['u33'] ] ]
    return math.fsum( u[i][j] * r[i] * r[j] * g[i][j]
                      for i in range(3) for j in range(3) ) / 3.0

def _rescale_uij( aniso, ueq, target ):
    """CIF uij values (dict) with Ueq=ueq rescaled to Ueq=target (exactly, for
    cross-platform reproducibility). Returned unchanged if ueq==target (or ueq
    is not positive, leaving the error for later)."""
    if not ueq > 0.0 or ueq == target:
        return aniso
    keys = list( aniso )
    return dict( zip( keys, _nc_stablemath.scale( [ aniso[k] for k in keys ],
                                                  target, ueq ) ) )

################################################################################
#Basis changes and displacement tensors. A transform is a pair (P,p) of numpy
#arrays acting on fractional coordinates as x' = P.x + p (as from spglib). A
#displacement tensor in fractional form, U_frac = N.U_cif.N with
#N=diag(a*,b*,c*), transforms as U_frac' = M.U_frac.M^T for both basis changes
#(M=P) and symmetry operations (M=R).

def _identity_transform():
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    return np.identity( 3 ), np.zeros( 3 )

def _compose_transforms( t1, t2 ):
    """The transform applying t1 and then t2."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    P1, p1 = ( np.asarray( e, dtype = float ) for e in t1 )
    P2, p2 = ( np.asarray( e, dtype = float ) for e in t2 )
    return np.dot( P2, P1 ), np.dot( P2, p1 ) + p2

def _transform_position( t, x ):
    """Fractional position x mapped with the transform t (not wrapped)."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    return np.dot( t[0], np.asarray( x, dtype = float ) ) + t[1]

def _cif_to_frac( aniso, cellsg ):
    """Fractional form (3x3 numpy array) of CIF uij values (dict with keys
    u11,...,u23) in the given cell."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    _, r = _cell_metric( cellsg )
    u = np.array( [ [ aniso['u11'], aniso['u12'], aniso['u13'] ],
                    [ aniso['u12'], aniso['u22'], aniso['u23'] ],
                    [ aniso['u13'], aniso['u23'], aniso['u33'] ] ],
                  dtype = float )
    return np.outer( r, r ) * u

def _frac_to_cif( ufrac, cellsg ):
    """CIF uij values (dict with keys u11,...,u23) of the fractional form of a
    tensor in the given cell."""
    _, r = _cell_metric( cellsg )
    def f( i, j ):
        return float( ufrac[i][j] ) / ( r[i] * r[j] )
    return { 'u11': f(0,0), 'u22': f(1,1), 'u33': f(2,2),
             'u12': f(0,1), 'u13': f(0,2), 'u23': f(1,2) }

def _transform_frac_tensor( m, ufrac ):
    """M.U_frac.M^T, for a basis change (M=P) or a symmetry operation (M=R).
    Computed exactly (cross-platform reproducible), with the entries of M
    snapped to simple rationals (spglib gives them as floats)."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    sm = _nc_stablemath
    m = [ [ sm.snap_rational( float( v ) ) for v in row ] for row in m ]
    u = [ [ float( v ) for v in row ] for row in ufrac ]
    return np.asarray( sm.transform_sym( m, u ), dtype = float )

def _cart_eigenvalues( ufrac, cellsg ):
    """Sorted eigenvalues of the Cartesian form of a tensor, i.e. of
    U_frac.G (G being the metric). These do not depend on the basis."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    g, _ = _cell_metric( cellsg )
    ev = np.linalg.eigvals( np.dot( np.asarray( ufrac, dtype = float ),
                                    np.asarray( g, dtype = float ) ) )
    return tuple( sorted( float( e.real ) for e in ev ) )

################################################################################
#Tolerances for symmetry detection and for merging site data (positions,
#displacements and occupancies). One relative tolerance, tol, applies to the
#site data and sets the position tolerance (spglib symprec) proportionally:

_TOL_DEFAULT = 1e-3
_TOL_SYMPREC_AT_DEFAULT = 0.01 #Aa (symprec for tol=_TOL_DEFAULT)
_TOL_MAX = 0.05 #max relative tolerance
_TOL_SYMPREC_MAX = 0.5 #Aa
#The largest tol must give the largest symprec:
assert abs( _TOL_SYMPREC_AT_DEFAULT * _TOL_MAX / _TOL_DEFAULT
            - _TOL_SYMPREC_MAX ) < 1e-12

def _resolve_tolerances( tol = None, *, tol_symprec = None, tol_uij = None,
                         tol_occu = None ):
    """Tolerances as dict with keys symprec (positions, in Aa), uij
    (displacements, uiso values included) and occu (occupancies). The latter
    two are relative. Values not given follow from tol (default _TOL_DEFAULT):
    uij=occu=tol, and symprec=0.01Aa*tol/_TOL_DEFAULT."""
    def check( name, v, vmax, unit ):
        try:
            v = float( v )
        except (TypeError, ValueError):
            pass
        if not ( isinstance( v, float ) and math.isfinite( v ) and v > 0.0 ):
            raise _nc_core.NCBadInput( f'Invalid {name} value: {v!r} (must be'
                                       ' a positive number)' )
        if not v <= vmax:
            raise _nc_core.NCBadInput( f'Invalid {name} value: {v!r} (must'
                                       f' not exceed {vmax:g}{unit})' )
        return v
    def rel( name, v, default ):
        return default if v is None else check( name, v, _TOL_MAX, '' )
    tol = rel( 'tol', tol, _TOL_DEFAULT )
    return { 'symprec': ( _TOL_SYMPREC_AT_DEFAULT * ( tol / _TOL_DEFAULT )
                          if tol_symprec is None
                          else check( 'tol_symprec', tol_symprec,
                                      _TOL_SYMPREC_MAX, 'Aa' ) ),
             'uij': rel( 'tol_uij', tol_uij, tol ),
             'occu': rel( 'tol_occu', tol_occu, tol ) }

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

################################################################################
#Orbits of sites with symmetry detection based on positions and site data.

def _periodic_distsq_fn( lattice ):
    """Function giving the squared Cartesian distance between two fractional
    positions, using the nearest periodic image along each axis (lattice: numpy
    array with the cell vectors as rows)."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    def distsq( x, y ):
        d = np.asarray( x, dtype = float ) - np.asarray( y, dtype = float )
        d = np.dot( lattice.T, d - np.round( d ) )
        return float( np.dot( d, d ) )
    return distsq

def _find_site_orbits( cellsg, sites, tols ):
    """Orbits of sites under the symmetry found by spglib (in the basis of the
    input cell, given as dict with keys a,b,c,alpha,beta,gamma). Each site is a
    dict with keys 'key' (exact, hashable, e.g. the element), 'pos'
    (fractional), and the optional 'occu', 'uiso' and 'aniso' (CIF uij dict).
    Sites are only equivalent if their keys are equal and their occupancies,
    uiso values and tensor eigenvalues agree within the tolerances (dict from
    _resolve_tolerances). Returns dict with the space group number, and the
    orbits as list of dicts with 'members' (site indices, the first being the
    representative) and 'aniso_ok' (False if a tensor is not the symmetry
    image of the representative tensor within tols['uij'], None if the orbit
    has no tensors)."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    spglib = _import_spglib()
    lattice = np.asarray( _cellparams_to_spglib_lattice( cellsg ),
                          dtype = float )
    ufrac = [ ( _cif_to_frac( s['aniso'], cellsg ) if s.get('aniso') else None )
              for s in sites ]
    eigs = [ ( _cart_eigenvalues( u, cellsg ) if u is not None else None )
             for u in ufrac ]
    tu, to = tols['uij'], tols['occu']
    def opt_same( a, b, cmp ):
        return ( a is None ) == ( b is None ) and ( a is None or cmp( a, b ) )
    def same( i, j ):
        si, sj = sites[i], sites[j]
        return ( si['key'] == sj['key']
                 and opt_same( si.get('occu'), sj.get('occu'),
                               lambda a, b : _relclose( a, b, to ) )
                 and opt_same( si.get('uiso'), sj.get('uiso'),
                               lambda a, b : _relclose( a, b, tu ) )
                 and opt_same( eigs[i], eigs[j],
                               lambda a, b : _eigclose( a, b, tu ) ) )
    types = _group_by_tol( list( range( len( sites ) ) ), same )
    positions = [ tuple( float( x ) for x in s['pos'] ) for s in sites ]
    ds = spglib.get_symmetry_dataset( ( lattice, positions, types ),
                                      symprec = tols['symprec'] )
    if not ds:
        raise _nc_core.NCBadInput( 'Could not determine symmetry of'
                                   ' structure with spglib' )
    def get( k ):
        return getattr( ds, k ) if hasattr( ds, k ) else ds[k]
    rots = np.asarray( get('rotations'), dtype = float )
    trans = np.asarray( get('translations'), dtype = float )
    orbits = {}
    for i, e in enumerate( get('equivalent_atoms') ):
        orbits.setdefault( int( e ), [] ).append( i )
    orbitlist = sorted( orbits.values(), key = lambda m : m[0] )
    distsq = _periodic_distsq_fn( lattice )
    def deviation( r, i ):
        #Smallest deviation of tensor i from the images of tensor r under the
        #operations mapping position r onto position i, as the largest
        #(absolute) Cartesian eigenvalue of the difference relative to the
        #largest eigenvalue of tensor r. Basis independent, and bounding the
        #eigenvalue differences.
        scale = max( abs( e ) for e in eigs[r] )
        def dev( R ):
            d = _transform_frac_tensor( R, ufrac[r] ) - ufrac[i]
            return max( abs( e ) for e in _cart_eigenvalues( d, cellsg ) )
        devs, dmin = [], None
        for R, t in zip( rots, trans ):
            dsq = distsq( np.dot( R, positions[r] ) + t, positions[i] )
            if dmin is None or dsq < dmin[0]:
                dmin = ( dsq, R )
            if dsq <= ( 2.0 * tols['symprec'] )**2:
                devs.append( dev( R ) )
        if not devs:
            #Should not happen, but use the nearest image just in case:
            devs.append( dev( dmin[1] ) )
        return min( devs ) / scale
    res = []
    for members in orbitlist:
        r = members[0]
        aniso_ok = None
        if ufrac[r] is not None:
            aniso_ok = all( deviation( r, i ) <= tu for i in members[1:] )
        res.append( { 'members': members, 'aniso_ok': aniso_ok } )
    return { 'sgno': int( get('number') ), 'orbits': res, 'types': types }

def _parse_xyz_op( op ):
    """Rotation (3x3 numpy array) and translation of a symmetry operation in
    the form produced by NCrystal (e.g. "-y+1/2,x-y,z+1/3")."""
    from fractions import Fraction

    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    R, t = np.zeros( ( 3, 3 ) ), np.zeros( 3 )
    comps = op.replace( ' ', '' ).split( ',' )
    if len( comps ) != 3:
        raise _nc_core.NCBadInput( f'Invalid symmetry operation: "{op}"' )
    for i, comp in enumerate( comps ):
        for term in comp.replace( '-', '+-' ).split( '+' ):
            if not term:
                continue
            sign = -1.0 if term.startswith( '-' ) else 1.0
            term = term.lstrip( '-' )
            if term in ( 'x', 'y', 'z' ):
                R[i][ 'xyz'.index( term ) ] += sign
            else:
                t[i] += sign * float( Fraction( term ) )
    return R, t

def _expand_sites( ops, cellsg, sites ):
    """All positions of sites, with their tensors rotated accordingly. The ops
    are strings as from _parse_xyz_op, and each site a dict with 'pos' and
    optionally 'aniso' (CIF uij dict). Returns list of dicts like the sites
    (other keys copied), one per position."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    ops = [ _parse_xyz_op( o ) for o in ops ]
    res = []
    for s in sites:
        x0 = np.asarray( s['pos'], dtype = float )
        uf = _cif_to_frac( s['aniso'], cellsg ) if s.get('aniso') else None
        found = []
        for R, t in ops:
            x = ( np.dot( R, x0 ) + t ) % 1.0
            if any( _unit_cell_point_distsq( x, f ) < 1e-12 for f in found ):
                continue
            found.append( x )
            e = dict( s )
            e['pos'] = tuple( _remap_fract_pos( v ) for v in x )
            if uf is not None:
                e['aniso'] = _frac_to_cif( _transform_frac_tensor( R, uf ),
                                           cellsg )
            res.append( e )
    return res

def _refine_sites( cellsg, sites, tols ):
    """Refine a structure with site data (sites as for _find_site_orbits, one
    per position) with spglib. Returns dict with the refined cell parameters
    ('cellsg', with 'spacegroup' number), spglib messages ('warnings',
    'msgs'), the space group symbol ('sgsymb_hm'), all refined positions
    ('allpos' as (site index of type representative, pos)), the site 'types',
    and one entry per
    orbit in 'sites': dicts with the representative 'index', the 'members'
    (site indices), its refined 'pos', and the 'aniso' tensor
    mapped to the refined cell (None if the orbit has no tensors, or if
    they failed the orientation check, as indicated by 'aniso_dropped').
    The 'sites' are None if spglib's basis change could not be tracked."""
    from ._numpy import _ensure_numpy
    np = _ensure_numpy()
    orb = _find_site_orbits( cellsg, sites, tols )
    types = orb['types']
    lattice = _cellparams_to_spglib_lattice( cellsg )
    cell = ( lattice, [ tuple( float( v ) for v in s['pos'] ) for s in sites ],
             types )
    d = _spglib_refine_cell( cell, symprec = tols['symprec'] )
    newcell = dict( d['cellparams_snapped'] )
    newcell['spacegroup'] = d['sgno']
    tr = d['transform']
    newlatt = np.asarray( _cellparams_to_spglib_lattice( newcell ),
                          dtype = float )
    refpos = list( zip( d['refined_cell'][1], d['refined_cell'][2] ) )
    distsq = _periodic_distsq_fn( newlatt )
    #Representative index of each type (orbits might share a type, but then
    #their positions differ, which is used below):
    type2rep = {}
    for o in orb['orbits']:
        type2rep.setdefault( types[ o['members'][0] ], o['members'][0] )
    allpos = [ ( type2rep[ typ ], tuple( _remap_fract_pos( v ) for v in x ) )
               for x, typ in refpos ]
    if tr is None:
        return { 'cellsg': newcell, 'sgsymb_hm': d['sgsymb_hm'],
                 'warnings': d['warnings'], 'msgs': d['msgs'],
                 'allpos': allpos, 'types': types, 'sites': None }
    res = []
    for o in orb['orbits']:
        r = o['members'][0]
        #The refined position of the representative is the one of its type
        #nearest to its transformed input position:
        xr = _transform_position( tr, sites[r]['pos'] )
        cands = [ x for x, typ in refpos if typ == types[r] ]
        pos = min( cands, key = lambda x : distsq( x, xr ) )
        e = { 'index': r, 'members': o['members'],
              'pos': tuple( _remap_fract_pos( v ) for v in pos ),
              'aniso': None, 'aniso_dropped': o['aniso_ok'] is False }
        if o['aniso_ok'] and _np_is_identity( tr[0] ):
            e['aniso'] = dict( sites[r]['aniso'] )
        elif o['aniso_ok']:
            uf = _transform_frac_tensor( tr[0], _cif_to_frac(
                sites[r]['aniso'], cellsg ) )
            u = _frac_to_cif( uf, newcell )
            #Remove numerical noise (from components which should be zero):
            umax = max( abs( v ) for v in u.values() )
            e['aniso'] = { k: ( 0.0 if abs( v ) < 1e-12 * umax else v )
                           for k, v in u.items() }
        res.append( e )
    return { 'cellsg': newcell, 'sgsymb_hm': d['sgsymb_hm'],
             'warnings': d['warnings'], 'msgs': d['msgs'],
             'allpos': allpos, 'types': types, 'sites': res }
