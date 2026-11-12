
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

Internal implementation of the decisions about which symmetry to use for
CIF data, for the cifutils module.

"""

__all__ = []

from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core
from ._cifsitesimpl import _uij_keys


class _SiteView:
    #Stand-in for a gemmi site at another position (and tensor):
    class _Aniso:
        def __init__( self, u ):
            self.nonzero = u is not None
            for k in _uij_keys:
                setattr( self, k, u[k] if u else 0.0 )
    def __init__( self, site, fract, aniso ):
        self.label, self.element = site.label, site.element
        self.occ, self.u_iso = site.occ, site.u_iso
        self.fract, self.aniso = fract, _SiteView._Aniso( aniso )

def _listed_ops( gemmi, cifblock ):
    """The symmetry operations listed in the CIF block: None if there are none,
    ('unparsable', why) if they can not be parsed, ('incomplete', ops) if they
    do not form a group, or ('group', ops, GroupOps)."""
    vals = []
    for tag in ( '_space_group_symop_operation_xyz',
                 '_symmetry_equiv_pos_as_xyz' ):
        vals = [ gemmi.cif.as_string( v ) for v in cifblock.find_values( tag ) ]
        if vals:
            break
    if not vals:
        return None
    try:
        ops = [ gemmi.Op( v ) for v in vals ]
    except (RuntimeError, ValueError) as e:
        return ( 'unparsable', str( e ).strip() )
    #Remove duplicates (as written by pymatgen for some files):
    seen, uops = set(), []
    for o in ops:
        if _op_key( o ) not in seen:
            seen.add( _op_key( o ) )
            uops.append( o )
    ops = uops
    go = gemmi.GroupOps( ops )
    go.add_missing_elements()
    #(GroupOps factorises into rotations and centring translations, so the
    #number of distinct listed operations must be compared with its order):
    if len( ops ) != len( go.sym_ops ) * len( go.cen_ops ):
        return ( 'incomplete', ops )
    return ( 'group', ops, go )

def _op_key( op ):
    #Symmetry operation with translations wrapped into [0,1), for comparisons:
    fs = op.float_seitz()
    return ( tuple( round( fs[i][j], 6 ) for i in range(3) for j in range(3) )
             + tuple( round( fs[i][3] % 1.0, 6 ) % 1.0 for i in range(3) ) )

def _symbol_spacegroup( gemmi, struct, cifblock ):
    """Space group from the symbol (H-M name, or else the number) alone."""
    hm = ' '.join( str( struct.spacegroup_hm ).split() )
    sg = ( gemmi.find_spacegroup_by_name( hm, struct.cell.alpha,
                                          struct.cell.gamma ) if hm else None )
    if not sg:
        for tag in ( '_space_group_IT_number', '_symmetry_Int_Tables_number' ):
            v = [ gemmi.cif.as_string( e ) for e in cifblock.find_values( tag ) ]
            if len( v ) == 1 and v[0].isdigit() and 1 <= int( v[0] ) <= 230:
                sg = gemmi.find_spacegroup_by_number( int( v[0] ) )
                break
    return sg or None

def _cif_multiplicities( gemmi, cifblock ):
    """Site label -> multiplicity, from the CIF block (empty if absent)."""
    lbls = [ gemmi.cif.as_string(v)
             for v in cifblock.find_values( '_atom_site_label' ) ]
    mults = list( cifblock.find_values( '_atom_site_symmetry_multiplicity' ) )
    res = {}
    if len( mults ) == len( lbls ):
        for lbl, m in zip( lbls, mults ):
            try:
                res[lbl] = int( gemmi.cif.as_string( m ) )
            except ValueError:
                pass
    return res

def _check_multiplicities( triples, nops, strict ):
    """Check the multiplicities m of CIF sites, as (label, m, n) with n the
    number of positions generated by the nops symmetry operations (in the
    cell). For the whole file, either m=n for all sites, or (unless strict)
    m*n=nops for all sites (SHELXL writes the order of the site symmetry)."""
    t = [ e for e in triples if e[1] is not None ]
    if not t or all( m == n for _, m, n in t ):
        return
    if not strict and all( m * n == nops for _, m, n in t ):
        return
    lbl, m, n = next( e for e in t if e[1] != e[2] )
    raise _nc_core.NCBadInput(
        f'CIF site "{lbl}" has multiplicity {m}, but its symmetry operations'
        f' give {n} positions' + ( '' if strict else ' (and the multiplicities'
                                   ' are not consistently the orders of the'
                                   ' site symmetries either)' )
        + ' (inconsistent symmetry information in the CIF data?)' )

def _expand_sites_with_ops( gemmi, struct, ops, mult, symprec, strict ):
    """All positions of the sites of struct under the symmetry operations ops,
    with tensors rotated accordingly (merging images within symprec), checked
    against the multiplicities (label -> int, cf. _check_multiplicities)."""
    cell = struct.cell
    cellsg = { k: getattr( cell, k )
               for k in ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' ) }
    sgs = _nc_sgsymimpl
    latt = sgs._cellparams_to_spglib_lattice( cellsg )
    def dist( x, y ):
        d = [ ( a - b + 0.5 ) % 1.0 - 0.5 for a, b in zip( x, y ) ]
        v = [ sum( d[i] * latt[i][k] for i in range(3) ) for k in range(3) ]
        return sum( e * e for e in v ) ** 0.5
    res, triples = [], []
    for site in struct.sites:
        x0 = ( site.fract.x, site.fract.y, site.fract.z )
        u = None
        if site.aniso.nonzero():
            u = { k: float( getattr( site.aniso, k ) ) for k in _uij_keys }
        uf = sgs._cif_to_frac( u, cellsg ) if u else None
        found = []
        for op in ops:
            fs = op.float_seitz()
            R = [ row[:3] for row in fs[:3] ]
            x = [ ( sum( R[i][j] * x0[j] for j in range(3) ) + fs[i][3] ) % 1.0
                  for i in range(3) ]
            if any( dist( x, f ) < symprec for f, _ in found ):
                continue
            ui = u
            if uf is not None and not sgs._np_is_identity( R ):
                ui = sgs._frac_to_cif( sgs._transform_frac_tensor( R, uf ),
                                       cellsg )
            found.append( ( x, ui ) )
        if float( site.occ ) > 0.0:
            triples.append( ( str( site.label ),
                              mult.get( str( site.label ) ), len( found ) ) )
        res += [ _SiteView( site, gemmi.Fractional( *x ), ui )
                 for x, ui in found ]
    _check_multiplicities( triples, len( ops ), strict )
    return res

def _symop_error( why, cifdata, permissive_hint ):
    msg = ( f'The symmetry operations listed in the CIF data {why}. Please use'
            ' a CIF file that either lists all atoms of the unit cell with'
            ' space group P 1, or lists the symmetry operations of a standard'
            ' space group setting consistent with its space group symbol.' )
    if permissive_hint:
        msg += ( ' Alternatively, try permissive mode (at your own risk).' )
    if 'pymatgen' in cifdata.lower():
        msg += ( ' This file was written by pymatgen, as are Materials Project'
                 ' downloads. For Materials Project entries, the "Conventional'
                 ' Standard" CIF download from the website works, or pass the'
                 ' ID as "mpid::<number>" instead of a file name to have'
                 ' NCrystal download that itself.' )
    return _nc_core.NCBadInput( msg )
