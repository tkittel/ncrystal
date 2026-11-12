
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

Internal implementation helpers for the cifutils module.

"""

__all__ = []

import re as _re
from fractions import Fraction as _Fraction

_re_cifnum = _re.compile( r'^([-+]?)(\d*)(?:(\.)(\d*))?'
                          r'(?:[eE]([-+]?\d+))?(\(\d+\))?$' )

def parse_cif_number( s ):
    """Value of a CIF numeric string and its rounding error: the largest
    possible error from the rounding of the printed number, i.e. half the
    step of its last digit (e.g. 0.00005 for "0.0283", 0.5 for "2344" and
    5e-7 for "3400e-6"). A zero without decimals ("0", "-0", "0e-3") is exact
    (rounding error 0.0), since it is in practice a constraint. The standard
    uncertainty in parentheses (e.g. "0.0283(4)") is ignored. Returns None if s
    is not a number (e.g. "?", ".", "1,5" or "inf")."""
    m = _re_cifnum.match( s.strip() )
    if not m:
        return None
    sign, ipart, point, fpart, exp, su = m.groups()
    fpart = fpart or ''
    if not ipart and not fpart:
        return None
    mant = ( ipart or '0' ) + ( '.' + fpart if point else '' )
    value = float( sign + mant + ( 'e' + exp if exp else '' ) )
    if value == 0.0:
        value = 0.0#(no negative zero)
    if not point and su is None and not ipart.strip( '0' ):
        return value, 0.0
    #(exact power of ten, correctly rounded, unlike 10.0**k via libm pow):
    k = ( int( exp ) if exp else 0 ) - len( fpart )
    return value, float( _Fraction( 10 ) ** k / 2 )

_uij_keys = ( 'u11', 'u22', 'u33', 'u12', 'u13', 'u23' )

def _cell_from_cif( gemmi, block ):
    #Cell parameters as (value,rounding_error), or None if not available:
    res = {}
    for k in ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' ):
        tag = ( '_cell_angle_' if k in ( 'alpha', 'beta', 'gamma' )
                else '_cell_length_' ) + k
        s = block.find_value( tag )
        v = parse_cif_number( gemmi.cif.as_string( s ) ) if s else None
        if v is None:
            return None
        res[k] = v
    return res

def uiso_from_cif( gemmi, block ):
    """Positive Uiso values of the sites of a CIF block (gemmi.cif.Block):
    label -> (uiso, rounding_error), cf. parse_cif_number."""
    res = {}
    for row in block.find( '_atom_site_', [ 'label', 'U_iso_or_equiv' ] ):
        v = parse_cif_number( gemmi.cif.as_string( row[1] ) )
        if v is not None and v[0] > 0.0:
            res[ gemmi.cif.as_string( row[0] ) ] = v
    return res

def ueq_from_cif( gemmi, block ):
    """Ueq of the anisotropic tensors of the sites of a CIF block
    (gemmi.cif.Block), computed from the printed tensor components and cell
    parameters: label -> (ueq, rounding_error), where the rounding error is the
    largest error in ueq which the rounding errors of these printed numbers can
    cause (cf. parse_cif_number). Sites or cells with values which can not be
    parsed or used (e.g. an impossible cell, rejected elsewhere) are
    omitted."""
    from . import _sgsymimpl as sgs
    cell = _cell_from_cif( gemmi, block )
    if cell is None:
        return {}
    cellsg = { k: v for k, ( v, _ ) in cell.items() }
    def ueq_or_none( aniso, c ):
        try:
            return sgs._cif_uij_ueq( aniso, c )
        except ( ValueError, ZeroDivisionError ):
            return None
    try:
        g, r = sgs._cell_metric( cellsg )
    except ( ValueError, ZeroDivisionError ):
        return {}
    idx = ( (0,0), (1,1), (2,2), (0,1), (0,2), (1,2) )
    res = {}
    for row in block.find( '_atom_site_aniso_',
                           [ 'label' ] + [ f'U_{k[1:]}' for k in _uij_keys ] ):
        vals = [ parse_cif_number( gemmi.cif.as_string( row[i] ) )
                 for i in range( 1, 7 ) ]
        if any( v is None for v in vals ):
            continue
        aniso = { k: v for k, ( v, _ ) in zip( _uij_keys, vals ) }
        ueq = ueq_or_none( aniso, cellsg )
        if ueq is None:
            continue
        #Rounding errors of the tensor components (Ueq is linear in them):
        err = 0.0
        for ( i, j ), ( _, e ) in zip( idx, vals ):
            f = r[i] * r[j] * g[i][j] * ( 1.0 if i == j else 2.0 ) / 3.0
            err += abs( f ) * e
        #Rounding errors of the cell parameters (estimated numerically):
        for k, ( v, e ) in cell.items():
            if e:
                c2 = dict( cellsg )
                c2[k] = v + e
                u2 = ueq_or_none( aniso, c2 )
                if u2 is None:
                    break
                err += abs( u2 - ueq )
        else:
            res[ gemmi.cif.as_string( row[0] ) ] = ( ueq, err )
    return res
