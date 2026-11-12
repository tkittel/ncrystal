
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

Internal implementation of reading CIF data with gemmi for the cifutils
module (including the numbers of the CIF data with their rounding errors).

"""

__all__ = []

import math as _math
import re as _re
from fractions import Fraction as _Fraction

from . import _common as _nc_common
from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core
from ._cifsrcimpl import _codid2url, _mpid2url, mp_origin_note


def _extract_descr_from_cif( raw_cif_dict, cifsrc, ciftextdata ):

    def normalise_result( _l ):
        #split any newlines and trim whitespace before returning:
        res = []
        for e in _l:
            for ee in e.splitlines():
                res.append( ee.rstrip() )
        return res

    if not raw_cif_dict:
        _nc_common.warn('Could not properly decode CIF metadata due to inability to access raw CIF data as dictionary.')
        return [], None, None, None

    #cif is cif section with publication info:
    ll=[]
    def extract(key,*altkeys,expectlist = False):
        if not any( key in s for _,s in sorted(raw_cif_dict.items())):
            if altkeys:
                return extract( altkeys[0], *altkeys[1:], expectlist = expectlist )
            return [] if expectlist else ''
        s = next( d[key] for _,d in sorted(raw_cif_dict.items()) if key in d )
        if not expectlist:
            if s is None:
                return ''
            if not hasattr(s,'strip') and hasattr(s,'__len__') and len(s)==1:
                s=s[0]#convert single element list to first element
            s=str(s).strip()
            return '' if s in ('?','.') else s
        else:
            _l=[ ('' if (e is None or (hasattr(e,'strip') and e.strip()=='?')) else str(e).strip()) for e in s ]
            return [ e for e in _l if e ]

    title = extract('_publ_section_title','_citation_title')
    authors = extract('_publ_author_name','_citation_author_name',expectlist=True)
    journalname = extract('_journal_name_full','_citation_journal_full')
    year = extract('_journal_year','_citation_journal_year','_citation_year')
    doi = extract('_journal_paper_doi')
    chemformsum = extract('_chemical_formula_sum','_cod_original_formula_sum')

    extracted_codid = extract('_cod_database_code')
    if extracted_codid:
        extracted_codid = int(extracted_codid) if extracted_codid.isdigit() else None

    extracted_mpid = None
    _sp = mp_origin_note
    for e in ciftextdata.splitlines():
        if _sp in e:
            _ = e.split(_sp,1)
            if len(_) >= 2:
                _ = _[1].split()[0].split('#',1)[0]
                if _.isdigit():
                    extracted_mpid = int(_)

    if extracted_codid is not None and cifsrc.codid is not None and cifsrc.codid != extracted_codid:
        _nc_common.warn(f'Embedded CODID {extracted_codid} is different than requested {cifsrc.codid}!')

    if extracted_mpid is not None and cifsrc.mpid is not None and cifsrc.mpid != extracted_mpid:
        _nc_common.warn(f'Embedded MPID {extracted_mpid} is different than requested {cifsrc.mpid}!')

    _thecodid = cifsrc.codid or extracted_codid
    _thempid = cifsrc.mpid or extracted_mpid

    audit_creation_method = extract('_audit_creation_method')
    audit_creation_date = extract('_audit_creation_date')
    audit_author_name = extract('_audit_author_name')

    nauthors = len(authors or [])
    if authors:
        if len(authors)==1:
            authors = authors[0]
        elif len(authors)==2:
            authors = f'{authors[0]} and {authors[1]}'
        else:
            authors = f'{authors[0]}, et al.'
    if title:
        ll.append(f'"{title}"')
    if  (journalname and year) and not authors:
        authors = '<unknown authors>'
    if journalname and year:
        _jy = f'{journalname}, {year}'
    elif journalname or year:
        _jy = journalname or year
    else:
        _jy = ''
    if authors:
        ll.append(f'{authors} [{_jy}]' if _jy else f'Author{"" if nauthors==1 else "s"}: {authors}')

    if doi:
        ll.append(f'DOI: https://dx.doi.org/{doi}')

    if audit_creation_method:
        ll.append(f'CIF creation method: {audit_creation_method}')
    if audit_creation_date:
        ll.append(f'CIF creation date: {audit_creation_date}')
    if audit_author_name:
        ll.append(f'CIF created by: {audit_author_name}')

    if _thecodid:
        ll += [ f'Crystallography Open Database entry {_thecodid}', _codid2url(_thecodid) ]
    if _thempid:
        ll += [ 'The Materials Project', _mpid2url(_thempid) ]

    return normalise_result( ll ), _thecodid, _thempid, chemformsum

_import_gemmi_cache = [None, None]

def _import_gemmi( *, sysexit = False ):
    if _import_gemmi_cache[0] is not None:
        return _import_gemmi_cache[0], _import_gemmi_cache[1]
    try:
        import gemmi  #both available on pypi and conda-forge
        import gemmi.cif
    except ImportError:
        m = ( 'Could not import gemmi modules needed to process CIF files.'
              +' The gemmi package is available on both PyPI ("python3 -mpip install'
              +' gemmi") and conda ("conda install -c conda-forge gemmi")' )
        if sysexit:
            raise SystemExit(m)
        else:
            raise ImportError(m)
    _import_gemmi_cache[0] = gemmi
    _import_gemmi_cache[1] = gemmi.cif
    return gemmi, gemmi.cif

_guessmap = [None]

def _guess_spacegroup_name( gemmi, s ):
    def _guess_keys( s ):
        s=''.join(s.split())
        return [ s, s.replace('/',''),s.replace('-',''),s.replace('-','').replace('/','' ) ]
    def _init_guess_map():
        guess = {}
        def _ag( i, s ):
            if s not in guess:
                guess[s] = {i}
            else:
                guess[s].add(i)
        def add_guess( i, s ):
            for k in _guess_keys(s):
                _ag( i, k )
        for i in range(1,230+1):
            sg = gemmi.find_spacegroup_by_number(i)
            add_guess( i, sg.xhm() )
            add_guess( i, sg.short_name() )
            add_guess( i, sg.hall )
            add_guess( i, sg.hm )
        return guess
    if not _guessmap[0]:
        _guessmap[0] = _init_guess_map()
    gm = _guessmap[0]
    possible_sgnos = set()
    for k in _guess_keys(s):
        for i in gm.get(k,[]):
            possible_sgnos.add( i )
    return sorted(possible_sgnos)

def _load_with_gemmi( cifblock, allow_fixup = True ):
    #Load cifblock into gemmi struct. We might perform in-place editing of the
    #cifblock, if gemmi does not immediately recognise the space group.
    gemmi, _gemmi_cif = _import_gemmi()

    struct = gemmi.make_small_structure_from_block( cifblock )

    assert struct
    if hasattr(struct,'find_spacegroup'):
        sg = struct.find_spacegroup()
    else:
        sg = struct.spacegroup

    if not sg:
        #Doing what the old struct.find_spacegroup() was doing:
        sg = gemmi.find_spacegroup_by_name( hm=struct.spacegroup_hm,
                                            alpha=struct.cell.alpha,
                                            gamma=struct.cell.gamma)

    if sg or not allow_fixup:
        return struct, sg

    #Let us see if we can find the spacegroup with a bit of manual
    #intervention. If we can, we update the cifblock and do a full reload
    #(otherwise the gemmi struct object won't have the spacegroup and images
    #properly initialised):

    _ = cifblock.find(['_space_group_IT_number'])
    if not _:
        _ = cifblock.find(['_symmetry_Int_Tables_number'])

    if _ and len(_)==1 and len(_[0])==1 and str(_[0][0]).isdigit() and (1<=int(str(_[0][0]))<=230):
        _sgnum = int(_[0][0])
        #NB: We can do this for some groups only (cubic?) Or perhaps we could use it as a double-check only?
        sg = gemmi.find_spacegroup_by_number(_sgnum)
        assert sg and sg.number == _sgnum

    sg_hm = str(struct.spacegroup_hm).strip()
    if sg_hm and not sg:
        def sg_searchname(x):
            return gemmi.find_spacegroup_by_name(x,
                                                 struct.cell.alpha,
                                                 struct.cell.gamma )
        attempts = [ sg_hm ]
        #TODO: We used to do this, but currently it does not make a difference:
        #for e in ('H','R'):#NB H+R is not enough, see table 6 at http://cci.lbl.gov/sginfo/hall_symbols.html
        #    if sg_hm.endswith(e):
        #        attempts.append( sg_hm[:-1] + ' : ' + sg_hm[-1] )
        for a in attempts:
            a = ' '.join(a.split())
            sg = sg_searchname( a )
            if sg:
                _nc_common.warn(f'Had to interpret spacegroup "{sg_hm}" as "{a}" before Gemmi could recognise it.')
                break

    if sg_hm and not sg:
        possible = _guess_spacegroup_name( gemmi, sg_hm )
        if len(possible) == 1:
            sgnum = possible[0]
            sg = gemmi.find_spacegroup_by_number(sgnum)
            _nc_common.warn(f'Had to interpret spacegroup "{sg_hm}" as "{sg.hm}" before Gemmi could recognise it.')
        elif len(possible) > 1:
            poshm = [ (sgnum,gemmi.find_spacegroup_by_number(sgnum).hm) for sgnum in possible ]
            poshm = ', '.join( f'"{sgstr}"(number {sgnum})' for sgnum,sgstr in poshm )
            _nc_common.warn(f'Failed to interpret spacegroup "{sg_hm}". Could be any of: {poshm}.')
    if not sg:
        #manual intervention did not help:
        return struct, None

    #patch up the cifblock with info about manually found spacegroup, and try again:
    _update_spacegroup_in_cifblock( cifblock, sg )

    return _load_with_gemmi( cifblock, allow_fixup = False )

def _update_spacegroup_in_cifblock_from_name( gemmi, cifblock, name ):
    sg = gemmi.find_spacegroup_by_name( name )
    if not sg:
        raise _nc_core.NCBadInput(f'Unknown spacegroup: "{name}"')
    _update_spacegroup_in_cifblock( cifblock, sg )

def _update_spacegroup_in_cifblock( cifblock, sg, use_xhm = True ):
    cifblock.set_pair('_space_group_name_H-M_alt', sg.xhm() if use_xhm else sg.hm )
    cifblock.set_pair('_space_group_IT_number', str(sg.number) )
    cifblock.set_pair('_space_group_name_Hall', sg.hall )
    cifblock.set_pair('_space_group_crystal_system', sg.crystal_system_str() )

_re_cifnum = _re.compile( r'^([-+]?)(\d*)(?:(\.)(\d*))?'
                          r'(?:[eE]([-+]?\d+))?(\(\d+\))?$' )

def parse_cif_number( s ):
    """Value of a CIF numeric string and its rounding error: the largest
    possible error from the rounding of the printed number, i.e. half the
    step of its last digit (e.g. 0.00005 for "0.0283", 0.5 for "2344" and
    5e-7 for "3400e-6"). A zero without decimals ("0", "-0", "0e-3") is exact
    (rounding error 0.0), since it is in practice a constraint. The standard
    uncertainty in parentheses (e.g. "0.0283(4)") is ignored. Returns None if s
    is not a number (e.g. "?", ".", "1,5" or "inf"), or if the value or rounding
    error is not a finite double (e.g. "1e999" or "0.0e310")."""
    m = _re_cifnum.match( s.strip() )
    if not m:
        return None
    sign, ipart, point, fpart, exp, su = m.groups()
    fpart = fpart or ''
    if not ipart and not fpart:
        return None
    mant = ( ipart or '0' ) + ( '.' + fpart if point else '' )
    value = float( sign + mant + ( 'e' + exp if exp else '' ) )
    if not _math.isfinite( value ):
        return None
    if value == 0.0:
        value = 0.0#(no negative zero)
    if not point and su is None and not ipart.strip( '0' ):
        return value, 0.0
    #(exact power of ten, correctly rounded, unlike 10.0**k via libm pow):
    e = 0
    if exp:
        #(huge exponents are clamped, avoiding the digit limit of int()):
        e = int( exp ) if len( exp ) < 12 else ( -1, 1 )[exp[0] != '-'] * 10**12
    k = e - len( fpart )
    if k < -400:
        return value, 0.0#(below the smallest double)
    if k > 308:
        return None#(e.g. "0.0e310": rounding error not a finite double)
    return value, float( _Fraction( 10 ) ** k / 2 )


def _cell_from_cif( gemmi, block ):
    #Cell parameters as (value,rounding_error), or None if not available:
    res = {}
    for k in ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' ):
        tag = ( '_cell_angle_' if k in ( 'alpha', 'beta', 'gamma' )
                else '_cell_length_' ) + k
        s = block.find_value( tag )
        v = parse_cif_number( gemmi.cif.as_string( s ) ) if s else None
        if v is None and k in ( 'alpha', 'beta', 'gamma' ) and s is None:
            v = ( 90.0, 0.0 )#(missing angles are 90 degrees, as in gemmi)
        if v is None:
            return None
        res[k] = v
    return res

def uiso_from_cif( gemmi, block ):
    """Positive Uiso values of the sites of a CIF block (gemmi.cif.Block):
    label -> (uiso, rounding_error), cf. parse_cif_number. Biso values
    (_atom_site_B_iso_or_equiv) are used for sites without a Uiso value, as
    Uiso = Biso/(8 pi^2), as gemmi does."""
    res = {}
    f = 1.0 / ( 8.0 * _math.pi * _math.pi )
    for row in block.find( '_atom_site_', [ 'label', '?U_iso_or_equiv',
                                            '?B_iso_or_equiv' ] ):
        v = ( parse_cif_number( gemmi.cif.as_string( row[1] ) )
              if row.has( 1 ) else None )
        if ( v is None or not v[0] > 0.0 ) and row.has( 2 ):
            b = parse_cif_number( gemmi.cif.as_string( row[2] ) )
            v = ( b[0] * f, b[1] * f ) if b is not None else None
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
    sgs = _nc_sgsymimpl
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
                           [ 'label' ] + [ f'U_{k[1:]}' for k in _nc_sgsymimpl._uij_keys ] ):
        vals = [ parse_cif_number( gemmi.cif.as_string( row[i] ) )
                 for i in range( 1, 7 ) ]
        if any( v is None for v in vals ):
            continue
        aniso = { k: v for k, ( v, _ ) in zip( _nc_sgsymimpl._uij_keys, vals ) }
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

def _plain_site( site ):
    #A gemmi small structure site as a plain dict:
    fract = ( site.fract.x, site.fract.y, site.fract.z )
    if site.occ > 0.0 and not all( _math.isfinite( x ) for x in fract ):
        #(e.g. a site with "?" for its coordinates; unoccupied sites are
        #ignored later):
        raise _nc_core.NCBadInput( f'CIF site "{site.label}" has missing or'
                                   ' invalid fractional coordinates' )
    aniso = { k: float( getattr( site.aniso, k ) ) for k in _nc_sgsymimpl._uij_keys }
    if all( v == 0.0 for v in aniso.values() ):
        aniso = None
    en, Z = str( site.element.name ).strip(), int( site.element.atomic_number )
    if en == 'X' and _re.sub( r'[^A-Za-z]', '', site.type_symbol ) == 'T':
        en, Z = 'T', 1#(tritium, unknown to gemmi)
    return { 'label': str( site.label ),
             'element_name': en,
             'element_Z': Z,
             'fract': fract,
             'occ': float( site.occ ),
             'u_iso': float( site.u_iso ),
             'aniso': aniso }

def read_cif( cifdata, *, quiet, override_spacegroup ):
    """Read CIF data (text) with gemmi, into the plain content needed for
    loading it. Returns a dict with:
      'cif_raw': the document as a dict (None if not decodable),
      'cell': dict with keys a,b,c,alpha,beta,gamma,
      'sites': the sites as dicts with keys label, element_name, element_Z,
               fract, occ, u_iso and aniso (CIF uij dict or None),
      'listed_ops': the listed symmetry operation strings (None if absent,
                    or if the space group is overridden),
      'symbol_hm', 'symbol_itnumbers': the space group symbol and the values
                    of the IT number fields (per tag),
      'orig_hm_alt': the original _space_group_name_H-M_alt value (or None),
      'multiplicities': site label -> multiplicity,
      'atom_type_number_in_cell': the raw values of that field,
      'uiso_printed', 'ueq_printed': cf. uiso_from_cif and ueq_from_cif.
    The gemmi module, structure and space group are also included (keys
    '_gemmi', '_struct' and '_sg') for the symmetry decisions."""
    from itertools import chain
    gemmi, gemmi_cif = _import_gemmi()
    if not quiet:
        _nc_common.print("Attempting to load CIF data with gemmi")
    try:
        cif_doc = gemmi_cif.read_string( cifdata )
    except ValueError as e:
        raise _nc_core.NCBadInput(f'CIF parsing error (from Gemmi): "{e}"')

    cif_doc_as_json = cif_doc.as_json()
    import json
    try:
        cif_doc_as_dict = json.loads( cif_doc_as_json )
    except json.decoder.JSONDecodeError as e:
        cif_doc_as_dict = None
        _nc_common.warn(f'Could not decode raw CIF data to dictionary (Gemmi bug?): JSONDecodeError("{e}")')

    if len(cif_doc) == 0:
        raise _nc_core.NCBadInput('CIF data had no blocks!')
    if len(cif_doc) == 1:
        blk = cif_doc.sole_block()
    else:
        blk = None
        for e in cif_doc:
            if not len(e.find(['_atom_site_fract_x'])):
                continue
            if blk:
                blk = None
                break
            blk = e
        if not blk:
            raise _nc_core.NCBadInput('Could not automatically determine block in CIF data with structure info')

    if override_spacegroup:
        assert isinstance(override_spacegroup,str)
        _nc_common.warn(f'Overriding spacegroup to "{override_spacegroup}" due to explicit request')
        _update_spacegroup_in_cifblock_from_name( gemmi, blk, override_spacegroup )

    #record original h-m-alt entry:
    orig_hm_alt = None
    _ = blk.find(['_space_group_name_H-M_alt'])
    if _ and len(_)==1 and len(_[0])==1:
        orig_hm_alt = gemmi.cif.as_string(_[0][0]).strip()#(without quotes)

    struct, sg = _load_with_gemmi( blk )

    if not struct:
        raise _nc_core.NCBadInput('Could not load structure with Gemmi')

    listed = None
    if not override_spacegroup:
        for tag in ( '_space_group_symop_operation_xyz',
                     '_symmetry_equiv_pos_as_xyz' ):
            listed = [ gemmi.cif.as_string( v )
                       for v in blk.find_values( tag ) ]
            if listed:
                break
        listed = listed or None
    _ = blk.find(['_atom_type_number_in_cell'])
    return { 'cif_raw': cif_doc_as_dict,
             'cell': { k: getattr( struct.cell, k )
                       for k in ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' ) },
             'sites': [ _plain_site( s ) for s in struct.sites ],
             'listed_ops': listed,
             'symbol_hm': ' '.join( str( struct.spacegroup_hm ).split() ),
             'symbol_itnumbers': [ [ gemmi.cif.as_string( e )
                                     for e in blk.find_values( tag ) ]
                                   for tag in ( '_space_group_IT_number',
                                                '_symmetry_Int_Tables_number' ) ],
             'orig_hm_alt': orig_hm_alt,
             'multiplicities': _cif_multiplicities( gemmi, blk ),
             'atom_type_number_in_cell': ( list(chain.from_iterable(_))
                                           if _ else [] ),
             'uiso_printed': uiso_from_cif( gemmi, blk ),
             'ueq_printed': ueq_from_cif( gemmi, blk ),
             '_gemmi': gemmi, '_struct': struct, '_sg': sg }
