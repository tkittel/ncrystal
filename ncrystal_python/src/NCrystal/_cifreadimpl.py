
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

Internal implementation of reading CIF data for the cifutils module: gemmi
only parses the CIF text, and the values needed are read from its tags and
loops as the printed strings (so their rounding errors are available).

"""

__all__ = []

from . import _common as _nc_common
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
    from ._imports import import_optional
    gemmi = import_optional( 'gemmi', purpose = 'to process CIF files',
                             sysexit = sysexit )
    import importlib
    _import_gemmi_cache[0] = gemmi
    _import_gemmi_cache[1] = importlib.import_module('gemmi.cif')
    return _import_gemmi_cache[0], _import_gemmi_cache[1]

def _decode_json( txt ):
    #gemmi <= 0.6.5 writes CIF numbers like -.065 as invalid JSON numbers, so
    #on failure add the missing zeros (outside JSON strings) and retry.
    import json
    try:
        return json.loads( txt )
    except json.decoder.JSONDecodeError:
        import re
        fixed = re.sub( r'("(?:[^"\\]|\\.)*")|(?<![0-9])(-?)\.(?=[0-9])',
                        lambda m : m.group(1) or m.group(2) + '0.', txt )
        if fixed == txt:
            raise
        return json.loads( fixed )

def _select_block( cif_doc ):
    #The block with the structure (the only one, or the only one with sites):
    if len( cif_doc ) == 0:
        raise _nc_core.NCBadInput('CIF data had no blocks!')
    if len( cif_doc ) == 1:
        return cif_doc.sole_block()
    blk = None
    for e in cif_doc:
        if not len( e.find( [ '_atom_site_fract_x' ] ) ):
            continue
        if blk:
            blk = None
            break
        blk = e
    if not blk:
        raise _nc_core.NCBadInput('Could not automatically determine block in'
                                  ' CIF data with structure info')
    return blk

_uij_tags = ( 'U_11', 'U_22', 'U_33', 'U_12', 'U_13', 'U_23' )

#Temperature fields (in order of preference):
cif_temperature_fields = ( '_diffrn_ambient_temperature',
                           '_cell_measurement_temperature' )

def read_cif( cifdata, *, quiet ):
    """Read CIF data (text). Values are the printed strings (without quotes),
    None if absent or null ("?" or "."). Returns a dict with:
      'cif_raw': the document as a dict (None if not decodable),
      'cell': a,b,c,alpha,beta,gamma (missing angles are "90"),
      'listed_ops': the listed symmetry operations (None if absent),
      'hall', 'hm', 'itnumber': the space group symbols and number,
      'formula': the chemical formula (_chemical_formula_sum),
      'sites': dicts with keys label, type_symbol, x, y, z, uiso, biso, occ and
               uij (the six aniso values, None if not all given),
      'atom_type_number_in_cell': the values of that field,
      'temperatures': the printed values of the temperature fields
                      (_diffrn_ambient_temperature and
                      _cell_measurement_temperature, None if absent)."""
    _, gemmi_cif = _import_gemmi()
    if not quiet:
        _nc_common.print("Attempting to load CIF data with gemmi")
    try:
        cif_doc = gemmi_cif.read_string( cifdata )
    except ( ValueError, RuntimeError ) as e:
        #(gemmi raises RuntimeError for e.g. duplicate tags)
        raise _nc_core.NCBadInput(f'CIF parsing error (from Gemmi): "{e}"')
    import json
    try:
        cif_raw = _decode_json( cif_doc.as_json() )
    except json.decoder.JSONDecodeError as e:
        cif_raw = None
        _nc_common.warn(f'Could not decode raw CIF data to dictionary (Gemmi bug?): JSONDecodeError("{e}")')
    blk = _select_block( cif_doc )

    def val( raw ):
        return None if raw is None or gemmi_cif.is_null( raw ) else gemmi_cif.as_string( raw )
    def first( *tags ):
        for tag in tags:
            v = blk.find_value( tag )
            if v is not None:
                return val( v )
    def column( tag ):
        return [ val( v ) for v in blk.find_values( tag ) ]

    cell = {}
    for k in ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' ):
        isangle = k in ( 'alpha', 'beta', 'gamma' )
        v = first( ( '_cell_angle_' if isangle else '_cell_length_' ) + k )
        if v is None and isangle and blk.find_value( '_cell_angle_' + k ) is None:
            v = '90'#(missing angles are 90 degrees)
        if v is None:
            raise _nc_core.NCBadInput(f'CIF data lacks the cell parameter {k}')
        cell[k] = v

    listed = None
    for tag in ( '_space_group_symop_operation_xyz',
                 '_symmetry_equiv_pos_as_xyz' ):
        listed = [ v for v in column( tag ) if v is not None ]
        if listed:
            break

    cols = [ 'label', '?type_symbol', '?fract_x', '?fract_y', '?fract_z',
             '?U_iso_or_equiv', '?B_iso_or_equiv', '?occupancy' ]
    keys = [ 'label', 'type_symbol', 'x', 'y', 'z', 'uiso', 'biso', 'occ' ]
    sites = []
    for row in blk.find( '_atom_site_', cols ):
        sites.append( { k: ( val( row[i] ) if row.has( i ) else None )
                        for i, k in enumerate( keys ) } )
    #Aniso rows go to the sites with their label, in order (some files reuse
    #labels):
    for s in sites:
        s['uij'] = None
    taken = set()
    for row in blk.find( '_atom_site_aniso_', [ 'label' ] + list( _uij_tags ) ):
        v = [ val( row[i] ) for i in range( 1, 7 ) ]
        i = next( ( i for i, s in enumerate( sites )
                    if s['label'] == val( row[0] ) and i not in taken ), None )
        if i is not None:
            taken.add( i )
            sites[i]['uij'] = None if any( e is None for e in v ) else v

    return { 'cif_raw': cif_raw,
             'cell': cell,
             'listed_ops': listed or None,
             'hall': first( '_space_group_name_Hall',
                            '_symmetry_space_group_name_Hall' ),
             'hm': first( '_space_group_name_H-M_alt',
                          '_symmetry_space_group_name_H-M' ),
             'itnumber': first( '_space_group_IT_number',
                                '_symmetry_Int_Tables_number' ),
             'formula': first( '_chemical_formula_sum' ),
             'sites': sites,
             'atom_type_number_in_cell': [ v for v in column( '_atom_type_number_in_cell' ) ],
             'temperatures': { t: first( t ) for t in cif_temperature_fields } }
