
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

Internal implementation of reading CIF data for the cifutils module, via the
["cif","uniquestructprov",...] query: the values needed are taken from its
data items as the printed strings (so their rounding errors are available).

"""

__all__ = []

from . import _common as _nc_common
from . import core as _nc_core
from ._cifsrcimpl import _codid2url, _mpid2url


def _extract_descr( prov, formula, cifsrc, *, warnfct ):
    #Description lines, COD and MP ids, and chemical formula, from the
    #provenance items (cf. read_cif, None if absent) and formula of CIF data:
    prov = prov or {}
    def get( *keys ):
        for k in keys:
            if k in prov:
                v = prov[k]
                return '' if v is None else v[0].strip()
        return ''
    def getlist( *keys ):
        for k in keys:
            if k in prov:
                return [ v[0].strip() for v in prov[k]
                         if v is not None and v[0].strip() ]
        return []
    def normalise_result( _l ):
        #split any newlines and trim whitespace before returning:
        return [ ee.rstrip() for e in _l for ee in e.splitlines() ]

    title = get('_publ_section_title','_citation_title')
    authors = getlist('_publ_author_name','_citation_author_name')
    journalname = get('_journal_name_full','_citation_journal_full')
    year = get('_journal_year','_citation_journal_year','_citation_year')
    doi = get('_journal_paper_doi')
    chemformsum = formula or get('_cod_original_formula_sum')

    extracted_codid = get('_cod_database_code')
    extracted_codid = ( int(extracted_codid) if extracted_codid.isdigit()
                        else None )
    _ = get('_[local]_ncrystal_mp_database_code')
    extracted_mpid = ( int(_[3:]) if _.startswith('mp-') and _[3:].isdigit()
                       else None )

    if extracted_codid is not None and cifsrc.codid is not None and cifsrc.codid != extracted_codid:
        warnfct(f'Embedded CODID {extracted_codid} is different than requested {cifsrc.codid}!')

    if extracted_mpid is not None and cifsrc.mpid is not None and cifsrc.mpid != extracted_mpid:
        warnfct(f'Embedded MPID {extracted_mpid} is different than requested {cifsrc.mpid}!')

    _thecodid = cifsrc.codid or extracted_codid
    _thempid = cifsrc.mpid or extracted_mpid

    audit_creation_method = get('_audit_creation_method')
    audit_creation_date = get('_audit_creation_date')
    audit_author_name = ', '.join( getlist('_audit_author_name') )

    ll = []
    nauthors = len(authors)
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

_uij_tags = ( 'U_11', 'U_22', 'U_33', 'U_12', 'U_13', 'U_23' )

#Temperature fields (in order of preference):
cif_temperature_fields = ( '_diffrn_ambient_temperature',
                           '_cell_measurement_temperature' )

class CIFContent:
    #Plain content of CIF data. Values are the printed strings (without
    #quotes), None if absent or null ("?" or ".").
    __slots__ = ( 'atom_type_number_in_cell', 'cell', 'formula', 'hall', 'hm',
                  'itnumber', 'listed_ops', 'provenance', 'sites',
                  'temperatures', 'text' )
    def __init__( self, *, text, provenance, cell, listed_ops, hall, hm,
                  itnumber, formula, sites, atom_type_number_in_cell,
                  temperatures ):
        self.text = text#the CIF text
        #the provenance items of the uniquestructprov query (or None):
        self.provenance = provenance
        self.cell = cell#a,b,c,alpha,beta,gamma (missing angles are "90")
        self.listed_ops = listed_ops#listed symmetry operations (or None)
        self.hall = hall#space group symbols and IT number
        self.hm = hm
        self.itnumber = itnumber
        self.formula = formula#_chemical_formula_sum
        #dicts with keys label, type_symbol, x, y, z, uiso, biso, occ and uij
        #(the six aniso U values, None if absent):
        self.sites = sites
        self.atom_type_number_in_cell = atom_type_number_in_cell
        #_diffrn_ambient_temperature and _cell_measurement_temperature:
        self.temperatures = temperatures

def read_cif( cifdata, *, name, permissive, quiet, warnfct ):
    #Returns the CIFContent of the CIF text cifdata (with the name for
    #messages), read via the ["cif","uniquestructprov",...] query in permissive
    #mode or not. Warnings are passed to warnfct.
    if not quiet:
        _nc_common.infomsg("Parsing CIF data")
    from .misc import evaluate_query
    q = evaluate_query( [ 'cif', 'uniquestructprov',
                          'permissive' if permissive else 'strict',
                          'RAWDATA', name or '', 'RAWDATA', cifdata ] )
    for w in q['warnings']:
        warnfct( w )
    st = q['structure']
    def val( v ):
        return None if v is None else v[0]
    def first( alternatives ):
        #(the first value which is not "?" or ".", current tags first):
        return next( ( val( a['value'] ) for a in alternatives
                       if a['value'] is not None ), None )

    cell = {}
    for k in ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' ):
        isangle = k in ( 'alpha', 'beta', 'gamma' )
        v = val( st[ ( '_cell_angle_' if isangle else '_cell_length_' ) + k ] )
        if v is None:
            raise _nc_core.NCBadInput(f'CIF data lacks the cell parameter {k}')
        cell[k] = v

    sym = st['symmetry']
    listed = None
    for o in sym['operations']:
        listed = [ v[0] for v in o['values'] if v is not None ]
        if listed:
            break

    cols = st['atom_site']
    keys = { 'label': '_atom_site_label',
             'type_symbol': '_atom_site_type_symbol',
             'x': '_atom_site_fract_x', 'y': '_atom_site_fract_y',
             'z': '_atom_site_fract_z',
             'uiso': '_atom_site_U_iso_or_equiv',
             'biso': '_atom_site_B_iso_or_equiv',
             'occ': '_atom_site_occupancy' }
    nsites = len( cols.get( '_atom_site_label', [] ) )
    sites = [ { k: ( val( cols[tag][i] ) if tag in cols else None )
                for k, tag in keys.items() } for i in range( nsites ) ]
    #Aniso rows (rows with "?" components are dropped by the query) go to the
    #sites with their label, in order (some files reuse labels):
    for s in sites:
        s['uij'] = None
    aniso = st.get( 'atom_site_aniso' )
    if aniso and aniso['form'] == 'B':
        #FIXME: Convert to U (with the rounding errors), preferably in C++:
        warnfct( 'Ignoring the anisotropic displacements, which are given as'
                 ' B (_atom_site_aniso_B_ij), which is not supported yet' )
    elif aniso:
        taken = set()
        comps = [ aniso[ '_atom_site_aniso_' + t ] for t in _uij_tags ]
        for r, lbl in enumerate( aniso['_atom_site_aniso_label'] ):
            i = next( ( i for i, s in enumerate( sites )
                        if s['label'] == val( lbl ) and i not in taken ),
                      None )
            if i is not None:
                taken.add( i )
                sites[i]['uij'] = [ val( c[r] ) for c in comps ]

    return CIFContent(
        text = cifdata,
        provenance = q['provenance'],
        cell = cell,
        listed_ops = listed or None,
        hall = first( sym['hall'] ),
        hm = first( sym['hm'] ),
        itnumber = first( sym['number'] ),
        formula = val( st.get( '_chemical_formula_sum' ) ),
        sites = sites,
        atom_type_number_in_cell = [ val( v ) for v in
                                     st.get( '_atom_type_number_in_cell', [] ) ],
        temperatures = { t: val( st.get( t ) )
                         for t in cif_temperature_fields } )
