
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
["cif","structure",...] query: the values needed are taken from its
data items as the printed strings (so their rounding errors are available).

"""

__all__ = []

from . import _common as _nc_common
from . import _sgsymimpl as _nc_sgsymimpl
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


#Temperature fields (in order of preference):
cif_temperature_fields = ( '_diffrn_ambient_temperature',
                           '_cell_measurement_temperature' )

def _atom( a ):
    #Atom of the query (cf. getAtoms in NCCIFSites.hh) as dict, with the
    #tensors as dicts (keys u11,...,u23):
    uk = _nc_sgsymimpl._uij_keys
    return { 'cif_labels': tuple( a['labels'] ),
             'composition': tuple( tuple( c ) for c in a['composition'] ),
             'occupancy': a['occupancy'],
             'uiso': a['uiso'],
             'aniso': dict( zip( uk, a['uij'] ) ) if a['uij'] else None,
             'positions': [ tuple( p ) for p in a['positions'] ],
             'uijs': ( [ dict( zip( uk, u ) ) for u in a['uijs'] ]
                       if a['uijs'] else None ),
             'split': a['split'] }

class CIFContent:
    #Content of CIF data, with printed strings (without quotes) for the cell,
    #formula and temperatures (None if absent or null, "?" or ".").
    __slots__ = ( 'atoms', 'cell', 'formula', 'provenance', 'setting',
                  'snapped_cell', 'temperatures' )
    def __init__( self, *, provenance, cell, snapped_cell, setting, atoms,
                  formula, temperatures ):
        #the provenance items of the structure query (or None):
        self.provenance = provenance
        self.cell = cell#a,b,c,alpha,beta,gamma (missing angles are "90")
        #the cell snapped to the setting (numbers), None if it does not fit
        #(and that is allowed):
        self.snapped_cell = snapped_cell
        self.setting = setting#the space group setting (e.g. "227:2")
        self.atoms = atoms#the atoms (cf. _atom)
        self.formula = formula#_chemical_formula_sum
        #_diffrn_ambient_temperature and _cell_measurement_temperature:
        self.temperatures = temperatures

def read_cif( cifdata, *, name, permissive, quiet, warnfct,
              override_spacegroup = None, tol_uij = None,
              allow_cell_misfit = False ):
    #Returns the CIFContent of the CIF text cifdata (with the name for
    #messages), read via the ["cif","structure",...] query in permissive mode
    #or not, with the space group setting overridden if requested, the
    #tolerance tol_uij (None for the default), and accepting a cell which does
    #not fit the setting if allow_cell_misfit. Warnings are passed to warnfct.
    if not quiet:
        _nc_common.infomsg("Parsing CIF data")
    from .misc import evaluate_query
    opts = ( [ f'override_spacegroup={override_spacegroup}' ]
             if override_spacegroup else [] )
    if tol_uij is not None:
        opts.append( f'tol_uij={float(tol_uij)!r}' )
    if allow_cell_misfit:
        opts.append( 'allow_cell_misfit=1' )
    q = evaluate_query( [ 'cif', 'structure',
                          'permissive' if permissive else 'strict', *opts,
                          'RAWDATA', name or '', 'RAWDATA', cifdata ] )
    for w in q['warnings']:
        warnfct( w )
    st = q['structure']
    def val( v ):
        return None if v is None else v[0]
    cell = {}
    for k in ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' ):
        isangle = k in ( 'alpha', 'beta', 'gamma' )
        cell[k] = val( st[ ( '_cell_angle_' if isangle else '_cell_length_' )
                           + k ] )
    return CIFContent(
        provenance = q['provenance'],
        cell = cell,
        snapped_cell = q['cell'],
        setting = q['spacegroup_setting'],
        atoms = [ _atom( a ) for a in q['atoms'] ],
        formula = val( st.get( '_chemical_formula_sum' ) ),
        temperatures = { t: val( st.get( t ) )
                         for t in cif_temperature_fields } )
