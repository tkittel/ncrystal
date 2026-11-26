
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
    #Content of CIF data, with printed strings (without quotes) for the cell
    #and temperatures (None if absent or null, "?" or ".").
    __slots__ = ( 'atoms', 'cell', 'description', 'setting', 'snapped_cell',
                  'temperatures' )
    def __init__( self, *, cell, snapped_cell, setting, atoms, description,
                  temperatures ):
        self.cell = cell#a,b,c,alpha,beta,gamma (missing angles are "90")
        #the cell snapped to the setting (numbers), None if it does not fit
        #(and that is allowed):
        self.snapped_cell = snapped_cell
        self.setting = setting#the space group setting (e.g. "227:2")
        self.atoms = atoms#the atoms (cf. _atom)
        #dict with 'lines' (for NCMAT comments), 'codid' and 'mpid' (int or
        #None) and 'formula' (cf. getDescription in NCCIFDescribe.hh):
        self.description = description
        #_diffrn_ambient_temperature and _cell_measurement_temperature:
        self.temperatures = temperatures

def read_cif( cifdata, *, name, permissive, quiet, warnfct,
              override_spacegroup = None, tol_uij = None,
              tol_symprec = None, allow_cell_misfit = False,
              requested_codid = None, requested_mpid = None ):
    #Returns the CIFContent of the CIF text cifdata (with the name for
    #messages), read via the ["cif","structure",...] query in permissive mode
    #or not, with the space group setting overridden if requested (else a
    #higher symmetry of the atoms in the same coordinates is used, unless
    #tol_symprec is 0), the tolerances tol_uij and tol_symprec (None for the
    #defaults), and accepting a cell which does not fit the setting if
    #allow_cell_misfit. The requested COD and MP ids (of the data source, None
    #if none) are used in the description. Warnings are passed to warnfct.
    if not quiet:
        _nc_common.infomsg("Parsing CIF data")
    from .misc import evaluate_query
    #(RAWDATA, as symbols can have any characters):
    opts = ( [ 'RAWDATA', f'override_spacegroup={override_spacegroup}' ]
             if override_spacegroup else [] )
    if tol_uij is not None:
        opts.append( f'tol_uij={float(tol_uij)!r}' )
    if tol_symprec is not None:
        opts.append( f'tol_symprec={float(tol_symprec)!r}' )
    if allow_cell_misfit:
        opts.append( 'allow_cell_misfit=1' )
    if requested_codid:
        opts.append( f'cod_id={int(requested_codid)}' )
    if requested_mpid:
        opts.append( f'mp_id={int(requested_mpid)}' )
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
        cell = cell,
        snapped_cell = q['cell'],
        setting = q['spacegroup_setting'],
        atoms = [ _atom( a ) for a in q['atoms'] ],
        description = { 'lines': q['description']['lines'],
                        'codid': q['description']['cod_id'],
                        'mpid': q['description']['mp_id'],
                        'formula': q['description']['formula'] },
        temperatures = { t: val( st.get( t ) )
                         for t in cif_temperature_fields } )
