
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

Internal implementation details for loading CIF data with the CIFLoader class
in cifutils.py.

"""

__all__ = []

from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core
from ._cifreadimpl import read_cif
from ._cifsitesimpl import _impl_find_sites, _impl_group_sites


def load_cif_data( cifsrc, cifdata, *, quiet, refine_with_spglib,
                   merge_equiv, tols, permissive, warnfct,
                   override_spacegroup = None ):
    #Returns the CIFLoader fields (a dict) for the CIF text cifdata, from the
    #CIFSource cifsrc: read with its space group setting and atoms (from C++,
    #cf. read_cif), and the final sites (refined with spglib if requested).
    #Warnings are passed to warnfct.
    from types import MappingProxyType

    from .misc import evaluate_query
    dict_ro, list_ro = MappingProxyType, tuple

    result = {}
    result['cifsrc'] = cifsrc
    result['cifdata'] = cifdata

    cif_content = read_cif( cifdata, name = cifsrc.name,
                            permissive = permissive, quiet = quiet,
                            warnfct = warnfct,
                            override_spacegroup = override_spacegroup,
                            tol_uij = tols['uij'],
                            #(spglib finds any higher symmetry itself):
                            tol_symprec = ( 0.0 if refine_with_spglib
                                            else tols['symprec'] ),
                            allow_cell_misfit = refine_with_spglib,
                            requested_codid = cifsrc.codid,
                            requested_mpid = cifsrc.mpid )
    result['cif_temperatures'] = cif_content.temperatures

    setting = cif_content.setting
    cell = cif_content.cell
    expect_sgno = None
    if cif_content.snapped_cell:
        cellsg = dict( cif_content.snapped_cell )
    else:
        #The cell does not quite have the symmetry of the space group. This is
        #only accepted if refining the structure finds the same space group:
        expect_sgno = int( setting.split(':')[0] )
        keys = ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' )
        vals = evaluate_query( [ 'cif', 'numbers',
                                 *_nc_sgsymimpl._rawlist( cell[k]
                                                          for k in keys ) ] )
        cellsg = { k: v[0] for k, v in zip( keys, vals ) }
    info = evaluate_query( [ 'sgsym', 'info', 'strict', setting ] )
    cellsg['spacegroup'] = dict_ro( { 'number': int( setting.split(':')[0] ),
                                      'hm': info['hm'] } )

    final_atoms = cif_content.atoms

    if not refine_with_spglib:
        warnfct('Structure refinement and verification with spglib was disabled due to explicit request.')

    cellsg, sites, setting = _impl_find_sites( cellsg, final_atoms, setting,
                                               tols = tols,
                                               refine = refine_with_spglib,
                                               permissive = permissive,
                                               quiet = quiet,
                                               warnfct = warnfct )
    if expect_sgno is not None:
        sgno = int( cellsg['spacegroup']['number'] )
        if sgno != expect_sgno:
            raise _nc_core.NCBadInput(
                'The unit cell parameters are not compatible with the space'
                f' group of the CIF data ({expect_sgno}), and refining the'
                f' structure gives space group {sgno}. If the declared'
                ' symmetry is trusted, try again with a larger tolerance'
                ' (tol).' )
        warnfct( 'The unit cell parameters do not exactly fit the'
                         f' space group of the CIF data ({expect_sgno}), so'
                         ' the structure was refined to exactly consistent'
                         ' cell parameters' )

    if int(cellsg['spacegroup']['number']) == 1:
        warnfct( 'The structure has no symmetry beyond space group 1'
                         ' (P1). This is rare for real crystals and might'
                         ' indicate a supercell or low-precision coordinates.'
                         ' Be wary!' if refine_with_spglib else
                         'The structure has space group 1 (P1), and was not'
                         ' checked for further symmetry. Be wary!' )

    atoms, site2atom = _impl_group_sites( sites, tols = tols,
                                          merge_equiv = merge_equiv )
    result['cellsg'] = dict_ro( cellsg )
    result['atoms'] = list_ro( dict_ro(a) for a in atoms )
    result['sites'] = list_ro( dict_ro(e) for e in sites )
    result['site2atom'] = tuple( site2atom )
    result['sgsetting'] = setting

    d = cif_content.description
    ( result['cifdescr'], result['actual_codid'], result['actual_mpid'],
      result['cif_chemformula'] ) = ( d['lines'], d['codid'], d['mpid'],
                                      d['formula'] )
    return result
