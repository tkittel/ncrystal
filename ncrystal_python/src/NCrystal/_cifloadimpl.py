
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

from . import _common as _nc_common
from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core
from ._cifreadimpl import _extract_descr_from_cif, read_cif
from ._cifsitesimpl import (
    _impl_find_sites,
    _impl_group_sites,
    check_split_sites,
    collect_atoms,
    merge_atoms,
)
from ._cifsymimpl import determine_setting, snapped_cell


def _expected_atom_count( cif_content ):
    #Total number of atoms in the cell, if stated in the CIF data:
    from .misc import evaluate_query
    v = cif_content.atom_type_number_in_cell
    if not v or any( e is None for e in v ):
        return None
    n = evaluate_query( [ 'cif', 'numbers', *_nc_sgsymimpl._rawargs( v ) ] )
    if any( e is None for e in n ):
        return None#(values which are not numbers)
    return sum( e[0] for e in n )

def load_cif_data( cifsrc, cifdata, *, quiet, refine_with_spglib,
                   merge_equiv, tols, permissive, override_spacegroup = None ):
    #Returns the CIFLoader fields (a dict) for the CIF text cifdata, from the
    #CIFSource cifsrc.
    from types import MappingProxyType

    from .misc import evaluate_query
    dict_ro, list_ro = MappingProxyType, tuple

    result = {}
    result['cifsrc'] = cifsrc
    result['cifdata'] = cifdata

    cif_content = read_cif( cifdata, quiet = quiet )
    result['cif_raw'] = cif_content.cif_raw

    setting = determine_setting( cif_content, permissive = permissive,
                                 override_spacegroup = override_spacegroup )
    cell = cif_content.cell
    expect_sgno = None
    snapped = snapped_cell( setting, cell )
    if snapped['cell'] is not None:
        cellsg = dict( snapped['cell'] )
    else:
        #The cell does not quite have the symmetry of the space group. This is
        #only accepted if refining the structure finds the same space group:
        if not refine_with_spglib:
            raise _nc_core.NCBadInput( 'Unit cell not compatible with the'
                                       ' space group of the CIF data:'
                                       f' {snapped["error_msg"]}' )
        expect_sgno = int( setting.split(':')[0] )
        keys = ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' )
        vals = evaluate_query( [ 'cif', 'numbers',
                                 *_nc_sgsymimpl._rawargs( cell[k]
                                                          for k in keys ) ] )
        cellsg = { k: v[0] for k, v in zip( keys, vals ) }
    info = evaluate_query( [ 'sgsym', 'info', 'strict', setting ] )
    cellsg['spacegroup'] = dict_ro( { 'number': int( setting.split(':')[0] ),
                                      'hm': info['hm'] } )

    collected_atoms = collect_atoms( cif_content, setting, cell, tols = tols,
                                     permissive = permissive )
    final_atoms = merge_atoms( collected_atoms, _expected_atom_count( cif_content ),
                               setting, cell )
    check_split_sites( final_atoms, permissive = permissive )

    if not refine_with_spglib:
        _nc_common.warn('Structure refinement and verification with spglib was disabled due to explicit request.')

    cellsg, sites, setting = _impl_find_sites( cellsg, final_atoms, setting,
                                               tols = tols,
                                               refine = refine_with_spglib,
                                               permissive = permissive,
                                               quiet = quiet )
    if expect_sgno is not None:
        sgno = int( cellsg['spacegroup']['number'] )
        if sgno != expect_sgno:
            raise _nc_core.NCBadInput(
                'The unit cell parameters are not compatible with the space'
                f' group of the CIF data ({expect_sgno}), and refining the'
                f' structure gives space group {sgno}. If the declared'
                ' symmetry is trusted, try again with a larger tolerance'
                ' (tol).' )
        _nc_common.warn( 'The unit cell parameters do not exactly fit the'
                         f' space group of the CIF data ({expect_sgno}), so'
                         ' the structure was refined to exactly consistent'
                         ' cell parameters' )

    if int(cellsg['spacegroup']['number']) == 1:
        #The following warning might be a bit obsolete now that we usually
        #refine cells, but we keep it here for now as an extra safeguard:
        _nc_common.warn('Space-group number 1 ("P1") detected in CIF input. A lot (but not all) of CIF data'
                       ' with this space group listed is in a non-conventional format not suitable for NCrystal. Be wary!')

    atoms, site2atom = _impl_group_sites( sites, tols = tols,
                                          merge_equiv = merge_equiv )
    result['cellsg'] = dict_ro( cellsg )
    result['atoms'] = list_ro( dict_ro(a) for a in atoms )
    result['sites'] = list_ro( dict_ro(e) for e in sites )
    result['site2atom'] = tuple( site2atom )
    result['sgsetting'] = setting

    result['cifdescr'], result['actual_codid'],result['actual_mpid'],result['cif_chemformula'] =_extract_descr_from_cif( result['cif_raw'], result['cifsrc'], result['cifdata'] )
    return result
