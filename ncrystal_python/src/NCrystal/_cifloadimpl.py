
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

Internal implementation of loading CIF data, for the CIFLoader class of the
cifutils module.

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
from ._cifsymimpl import choose_setting, snapped_cell


def _expected_atom_count( content ):
    #Total number of atoms in the cell, if stated in the CIF data:
    from .misc import evaluate_query
    v = content['atom_type_number_in_cell']
    if not v or any( e is None for e in v ):
        return None
    n = evaluate_query( [ 'cif', 'numbers', *_nc_sgsymimpl._rawargs( v ) ] )
    if any( e is None for e in n ):
        return None#(values which are not numbers)
    return sum( e[0] for e in n )

def load_cif_data( cifsrc, cifdata, *, quiet, refine_with_spglib,
                   merge_equiv, tols, permissive, override_spacegroup = None ):
    """Load CIF data (text, from the CIFSource cifsrc): read it, choose the
    space group setting, collect and merge the atoms of the sites, and find
    the final sites (refined with spglib if requested). Returns the results
    for CIFLoader (a dict)."""
    from types import MappingProxyType

    from .misc import evaluate_query
    dict_ro, list_ro = MappingProxyType, tuple

    result = {}
    result['cifsrc'] = cifsrc
    result['cifdata'] = cifdata

    content = read_cif( cifdata, quiet = quiet )
    result['cif_raw'] = content['cif_raw']
    result['cif_temperatures'] = content['temperatures']

    setting = choose_setting( content, permissive = permissive,
                              override_spacegroup = override_spacegroup,
                              cifdata = cifdata )
    cell = content['cell']
    expect_sgno = None
    try:
        cellsg = dict( snapped_cell( setting, cell ) )
    except _nc_core.NCBadInput as e:
        #The cell does not quite have the symmetry of the space group. This is
        #only accepted if refining the structure finds the same space group:
        if not refine_with_spglib:
            raise _nc_core.NCBadInput( f'Unit cell not compatible with the'
                                       f' space group of the CIF data: {e}' )
        expect_sgno = int( setting.split(':')[0] )
        keys = ( 'a', 'b', 'c', 'alpha', 'beta', 'gamma' )
        vals = evaluate_query( [ 'cif', 'numbers',
                                 *_nc_sgsymimpl._rawargs( cell[k]
                                                          for k in keys ) ] )
        cellsg = { k: v[0] for k, v in zip( keys, vals ) }
    info = evaluate_query( [ 'sgsym', 'info', 'strict', setting ] )
    cellsg['spacegroup'] = dict_ro( { 'number': int( setting.split(':')[0] ),
                                      'hm': info['hm'] } )

    collected_atoms = collect_atoms( content, setting, cell, tols = tols,
                                     permissive = permissive )
    final_atoms = merge_atoms( collected_atoms, _expected_atom_count( content ),
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
