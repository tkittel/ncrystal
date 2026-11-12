
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
from . import core as _nc_core
from ._cifreadimpl import _extract_descr_from_cif, parse_cif_number, read_cif
from ._cifsitesimpl import (
    _impl_find_sites,
    _impl_group_sites,
    collect_atoms,
    merge_atoms,
)
from ._cifsymimpl import decide_symmetry


def _expected_atom_count( content ):
    #Total number of atoms in the cell, if stated in the CIF data:
    v = [ ( parse_cif_number( e ) if e is not None else None )
          for e in content['atom_type_number_in_cell'] ]
    if v and not any( e is None for e in v ):
        return sum( e[0] for e in v )
    return None#(absent, or with values which are not numbers like "?")

def load_cif_data( cifsrc, cifdata, *, quiet, refine_with_spglib,
                   merge_equiv, tols, permissive, override_spacegroup = None ):
    """Load CIF data (text, from the CIFSource cifsrc): read it (gemmi),
    decide which symmetry to use, collect and merge the atoms of the sites,
    and find the final sites (spglib). Returns the results for CIFLoader (a
    dict), and the gemmi module and structure."""
    from types import MappingProxyType
    dict_ro, list_ro = MappingProxyType, tuple

    result = {}
    result['cifsrc'] = cifsrc
    result['cifdata'] = cifdata

    content = read_cif( cifdata, quiet = quiet,
                        override_spacegroup = override_spacegroup )
    result['cif_raw'] = content['cif_raw']

    choice = decide_symmetry( content, permissive = permissive, tols = tols,
                              quiet = quiet,
                              refine_with_spglib = refine_with_spglib,
                              override_spacegroup = override_spacegroup,
                              cifdata = cifdata )
    expected_tot_atom_in_orig_cell = _expected_atom_count( content )
    expect_sgno = choice['expect_sgno']

    cellsg = dict( content['cell'] )
    #Although the sg object also provides sg.ccp4 and sg.hall, we record here
    #just sgnumber and xhm, since that is what spglib refinement also provides:
    cellsg['spacegroup'] = dict_ro( { 'number': choice['sgnumber'],
                                      'hm': choice['hm'] } )

    collected_atoms = collect_atoms( content, choice, tols = tols,
                                     permissive = permissive )
    final_atoms = merge_atoms( collected_atoms, expected_tot_atom_in_orig_cell )

    if not refine_with_spglib:
        _nc_common.warn('Structure refinement and verification with spglib was disabled due to explicit request.')

    cellsg, sites, setting = _impl_find_sites( cellsg, final_atoms,
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
    return result, content['_gemmi'], content['_struct']
