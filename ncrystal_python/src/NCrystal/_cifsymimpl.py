
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

Internal implementation details for determining the space group setting of
CIF data in cifutils.py.

"""

__all__ = []

from . import _common as _nc_common
from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core
from ._cifreadimpl import CIFContent

_rawargs = _nc_sgsymimpl._rawargs


def _query( *args ):
    from .misc import evaluate_query
    return evaluate_query( list( args ) )

def cell_fits( setting, cell ):
    #Whether the cell (dict of printed values) fits the setting (with lenient
    #tolerances, cf. the sgsym cell query).
    return snapped_cell( setting, cell )['cell'] is not None

def snapped_cell( setting, cell ):
    #Returns the sgsym cell query result ({'cell':..,'error_msg':..}) for the
    #cell (dict of printed values) in the setting, with lenient tolerances.
    return _query( 'sgsym', 'cell', 'lenient', '0', setting,
                   *_rawargs( cell[k] for k in ( 'a', 'b', 'c',
                                                 'alpha', 'beta', 'gamma' ) ) )

def _symbol_settings( cif_content ):
    #Settings from the space group symbol (Hall, else H-M, else IT number), and
    #a description of the symbol used (or of the last one not recognised, or
    #None if there is none):
    descr = None
    if cif_content.hall:
        descr = f'Hall symbol "{cif_content.hall}"'
        s = _query( 'sgsym', 'fromhall', *_rawargs( [ cif_content.hall ] ) )
        if s:
            return s, descr
    if cif_content.hm:
        descr = f'space group symbol "{cif_content.hm}"'
        s = _query( 'sgsym', 'fromhm', *_rawargs( [ cif_content.hm ] ) )
        if s:
            return s, descr
    n = cif_content.itnumber
    if n:
        descr = f'space group number "{n}"'
        if n.strip().isdigit() and 1 <= int( n ) <= 230:
            return ( _query( 'sgsym', 'fromnumber', n.strip() ),
                     f'space group number {int(n)}' )
    return [], descr

def _number( setting ):
    return int( setting.split( ':' )[0] )

def _resolve( settings, cell, what ):
    #The one of the candidate settings which the cell fits:
    if not settings:
        raise _nc_core.NCBadInput( 'Could not determine space group from CIF'
                                   + ( f' data ({what} not recognised)'
                                       if what else ' data' ) )
    if len( settings ) == 1:
        return settings[0]
    fit = [ s for s in settings if cell_fits( s, cell ) ]
    if len( fit ) == 1:
        return fit[0]
    if not fit:
        raise _nc_core.NCBadInput( f'The unit cell does not fit any of the'
                                   f' settings of the {what}'
                                   f' ({", ".join(settings)})' )
    raise _nc_core.NCBadInput( f'The {what} has several settings fitting the'
                               f' unit cell ({", ".join(fit)}). Please specify'
                               ' which one to use (e.g. override_spacegroup='
                               f'"{fit[0]}", or --spacegroup={fit[0]} for'
                               ' ncrystal_cif2ncmat)' )

def determine_setting( cif_content, *, permissive, override_spacegroup ):
    #Returns the space group setting (e.g. "227:2") of cif_content:
    #override_spacegroup (symbol or setting) if given, else from the listed
    #symmetry operations checked against the symbol, else from the symbol and
    #the cell. Inconsistencies raise, or warn if permissive.
    assert isinstance( cif_content, CIFContent )
    cell = cif_content.cell
    if override_spacegroup:
        o = str( override_spacegroup )
        _nc_common.warn(f'Overriding spacegroup to "{o}" due to explicit request')
        s = _query( 'sgsym', 'fromhm', *_rawargs( [ o ] ) )
        if not s:
            r = _query( 'sgsym', 'setting', '0', *_rawargs( [ o ] ) )
            if r['status'] == 'syntax':
                raise _nc_core.NCBadInput( f'Invalid override_spacegroup "{o}":'
                                           ' neither the H-M symbol of a'
                                           ' tabulated setting nor a setting'
                                           f' ({r["error_msg"]})' )
            if r['status'] != 'ok':
                raise _nc_core.NCBadInput( f'Invalid override_spacegroup "{o}":'
                                           f' {r["error_msg"]}' )
            s = [ r['setting'] ]
        return _resolve( s, cell, f'space group "{o}"' )
    ops = cif_content.listed_ops
    symsettings, symdescr = _symbol_settings( cif_content )
    if ops:
        try:
            s = _query( 'sgsym', 'fromops', *_rawargs( ops ) )
            why = ( None if s else 'are not those of a tabulated space group'
                    ' setting' )
        except _nc_core.NCBadInput as e:
            s, why = [], f'could not be parsed ({e})'
        if s:
            setting = s[0]#(a few settings share operations, then any is fine)
            if symsettings and _number( setting ) not in { _number( e ) for e
                                                           in symsettings }:
                why = ( f'are those of space group {_number(setting)}, but the'
                        f' {symdescr} is that of space group'
                        f' {_number(symsettings[0])}' )
                if not permissive:
                    raise _symop_error( why, cif_content.text, True )
                _nc_common.warn( f'The symmetry operations listed in the CIF'
                                 f' data {why}. The listed operations are used'
                                 ' (accepted due to permissive mode)' )
            return setting
        if not permissive:
            raise _symop_error( why, cif_content.text, bool( symsettings ) )
        setting = _resolve( symsettings, cell, symdescr )
        _nc_common.warn( f'The symmetry operations listed in the CIF data {why}'
                         f', so the {symdescr} is used (accepted due to'
                         ' permissive mode)' )
        return setting
    return _resolve( symsettings, cell, symdescr )

def _symop_error( why, cifdata, permissive_hint ):
    msg = ( f'The symmetry operations listed in the CIF data {why}. Please use'
            ' a CIF file that either lists all atoms of the unit cell with'
            ' space group P 1, or lists the symmetry operations of a standard'
            ' space group setting consistent with its space group symbol.' )
    if permissive_hint:
        msg += ( ' Alternatively, try permissive mode (at your own risk).' )
    if 'pymatgen' in cifdata.lower():
        msg += ( ' This file was written by pymatgen, as are Materials Project'
                 ' downloads. For Materials Project entries, the "Conventional'
                 ' Standard" CIF download from the website works, or pass the'
                 ' ID as "mpid::<number>" instead of a file name to have'
                 ' NCrystal download that itself.' )
    return _nc_core.NCBadInput( msg )
