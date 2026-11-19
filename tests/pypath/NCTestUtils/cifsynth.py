
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


# Tests of CIF loading with the synthetic CIF files in tests/data/cifsynth.
#
# The synthesised files are written so as to replicate real issues found in a
# large survey of CIF files from various sources (COD, Materials Project, McStas
# data dir, ...), covering all space group settings, symmetry information,
# number formats, site and displacement data, warnings and errors. Files are
# written with code coverage in mind, and are unsuitable for physics usage.
#
# Each file is loaded (in permissive mode if the default mode fails), and a
# summary is printed. For some files, the material is also created.

import pathlib
import warnings

import NCrystalDev as NC
import NCrystalDev.cifutils as nc_cifutils

from .dirs import get_named_test_data_dir

#The files of the normal test (selected to cover all errors, warnings, symmetry
#information and number formats of the full set, but not all space group
#settings). The full set is tested in the long tests:
_short_numbers = """
    0005 0006 0007 0009 0010 0013 0016 0017 0018 0019 0022 0023 0024 0025
    0032 0035 0039 0046 0049 0052 0053 0054 0059 0066 0068 0069 0070 0074
    0078 0079 0082 0095 0103 0108 0110 0117 0122 0131 0132 0133 0139 0179
    0237 0246 0248 0253 0255 0259 0263 0266 0267 0271 0287 0291 0294 0295
    0301 0302 0304 0314 0329 0343 0357 0358 0361 0383 0385 0394 0395 0416
    0424 0442 0450 0481 0495 0500 0509 0519 0520 0521 0522
"""
short_test_files = tuple( f'cs{n}.cif' for n in _short_numbers.split() )

def all_files():
    d = pathlib.Path( get_named_test_data_dir( 'cifsynth' ) )
    return sorted( d.glob( 'cs*.cif' ) )

def _load( f, permissive ):
    with warnings.catch_warnings( record = True ) as wl:
        warnings.simplefilter( 'always' )
        try:
            lc = nc_cifutils.CIFLoader( f, quiet = True,
                                        permissive = permissive )
        except NC.NCException as e:
            msg = f'{e.__class__.__name__}: {e}'
            #(details of gemmi's parse errors depend on its version)
            p = 'CIF parsing error (from Gemmi)'
            return None, ( msg[:msg.index(p)+len(p)] if p in msg else msg ), []
    return lc, None, [ str( w.message ) for w in wl ]

def _summary( lc ):
    sg = lc.cellsg.get( 'spacegroup', {} )
    naniso = sum( 1 for s in lc.sites if s.get( 'aniso' ) )
    return ( f'SG-{sg.get("number")} ({lc._sgsetting}),'
             f' {len(lc.sites)} sites, {len(lc.atoms)} atoms,'
             f' {naniso} with uij' )

def _material( lc ):
    try:
        txt = lc.create_ncmat( quiet = True, no_formula_check = True )
        NC.registerInMemoryFileData( 'cifsynth.ncmat', txt )
        cfg = 'cifsynth.ncmat;dcutoff=0.5'
        info = NC.createInfo( cfg )
        NC.createScatter( cfg + ';comp=bragg' )
    except NC.NCException as e:
        return f'material: {e.__class__.__name__}: {e}'
    hkl = list( info.hklObjects() )
    dmax = max( h.d for h in hkl ) if hkl else None
    return ( f'material: {len(hkl)} hkl families'
             + ( f', dmax={dmax:.5g}' if dmax else '' ) )

def check_file( f, with_material ):
    print( f'{f.name}:' )
    lc = None
    for permissive in ( False, True ):
        lc, err, wl = _load( f, permissive )
        mode = 'permissive' if permissive else 'default'
        for w in wl:
            print( f'  {mode} warning: {w}' )
        if lc is None:
            print( f'  {mode} error: {err}' )
            continue
        print( f'  {mode}: {_summary( lc )}' )
        break
    if lc is not None and with_material:
        print( f'  {_material( lc )}' )

def run( files ):
    """Test the files. Every fourth file also creates the material (if loading
    succeeds)."""
    for f in files:
        check_file( f, with_material = ( int( f.stem[2:] ) % 4 == 0 ) )

def run_short():
    run( [ f for f in all_files() if f.name in short_test_files ] )

def run_all( part, nparts ):
    """Test part number part (0-based) of nparts of all files, with parts of
    similar run time (balanced by file sizes)."""
    parts, sizes = [ [] for _ in range( nparts ) ], [ 0 ] * nparts
    for f in sorted( all_files(), key = lambda f : ( -f.stat().st_size, f.name ) ):
        i = sizes.index( min( sizes ) )
        parts[i].append( f )
        sizes[i] += f.stat().st_size
    run( sorted( parts[part] ) )
