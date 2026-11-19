
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

#The files of each long test (balanced by run time):
_long_numbers = (
    """
    0004 0012 0013 0014 0019 0020 0029 0037 0042 0045 0049 0050 0051 0053
    0057 0063 0066 0080 0082 0083 0086 0087 0090 0096 0097 0103 0110 0115
    0118 0119 0126 0130 0132 0133 0135 0143 0148 0159 0170 0186 0189 0193
    0202 0208 0209 0213 0215 0216 0218 0225 0226 0228 0231 0233 0238 0240
    0241 0243 0248 0254 0256 0257 0259 0265 0268 0275 0276 0277 0288 0289
    0304 0305 0306 0310 0321 0324 0328 0329 0331 0335 0337 0338 0347 0350
    0352 0353 0355 0356 0358 0360 0365 0373 0375 0386 0391 0396 0397 0399
    0401 0402 0403 0412 0413 0415 0427 0431 0437 0438 0441 0442 0445 0453
    0460 0462 0463 0467 0469 0472 0473 0476 0482 0485 0487 0493 0500 0501
    0507 0523
    """,
    """
    0001 0002 0007 0016 0017 0018 0021 0022 0027 0028 0031 0034 0035 0038
    0039 0041 0046 0052 0059 0062 0064 0065 0072 0073 0074 0081 0088 0089
    0092 0095 0102 0106 0107 0109 0122 0123 0128 0131 0136 0140 0142 0144
    0145 0149 0156 0160 0165 0169 0172 0175 0176 0178 0182 0183 0184 0187
    0190 0191 0192 0203 0205 0207 0210 0221 0222 0227 0245 0253 0258 0263
    0270 0278 0279 0280 0281 0293 0296 0300 0301 0303 0308 0312 0316 0319
    0320 0334 0336 0339 0343 0344 0357 0359 0361 0368 0369 0371 0374 0376
    0378 0381 0382 0388 0394 0400 0406 0414 0421 0423 0426 0436 0439 0443
    0444 0446 0449 0456 0465 0466 0470 0471 0479 0480 0481 0484 0494 0498
    0504 0514 0517 0522
    """,
    """
    0006 0008 0009 0011 0023 0025 0026 0032 0033 0040 0043 0044 0055 0056
    0060 0067 0068 0069 0070 0078 0079 0091 0094 0100 0105 0108 0112 0113
    0116 0125 0129 0141 0147 0150 0151 0152 0153 0158 0161 0162 0166 0167
    0171 0179 0180 0188 0195 0197 0201 0206 0214 0217 0219 0223 0230 0232
    0234 0235 0236 0237 0242 0244 0249 0251 0252 0260 0261 0262 0264 0267
    0272 0274 0284 0285 0286 0291 0294 0297 0299 0313 0314 0315 0317 0322
    0325 0326 0327 0330 0351 0354 0362 0363 0366 0372 0377 0380 0390 0393
    0395 0409 0417 0418 0420 0424 0428 0430 0432 0433 0434 0451 0452 0455
    0458 0461 0474 0483 0486 0488 0489 0495 0496 0497 0502 0505 0506 0508
    0509 0513 0516 0518 0519 0520 0521
    """,
    """
    0003 0005 0010 0015 0024 0030 0036 0047 0048 0054 0058 0061 0071 0075
    0076 0077 0084 0085 0093 0098 0099 0101 0104 0111 0114 0117 0120 0121
    0124 0127 0134 0137 0138 0139 0146 0154 0155 0157 0163 0164 0168 0173
    0174 0177 0181 0185 0194 0196 0198 0199 0200 0204 0211 0212 0220 0224
    0229 0239 0246 0247 0250 0255 0266 0269 0271 0273 0282 0283 0287 0290
    0292 0295 0298 0302 0307 0309 0311 0318 0323 0332 0333 0340 0341 0342
    0345 0346 0348 0349 0364 0367 0370 0379 0383 0384 0385 0387 0389 0392
    0398 0404 0405 0407 0408 0410 0411 0416 0419 0422 0425 0429 0435 0440
    0447 0448 0450 0454 0457 0459 0464 0468 0475 0477 0478 0490 0491 0492
    0499 0503 0510 0511 0512 0515
    """ )
long_test_files = tuple( tuple( f'cs{n}.cif' for n in e.split() )
                         for e in _long_numbers )

def all_files():
    d = pathlib.Path( get_named_test_data_dir( 'cifsynth' ) )
    return sorted( d.glob( 'cs*.cif' ) )

def _errmsg( e ):
    msg = f'{e.__class__.__name__}: {e}'
    #Messages printed without details which depend on versions. Remove the
    #first when requiring gemmi>=0.6.6 (changed format of parse errors), the
    #second when requiring spglib>=2.7 (older versions give no reason):
    for p in ( 'CIF parsing error (from Gemmi)',
               'spglib could not analyse the crystal structure' ):
        if p in msg:
            msg = msg[:msg.index(p)+len(p)]
    return msg

def _load( f, permissive ):
    with warnings.catch_warnings( record = True ) as wl:
        warnings.simplefilter( 'always' )
        try:
            lc = nc_cifutils.CIFLoader( f, quiet = True,
                                        permissive = permissive )
        except NC.NCException as e:
            return None, _errmsg( e ), []
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
        return f'material: {_errmsg( e )}'
    hkl = list( info.hklObjects() )
    dmax = max( h.d for h in hkl ) if hkl else None
    return ( f'material: {len(hkl)} hkl families'
             + ( f', dmax={dmax:.5g}' if dmax else '' ) )

def check_file( f, with_material ):
    #spglib<2.5 (e.g. 2.3.1 in Ubuntu 24.04) prints warnings from C code to
    #stderr, which are filtered out here (remove when requiring spglib>=2.5):
    import os
    import sys

    import spglib
    old_spglib = tuple( int(e) for e in spglib.__version__.split('.')[:2] ) < (2,5)
    if not old_spglib or os.name != 'posix':
        return _check_file( f, with_material )
    import contextlib
    import ctypes
    import tempfile
    libc = ctypes.CDLL( None )
    streams = ( ( 1, sys.stdout ), ( 2, sys.stderr ) )
    def flush():
        for _, s in streams:
            s.flush()
        libc.fflush( None )
    flush()
    with contextlib.ExitStack() as stack:
        saved = [ ( fd, os.dup( fd ), s,
                    stack.enter_context( tempfile.TemporaryFile() ) )
                  for fd, s in streams ]
        for fd, _, _, tmp in saved:
            os.dup2( tmp.fileno(), fd )
        try:
            _check_file( f, with_material )
        finally:
            flush()
            for fd, orig, s, tmp in saved:
                os.dup2( orig, fd )
                os.close( orig )
                tmp.seek( 0 )
                for line in tmp.read().decode().splitlines( keepends = True ):
                    if not line.startswith( 'spglib:' ):
                        s.write( line )
            flush()

def _check_file( f, with_material ):
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

def run_long( part ):
    """Test the files of long test number part (0-based)."""
    allnames = sorted( f.name for f in all_files() )
    listed = sorted( n for e in long_test_files for n in e )
    assert listed == allnames, 'long_test_files must list each file once'
    d = pathlib.Path( get_named_test_data_dir( 'cifsynth' ) )
    run( [ d / n for n in sorted( long_test_files[part] ) ] )
