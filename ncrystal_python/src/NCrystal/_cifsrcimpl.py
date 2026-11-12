
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

Internal implementation of fetching CIF data for the cifutils module
(COD, Materials Project and the local cache).

"""

__all__ = []

from . import _common as _nc_common
from . import core as _nc_core


def _codid2url( codid ):
    return f'https://www.crystallography.net/cod/{codid}.html'

def _mpid2url( mpid ):
    return f'https://www.materialsproject.org/materials/mp-{mpid}'

def _getMaterialsProjectAPIKEY():
    import os
    apikey = os.environ.get('MATERIALSPROJECT_USER_API_KEY',None)
    if not apikey or not apikey.strip():
        raise _nc_core.NCException('ERROR: Missing API key for materialsproject.org access. To fix, '
                           +'please make sure the environment variable MATERIALSPROJECT_USER_API_KEY'
                           +' contains your personal access key that you see on'
                           +' https://www.materialsproject.org/dashboard (after logging in).')
    if len(apikey.strip()) < 25:
        raise _nc_core.NCException('ERROR: Legacy API key found in MATERIALSPROJECT_USER_API_KEY. Please instead provide a new API key from https://materialsproject.org/api (Keys for the new API are ~32 characters, whereas keys for the legacy API are ~16 characters).')
    return apikey

def _use_local_cif_cache( fn, text_data = None, quiet = False ):
    #NB: This simple implementation use no locking to guard against race
    #conditions!!! But we do perform write+move instead of simply write, which
    #is a bit more "atomic".
    if _nc_common.ncgetenv_bool('ONLINEDB_FORBID_NETWORK'):
        def notfound():
            n = _nc_common.expand_envname('ONLINEDB_FORBID_NETWORK')
            raise RuntimeError('Error: Trying to access remote DB but'
                               f' {n} is set')
        time_limit_hours = 24*7*365*1000#revisit this in 3023
    else:
        def notfound():
            pass
        time_limit_hours = 24*7
    d = _nc_common.ncgetenv('ONLINEDB_CACHEDIR')
    if not d:
        return notfound()
    import os.path
    import pathlib
    p = pathlib.Path(os.path.expanduser(d)).resolve().absolute() if d else None
    if not p:
        return notfound()
    pfn = ( p / fn )
    if text_data:
        #store data into local cache:
        if not quiet:
            _nc_common.print(f"Adding {fn} to local file cache in $NCRYSTAL_ONLINEDB_CACHEDIR")
        p.mkdir(parents=True, exist_ok=True)
        pfn_tmp = p / f'{fn}_tmp_{os.getpid()}'
        _nc_common.write_text(pfn_tmp,text_data)
        pfn_tmp.replace(pfn)
        return notfound()
    #retrieve data from local cache:
    if not pfn.exists():
        return notfound()
    #Since online data can change, entries expire eventually:
    import datetime
    timelim = datetime.timedelta( hours = time_limit_hours )
    now = datetime.datetime.now(tz=datetime.timezone.utc)
    mtime = datetime.datetime.fromtimestamp(pfn.stat().st_mtime, tz=datetime.timezone.utc)
    if ( now-mtime ) > timelim:
        #entry expired:
        if not quiet:
            _nc_common.print(f"Ignoring (and removing) expired {fn} from local file cache in $NCRYSTAL_ONLINEDB_CACHEDIR")
        pfn.unlink()
        return notfound()
    if not quiet:
        _nc_common.print(f"Getting {fn} from local file cache in $NCRYSTAL_ONLINEDB_CACHEDIR")
    return pfn.read_text()

_cod_cache = []

def _cod_get_cifdata( codid, quiet = False ):
    for _codid, _result in _cod_cache:
        if codid == _codid:
            if not quiet:
                _nc_common.print(f"Using cached Crystallography Open Database result for entry {codid}")
            return _result
    cache_fn = f'cod_{int(codid)}.cif'
    #check file cache:
    c = _use_local_cif_cache( cache_fn, quiet = quiet )
    if c:
        return c
    if not quiet:
        _nc_common.print(f"Querying the Crystallography Open Database for entry {codid}")

    mirror_urls = [
        "https://www.crystallography.net/cod/%i.cif",#canonical first
        'https://qiserver.ugr.es/cod/%i.cif',
        'http://cod.ibt.lt/cod/%i.cif',
    ]

    result = None
    while result is None:
        #Fail gently and try next mirror - except for the last attempt.
        url = mirror_urls.pop(0)%codid
        result = _nc_common.download_url( url,
                                          timeout = 10.0,
                                          quiet_network_fail = bool(mirror_urls) )
        if result or not mirror_urls:
            break
        nextdescr = '/'.join(mirror_urls[0].split('/')[0:-1])+'/'
        _nc_common.print(f'Retrival failed. Trying mirror at: {nextdescr}')

    assert result is not None

    if len(_cod_cache)==10:
        _cod_cache.pop(0)
    _cod_cache.append( (codid, result ) )
    #update file cache:
    _use_local_cif_cache( cache_fn, quiet = quiet, text_data = result )
    return result

_mp_cache = []

def _mp_get_cifdata( mpid, quiet = False, apikey = None ):
    if hasattr(mpid,'startswith') and mpid.startswith('mp-'):
        mpid = mpid[3:]
    if hasattr(mpid,'startswith') and mpid.startswith('mpid::'):
        mpid = mpid[6:]
    mpid = int(mpid)
    for _mpid, _result in _mp_cache:
        if mpid == _mpid:
            if not quiet:
                _nc_common.print(f"Using cached materialsproject.org result for entry mp-{mpid}")
            return _result
    #check file cache:
    cache_fn = f'mp_{mpid}.cif'
    c = _use_local_cif_cache( cache_fn, quiet = quiet )
    if c:
        return c

    if not quiet:
        _nc_common.print(f"Querying materialsproject.org for entry mp-{mpid}")

    if apikey is None:
        apikey = _getMaterialsProjectAPIKEY()

    with _nc_common.WarningSpy(blockfct = lambda msg, cat : cat in ('PendingDeprecationWarning','DeprecationWarning') ):
        try:
            import mp_api.client  #NB: This might trigger a spurious FPE
        except ImportError:
            raise ImportError('Could not import mp_api.client. Installing the mp-api package will most likely'
                              ' fix this (perhaps with a command like "conda install -c conda-forge mp-api"'
                              ' or "python3 -mpip install mp-api").')

    with _nc_common.WarningSpy(
            blockfct= ( lambda msg, cat : msg.lower()
                        .startswith('mpcontribs-client not installed') ) ), mp_api.client.MPRester(apikey) as mpr:
        s = mpr.get_structure_by_material_id( f'mp-{mpid}', conventional_unit_cell=True )
        #If we do no use the symprec argument in the next line, we will get a an unrefined P1 structure:
        result = s.to(fmt='cif',symprec=1e-4, significant_figures=15, angle_tolerance=5.0, refine_struct=True)
        #But we want to sanity check that this refinement gives the same spacegroup result as listed in the MP database:
        mp_expected_sg_number = s.get_space_group_info()[1]
    sg_checked = False
    errmsg = f'Unable to reliably determine spacegroup when trying to retrieve structure for mp-{mpid} from materialsproject.org'
    for ll in result.splitlines():
        p = ll.split('#',1)[0].split()
        if p and p[0]=='_symmetry_Int_Tables_number':
            if not len(p)>=2 or not p[1].isdigit():
                raise _nc_core.NCBadInput(errmsg+f' (unexpected format of line: "{ll}")')
            _sgnum = int(p[1])
            if _sgnum == mp_expected_sg_number:
                sg_checked = True
            else:
                raise _nc_core.NCBadInput(errmsg+f' (expected SG-{mp_expected_sg_number} but got SG-{_sgnum})')

    if not sg_checked:
        raise _nc_core.NCBadInput(errmsg+' (no line with "_symmetry_Int_Tables_number" produced)')

    #Finally embed origin as a note (so we can extract it later for the ncmat
    #header comments). Keep the format below synchronised with the reader
    #elsewhere in this file!!
    result += f'\n# Note from NCrystal.cifutils: Data from materialsproject.org / mp-{mpid}\n'

    #Update cache and return:
    if len( _mp_cache ) == 10:
        _mp_cache.pop(0)
    _mp_cache.append( (mpid, result ) )

    #update file cache:
    _use_local_cif_cache( cache_fn, quiet = quiet, text_data = result )

    return result
