
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

Internal implementation of creating NCMATComposer objects from loaded CIF
data, for the cifutils module.

"""

__all__ = []

from . import _common as _nc_common
from . import _fmt as _nc_fmt
from . import _ncmatimpl as _nc_ncmatimpl
from . import core as _nc_core
from . import ncmat as _nc_ncmat
from ._cifsitesimpl import _uij_keys
from ._cifsrcimpl import _codid2url, _mpid2url


def _impl_create_ncmat_composer( cifloader, *,
                                 uiso_temperature,
                                 skip_dyninfo,
                                 top_comments,
                                 fallback_debye_temp,
                                 #merge_equiv,
                                 quiet,
                                 remap,
                                 no_formula_check ):
        with _nc_common.WarningSpy( block = quiet ) as extra_ana_warnings:
            composer = _impl_create_ncmat_composer_internal( cifloader = cifloader,
                                                             #merge_equiv = merge_equiv,
                                                             uiso_temperature = uiso_temperature,
                                                             skip_dyninfo = skip_dyninfo,
                                                             top_comments = top_comments,
                                                             remap = remap,
                                                             no_formula_check = no_formula_check )
            if not skip_dyninfo and fallback_debye_temp and fallback_debye_temp > 0.0:
                composer.allow_fallback_dyninfo( fallback_debye_temp )

        def fmtwarnings( warnings, descrtxt, nwmax ):
            ll=[]
            if not warnings:
                return ll
            ll.append(f'Notice: The following WARNINGS were emitted {descrtxt}:')
            ll.append('')

            wleft=warnings[::-1]
            while wleft:
                w = wleft.pop()
                ncount = 1 + wleft.count(w)
                if ncount > 1:
                    wleft = [e for e in wleft if e!=w]
                wmsg,wcat = w
                wmsg = [e.strip() for e in wmsg.splitlines() if e.strip()]
                s0 = f'  {wcat} :' if ncount==1 else f'  ({ncount}x) {wcat} : '
                for i,m in enumerate(wmsg):
                    if i==nwmax:
                        break
                    pref=s0 if i==0 else ' '*(len(s0))
                    if i+1==nwmax and len(wmsg)>nwmax:
                        m=f'<{len(wmsg)-nwmax} lines of output hidden>'
                    ll.append(f'{pref} {m}')
            return ll

        composer.add_comments( fmtwarnings( list(cifloader.warnings) + extra_ana_warnings,
                                            'when loading the CIF data',nwmax=10),
                               add_empty_line_divider = True )
        return composer

def _impl_create_ncmat_composer_internal( cifloader, *, uiso_temperature, skip_dyninfo, top_comments, remap, no_formula_check ):

    src = cifloader

    src_atoms = src.atoms

    #Apply remap:
    def _stdatomname( name ):
        return {'D':'H2','T':'H3'}.get(name,name)

    remap_decoded = {}
    for elemiso, compos in (remap or []):
        elemiso, compos = _nc_ncmatimpl._decode_composition(elemiso,compos)
        elemiso = _stdatomname( elemiso )
        remap_decoded[ elemiso ] = compos
    remap = remap_decoded
    remap_str = ''

    if remap:
        new_src_atoms = []
        for atom in src_atoms:
            if not any( e in remap for fr,e in atom['composition'] ):
                new_src_atoms.append( atom )
                continue
            newcompos = []
            for fr,e in atom['composition']:
                remap_compos = remap.get(e,None)
                if remap_compos is None:
                    newcompos.append( (fr,e) )
                else:
                    for remap_fr, remap_elem in remap_compos:
                        newcompos.append( ( fr*remap_fr, remap_elem ) )
                d = { k: v for k,v in atom.items() if k!='composition' }
                d['composition'] = newcompos
                new_src_atoms.append( d )
        src_atoms = tuple( new_src_atoms )
        #For the description:
        remap_strs = []
        for elemiso, compos in sorted(remap.items()):
            compos_str = compos[0][1] if len(compos)==1 else ' '.join( f'{_f:g} {_n}' for _f,_n in compos)
            remap_strs.append(f'{elemiso} -> {compos_str}')
        remap_str = ', '.join(remap_strs)

    #DYNINFO analysis can only be done now, where we know the skip_dyninfo and uiso_temperature values:

    ncmat = _nc_ncmat.NCMATComposer()
    for idx, a in enumerate( src_atoms ):
        lbl = f'cif_species_{idx}'
        compos = [ (fr,e) for fr,e in a['composition'] ]
        if len(compos)==1 and compos[0][0]==1.0:
            composstr = compos[0][1]
        else:
            composstr = ' '.join(f'{f:g} {n}' for f,n in compos)
        ncmat.set_composition( lbl, compos )
        if not skip_dyninfo:
            uiso = a.get('uiso',None)
            if uiso is not None:
                assert uiso>0.0
                if uiso_temperature is None:
                    _nc_common.warn(f'ignoring uiso info present in CIF input for "{composstr}" since uiso_temperature parameter value is not provided')
                else:
                    ncmat.set_dyninfo_msd( lbl, msd=uiso, temperature=uiso_temperature )


    #Sites (uij must be given for all or none of the sites of a label):
    site2atom = src._site2atom
    lbl_has_uij = {}
    for s, iatom in zip( src.sites, site2atom ):
        lbl_has_uij.setdefault( iatom, set() ).add( s['aniso'] is not None )
    sites = []
    for s, iatom in zip( src.sites, site2atom ):
        e = ( f'cif_species_{iatom}', ) + tuple( s['position'] )
        occu = src_atoms[iatom].get('occupancy',1.0)
        e += ( occu, ) if occu != 1.0 else ()
        if s['aniso'] and lbl_has_uij[iatom] == { True }:
            e += ( tuple( s['aniso'][k] for k in _uij_keys ), )
        sites.append( e )
    for iatom, has in sorted( lbl_has_uij.items() ):
        if len( has ) > 1:
            _nc_common.warn( 'Anisotropic displacements ignored for CIF sites'
                             f' {", ".join( src_atoms[iatom]["cif_labels"] )},'
                             ' since they are not available for all sites'
                             ' with the same NCMAT label' )

    st = src.cellsg
    ncmat.set_cellsg( a=st['a'], b=st['b'], c=st['c'],
                      alpha=st['alpha'], beta=st['beta'], gamma=st['gamma'],
                      spacegroup=src._sgsetting )
    ncmat.set_atomsites( sites )

    if top_comments is not None:
        ncmat.add_comments(top_comments)

    _extracted_description = src.extracted_description
    if not _extracted_description:
        if src.cifsrc.name:
            _ds = 'CIF data' if src.cifsrc.is_remote else 'CIF file'
            _extracted_description = [f'{_ds}: {src.cifsrc.name}']
        else:
            _extracted_description = ['Anonymous CIF data']
    if _extracted_description:
        ncmat.add_comments([('Structure converted (with NCrystal'
                             '.cifutils module) from:'),''])

        _thedescr = _extracted_description
        if remap_str:
            _codidurl = _codid2url(src.actual_codid) if src.actual_codid else None
            _mpidurl = _mpid2url(src.actual_mpid) if src.actual_mpid else None
            if _codidurl and _codidurl in _thedescr:
                _=f'{_codidurl} [with {remap_str}]'
                _thedescr = [e.replace(_codidurl,_) for e in _thedescr]
            elif _mpidurl and _mpidurl in _thedescr:
                _=f'{_mpidurl} [with {remap_str}]'
                _thedescr = [e.replace(_mpidurl,_) for e in _thedescr]
            else:
                _thedescr+=['',f'Note: With custom remappings: {remap_str}']

        ncmat.add_comments(['  '+e for e in _thedescr])
        ncmat.add_comments([''])

    ncmat.add_comments('IMPORTANT NOTICE: This is a mostly automatic conversion which has not been',add_empty_line_divider=True)
    ncmat.add_comments('                  verified!  In particular the @DYNINFO sections might need')
    ncmat.add_comments('                  post-editing. Before distributing this file to other people,')
    ncmat.add_comments('                  please review this, amend the comments here to document,')
    ncmat.add_comments('                  anything done, and remove this notice.')

   #formula
    def _extract_atomcount( atom ):
        occu, npos = atom['occupancy'], len(atom['equivalent_positions'])
        return npos * ( int(occu) if occu == int(occu) else float(occu) )
    total_composition = []
    for a in src_atoms:
        atomcount = _extract_atomcount( a )
        for elemfrac, eleminfo in a['composition']:
            total_composition.append( ( eleminfo, elemfrac * atomcount ) )
    formula = _nc_fmt.format_chemform( total_composition )


    expected_formula = cifloader.raw_cif_chemformula

    if expected_formula:
        def formula_to_dict( f ):
            if not f:
                return None
            if isinstance(f,str):
                return formula_to_dict( _decode_formula(f) )
            if isinstance(f,dict):
                return f
            res = {}
            for lbl,count in f:
                assert count>=0.0
                if count==0:
                    continue
                if lbl in res:
                    res[lbl] += count
                else:
                    res[lbl] = count
            return res or None
        from .atomdata import allElementNames
        _all_elements = allElementNames()
        def _decode_formula( s ):
            l1 = [e for e in _all_elements if len(e)==1]+['D','T']
            l2 = [e for e in _all_elements if len(e)==2]
            res = []
            while s.strip():
                s=s.strip()
                c = [ e for e in l2 if s.startswith(e) ]
                if not c:
                    c = [ e for e in l1 if s.startswith(e) ]
                if len(c)!=1:
                    return None
                c = c[0]
                s=s[len(c):].strip()
                leaddigits = ''
                while s and s[0].isdigit():
                    leaddigits += s[0]
                    s=s[1:]#nb: no strip here, we want to stop on a space!
                n = int(leaddigits) if leaddigits else 1
                res.append( (c,n) )
            return formula_to_dict(res)
        def formulas_incompatible( formula1, formula2 ):
            f1 = formula_to_dict(formula1)
            f2 = formula_to_dict(formula2)
            if not f1 or not f2:
                return False
            #map D,T,H2,H3 -> H to ensure fewer false positives:
            def _tmp(d):
                return { ('H' if k in ('D','T','H2','H3') else k): v
                         for k,v in d.items() }
            f1,f2 = _tmp(f1),_tmp(f2)
            if f1 == f2:
                return False
            if set(f1.keys()) != set(f2.keys()):
                return True
            kref = next(iter(f1.keys()))
            assert f1[kref] > 0.0 and f2[kref]>0.0
            f2scale = f1[kref]/f2[kref]
            tol = 1e-2#not tighter to avoid false positives!
            for k,v1 in f1.items():
                v2 = f2[k]*f2scale
                if abs(v1-v2)>tol and abs(v1-v2)/(1e-300+abs(v1)+abs(v2)) > tol:
                    return True
            return False

        if not no_formula_check and expected_formula:
            actual_expected_formula_dict = formula_to_dict( expected_formula )
            if remap:
                newf = {}
                def _collect( _name, _frac ):
                    _name = _stdatomname(_name)
                    if _name not in newf:
                        newf[_name] = _frac
                    else:
                        newf[_name] += _frac
                for name, frac in actual_expected_formula_dict.items():
                    _remapped = remap.get(_stdatomname(name),None)
                    if _remapped is not None:
                        for remap_frac, remap_name in _remapped:
                            _collect( remap_name, remap_frac * frac )
                    else:
                        _collect( name, frac )
                actual_expected_formula_dict = newf
            if formulas_incompatible(total_composition,actual_expected_formula_dict):
                s = f'"{expected_formula}"'
                if remap:
                    _ = _nc_fmt.format_chemform( sorted(actual_expected_formula_dict.items() ) )
                    s += f' remapped to "{_}"'
                raise _nc_core.NCBadInput(f'Formula encoded in CIF data ({s}) is not compatible with formula of loaded structure ("{formula}")')

    return ncmat

def _suggest_filename( ncmat_metadata, cifloader ):
    ll = ['autogen']
    ll.append( ncmat_metadata['chemform'] )
    sgnum = ncmat_metadata.get('cellsg',{}).get('spacegroup',None)
    if sgnum:
        ll.append( f'sg{int(sgnum)}' )
    if cifloader.actual_codid:
        ll.append( f'cod{cifloader.actual_codid}' )
    if cifloader.actual_mpid:
        ll.append( f'mp{cifloader.actual_mpid}' )
    return '_'.join(ll) + '.ncmat'
