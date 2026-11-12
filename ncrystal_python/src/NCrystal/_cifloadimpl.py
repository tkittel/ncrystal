
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
from ._cifreadimpl import (
    _extract_descr_from_cif,
    _gemmi_wrap_to_unit,
    _import_gemmi,
    _load_with_gemmi,
    _update_spacegroup_in_cifblock_from_name,
    ueq_from_cif,
    uiso_from_cif,
)
from ._cifsitesimpl import _impl_find_sites, _impl_group_sites
from ._cifsymimpl import (
    _check_multiplicities,
    _cif_multiplicities,
    _expand_sites_with_ops,
    _listed_ops,
    _op_key,
    _symbol_spacegroup,
    _symop_error,
)


def _actual_init_gemmicif( cifsrc, *, quiet, mp_apikey, refine_with_spglib,
                           merge_equiv, tols, permissive,
                           override_spacegroup = None ):

    from itertools import chain
    from math import fsum as _math_fsum

    from .cifutils import CIFSource

    result = {}
    pos_tolerance = 0.0001#NB: Best if value matches the one in NCInfoBuilder.cc

    from types import MappingProxyType
    dict_ro, list_ro = MappingProxyType, tuple

    gemmi, gemmi_cif = _import_gemmi()

    if not isinstance( cifsrc, CIFSource ):
        cifsrc = CIFSource( cifsrc )

    result['cifsrc'] = cifsrc

    _cifdata = cifsrc.load_data( quiet = quiet, mp_apikey = mp_apikey )
    if not quiet:
        _nc_common.print("Attempting to load CIF data with gemmi")

    result['cifdata'] = _cifdata
    try:
        cif_doc = gemmi_cif.read_string( _cifdata )
    except ValueError as e:
        errmsg = f'CIF parsing error (from Gemmi): "{e}"'
        cif_doc = None
    if cif_doc is None:
        raise _nc_core.NCBadInput(errmsg or 'Unknown CIF parsing error from Gemmi')

    cif_doc_as_json = cif_doc.as_json()
    import json
    try:
        cif_doc_as_dict = json.loads( cif_doc_as_json )
    except json.decoder.JSONDecodeError as e:
        cif_doc_as_dict = None
        _nc_common.warn(f'Could not decode raw CIF data to dictionary (Gemmi bug?): JSONDecodeError("{e}")')

    result['cif_raw'] = cif_doc_as_dict
    if len(cif_doc) == 0:
        raise _nc_core.NCBadInput('CIF data had no blocks!')
    if len(cif_doc) == 1:
        cif_block_with_structure = cif_doc.sole_block()
    else:
        cif_block_with_structure = None
        for e in cif_doc:
            if not len(e.find(['_atom_site_fract_x'])):
                continue
            if cif_block_with_structure:
                cif_block_with_structure = None
                break
            cif_block_with_structure = e
        if not cif_block_with_structure:
            raise _nc_core.NCBadInput('Could not automatically determine block in CIF data with structure info')

    assert cif_block_with_structure

    if override_spacegroup:
        assert isinstance(override_spacegroup,str)
        _nc_common.warn(f'Overriding spacegroup to "{override_spacegroup}" due to explicit request')
        _update_spacegroup_in_cifblock_from_name( gemmi, cif_block_with_structure, override_spacegroup )

    #record original h-m-alt entry:
    orig_hm_alt = None
    _ = cif_block_with_structure.find(['_space_group_name_H-M_alt'])
    if _ and len(_)==1 and len(_[0])==1:
        orig_hm_alt = str(_[0][0]).strip()

    struct, sg = _load_with_gemmi( cif_block_with_structure )

    if not struct:
        raise _nc_core.NCBadInput('Could not load structure with Gemmi')

    #Decide which symmetry to use, from the listed symmetry operations and the
    #space group symbol. Inconsistent information is an error, unless
    #permissive (and even then some evidence is required):
    blk = cif_block_with_structure
    mult = _cif_multiplicities( gemmi, blk )
    p1_sites, require_mult, expect_sgno = None, False, None
    listed = None if override_spacegroup else _listed_ops( gemmi, blk )
    symsg = _symbol_spacegroup( gemmi, struct, blk ) if listed else None
    def permissive_warn( what ):
        _nc_common.warn( f'{what} (accepted due to permissive mode)' )
    if listed is None:
        pass#symbol only (or override), as determined by _load_with_gemmi
    elif listed[0] == 'unparsable':
        why = f'could not be parsed ({listed[1]})'
        if not permissive or not symsg:
            raise _symop_error( why, _cifdata, symsg and not permissive )
        permissive_warn( f'The symmetry operations listed in the CIF data {why}'
                         f', so the space group symbol ({symsg.xhm()}) is used' )
        sg = symsg
    elif listed[0] == 'incomplete':
        why = 'are incomplete (they do not form a group)'
        ok = symsg and { _op_key( o ) for o in listed[1] } <= {
            _op_key( o ) for o in symsg.operations() }
        if not permissive or not ok:
            raise _symop_error( why, _cifdata, ok and not permissive )
        permissive_warn( f'The symmetry operations listed in the CIF data {why}'
                         ' but belong to the group of the space group symbol'
                         f' ({symsg.xhm()}), which is used' )
        sg = symsg
    else:
        _, ops, go = listed
        tab = gemmi.find_spacegroup_by_ops( go )
        if tab:
            sg = tab
            if symsg and symsg.number != tab.number:
                why = ( f'are those of space group {tab.xhm()}, but the space'
                        f' group symbol is {symsg.xhm()}' )
                if not permissive:
                    raise _symop_error( why, _cifdata, True )
                permissive_warn( f'The symmetry operations listed in the CIF'
                                 f' data {why}. The listed operations are used' )
                require_mult = True
        else:
            nrot = len( symsg.operations().sym_ops ) if symsg else None
            if nrot != len( go.sym_ops ):
                why = ( f'have {len(go.sym_ops)} rotations, but space group'
                        f' {symsg.xhm()} has {nrot}' if symsg else
                        'are not those of a tabulated space group setting, and'
                        ' can not be checked against a space group symbol' )
                if not permissive:
                    raise _symop_error( why, _cifdata, True )
                permissive_warn( f'The symmetry operations listed in the CIF'
                                 f' data {why}. The listed operations are used' )
                require_mult = True
            #Not a tabulated setting (e.g. a primitive cell of a centred
            #lattice, as written by pymatgen), so expand the sites with them:
            if require_mult and not mult:
                raise _nc_core.NCBadInput( 'Using the listed symmetry'
                                           ' operations requires site'
                                           ' multiplicities in the CIF data' )
            #(the multiplicities must be true ones, as written by pymatgen):
            p1_sites = _expand_sites_with_ops( gemmi, struct, ops, mult,
                                               tols['symprec'], strict = True )
            sg = gemmi.find_spacegroup_by_number( 1 )
            if not quiet:
                _nc_common.print('Expanded sites with the symmetry operations'
                                 ' listed in the CIF data, since they do not'
                                 ' correspond to a tabulated space group'
                                 ' setting')
    if require_mult and not mult:
        raise _nc_core.NCBadInput( 'Using the listed symmetry operations'
                                   ' requires site multiplicities in the CIF'
                                   ' data' )
    if sg and p1_sites is None and struct.spacegroup != sg:
        struct.spacegroup = sg
        struct.setup_cell_images()

    if not sg:
        raise _nc_core.NCBadInput('Could not determine space group from CIF data')

    #warn if spacegroup setting might be ambiguous:
    if sg.number != 1:
        _ = sorted( (e.number,e.xhm(),e.is_reference_setting()) for e in gemmi.spacegroup_table() if ( e.number==sg.number and ':' in e.xhm()))
        #_ = [ (no,xhm,isref) for no,xhm,isref in _ if isref ]
        if len(_)>1:
            _str = '", "'.join( xhm for no,xhm,isref in _)
            if not orig_hm_alt or not any( (orig_hm_alt==e[1] or orig_hm_alt.replace(' ','')==e[1].replace(' ','')) for e in _ ):
                _nc_common.warn(f'SG-{sg.number} available in multiple'
                                f' choices ("{_str}") and which one was not'
                                ' encoded explicitly in the '
                                '_space_group_name_H-M_alt CIF field. Consider'
                                ' overriding the space group explicitly when'
                                ' loading this file.')

    _ = cif_block_with_structure.find(['_atom_type_number_in_cell'])
    #was: _ = sum((sum(([e] for e in ll),[]) for ll in _),[]) if _ else []
    #but replaced with:
    _ = list(chain.from_iterable(_)) if _ else []
    if _ and not any( e is None for e in _ ):
        expected_tot_atom_in_orig_cell = sum(float(e) for e in _)
    else:
        expected_tot_atom_in_orig_cell = None
    if not ( isinstance(sg.number,int) or sg.number.isdigit() ) or not ( 1<=int(sg.number)<=230 ):
        raise _nc_core.NCBadInput(f'Could not determine space group from CIF data (it loaded with invalid SG number: {sg.number}')
    if not struct.cell.is_crystal():
        raise _nc_core.NCBadInput('Loaded structure is non-crystalline (according to Gemmi).')

    sgnumber = int(sg.number)

    if not struct.cell.is_compatible_with_spacegroup( sg ):
        #The cell does not quite have the symmetry of the space group (e.g. due
        #to rounding). Expand the sites with its operations and refine. This is
        #only accepted if spglib then finds the same space group:
        if not refine_with_spglib:
            raise _nc_core.NCBadInput('Loaded unit cell is not compatible with deduced space group (according to Gemmi).')
        expect_sgno = int( sg.number )
        p1_sites = _expand_sites_with_ops( gemmi, struct, list( sg.operations() ),
                                           mult, tols['symprec'],
                                           strict = require_mult )
        sg = gemmi.find_spacegroup_by_number( 1 )
        sgnumber = 1

    cellsg = { 'a': struct.cell.a, 'b': struct.cell.b, 'c': struct.cell.c,
               'alpha': struct.cell.alpha,
               'beta': struct.cell.beta,
               'gamma': struct.cell.gamma }

    fractcoord_approx_1angstrom = 1.0 / ( (struct.cell.a+struct.cell.b+struct.cell.c)/3.0 )


    #Although the sg object also provides sg.ccp4 and sg.hall, we record here
    #just sgnumber and xhm, since that is what spglib refinement also provides:
    cellsg['spacegroup'] = dict_ro( { 'number': sgnumber,
                                      'hm': sg.xhm() } )

    collected_atoms = []
    cell_images = [] if p1_sites is not None else struct.cell.images

    def _expand_coord_to_all_other_images( coord, allow_finetune = 3 ):
        pos0 = _gemmi_wrap_to_unit(coord)
        #First find all brute-force expanded coords (ignore those within machine
        #precision of each other!):
        ll = [ pos0 ]
        for candidate in ( _gemmi_wrap_to_unit(img.apply( pos0 )) for img in cell_images):
            use = True
            for c in ll:
                if _nc_sgsymimpl._unit_cell_point_dist(candidate,c) < 1e-10:
                    use = False
                    break
            if use:
                ll.append( candidate )
        if not allow_finetune:
            return ll
        #Now, check how many of these are very close to the initial point:
        _ucpdist = _nc_sgsymimpl._unit_cell_point_dist
        lclose = [ _nc_sgsymimpl._remap_fract_pos_pt(e) for e in ll
                   if _ucpdist(pos0,e) < 0.01*fractcoord_approx_1angstrom ]
        assert len(lclose) > 0
        if len(lclose) == 1:
            #no issues, just return:
            return ll
        #Input might have had inexact coordinates for atoms at special
        #positions. Fine-tune pos0 as average over the close points, and rerun:
        def _fract_delta( x1, x0 ):
            #remember that distance of (0,0,eps) and (0,0,1-eps) is 2eps, not 1-2eps.
            dx = x1 - x0
            if dx > 0.5:
                dx -= 1.0
            if dx < -0.5:
                dx += 1.0
            return dx
        def _fract_delta_pt( xyz1, xyz0 ):
            return tuple( _fract_delta( xyz1[i], xyz0[i] ) for i in range(3) )

        pos0_new = [ _math_fsum( [lclose[0][i]] + [ _fract_delta(c[i],lclose[0][i])/len(lclose) for c in lclose ] ) for i in range(3) ]
        #snap to (0,0,0):
        for i in range(3):
            if abs(pos0_new[i])< 1e-16:
                pos0_new[i] = 0.0
        def fmt( c ):
            return f'({c[0]:.15g},{c[1]:.15g},{c[2]:.15g})'
        _nc_common.warn(f'Fractional coordinate {fmt(coord)} interpreted as special position {fmt(pos0_new)} to avoid numerical precision issues')
        return _expand_coord_to_all_other_images( gemmi.Fractional(pos0_new[0],pos0_new[1],pos0_new[2]), allow_finetune = (allow_finetune-1) )


    def _image_rotations( coords ):
        #Rotation parts of the operations mapping coords[0] to each coord:
        ops = [ None ] + list( cell_images )
        res = []
        for x in coords:
            for img in ops:
                y = coords[0] if img is None else _gemmi_wrap_to_unit(
                    img.apply( coords[0] ) )
                if _nc_sgsymimpl._unit_cell_point_dist( y, x ) < 1e-6:
                    res.append( ( ( 1, 0, 0 ), ( 0, 1, 0 ), ( 0, 0, 1 ) )
                                if img is None else img.mat.tolist() )
                    break
            else:
                raise _nc_core.NCLogicError('Could not find symmetry'
                                            ' operation of image')
        return res

    #Uiso and Ueq from the printed values, with their rounding errors:
    uiso_printed = uiso_from_cif( gemmi, blk )
    ueq_printed = ueq_from_cif( gemmi, blk )

    for site in ( p1_sites if p1_sites is not None else struct.sites ):
        if not float( site.occ ) > 0.0:
            continue
        pos0 = _gemmi_wrap_to_unit( site.fract )
        expanded_coords = _expand_coord_to_all_other_images( pos0 )

        aniso = None
        if site.aniso.nonzero:
            aniso = { k: float(getattr(site.aniso,k)) for k in ('u11','u22','u33','u12','u13','u23') }
            if all( v==0.0 for v in aniso.values() ):
                aniso = None
        u_iso = float( site.u_iso )
        if not u_iso > 0.0:
            u_iso = None

        if aniso:
            #Resolve Uiso versus Ueq here, so each site carries one magnitude
            #(the tensor is rescaled to it). Accepted differences are those
            #explained by the rounding errors of the printed values, plus
            #tol_uij:
            sgs = _nc_sgsymimpl
            ueq = sgs._cif_uij_ueq( aniso, cellsg )
            lbl = str( site.label )
            allow = ( uiso_printed.get( lbl, ( None, 0.0 ) )[1]
                      + ueq_printed.get( lbl, ( None, 0.0 ) )[1]
                      + tols['uij'] * max( ueq, u_iso or 0.0 ) )
            if not u_iso:
                u_iso = ueq
            elif abs( ueq - u_iso ) > allow:
                msg = ( f'CIF site "{site.label}" has Uiso={u_iso:g} which'
                        f' differs from the Ueq={ueq:g} of its anisotropic'
                        ' displacements by more than the rounding errors of'
                        ' the printed values explain' )
                if not permissive:
                    raise _nc_core.NCBadInput( msg + ' (try permissive mode'
                                               ' to use the Uiso value)' )
                _nc_common.warn( msg + ' (using the Uiso value due to'
                                 ' permissive mode)' )
            aniso = sgs._rescale_uij( aniso, ueq, u_iso )

        site_en = str(site.element.name).strip()
        if site_en == 'X' or not site_en:
            #Issue seen in codid::9005777, which had no atomic symbols but
            #labels like "CaM1", "CaM2", "OA1".."OA3", "OB1", "OB2", "OC1",
            #"OC3",... Correct if beginning of string starts with exactly one
            #element name, followed by either a digit or an upper case string
            #followed by digits.
            from .atomdata import allElementNames
            def guess_elem_name( site_label ):
                sl, traildigit = _nc_common._split_trailing_digit( site_label )
                if traildigit is None:
                    return
                #first check element names with 2 chars:
                candidates = [ e for e in allElementNames() if ( sl.startswith(e) and len(e)==2 ) ]
                if len(candidates) > 1:
                    return
                #Then those with 1 char if no hit already (also check for 'D' and 'T'):
                if not candidates:
                    candidates = [ e for e in allElementNames() if ( sl.startswith(e) and len(e)==1 ) ]
                    candidates += [ e for e in ('D','T') if ( sl.startswith(e) and len(e)==1 ) ]
                if len(candidates)!=1:
                    return
                c = candidates[0]
                assert site_label.startswith(c)
                leftover = site_label[len(c):].strip()
                while leftover and leftover[0].isupper()  and leftover[0].isalpha():
                    leftover = leftover[1:].strip()
                if not leftover or leftover.isdigit():
                    return c
            elem_name = guess_elem_name( site.label )
            if elem_name:
                _nc_common.warn(f'Assuming atomic symbol "{elem_name}" based on CIF label "{site.label}"')
            else:
                raise _nc_core.NCBadInput(f'Neither Gemmi nor NCrystal\'s custom'
                                          ' code could deduce an atomic symbol'
                                          f' for the CIF label "{site.label}"')
            elem_Z = None
        else:
            elem_name, elem_Z = site_en, int(site.element.atomic_number)

        elem_A = None
        if elem_name == 'D':
            elem_name, elem_Z, elem_A = 'H',1,2
        if elem_name == 'T':
            elem_name, elem_Z, elem_A = 'T',1,3
        elem_marker = elem_name + ( str(elem_A) if elem_A else '' )

        from .atomdata import elementNameToZValue
        z2name = elementNameToZValue( elem_name, allow_isotopes = False )
        if z2name is None:
            raise _nc_core.NCBadInput(f'Not an element name: "{elem_name}"')

        if elem_Z is None:
            elem_Z = z2name
        if elem_Z != z2name:
            raise _nc_core.NCBadInput(f'Wrong atomic number ({elem_Z}) provided'
                                      f' for element "{elem_name}"')

        collected_atoms.append( { 'expanded_coords': expanded_coords,
                                  'rotations':
                                  _image_rotations( expanded_coords ),
                                  'element_marker': elem_marker,
                                  'occupancy': float(site.occ),
                                  'cif_label': str(site.label),
                                  'uiso': u_iso,
                                  'aniso': aniso } )

    if not collected_atoms:
        raise _nc_core.NCBadInput('No atom sites found in CIF data.')

    if p1_sites is None:
        #(already checked when expanding p1_sites)
        _check_multiplicities( [ ( a['cif_label'], mult.get( a['cif_label'] ),
                                   len( a['expanded_coords'] ) )
                                 for a in collected_atoms ],
                               len( list( sg.operations() ) ),
                               strict = require_mult )

    #Now merge collected atoms which occupy the same sites:
    def has_same_sites( atom1, atom2 ):
        l1,l2 = atom1['expanded_coords'],atom2['expanded_coords']
        assert len(l1)>0 and len(l2)>0
        return len(l1) == len(l2) and any( ( _nc_sgsymimpl._unit_cell_point_dist(c2,l1[0])<pos_tolerance) for c2 in l2 )

    atoms_to_process, final_atoms = collected_atoms, []
    while atoms_to_process:
        atom1 = atoms_to_process[0]
        atoms_with_same_pos = [ atom1 ]
        atoms_with_different_pos = []
        for atom2 in atoms_to_process[1:]:
            if has_same_sites(atom1,atom2):
                if atom1['uiso'] != atom2['uiso'] or atom1['aniso'] != atom2['aniso']:
                    raise _nc_core.NCBadInput('atoms labelled "{}" and "{}" have same cell '.format(atom1['cif_label'],atom2['cif_label'])
                                              +'positions but different uiso information. This is not supported.')
                atoms_with_same_pos.append( atom2 )
            else:
                atoms_with_different_pos.append( atom2 )

        cif_labels = [ a['cif_label'] for a in atoms_with_same_pos ]
        cif_labels_str = '"{}"'.format('", "'.join(cif_labels))

        composition = [ ( a['occupancy'], a['element_marker'] ) for a in atoms_with_same_pos ]
        occupancy = _math_fsum( fr for fr, elem in composition )
        assert occupancy > 0.0
        if occupancy != 1.0:
            composition = [ (fr/occupancy, elem) for fr, elem in composition ]
            if occupancy > (1.0 + 1e-6) :
                m=f'Error: Too high total occupancy ({occupancy} which is >1) for atoms with CIF labels: {cif_labels_str}'
                if abs(int(occupancy+0.5)-occupancy)<1e-10:
                    m+='. A possible cause could be multiple entries in the CIF file using atomic positions that are actually symmetrically equivalent!'
                raise _nc_core.NCBadInput(m)
            occupancy = min( 1.0, occupancy )

        final_atoms.append( { 'cif_labels': list_ro( cif_labels ),
                              'positions': [ (c.x,c.y,c.z) for c in
                                             atoms_with_same_pos[0]['expanded_coords'] ],
                              'rotations': atoms_with_same_pos[0]['rotations'],
                              'composition': list_ro( sorted( composition ) ),
                              'occupancy': occupancy,
                              'uiso': atoms_with_same_pos[0]['uiso'],
                              'aniso': atoms_with_same_pos[0]['aniso'] } )
        atoms_to_process = atoms_with_different_pos

    if expected_tot_atom_in_orig_cell:
        n_actual = sum( len(a['positions']) for a in final_atoms )
        def rd( x,y ):
            return abs(x-y)/(max(1e-300,abs(x)+abs(y)))
        if rd( n_actual, expected_tot_atom_in_orig_cell ) > 0.01:
            raise _nc_core.NCBadInput(f'Expanded number of atoms per cell ({n_actual}) is different '
                                      'from what is stated explicitly in the CIF data'
                                      f' ({expected_tot_atom_in_orig_cell}). Wrong spacegroup (or'
                                      ' choice/setting improperly specified)?')

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
    return result, gemmi, struct
