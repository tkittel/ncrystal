#!/usr/bin/env python3

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

# NEEDS: numpy

# AtomInfo.sites: one site per orbit of the space group setting, consistent
# with the unit cell positions (checked by expanding each site with the
# symmetry operations of the setting from the sgsym query).

import NCTestUtils.enable_fpe # noqa F401
import numpy as np
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query

def ops_of( setting ):
    #Parse the CIF-like operation strings ("-y,x-y,z+1/3") of the sgsym query:
    ops = []
    for s in evaluate_query( ['util','sgsym','strict',setting] )['ops']:
        rot = np.zeros((3,3)); tr = np.zeros(3)
        for i, comp in enumerate( s.replace(' ','').split(',') ):
            for term in comp.replace('-','+-').split('+'):
                if not term:
                    continue
                sign = -1.0 if term.startswith('-') else 1.0
                t = term.lstrip('-')
                if t in 'xyz':
                    rot[i,'xyz'.index(t)] += sign
                else:
                    num, _, den = t.partition('/')
                    tr[i] += sign * float(num) / ( float(den) if den else 1.0 )
        ops.append( (rot,tr) )
    return ops

def pdist( a, b ):
    #largest periodic coordinate difference
    d = np.abs( np.asarray(a) - np.asarray(b) ) % 1.0
    return np.max( np.minimum( d, 1.0 - d ) )

def expand( ops, site ):
    orbit = []
    for rot, tr in ops:
        p = ( rot @ np.asarray(site) + tr ) % 1.0
        if not any( pdist(p,q) < 1e-6 for q in orbit ):
            orbit.append( p )
    return orbit

def main():
    nfiles = nsites = npos = 0
    for f in sorted( NC.browseFiles( factory = 'stdlib' ), key = lambda x: x.name ):
        try:
            info = NC.createInfo( f'{f.name};inelas=0' )
        except NC.NCException:
            continue
        if not info.hasAtomInfo() or not info.structure_info['spacegroup']:
            continue
        ops = ops_of( info.structure_info['spacegroup_setting'] )
        nfiles += 1
        for ai in info.atominfos:
            sites = ai.sites
            assert sites is not None, f'{f.name}: sites unknown'
            assert sum( s['multiplicity'] for s in sites ) == len(ai.positions)
            expanded = []
            for s in sites:
                orbit = expand( ops, s['position'] )
                assert len(orbit) == s['multiplicity'], ( f.name, s, len(orbit) )
                expanded += orbit
            #Same set of positions:
            for p in ai.positions:
                assert any( pdist(p,q) < 1e-6 for q in expanded ), ( f.name, p )
            nsites += len(sites); npos += len(ai.positions)
    print( f'Checked {nfiles} materials: {nsites} sites expanding to {npos} positions' )
    for name in ( 'Al_sg225.ncmat', 'Al2O3_sg167_Corundum.ncmat',
                  'PbO-alpha_sg129_Litharge.ncmat' ):
        info = NC.createInfo( f'{name};inelas=0' )
        print( f'{name} (setting {info.structure_info["spacegroup_setting"]}):' )
        for ai in info.atominfos:
            for s in ai.sites:
                p = s['position']
                print( f'  {ai.atomData.displayLabel():4s} ({p[0]:.6g}, {p[1]:.6g},'
                       f' {p[2]:.6g}) multiplicity {s["multiplicity"]}' )
    NC.registerInMemoryFileData( 'nosg.ncmat', """NCMAT v7
@CELL
  cubic 4.05
@ATOMPOSITIONS
  Al 0 0 0
  Al 0 1/2 1/2
  Al 1/2 0 1/2
  Al 1/2 1/2 0
@DYNINFO
  element Al
  fraction 1
  type vdosdebye
  debye_temp 400
""" )
    info = NC.createInfo( 'nosg.ncmat;inelas=0' )
    assert info.atominfos[0].sites is None
    print( 'Without space group: sites unknown' )
    #Single-phase information must fail clearly for multiphase objects (not
    #just look absent), whichever way it is requested:
    multi = NC.createInfo( 'phases<0.5*Al_sg225.ncmat&0.5*Cu_sg225.ncmat>' )
    for descr, fct in [ ('structure query', lambda : multi._query('structure')),
                        ('sites query', lambda : multi._query('sites')),
                        ('hasStructureInfo()', multi.hasStructureInfo),
                        ('structure_info', lambda : multi.structure_info),
                        ('atominfos', lambda : multi.atominfos) ]:
        try:
            fct()
        except NC.NCLogicError as e:
            print( f'Multiphase, {descr}: {e}' )
        else:
            raise SystemExit(f'Multiphase, {descr}: did not raise')

if __name__ == '__main__':
    main()
