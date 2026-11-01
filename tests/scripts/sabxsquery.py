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

# Test of the ["sab","xs",...] query, against the cross sections of the
# materials (which always use the SCT extender when an effective temperature is
# available, so the free-gas extender results are only monitored via the log).

import NCrystalDev as NC
from NCrystalDev._numpy import _np_geomspace
from NCrystalDev.misc import evaluate_query

energies = list( _np_geomspace( 1e-6, 1e3, 10 ) )

def query( cfgstr, lbl, extender ):
    r = evaluate_query( [ 'sab', 'xs', cfgstr, lbl, extender,
                          *[ str(e) for e in energies ] ] )
    assert len(r['xs']) == len(energies)
    return r

def material_xs( cfgstr ):
    sc = NC.createScatter( cfgstr + ';comp=inelas' )
    return [ sc.xsect( ekin = e ) for e in energies ]

def test( cfgstr, extender, lbls = ( '', ) ):
    rs = [ query( cfgstr, lbl, extender ) for lbl in lbls ]
    xs = [ sum( r['xs'][i] for r in rs ) for i in range(len(energies)) ]
    atoms = '' if lbls == ( '', ) else f', atoms {",".join(lbls)}'
    print(f'{cfgstr} ({extender}{atoms}):'
          f' extension_method={rs[0]["extension_method"]}')
    for e, x in zip( energies, xs ):
        print(f'  E={e:<10.4g} xs={x:.5g} barn')
    if extender == 'sct':
        for e, x, xref in zip( energies, xs, material_xs( cfgstr ) ):
            if len(lbls) == 1:
                assert x == xref, (e, x, xref)
            else:
                assert abs( x - xref ) <= 1e-14 * xref, (e, x, xref)
    assert all( r['extension_method'] == extender for r in rs )

def main():
    test( 'stdlib::Al_sg225.ncmat;knllux=3', 'sct' )
    test( 'stdlib::Al_sg225.ncmat;temp=20K;knllux=1', 'sct' )
    test( 'stdlib::Al_sg225.ncmat;knllux=3', 'freegas' )
    test( 'stdlib::Polyethylene_CH2.ncmat;knllux=2', 'sct', lbls = ( 'C', 'H' ) )
    test( 'stdlib::Polyethylene_CH2.ncmat;knllux=2', 'freegas',
          lbls = ( 'C', 'H' ) )
    for bad in ( [ 'stdlib::Al_sg225.ncmat', '', 'foo', '0.025' ],
                 [ 'stdlib::Al_sg225.ncmat', '', 'sct', 'x' ],
                 [ 'stdlib::Al_sg225.ncmat', '', 'sct' ],
                 [ 'stdlib::Polyethylene_CH2.ncmat', '', 'sct', '0.025' ] ):
        try:
            evaluate_query( [ 'sab', 'xs', *bad ] )
        except NC.NCBadInput as e:
            print(f'Expected error: {e}')
        else:
            raise SystemExit(f'Missing error for: {bad}')

if __name__ == '__main__':
    main()
