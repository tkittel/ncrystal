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


# NEEDS: gemmi

# Tests of _cifsymimpl.choose_symmetry: which symmetry to use for CIF data,
# from the listed symmetry operations and the space group symbol, in default
# and permissive mode.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev._cifsymimpl import choose_symmetry
import gemmi
import re

def no_gemmi_details( msg ):
    #(without details from gemmi, which might change between versions):
    return re.sub( r'could not be parsed \([^)]*\)',
                   'could not be parsed (...)', msg )

ops_p21c = [ 'x,y,z', '-x,y+1/2,-z+1/2', '-x,-y,-z', 'x,-y+1/2,z+1/2' ]
ops_p1bar = [ 'x,y,z', '-x,-y,-z' ]
#P-1 with the inversion centre at (1/4,1/4,1/4), not a tabulated setting:
ops_p1bar_shifted = [ 'x,y,z', '-x+1/2,-y+1/2,-z+1/2' ]
sg_p21c = gemmi.find_spacegroup_by_name( 'P 1 21/c 1' )
sg_p1bar = gemmi.find_spacegroup_by_name( 'P -1' )
sg_p1 = gemmi.find_spacegroup_by_number( 1 )

def run( descr, listed, symsg, *, load_sg = None, has_mult = True,
         permissive = False, cifdata = 'data_x', expect = None ):
    warnings = []
    try:
        c = choose_symmetry( gemmi, listed, symsg, load_sg, has_mult = has_mult,
                             permissive = permissive, cifdata = cifdata,
                             warn = warnings.append )
    except NC.NCBadInput as e:
        print( f'{descr}: rejected: {no_gemmi_details( str( e ) )}' )
        for w in warnings:
            print( '     (warning before rejection:', no_gemmi_details( w ), ')' )
        assert expect == 'error', descr
        return None
    sg = c['sg']
    print( f'{descr}: space group {sg.xhm() if sg else None},'
           f' expand with {len(c["expand_ops"]) if c["expand_ops"] else 0}'
           f' listed operations, require_mult={c["require_mult"]}' )
    for w in warnings:
        print( '     warning:', no_gemmi_details( w ) )
    assert expect != 'error', descr
    if expect is not None:
        assert ( sg.number if sg else None, c['expand_ops'] is not None,
                 c['require_mult'], len( warnings ) ) == expect, ( descr, expect )
    return c

def main():
    #(expect: space group number, expanded, require_mult, number of warnings)
    run( 'no listed operations', None, None, load_sg = sg_p21c,
         expect = ( 14, False, False, 0 ) )
    run( 'listed = symbol', ops_p21c, sg_p21c, expect = ( 14, False, False, 0 ) )
    run( 'listed twice (pymatgen duplicates)', ops_p21c * 2, sg_p21c,
         expect = ( 14, False, False, 0 ) )
    run( 'listed, no symbol', ops_p21c, None, expect = ( 14, False, False, 0 ) )
    #Tabulated, but another space group than the symbol:
    run( 'listed P-1, symbol P21/c', ops_p1bar, sg_p21c, expect = 'error' )
    run( 'listed P-1, symbol P21/c, permissive', ops_p1bar, sg_p21c,
         permissive = True, expect = ( 2, False, True, 1 ) )
    run( 'listed P-1, symbol P21/c, permissive, no multiplicities', ops_p1bar,
         sg_p21c, permissive = True, has_mult = False, expect = 'error' )
    #Not tabulated, same number of rotations as the symbol: expanded.
    run( 'shifted P-1, symbol P-1', ops_p1bar_shifted, sg_p1bar,
         expect = ( 1, True, False, 0 ) )
    #Not tabulated, other number of rotations:
    run( 'shifted P-1, symbol P21/c', ops_p1bar_shifted, sg_p21c,
         expect = 'error' )
    run( 'shifted P-1, symbol P21/c, permissive', ops_p1bar_shifted, sg_p21c,
         permissive = True, expect = ( 1, True, True, 1 ) )
    run( 'shifted P-1, symbol P21/c, permissive, no multiplicities',
         ops_p1bar_shifted, sg_p21c, permissive = True, has_mult = False,
         expect = 'error' )
    run( 'shifted P-1, no symbol', ops_p1bar_shifted, None, expect = 'error' )
    run( 'shifted P-1, no symbol, permissive', ops_p1bar_shifted, None,
         permissive = True, expect = ( 1, True, True, 1 ) )
    #Unparsable:
    bad = [ 'x,y,z', 'x+500/999+q,y,z' ]
    run( 'unparsable', bad, sg_p21c, expect = 'error' )
    run( 'unparsable, permissive', bad, sg_p21c, permissive = True,
         expect = ( 14, False, False, 1 ) )
    run( 'unparsable, no symbol, permissive', bad, None, permissive = True,
         expect = 'error' )
    #Incomplete (three of the four operations of P21/c):
    inc = ops_p21c[:3]
    run( 'incomplete', inc, sg_p21c, expect = 'error' )
    run( 'incomplete, permissive', inc, sg_p21c, permissive = True,
         expect = ( 14, False, False, 1 ) )
    run( 'incomplete, permissive, not in the symbol group', inc, sg_p1bar,
         permissive = True, expect = 'error' )
    run( 'incomplete, permissive, no symbol', inc, None, permissive = True,
         expect = 'error' )
    #Hint for files written by pymatgen:
    run( 'pymatgen file, other rotations', ops_p1bar_shifted, sg_p21c,
         cifdata = '# generated using pymatgen\ndata_x', expect = 'error' )
    assert sg_p1.number == 1

if __name__ == '__main__':
    main()
