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

# NEEDS: numpy spglib gemmi

# Tests of the sgsym component (space groups and their symmetries),
# mostly by comparing with spglib and gemmi.

from fractions import Fraction

import gemmi
import NCTestUtils.enable_fpe  # noqa: F401
import numpy as np
import spglib
from NCTestUtils.loadlib import Lib

lib = Lib('testsgsym')

# Verifies NCrystal's hardwired table of the 530 space group settings (Hall
# numbers, space group numbers, setting choice codes, Hermann-Mauguin and Hall
# symbols) against spglib and gemmi, including that the Hall symbols give the
# same symmetry operations as those in spglib's database. NCrystal follows the
# printed ITVB table in spelling the trigonal axes codes "h"/"r" (spglib and
# gemmi: "H"/"R"), and writes the "abc" code of the default setting of
# orthorhombic numbers with axis permutations explicitly (left out by other
# codes).
def test_table():
    def nc_choice( number, choice ):
        choice = { 'H' : 'h', 'R' : 'r' }.get( choice, choice )
        if not choice and number in ortho_abc:
            choice = 'abc'
        return choice

    def itvb_hall( h ):
        #ITVB's "Hall entry" column (and NCrystal) spell the lattice letter in
        #lower case, spglib and gemmi in upper case:
        h = ' '.join( h.split() )
        i = next( i for i, c in enumerate( h ) if c.isalpha() )
        return h[:i] + h[i].lower() + h[i+1:]

    def spglib_type( hn ):
        t = spglib.get_spacegroup_type( hn )
        get = ( lambda k : getattr( t, k ) ) if hasattr( t, 'number' ) else t.get
        return ( get('number'), nc_choice( get('number'), get('choice') ),
                 itvb_hall( get('hall_symbol') ) )

    def ops_set( rots, trans ):
        return frozenset( ( tuple( int(x) for x in np.asarray(r).flatten() ),
                            tuple( Fraction( round( float(x)*24 ) % 24, 24 )
                                   for x in t ) )
                          for r, t in zip( rots, trans ) )

    def spglib_ops( hn ):
        d = spglib.get_symmetry_from_database( hn )
        return ops_set( d['rotations'], d['translations'] )

    def gemmi_ops( hall ):
        #(gemmi wants the upper-case lattice letter)
        i = next( i for i, c in enumerate( hall ) if c.isalpha() )
        ops = gemmi.symops_from_hall( hall[:i] + hall[i].upper() + hall[i+1:] )
        return ops_set( [ np.array(o.rot)//gemmi.Op.DEN for o in ops ],
                        [ np.array(o.tran)/gemmi.Op.DEN for o in ops ] )

    gemmi_itb = list( gemmi.spacegroup_table_itb() )
    assert len( gemmi_itb ) == 530
    #The 46 orthorhombic numbers whose settings are only axis permutations
    #(no origin choices), so their default setting has the implicit code abc:
    ortho_abc = { g.number for g in gemmi_itb if 16 <= g.number <= 74
                  and g.qualifier and not g.ext.strip('\x00 ') }
    assert len( ortho_abc ) == 46, len( ortho_abc )
    first_hn = {}
    nsettings = {}
    for hn in range( 1, 531 ):
        number = lib.nctest_sgsym_number( hn )
        choice = lib.nctest_sgsym_choice( hn )
        hall = lib.nctest_sgsym_hallsymbol( hn )
        #spglib:
        assert ( number, choice, hall ) == spglib_type( hn ), ( hn, number, choice,
                                                                hall )
        assert gemmi_ops( hall ) == spglib_ops( hn ), hn
        #gemmi:
        g = gemmi_itb[hn-1]
        gchoice = nc_choice( g.number, ( g.ext + g.qualifier ).strip('\x00 ') )
        assert ( g.number, gchoice, itvb_hall( g.hall ) ) == ( number, choice,
                                                              hall ), hn
        #Hermann-Mauguin symbol (gemmi's extended symbol, with the axes codes
        #in lower case as in ITVB):
        hm = lib.nctest_sgsym_hmsymbol( hn )
        assert hm == g.xhm().replace( ':H', ':h' ).replace( ':R', ':r' ), ( hn, hm )
        assert hm.split(':')[0] == g.hm, ( hn, hm )
        #String round trips:
        s = lib.nctest_sgsym_tostring( hn )
        assert s == ( f'{number}:{choice}' if choice else str(number) )
        assert lib.nctest_sgsym_parse( s ) == hn
        assert lib.nctest_sgsym_parse_strict( s ) == hn
        nsettings[number] = nsettings.get( number, 0 ) + 1
        first_hn.setdefault( number, hn )
        assert lib.nctest_sgsym_parse( str(number) ) == first_hn[number]
    assert sorted( first_hn ) == list( range( 1, 231 ) )
    assert lib.nctest_sgsym_parse( '227:3' ) == 0
    #Strict parsing only accepts a bare number for the 140 single-setting
    #numbers (relaxed parsing takes the default setting), and the messages
    #name the settings:
    for number, hn in first_hn.items():
        assert lib.nctest_sgsym_parse_strict( str(number) ) == ( hn if nsettings[number] == 1 else 0 )
    assert sum( 1 for n in nsettings.values() if n == 1 ) == 140
    err = lib.nctest_sgsym_parse_error
    assert 'has several settings, so one of them must be specified: 227:1 (default), 227:2' in err('227')
    assert 'the settings are: 227:1 (default), 227:2' in err('227:3')
    assert 'change it to "166:h" if that is the intention' in err('166:H')
    assert 'has only one setting, which is written as just "225"' in err('225:1')
    assert err('227:2') == ''
    print('All 530 space group settings (numbers, choice codes, Hermann-Mauguin'
          ' and Hall symbols) consistent with spglib and gemmi')

# Verifies NCrystal's parsing and formatting of symmetry operations (like
# "-x,y+1/2,-z+1/2") against gemmi, using all operations of the 530 space
# group settings in gemmi's table.
def test_ops():
    def ncrystal_symop( s ):
        res = lib.nctest_sgsym_symop( s )
        assert res, f'NCrystal failed to parse "{s}"'
        parts = res.split()
        return parts[0], [ int(x) for x in parts[1:10] ], [ int(x) for x in parts[10:] ]

    def gemmi_rot_trans( op ):
        #gemmi.Op.DEN is 24, like NCrystal's translation unit:
        assert gemmi.Op.DEN == 24
        rot = [ v // gemmi.Op.DEN for row in op.rot for v in row ]
        trans = [ v % gemmi.Op.DEN for v in op.tran ]
        return rot, trans

    ops = {}
    for sg in gemmi.spacegroup_table_itb():
        for op in sg.operations():
            ops[ op.triplet() ] = op
    print(f'Distinct operations in the 530 settings of gemmi: {len(ops)}')
    nsame_str = 0
    for triplet, op in sorted( ops.items() ):
        s, rot, trans = ncrystal_symop( triplet )
        assert ( rot, trans ) == gemmi_rot_trans( op ), (triplet, s)
        #gemmi must parse our format to the same operation:
        assert gemmi_rot_trans( gemmi.Op( s ) ) == ( rot, trans ), (triplet, s)
        #And we must parse our own format to the same operation:
        assert ncrystal_symop( s ) == ( s, rot, trans )
        if s == triplet:
            nsame_str += 1
    print(f'All parsed identically by NCrystal and gemmi (string representations'
          f' identical for {nsame_str} of them)')

# Verifies the symmetry operations, centring vectors, coset representatives
# and centrosymmetry which NCrystal derives from the Hall symbols of the 530
# space group settings, against spglib and gemmi. The Hall symbol parser is
# additionally tested on gemmi's non-ITVB settings.
def test_symmetry():
    def triplet( rot, trans ):
        op = gemmi.Op()
        op.rot = [ [ int(v) * gemmi.Op.DEN for v in row ] for row in rot ]
        op.tran = [ round( float(x) * gemmi.Op.DEN ) % gemmi.Op.DEN
                    for x in trans ]
        return op.triplet()

    def spglib_ops( hn ):
        d = spglib.get_symmetry_from_database( hn )
        return { triplet( np.asarray(r), t )
                 for r, t in zip( d['rotations'], d['translations'] ) }

    def gemmi_ops( gops ):
        return { op.triplet() for op in gops }

    def gemmi_centring( gops ):
        return { triplet( np.eye( 3, dtype = int ), np.asarray(c) / gemmi.Op.DEN )
                 for c in gops.cen_ops }

    def ncrystal_symmetry( hn ):
        head, ops, reps, cent = lib.nctest_sgsym_symmetry( hn ).split('|')
        lattice, centro, order = head.split()
        return ( lattice, centro == '1', int(order), ops.split(';'),
                 reps.split(';'), cent.split(';') )

    gemmi_itb = list( gemmi.spacegroup_table_itb() )
    assert len( gemmi_itb ) == 530
    for hn in range( 1, 531 ):
        lattice, centro, order, ops, reps, cent = ncrystal_symmetry( hn )
        g = gemmi_itb[hn-1]
        gops = g.operations()
        assert order == len( ops ) == len( set( ops ) )
        assert set( ops ) == spglib_ops( hn ) == gemmi_ops( gops ), hn
        assert ops[0] == reps[0] == cent[0] == 'x,y,z'
        assert set( cent ) == gemmi_centring( gops ), hn
        assert len( reps ) == len( gops.sym_ops ), hn
        assert centro == g.is_centrosymmetric(), hn
        hall = lib.nctest_sgsym_hallsymbol( hn )
        assert lattice == hall.lstrip('-')[0], hn
        #Raw parser gives the same (canonically ordered) operations:
        assert lib.nctest_sgsym_rawhall( hall ).split(';') == ops, hn
    print('All 530 space group settings: operations, centring vectors,'
          ' representatives and centrosymmetry consistent with spglib and gemmi')

    #Hall symbol parser on gemmi's settings which are not among the 530 of ITVB
    #(origin shifts and non-conventional centred cells):
    itb_halls = { ' '.join( g.hall.split() ) for g in gemmi_itb }
    nextra = 0
    for g in gemmi.spacegroup_table():
        hall = ' '.join( g.hall.split() )
        if hall in itb_halls:
            continue
        #(the parser only accepts ITVB's spelling, i.e. a lower-case lattice
        #letter, where gemmi uses upper case):
        i = next( i for i, c in enumerate( hall ) if c.isalpha() )
        res = lib.nctest_sgsym_rawhall( hall[:i] + hall[i].lower() + hall[i+1:] )
        assert not res.startswith('ERROR'), (hall,res)
        assert set( res.split(';') ) == gemmi_ops( g.operations() ), hall
        nextra += 1
    print(f'Hall symbols of {nextra} additional (non-ITVB) gemmi settings'
          ' consistent with gemmi')

def main():
    test_table()
    test_ops()
    test_symmetry()

if __name__ == '__main__':
    main()
