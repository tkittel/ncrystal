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

# Test of the grouped treatment of VDOSDebye entries (VDOSDebyeGroupScatter,
# at vdoslux 0-5) against separate expansions for all entries (vdoslux=103),
# for five oxygen labels: cross sections agree to ~1e-3 (limited by the
# precision of the phonon expansions), and sampled distributions agree within
# the clustering tolerances. With hydrogen added (where the expansions are
# ~2e-3 low even at vdoslux=6), checks the dumps and JSON (also of SABScatter
# at vdoslux=103), the exact free-atom limit at high energies, and which
# vdoslux values use the grouping.

import json

import NCrystalDev as NC


def material( with_h ):
    #Five oxygen labels with different Debye temperatures, and hydrogen:
    tds = ( 250.0, 320.0, 400.0, 520.0, 700.0 )
    lines = [ 'NCMAT v7', '@STATEOFMATTER', '  solid',
              '@DENSITY', '  1.0 g_per_cm3', '@ATOMDB' ]
    lines += [ f'  X{i+1} is O' for i in range( len( tds ) ) ]
    for i, td in enumerate( tds ):
        lines += [ '@DYNINFO', f'  element X{i+1}',
                   '  fraction %g' % ( 0.16 if with_h else 0.2 ),
                   '  type vdosdebye', f'  debye_temp {td}' ]
    if with_h:
        lines += [ '@DYNINFO', '  element H', '  fraction 0.2',
                   '  type vdosdebye', '  debye_temp 600' ]
    fn = 'vdosdebyegroup%s.ncmat' % ( '_h' if with_h else '' )
    NC.registerInMemoryFileData( fn, '\n'.join( lines ) )
    return fn + ';comp=inelas'

def load( vdoslux, with_h = False ):
    return NC.createScatter( material( with_h ) + f';vdoslux={vdoslux}' )

def moments( sc, e, n = 100000 ):
    ef, mu = sc.sampleScatterIsotropic( e, repeat = n )
    return ef.mean() / e, ef.std() / e, mu.mean(), mu.std()

def main():
    sa, sn = load( 3 ), load( 103 )
    energies = [ 1e-5, 1e-4, 1e-3, 0.0253, 0.1, 0.5, 2.0, 10.0 ]
    xa = [ sa.xsect( wl = NC.ekin2wl( e ) ) for e in energies ]
    xn = [ sn.xsect( wl = NC.ekin2wl( e ) ) for e in energies ]
    for e, a, n in zip( energies, xa, xn ):
        assert abs( a / n - 1.0 ) < 2e-3, ( e, a, n )
    assert xa != xn, 'grouped treatment not in effect'
    print( 'Cross sections agree within 2e-3' )
    #Sampled distributions: means and spreads of E'/E and mu:
    for e in ( 0.003, 0.025, 0.2 ):
        ma, mn = moments( sa, e ), moments( sn, e )
        for a, n in zip( ma, mn ):
            assert abs( a - n ) < 0.02 * max( 1.0, abs( n ) ), ( e, ma, mn )
    print( 'Sampled distributions agree' )
    #High energy limit: sigma_free per atom:
    def sigfree( el ):
        ad = NC.atomDB( el )
        return ad.scatteringXS() * ( ad.averageMassAMU() / ( ad.averageMassAMU() + 1.00866491606 ) )**2
    sh = load( 3, with_h = True )
    print( sh.dump_str(), end = '' )
    #JSON of the VDOSDebyeGroupScatter processes (without uids, and floats
    #with 6 significant digits):
    def rnd( o ):
        if isinstance( o, dict ):
            return { k : rnd( v ) for k, v in o.items() if not k.endswith( 'uid' ) }
        if isinstance( o, list ):
            return [ rnd( v ) for v in o ]
        return float( f'{o:.6g}' ) if isinstance( o, float ) else o
    for scale, p in sh.getSummary( short = False )['specific']['components']:
        print( scale, json.dumps( rnd( p ), indent = 1 ) )
    #The same without grouping (SABScatter for each entry, H only for brevity):
    sh103 = load( 103, with_h = True )
    print( sh103.dump_str(), end = '' )
    for scale, p in sh103.getSummary( short = False )['specific']['components']:
        if p['specific']['summarystr'].endswith( 'sigma_free=20.4902barn' ):
            print( scale, json.dumps( rnd( p ), indent = 1 ) )
    ref = 0.8 * sigfree( 'O' ) + 0.2 * sigfree( 'H' )
    xhi = sh.xsect( wl = NC.ekin2wl( 1e3 ) )
    assert abs( xhi / ref - 1.0 ) < 1e-3, ( xhi, ref )
    print( 'High energy limit is the free-atom cross section' )
    #Grouping at vdoslux 0-5 only (a single oxygen entry):
    NC.registerInMemoryFileData( 'vdosdebyegroup_one.ncmat',
                                 'NCMAT v7\n@STATEOFMATTER\n  solid\n'
                                 '@DENSITY\n  1.0 g_per_cm3\n@DYNINFO\n'
                                 '  element O\n  fraction 1\n'
                                 '  type vdosdebye\n  debye_temp 400\n' )
    for lux in ( 0, 5, 6, 100, 105, 1003 ):
        sc = NC.createScatter( f'vdosdebyegroup_one.ncmat;vdoslux={lux};comp=inelas' )
        print( f'vdoslux={lux}:', sc.getSummary( short = False )['name'] )

def test_lowering():
    #Ten oxygen entries with similar Debye temperatures, and two light
    #outliers: neff>3, so all clusters are sampled with vdoslux and knllux
    #two levels lower (never below 1), visible as smaller kernels:
    tds = [ 380.0 + 4.0 * i for i in range( 10 ) ] + [ 150.0, 1000.0 ]
    fracs = [ '9/100' ] * 10 + [ '1/20' ] * 2
    lines = [ 'NCMAT v7', '@STATEOFMATTER', '  solid',
              '@DENSITY', '  1.0 g_per_cm3', '@ATOMDB' ]
    lines += [ f'  X{i+1} is O' for i in range( len( tds ) ) ]
    for i, ( td, f ) in enumerate( zip( tds, fracs ) ):
        lines += [ '@DYNINFO', f'  element X{i+1}', f'  fraction {f}',
                   '  type vdosdebye', f'  debye_temp {td}' ]
    NC.registerInMemoryFileData( 'vdosdebyegroup_low.ncmat', '\n'.join( lines ) )
    for lux in ( 1, 2, 3, 5 ):
        sc = NC.createScatter( f'vdosdebyegroup_low.ncmat;comp=inelas;vdoslux={lux}' )
        print( f'vdoslux={lux}:' )
        for c in sc.getSummary( short = False )['specific']['clusters']:
            print( '   ', c )

if __name__ == '__main__':
    main()
    test_lowering()
