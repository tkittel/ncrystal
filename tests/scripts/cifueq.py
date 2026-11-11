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

# Ueq of CIF anisotropic displacements (_sgsymimpl._cif_uij_ueq), against an
# independent calculation via the Cartesian tensor.

import NCTestUtils.enable_fpe # noqa F401
import numpy as np
from NCrystalDev._sgsymimpl import _cif_uij_ueq

def lattice( a, b, c, alpha, beta, gamma ):
    #Columns are the lattice vectors (a along x, b in the xy plane):
    ca, cb, cg = ( np.cos( np.radians( x ) ) for x in ( alpha, beta, gamma ) )
    sg = np.sin( np.radians( gamma ) )
    cx = c * cb
    cy = c * ( ca - cb * cg ) / sg
    cz = np.sqrt( c*c - cx*cx - cy*cy )
    return np.array( [ [ a, b*cg, cx ], [ 0.0, b*sg, cy ], [ 0.0, 0.0, cz ] ] )

def ref_ueq( u, cell ):
    lat = lattice( *cell )
    rlen = np.linalg.norm( np.linalg.inv( lat ), axis = 1 )#|a*_i|
    d = np.diag( rlen )
    ucart = lat @ d @ u @ d @ lat.T
    return np.trace( ucart ) / 3.0

def cif_from_cart( ucart, cell ):
    lat = lattice( *cell )
    linv = np.linalg.inv( lat )
    dinv = np.diag( 1.0 / np.linalg.norm( linv, axis = 1 ) )
    return dinv @ linv @ ucart @ linv.T @ dinv

def as_dict( u ):
    return { 'u11': u[0,0], 'u22': u[1,1], 'u33': u[2,2],
             'u12': u[0,1], 'u13': u[0,2], 'u23': u[1,2] }

def cellsg( cell ):
    return dict( zip( ('a','b','c','alpha','beta','gamma'), cell ) )

def main():
    cells = [ ( 4.05, 4.05, 4.05, 90, 90, 90 ),
              ( 5.1, 6.3, 7.2, 90, 90, 90 ),
              ( 3.21, 3.21, 5.21, 90, 90, 120 ),
              ( 5.1, 6.3, 7.2, 90, 105, 90 ),
              ( 5.1, 6.3, 7.2, 81, 97, 103 ) ]
    rng = np.random.default_rng( 1234 )
    for cell in cells:
        #Isotropic displacements, for which Ueq is the same in any cell:
        uiso = cif_from_cart( 0.01 * np.eye(3), cell )
        ueq_iso = _cif_uij_ueq( as_dict( uiso ), cellsg( cell ) )
        #Random positive definite Cartesian tensors:
        maxdiff = 0.0
        for _ in range( 20 ):
            m = rng.normal( size = (3,3) ) * 0.05
            u = cif_from_cart( m @ m.T + 0.001 * np.eye(3), cell )
            ueq = _cif_uij_ueq( as_dict( u ), cellsg( cell ) )
            ref = ref_ueq( u, cell )
            maxdiff = max( maxdiff, abs( ueq - ref ) / ref )
        assert maxdiff < 1e-13, maxdiff
        assert abs( ueq_iso - 0.01 ) < 1e-15
        #A fixed anisotropic tensor, where the mean of the diagonal is only
        #correct for orthogonal axes:
        u = { 'u11': 0.01, 'u22': 0.012, 'u33': 0.02,
              'u12': 0.0, 'u13': 0.004, 'u23': 0.0 }
        ueq = _cif_uij_ueq( u, cellsg( cell ) )
        diagmean = ( u['u11'] + u['u22'] + u['u33'] ) / 3.0
        print( f'cell {cell}: Ueq={ueq:.6g} (mean of the diagonal'
               f' {diagmean:.6g}), isotropic and random tensors OK' )

if __name__ == '__main__':
    main()
