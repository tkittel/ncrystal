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

#SABAnalyser on sparse beta grids: Teff estimates must not depend on how
#densely the (linearly interpolated) kernel is sampled in beta, and msd must
#be refused for kernels too coarse to resolve it, but not for dense ones.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query as q
import numpy as np

def analyse( cfgstr, lbl = '', diag = False ):
    return q( ['sab','analyse',cfgstr,lbl] + ( ['diag'] if diag else [] ) )[
        'sabanalyse']

def scatknl_section( el, fraction, temp, alpha, beta, sab ):
    fmt = lambda v : ' '.join( repr(float(x)) for x in v )
    return ( f'@DYNINFO\n  element {el}\n  fraction {fraction}\n'
             f'  type scatknl\n  temperature {temp}\n'
             f'  alphagrid {fmt(alpha)}\n  betagrid {fmt(beta)}\n'
             f'  sab {fmt(sab)}\n' )

def linearly_thickened( beta, sab, nextra ):
    #nextra extra beta rows, linear in S, between each pair of rows:
    nb = len(beta)
    s = np.asarray(sab).reshape( nb, -1 )
    b_out, s_out = [], []
    for i in range( nb - 1 ):
        for k in range( nextra + 1 ):
            f = k / ( nextra + 1 )
            b_out.append( beta[i] * ( 1 - f ) + beta[i+1] * f )
            s_out.append( s[i] * ( 1 - f ) + s[i+1] * f )
    b_out.append( beta[-1] )
    s_out.append( s[-1] )
    return np.array( b_out ), np.concatenate( s_out )

def grid_invariance():
    #Per-row Teff estimates of rows accepted for both the sparse test kernel
    #and a linearly thickened copy of it must agree:
    f = 'C_from_benzene_solid_20K_sabsmall.ncmat'
    info = NC.createInfo( f )
    k = info.dyninfos[0].loadKernel()
    bt, st = linearly_thickened( k['beta'], k['sab'], 3 )
    NC.registerInMemoryFileData( 'thick_c.ncmat', 'NCMAT v7\n'
                                 '@DENSITY\n  1 g_per_cm3\n'
                                 + scatknl_section( 'C', 1, 20, k['alpha'],
                                                    bt, st ) )
    d0 = analyse( f, diag = True )['diagnostics']
    d1 = analyse( 'virtual::thick_c.ncmat', diag = True )['diagnostics']
    assert d0['alpha'] == d1['alpha']
    both = [ (a,b) for a,b,sa,sb in zip( d0['teff_row'], d1['teff_row'],
                                         d0['teff_row_status'],
                                         d1['teff_row_status'] )
             if sa == 0 and sb == 0 ]
    assert len( both ) >= 10
    maxdiff = max( abs( b / a - 1 ) for a, b in both )

    assert maxdiff < 1e-9
    print( f'==> per-row Teff of {f} unchanged by thickening the beta grid'
           f' ({len(both)} rows): OK' )

def msd_resolution():
    #The Li2O test kernel keeps only every 8th beta point, too coarse for msd,
    #while kernels expanded from the same VDOS are dense enough:
    ref = 'stdlib::Li2O_sg225_LithiumOxide.ncmat;temp=10K'
    info = NC.createInfo( ref )
    secs = []
    for di in info.dyninfos:
        el = di.atomData.displayLabel()
        k = di.loadKernel( vdoslux = 1 )
        secs.append( scatknl_section( el, '2/3' if el == 'Li' else '1/3', 10,
                                      k['alpha'], k['beta'], k['sab'] ) )
    NC.registerInMemoryFileData( 'dense_li2o.ncmat',
                                 'NCMAT v7\n@DENSITY\n  2.013 g_per_cm3\n'
                                 + ''.join( secs ) )
    for el in ( 'Li', 'O' ):
        coarse = analyse( 'Li2O_sg225_LithiumOxide_sabsmall_temp10K.ncmat', el )
        assert coarse['auto']['msd'] is None
        dense = analyse( 'virtual::dense_li2o.ncmat', el )
        truth = analyse( ref + ';vdoslux=1', el )['msd_vdos']
        rel = dense['auto']['msd'] / truth - 1
        assert abs( rel ) < 0.005
        print( f'==> Li2O {el} msd: refused for the coarse test kernel,'
               f' within 0.5% of the VDOS value for the dense one: OK' )

grid_invariance()
msd_resolution()
