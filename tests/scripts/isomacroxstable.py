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

# Test the "isomacroxstable" query, which provides piecewise cubic tables of
# total cross sections vs. wavelength.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query as ncquery
from NCrystalDev.exceptions import NCBadInput, NCCalcError
from NCTestUtils.common import ensure_error
import numpy as np

def exact_xs( cfg, wl ):
    #Macroscopic cross section (1/cm):
    return ( ( NC.createScatter(cfg).xsect(wl=wl)
               + NC.createAbsorption(cfg).xsect(wl=wl) )
             * NC.createInfo(cfg).numberdensity )

def bragg_edges( cfg ):
    res = []
    def collect( info ):
        if info.isMultiPhase():
            for _,ph in info.phases:
                collect(ph)
        elif info.hasHKLInfo():
            res.extend( 2.0*h.d for h in info.hklObjects() )
    collect( NC.createInfo(cfg) )
    return np.unique(res)

def table_eval( wl_t, xs_t, s_t, wl ):
    #Cubic Hermite interpolation, where a pair of identical wavelengths (a
    #Bragg edge) uses the second point (above the edge) for wl >= edge:
    i = np.clip( np.searchsorted(wl_t, wl, side='right') - 1, 0, len(wl_t)-2 )
    x0, x1 = wl_t[i], wl_t[i+1]
    h = np.where( x1 > x0, x1 - x0, 1.0 )
    u = np.where( x1 > x0, (wl-x0)/h, 1.0 )
    v = 1.0 - u
    return ( v*v*( (1.0+2.0*u)*xs_t[i] + u*h*s_t[i] )
             + u*u*( (3.0-2.0*u)*xs_t[i+1] - v*h*s_t[i+1] ) )

def decode( r ):
    #Wavelengths, macroscopic cross sections and slopes, from the table in
    #r['raw_result'] (see examples/ncrystal_example_filter.c for the format):
    raw = r['raw_result']
    n = int( raw[2] )
    assert raw[0] == raw[-1] == 4631638934963022.0 and raw[1] == 1.0
    assert len(raw) == 3 * n + 4 and n >= 2
    return raw[3:3+n], raw[3+n:3+2*n], raw[3+2*n:3+3*n]

def test( cfg, *opts ):
    r = ncquery( ['isomacroxstable',cfg] + list(opts), huge_arrays = True )
    wl, xs, sl = decode( r )
    tol, wlmax = r['tol'], r['wlmax']
    assert wl[0] == 0.0 and wl[-1] == wlmax
    assert np.all( np.diff(wl) >= 0.0 )
    assert np.all( xs >= 0.0 )
    #Identical wavelengths only come in pairs. Pairs with different values are
    #discontinuities (for the supported processes, these are the Bragg edges),
    #while pairs with identical values only change the slope:
    idup = np.where( np.diff(wl) == 0.0 )[0]
    assert np.all( np.diff(idup) > 1 )
    idup = idup[ xs[idup] != xs[idup+1] ]
    edges = bragg_edges( cfg )
    for e in wl[idup]:
        assert np.abs(edges-e).min() <= 1e-9 * e
    #The first point is the limit for wavelength -> 0 (which might be 0),
    #compared on the scale of the cross section at 1e-5 Aa:
    xs_lim, xs_scale = exact_xs( cfg, np.asarray([1e-14,1e-5]) )
    assert abs( xs[0] - xs_lim ) <= 1e-6 * max( xs_lim, xs_scale )
    #The other table points (except at edges, which have the values just
    #below or above the edge) are exact cross sections:
    if len(edges):
        j = np.clip( np.searchsorted( edges, wl ), 1, len(edges)-1 )
        dist = np.minimum( np.abs(wl-edges[j-1]), np.abs(wl-edges[j]) )
        iexact = np.where( dist > 1e-9 * wl )[0]
    else:
        iexact = np.arange(len(wl))
    iexact = iexact[ iexact > 0 ]
    xsx = exact_xs( cfg, wl[iexact] )
    assert np.allclose( xs[iexact], xsx, rtol = 1e-12, atol = 0.0 )
    #Interpolated values are within the tolerance (which is not required below
    #1e-5 Aa), at random wavelengths (log-uniform, and uniform at the shortest
    #wavelengths) and close to the Bragg edges:
    rng = np.random.default_rng( 12345 )
    wlmin = 1e-5
    wlt = np.concatenate( [ np.exp( rng.uniform( np.log(wlmin), np.log(wlmax),
                                                 20000 ) ),
                            rng.uniform( wlmin, 1e-3, 2000 ),
                            rng.uniform( wlmin, 1e-4, 2000 ) ] )
    e = edges[ (edges > 2*wlmin) & (edges < wlmax/1.001) ]
    wlt = np.concatenate( [wlt] + [ e*f for f in (1-1e-8,1+1e-8,1-1e-5,1+1e-5) ] )
    xst = exact_xs( cfg, wlt )
    xsi = table_eval( wl, xs, sl, wlt )
    assert np.all( xsi >= 0.0 )
    #Where the cross section is 0, the table must be exactly 0:
    assert np.all( xsi[xst == 0.0] == 0.0 )
    pos = xst > 0.0
    maxrelerr = ( ( np.abs( xsi[pos]-xst[pos] ) / xst[pos] ).max()
                  if pos.any() else 0.0 )
    assert maxrelerr <= tol, f'max relative error {maxrelerr} > tol={tol}'
    print(f'{cfg} {" ".join(opts)}'.strip())
    print(f'   ndiscontinuities={r["ndiscontinuities"]}'
          f' tol={tol:g} wlmax={wlmax:g}'
          f' maxrelerr<tol: OK')
    print(f'   processes: {" ".join(r["processes"])}')
    print(f'   first: ({wl[0]:.6g} Aa, {xs[0]:.6g}/cm)'
          f' last: ({wl[-1]:.6g} Aa, {xs[-1]:.6g}/cm)')

def main():
    test('stdlib::Al_sg225.ncmat')
    test('stdlib::Al_sg225.ncmat','tol=1e-2')
    test('stdlib::Al_sg225.ncmat','tol=1e-3')
    test('stdlib::Al_sg225.ncmat','wlmax=10')
    test('stdlib::Al_sg225.ncmat','wlmax=4.7')
    test('stdlib::Al_sg225.ncmat;density=0.5x')
    test('stdlib::Al_sg225.ncmat;inelas=0')
    test('stdlib::Be_sg194.ncmat;temp=80K')
    test('stdlib::Polyethylene_CH2.ncmat')
    #Next-gen kernels, so all processes have tables:
    test('stdlib::Al_sg225.ncmat;knllux=3')
    test('stdlib::Be_sg194.ncmat;temp=80K;knllux=3')
    test('stdlib::Polyethylene_CH2.ncmat;knllux=3')
    test('stdlib::Y2SiO5_sg15_YSO.ncmat')
    test('stdlib::B4C_sg166_BoronCarbide.ncmat')
    test('gasmix::air')
    test('phases<0.3*stdlib::Al_sg225.ncmat&0.7*stdlib::Cu_sg225.ncmat>')
    test('stdlib::void.ncmat')

    base = ['isomacroxstable','stdlib::Al_sg225.ncmat']
    oriented = ('stdlib::Al_sg225.ncmat;mos=0.3deg;dir1=@crys_hkl:0,0,1@lab:0,0,1'
                ';dir2=@crys_hkl:0,1,0@lab:0,1,0')
    with ensure_error(NCBadInput,'IsoMacroXSTable: only isotropic materials are'
                      ' supported'):
        ncquery(['isomacroxstable',oriented])
    #Processes which have not been validated for IsoMacroXSTable give a warning:
    test('stdlib::Polyethylene_CH2.ncmat;ucnmode=remove')
    #Discontinuities which are not known in advance are detected (here by
    #disabling the known ones):
    with ensure_error(NCCalcError,'IsoMacroXSTable: the cross section has an'
                      ' unexpected discontinuity (or an extremely sharp'
                      ' feature) at wavelength 1.65323 Aa'):
        ncquery(base+['discontinuities=0'])
    with ensure_error(NCBadInput,'IsoMacroXSTable: invalid tolerance 0 (must be'
                      ' in the range (0,1))'):
        ncquery(base+['tol=0'])
    with ensure_error(NCBadInput,'IsoMacroXSTable: invalid wlmax 1e-05 (must be'
                      ' larger than 0.0001 Aa)'):
        ncquery(base+['wlmax=1e-5'])
    with ensure_error(NCBadInput,'Invalid isomacroxstable query:'
                      ' ["isomacroxstable","stdlib::Al_sg225.ncmat","wlmin=1"]'
                      ' (unknown option wlmin)'):
        ncquery(base+['wlmin=1'])
    with ensure_error(NCBadInput,'Invalid isomacroxstable query:'
                      ' ["isomacroxstable","stdlib::Al_sg225.ncmat","tol"]'
                      ' (options must have the form KEY=VALUE)'):
        ncquery(base+['tol'])

if __name__ == '__main__':
    main()
