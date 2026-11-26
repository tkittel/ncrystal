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

"""Generates the reference values in tests/data/mosaicrefs/ for tests of the
mosaic crystal numerics, with high precision mpmath calculations which are
independent of NCrystal. The model is a Gaussian mosaicity distribution on the
sphere, truncated at ntrunc*sigma and normalised. Usage:

  python3 gen_mosaicrefs.py <outdir>

Needs mpmath and takes some minutes (all cores are used)."""

import multiprocessing
import pathlib
import random
import sys

import mpmath as mp

mp.mp.dps = 25

def fwhm2sigma( fwhm ):
    return mp.mpf(fwhm) / ( 2 * mp.sqrt( 2 * mp.log(2) ) )

class Mosaic:
    """Truncated Gaussian density on the unit sphere."""

    def __init__( self, fwhm, ntrunc ):
        self.sigma = fwhm2sigma( fwhm )
        self.ta = ntrunc * self.sigma
        self.cta = mp.cos( self.ta )
        self.k = -1 / ( 2 * self.sigma**2 )
        integral = 2 * mp.pi * mp.quad( lambda x : mp.sin(x)
                                        * mp.exp( self.k * x * x ),
                                        [ 0, self.ta ] )
        self.norm = 1 / integral

    def density( self, x ):
        return self.norm * mp.exp( self.k * x * x ) if x < self.ta else 0

    def circle( self, cg, ca ):
        """Density integrated along a circle on the sphere, with opening angle
        alpha (cos(alpha)=ca) and at an angle gamma (cos(gamma)=cg) from the
        center of the density (cf. GaussOnSphere::circleIntegral)."""
        sg, sa = mp.sqrt( 1 - cg * cg ), mp.sqrt( 1 - ca * ca )
        sasg, cacg = sa * sg, ca * cg
        if cacg + sasg <= self.cta:
            return mp.mpf(0)
        if sasg == 0:
            return 2 * mp.pi * sa * self.density( mp.acos(ca) )
        cos_tmax = ( self.cta - cacg ) / sasg
        tmax = mp.pi if cos_tmax <= -1 else mp.acos( cos_tmax )
        def f( t ):
            cb = sasg * mp.cos(t) + cacg
            return self.density( mp.acos( mp.mpf(1) if cb > 1 else cb ) )
        return 2 * sa * mp.quad( f, [ 0, tmax ] )

def lc_xs_integral( mosaic, dsp, alpha, fsq, wl, c3, sign ):
    """Cross section of a planeset in a layered crystal, integrated over the
    rotation phi in [0,pi] of the crystallites around the lcaxis (cf.
    LCStdFrame::calcXSIntegral)."""
    sinth = mp.mpf(wl) / ( 2 * mp.mpf(dsp) )
    costh = mp.sqrt( 1 - sinth * sinth )
    q = mp.mpf(wl)**3 / ( sinth * costh ) * mp.mpf(fsq) * mp.mpf('0.5')
    s3 = mp.sqrt( 1 - mp.mpf(c3)**2 )
    A = mp.sin(alpha) * s3 * sign
    B = mp.cos(alpha) * mp.mpf(c3) * sign
    #Break points at the phi values where the circle meets the truncation:
    oa = mp.acos( sinth )#opening angle
    pts = [ mp.mpf(0), mp.pi ]
    for g in ( oa - mosaic.ta, oa + mosaic.ta, mosaic.ta - oa ):
        if A != 0 and 0 <= g <= mp.pi:
            x = ( mp.cos(g) - B ) / A
            if -1 < x < 1:
                pts.append( mp.acos(x) )
    pts = sorted( pts )
    def f( phi ):
        cg = A * mp.cos(phi) + B
        return mosaic.circle( max( -1, min( 1, cg ) ), sinth )
    return q * mp.quad( f, pts )

def gen_lcinteg_case( args ):
    ( mosdeg, ntrunc, _prec, dsp, alpha, fsq, wl, c3, sign ) = args
    mosaic = Mosaic( mp.mpf(mosdeg) * mp.pi / 180, ntrunc )
    ref = lc_xs_integral( mosaic, dsp, alpha, fsq, wl, c3, sign )
    return ( args, ref )

def lcinteg_cases():
    """Planesets and neutrons with nonzero integrals, for a range of
    mosaicities, truncations and precision parameters (passed on to the
    tests). Both near and far from back-scattering, and with the window of
    the mosaicity distribution near or far from the edges of the range of
    angles between the neutron and the normals."""
    rng = random.Random( 12345 )
    res = []
    for mosdeg, ntrunc, prec in ( ( 0.01, 3.7, 1e-3 ), ( 0.5, 3.7, 1e-3 ),
                                  ( 0.5, 5.6, 1e-5 ), ( 5.0, 3.7, 1e-3 ),
                                  ( 5.0, 5.6, 1e-5 ), ( 20.0, 3.0, 1e-2 ) ):
        ta = float( fwhm2sigma( mosdeg * mp.pi / 180 ) ) * ntrunc
        n = 0
        while n < 20:
            dsp = rng.uniform( 0.5, 3.5 )
            alpha = rng.choice( ( mp.pi / 2, rng.uniform( 0.02, 1.0 )
                                  * mp.pi / 2 ) )
            wl = 2 * dsp * rng.choice( ( rng.uniform( 0.05, 0.95 ),
                                         rng.uniform( 0.95, 1.0 ) ) )
            c3 = rng.choice( ( rng.uniform( 0.0, 1.0 ),
                               rng.uniform( 0.0, 0.001 ) ) )
            sign = rng.choice( ( -1, 1 ) )
            #Skip vanishing integrals:
            oa = float( mp.acos( wl / ( 2 * dsp ) ) )
            s3 = ( 1 - c3 * c3 )**0.5
            A = abs( float( mp.sin(alpha) ) * s3 )
            B = float( mp.cos(alpha) ) * c3 * sign
            gmin = float( mp.acos( min( 1.0, B + A ) ) )
            gmax = float( mp.acos( max( -1.0, B - A ) ) )
            if gmax <= oa - ta or gmin >= oa + ta:
                continue
            res.append( ( mosdeg, ntrunc, prec, dsp, float(alpha),
                          rng.uniform( 0.1, 1.1 ), wl, c3, sign ) )
            n += 1
    return res

def default_ntrunc( prec ):
    """As GaussOnSphere::estimateNTruncFromPrec."""
    return min( 8.0, max( 3.0, 1.1 * float( mp.sqrt( -2 * mp.log(prec) ) ) ) )

def gen_circle_case( args ):
    ( mosdeg, ntrunc, _prec, alpha, gamma ) = args
    mosaic = Mosaic( mp.mpf(mosdeg) * mp.pi / 180, ntrunc )
    ref = mosaic.circle( mp.cos( mp.mpf(gamma) ), mp.cos( mp.mpf(alpha) ) )
    return ( args, ref )

def circle_cases():
    """Circles for mosaicities from sub-arcsec to 40 degrees and precision
    parameters from 1e-2 to 1e-7, so all methods in GaussOnSphere are used
    (closed form, all table sizes, and numerical integration), with angles
    alpha and gamma-alpha on a grid including the edges of the truncation."""
    rng = random.Random( 54321 )
    res = []
    for mosdeg, prec, ntrunc in ( ( 0.0001, 1e-3, None ), ( 0.0001, 1e-7, None ),
                                  ( 0.001, 1e-5, None ), ( 0.01, 1e-3, None ),
                                  ( 0.1, 1e-7, None ), ( 1.0, 1e-3, None ),
                                  ( 1.0, 1e-7, None ), ( 1.0, 1e-5, 40.0 ),
                                  ( 5.0, 1e-2, None ), ( 5.0, 1e-4, None ),
                                  ( 10.0, 1e-4, None ), ( 20.0, 1e-4, None ),
                                  ( 20.0, 1e-7, None ), ( 40.0, 1e-3, None ),
                                  ( 40.0, 1e-4, None ) ):
        if ntrunc is None:
            ntrunc = default_ntrunc( prec )
        ta = float( fwhm2sigma( mosdeg * mp.pi / 180 ) ) * ntrunc
        alphas = [ a for a in ( 0.3*ta, 0.95*ta, 1.1*ta, 1.25*ta, 2*ta, 0.2,
                                0.7, 1.2, 1.5 ) if a < float( mp.pi / 2 ) ]
        deltas = [ f * ta for f in ( -0.999, -0.8, -0.5, -0.1, 0.0, 0.05, 0.3,
                                     0.7, 0.95, 0.999 ) ]
        combos = [ ( a, a + d ) for a in alphas for d in deltas
                   if 0.0 < a + d < float( mp.pi ) and a + d > 1e-3 * ta ]
        rng.shuffle( combos )
        for a, g in sorted( combos[:16] ):
            res.append( ( mosdeg, ntrunc, prec, a, g ) )
    return res

def fmt( x ):
    return repr( float(x) )

def main():
    outdir = pathlib.Path( sys.argv[1] )
    outdir.mkdir( parents = True, exist_ok = True )
    def write( fn, columns, results ):
        lines = [ '# Generated by devel/misc/gen_mosaicrefs.py. Columns:'
                  + columns ]
        for args, ref in results:
            lines.append( ' '.join( fmt(e) for e in args ) + ' ' + fmt(ref) )
        ( outdir / fn ).write_text( '\n'.join( lines ) + '\n' )
    with multiprocessing.Pool() as pool:
        write( 'circleint.txt',
               ' mosaicity FWHM (deg), ntrunc, prec, alpha, gamma, integral',
               pool.map( gen_circle_case, circle_cases() ) )
        write( 'lcinteg.txt',
               ' mosaicity FWHM (deg), ntrunc, prec, dspacing, alpha, fsq,'
               ' wavelength, c3, sign, integral',
               pool.map( gen_lcinteg_case, lcinteg_cases() ) )

if __name__ == '__main__':
    main()
