
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

Utilities for tabulating macroscopic cross section of isotropic materials.

"""

__all__ = [ 'IsoMacroXSTable' ]

class IsoMacroXSTable:

    """Macroscopic total cross section (scattering plus absorption) of an
    isotropic material, based on an automated tabulation of the normal physics
    processes provided by NCrystal. This could for instance be useful for very
    fast calculation of transmission probabilities for a neutron beam passing
    through a material unscattered and unabsorbed (e.g. for a filter or a
    window). A neutron travelling a distance L (in cm) through the material is
    transmitted with probability exp(-xsect*L). This completely ignores multiple
    scattering effects.

    The cross sections are evaluated from a piecewise cubic table vs. the
    neutron wavelength, reproducing the cross section of the material to a
    relative precision of 0.01%. An exception is raised if the cross section can
    not be tabulated reliably.

    """

    def __init__( self, cfgstr ):
        from ._chooks import _get_raw_cfcts
        self.__cfgstr = str(cfgstr)
        block = _get_raw_cfcts()['isomacroxstable']( self.__cfgstr )
        block.setflags( write = False )
        n = int( block[2] )
        self.__block = block
        self.__wl = block[3:3+n]
        self.__macroxs = block[3+n:3+2*n]
        self.__slope = block[3+2*n:3+3*n]

    @property
    def cfgstr( self ):
        """The cfg-string of the material."""
        return self.__cfgstr

    @property
    def block( self ):
        """The table as a single read-only numpy array, exactly as returned by
        the C API function ncrystal_isomacroxstable (see the format description
        in the ncrystal_example_filter.c file in the NCrystal repository)."""
        return self.__block

    @property
    def table( self ):
        """The table, as three read-only numpy arrays (views into the block):
        the wavelengths (in Aa), the macroscopic cross sections (in 1/cm) and
        their slopes (in 1/(cm*Aa)). Between neighbouring points, the cross
        section is the cubic Hermite interpolation of these. The first point is
        the limit for wavelength -> 0. Discontinuities (e.g. Bragg edges) are
        represented by two points with the same wavelength, where the first has
        the value just below the discontinuity, and the second the value just
        above it."""
        return self.__wl, self.__macroxs, self.__slope

    def xsect( self, ekin = None, wl = None ):
        """Macroscopic cross section (in 1/cm) at the given neutron kinetic
        energy (in eV) or wavelength (in Aa). These can be numbers (returning a
        float), or arrays (returning a numpy array)."""
        from ._numpy import _ensure_numpy, _np
        from .constants import ekin2wl
        from .exceptions import NCBadInput
        if ( ekin is None ) == ( wl is None ):
            raise NCBadInput('Please provide exactly one of the "ekin" or'
                             ' "wl" parameters.')
        _ensure_numpy()
        w = _np.asarray( ekin2wl( ekin ) if wl is None else wl, dtype = float )
        x, y, s = self.__wl, self.__macroxs, self.__slope
        #Cubic Hermite interpolation in the interval with x[i] <= w < x[i+1]
        #(at a discontinuity, this is the interval above it):
        wi = _np.clip( w, x[0], x[-1] )
        i = _np.clip( _np.searchsorted( x, wi, side = 'right' ) - 1,
                      0, len(x) - 2 )
        x0, x1 = x[i], x[i+1]
        h = _np.where( x1 > x0, x1 - x0, 1.0 )
        u = ( wi - x0 ) / h
        v = 1.0 - u
        res = ( v * v * ( ( 1.0 + 2.0 * u ) * y[i] + u * h * s[i] )
                + u * u * ( ( 3.0 - 2.0 * u ) * y[i+1] - v * h * s[i+1] ) )
        #Linear extrapolation of the value and slope of the last point beyond
        #the table (clamped at 0):
        if s[-1] != 0.0:
            extrap = _np.maximum( y[-1] + ( w - x[-1] ) * s[-1], 0.0 )
        else:
            extrap = _np.full_like( w, y[-1] )
        res = _np.where( w >= x[-1], extrap, res )
        res = _np.where( w <= x[0], y[0], res )
        return float( res ) if res.ndim == 0 else res
