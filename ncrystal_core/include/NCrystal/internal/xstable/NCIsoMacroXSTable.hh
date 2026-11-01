#ifndef NCrystal_IsoMacroXSTable_hh
#define NCrystal_IsoMacroXSTable_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/internal/utils/NCSpan.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace XSTable {

    ////////////////////////////////////////////////////////////////////////////
    //
    // Piecewise cubic tables of the macroscopic total cross section
    // (scattering + absorption, in 1/cm) of an isotropic material, as a
    // function of the neutron wavelength. This allows generic production of
    // "just the cross section curve as a table", decoupled from our object
    // instances. One intended use-case is for simple attenuation-only
    // GPU-compatible filter components in McStas.
    //
    // Each table point has a wavelength, a cross section and its slope, and
    // between neighbouring points the cross section is the cubic Hermite
    // interpolation of these. To support discontinuities like Bragg edges, the
    // tables can contain two entries at the same point. It is therefore
    // important that consumers evaluating the code use an algorithm supporting
    // that, and the implementation of the evalIsoMacroXSTable function below
    // serve as a reference implementation.
    //
    // It is also very much the intention that this should be a non-fiddly
    // turn-key feature, so the IsoMacroXSTableParams options are only intended
    // for expert investigations, and should normally be left at default values.
    //
    ////////////////////////////////////////////////////////////////////////////

    // Technical details about the tables:
    //
    // The table covers wavelengths from 0 to wlmax. Beyond wlmax, the value
    // and slope of the last point can be extrapolated linearly (clamped at
    // 0). The value at 0 is the limit of the cross section for wavelength ->
    // 0.
    //
    // The table points are a subset of points where the cross section is
    // evaluated, and the slopes at the table points are those of the cross
    // section, limited where needed such that the interpolation can not
    // become negative (J. W. Schmidt and W. Hess 1988, DOI:
    // 10.1007/BF01934097). The table reproduces the cross section within the
    // tolerance at all wavelengths (relative to the cross section, except
    // below denseWavelengthMin(), see below).
    //
    // The tolerance is relative to the cross section, except at wavelengths
    // below denseWavelengthMin(), where the table is simply the interpolation
    // between the limit at wavelength 0 and the first sample point. This is
    // needed for cross sections vanishing faster than linearly at wavelength 0
    // (e.g. Bragg diffraction and incoherent elastic scattering, going like
    // wl^2, without absorption): the relative error of the interpolation in
    // the first interval is then the same, no matter how close to 0 the first
    // sample point is, so no table could meet the tolerance there.
    //
    // Known discontinuities (the breakpts of the isoXSDescr() of the
    // processes, e.g. Bragg edges, or the boundaries of their energy domains)
    // are represented by two table points with the same wavelength: the first
    // with the cross section just below the discontinuity, the second with the
    // cross section just above it. A pair can also have the same cross section
    // and different slopes (e.g. at a kink).
    //
    // Features narrower than the dense sampling, which are not known in
    // advance, can not be found reliably by sampling. Therefore a warning is
    // emitted for physics processes not known to be free of such features
    // (e.g. a future process with absorption resonances). Exceptions are thrown
    // if the cross section is not finite or negative, if it does not converge
    // for wavelength -> 0, or if an unexpected discontinuity is found.

    struct IsoMacroXSTableParams {
      //NB: The defaults here are what will be used by the C API
      double wlmax = 500.0;//Angstrom
      double tol = 1e-3;//Max relative error of the cross section (see above)
      unsigned ndense = 0;//Log-spaced dense sample points (plus
                          //discontinuities). 0 means auto (see autoNDense)
      bool discontinuities = true;//Include points at the known discontinuities
    };

    //Information about the creation of a table, for diagnostics:
    struct IsoMacroXSTableDiagnostics {
      std::size_t nDiscontinuities = 0;//Known discontinuities (clusters) in
                                       //range, excluding kinks
      std::size_t nDense = 0;//Number of dense sample points (after refinement)
      std::size_t nRefine = 0;//Number of refinement iterations
      std::vector<std::string> processes;//Names of the (leaf) physics processes
    };

    class IsoMacroXSTable;

    IsoMacroXSTable
    createIsoMacroXSTable( const MatCfg&, const IsoMacroXSTableParams&,
                           IsoMacroXSTableDiagnostics* = nullptr );

    //The same, for a general function of the wavelength, whose values are
    //taken to be macroscopic cross sections (1/cm), with its derivative
    //d(xs)/d(wl) and the discontinuities given explicitly (processes in the
    //diagnostics will be empty). The functions are only evaluated at
    //wavelengths away from the discontinuities, where the derivative is
    //well-defined:
    IsoMacroXSTable
    createIsoMacroXSTable( const std::function<double(double)>& xs_of_wl,
                           const std::function<double(double)>& dxs_of_wl,
                           VectD discontinuities,
                           const IsoMacroXSTableParams&,
                           IsoMacroXSTableDiagnostics* = nullptr );

    //The table as a single block of 3*n+4 doubles, which is also what the C
    //API function ncrystal_isomacroxstable returns. The format is documented
    //in examples/ncrystal_example_filter.c:
    class IsoMacroXSTable final : private MoveOnly {
    public:
      static constexpr double magic() { return 4631638934963022.0; }
      static constexpr unsigned formatVersion() { return 1; }

      //Number of points (0 after release()):
      std::size_t n() const;

      Span<const double> raw() const;//all 3*n+4 values
      Span<const double> wavelengths() const;//Angstrom
      Span<const double> macroXS() const;//1/cm
      Span<const double> slopes() const;//d(macroxs)/d(wl), 1/(cm*Angstrom)

      //Transfer the block (allocated with new double[]) to the caller,
      //leaving this object empty:
      std::unique_ptr<double[]> release();

    private:
      explicit IsoMacroXSTable( std::unique_ptr<double[]> );
      friend IsoMacroXSTable
      createIsoMacroXSTable( const std::function<double(double)>&,
                             const std::function<double(double)>&, VectD,
                             const IsoMacroXSTableParams&,
                             IsoMacroXSTableDiagnostics* );
      std::unique_ptr<double[]> m_data;
    };

    //The lowest wavelength of the log-spaced dense sample (Angstrom). Below
    //it (i.e. above neutron energies of about 0.8 GeV), the table is the
    //interpolation between the limit at wavelength 0 and the first sample
    //point, and the tolerance is not required:
    constexpr double denseWavelengthMin() { return 1e-5; }

    //The number of log-spaced dense sample points used for a given tolerance
    //and wlmax, when IsoMacroXSTableParams::ndense is 0:
    unsigned autoNDense( double tol, double wlmax );

    //Reduce a dense sample (x must be sorted, but may contain pairs of
    //identical values to represent discontinuities, and y must be
    //non-negative), with the slopes of the cross section at the sample points,
    //to a subset of points, such that the cubic Hermite interpolation with
    //these slopes reproduces all y values within the tolerance tol (relative
    //to y). The slopes are only modified where needed to ensure that the
    //interpolation is never negative:
    struct ReducedHermitePoints {
      std::vector<std::size_t> idx;//indices of the selected points
      VectD slope;//slopes at the selected points
    };
    ReducedHermitePoints reduceHermitePoints( const VectD& x, const VectD& y,
                                              const VectD& slopes, double tol );

    //Evaluate a table at a given wavelength: cubic Hermite interpolation,
    //where a pair of identical wavelengths (a discontinuity) uses the second
    //point for wavelengths >= the discontinuity, and linear extrapolation of
    //the value and slope of the last point (clamped at 0) beyond the table:
    double evalIsoMacroXSTable( const double* wl, const double* xs,
                                const double* slope, std::size_t n,
                                double wavelength );

  }

}

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace XSTable {
    inline std::size_t IsoMacroXSTable::n() const
    {
      return m_data ? static_cast<std::size_t>( m_data[2] ) : 0;
    }
    inline Span<const double> IsoMacroXSTable::raw() const
    {
      const std::size_t nn = n();
      return nn ? Span<const double>( m_data.get(), m_data.get() + 3 * nn + 4 )
                : Span<const double>();
    }
    inline Span<const double> IsoMacroXSTable::wavelengths() const
    {
      const std::size_t nn = n();
      return nn ? Span<const double>( m_data.get() + 3, m_data.get() + 3 + nn )
                : Span<const double>();
    }
    inline Span<const double> IsoMacroXSTable::macroXS() const
    {
      const std::size_t nn = n();
      return nn ? Span<const double>( m_data.get() + 3 + nn,
                                      m_data.get() + 3 + 2 * nn )
                : Span<const double>();
    }
    inline Span<const double> IsoMacroXSTable::slopes() const
    {
      const std::size_t nn = n();
      return nn ? Span<const double>( m_data.get() + 3 + 2 * nn,
                                      m_data.get() + 3 + 3 * nn )
                : Span<const double>();
    }
    inline std::unique_ptr<double[]> IsoMacroXSTable::release()
    {
      return std::move( m_data );
    }
  }
}

#endif
