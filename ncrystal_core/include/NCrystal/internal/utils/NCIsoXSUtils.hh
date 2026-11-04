#ifndef NCrystal_IsoXSUtils_hh
#define NCrystal_IsoXSUtils_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCFastSearch.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCFMA.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NCRYSTAL_NAMESPACE {

  //Utilities for IsoXSTable (NCTypes.hh), which describes a cross section
  //xs(e) (energy e in eV, xs in barn) as:
  //
  //  e < E[0]:           xs = tail[0]/sqrt(e) + tail[1]
  //  E[i] <= e < E[i+1]: e*xs = P_i(u) + sqrt(e)*Q_i(u) + r[i]/e
  //  e >= E.back():      xs = tail[3] + tail[2]/e + tail[4]/e^2
  //
  //with u=(e-E[i])/(E[i+1]-E[i]), and P_i and Q_i cubics in u with
  //coefficients p[4i..4i+3] and q[4i..4i+3] (lowest order first). Empty p, q
  //or r means zero. E must be positive and strictly increasing, with at least 2
  //entries, and e*xs must be finite within [E[0],E.back()] (no problem at
  //physical energies). The r[i] terms allow a sum of tabulations to be a
  //tabulation. At e=0 the first form uses the smallest positive double
  //instead (keeping xs finite).

  class IsoXSEvaluator final : private MoveOnly {
  public:
    //Provides IsoXSTable evaluations. Holds an IsoXSTable, and adds relevant
    //cached data allowing for subsequent fast repeated evaluations.
    //Alternatively, the oneShotEval method provides bit-identical results
    //without initialisation, at the cost of a binary search.

    explicit IsoXSEvaluator( IsoXSTable );
    const IsoXSTable& table() const noexcept { return m_table; }

    //Evaluations accept ekin in eV and return xs in barn:
    double eval( double ekin ) const ncnoexceptndebug;
    static double oneShotEval( const IsoXSTable&, double ekin );

    //Scaled evaluation (cf. scaleIsoXSTable):
    double evalScaled( double ekin,
                       double escale,
                       double xsscale ) const ncnoexceptndebug;


  private:
    IsoXSTable m_table;
    VectD m_invw;
    LogGuideTable m_guide;
  };

  //Form of a cross section at high energies: xs(E) = B + A/E for E >= emin.
  //Note: emin must be large enough that the form is exact to double precision
  //      (within a few ulp). It becomes the high-energy tail of an
  //      IsoXSTable (tail[3]=B, tail[2]=A), so only cross sections
  //      with such a form for some finite emin can currently be tabulated
  //      this way.
  struct IsoXSHighEForm final {
    double B, A;
    NeutronEnergy emin;
  };

  //Merge the tables of several components, each with a scale factor, into a
  //table of their summed cross section: the exact sum of the scaled tables
  //(agreeing with the sum of the scaled inputs to within a few ulp), with E the
  //sorted union of all E. Scale factors must be non-negative and finite, and
  //inputs with scale factor 0 are ignored, but at least one input must
  //remain:
  IsoXSTable mergeIsoXSTables( Span<const std::pair<double,const IsoXSTable*>> );

  //Explicit copy:
  IsoXSTable cloneIsoXSTable( const IsoXSTable& );

  //Sort energies of kinks and discontinuities, removing duplicates:
  void sortBreakpoints( VectD& );

  //The derivative d(xs)/de of the tabulated cross section at e > 0, in the
  //same piece of the table as evalIsoXS (so at the energies E[i] it is the
  //derivative from above):
  double evalIsoXSDerivative( const IsoXSTable&, double e );

  //Fixme: Remove this standalone evalIsoXS and use two functions above. Also
  //       make all 3 evaluations methods in IsoXSEvaluator return CrossSect
  //       type:
  double evalIsoXS( const IsoXSTable&, double e );

  //The table of xs'(e) = xsscale*xs(escale*e), where xs is the cross section
  //of the given table (escale>0 and xsscale>=0, both finite).
  //For instance, all free-gas cross sections are such scaled versions of a
  //single universal one. It agrees with IsoXSEvaluator::evalScaled to within
  //rounding (at discontinuities, energies within rounding of the jump might
  //evaluate on either side):
  IsoXSTable scaleIsoXSTable( const IsoXSTable&, double escale, double xsscale );

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::IsoXSEvaluator::eval( double e ) const ncnoexceptndebug
{
  nc_assert( e >= 0.0 );
  const auto& t = m_table;
  const VectD& E = t.E;
  if ( e < E.front() ) {
    if ( !t.tail[0] )
      return t.tail[1];
    const double ee = ncmax( std::numeric_limits<double>::denorm_min(), e );
    return t.tail[0] / std::sqrt( ee ) + t.tail[1];
  }
  if ( !( e < E.back() ) )
    return t.tail[3] + ( t.tail[2] + t.tail[4] / e ) / e;
  const std::size_t i = m_guide.upperBoundIdxScan( E.data(), E.size(), e ) - 1;
  const double u = ( e - E[i] ) * m_invw[i];
  double exs = 0.0;
  if ( !t.p.empty() ) {
    const double * c = t.p.data() + 4 * i;
    exs = ncfma( u, ncfma( u, ncfma( u, c[3], c[2] ), c[1] ), c[0] );
  }
  if ( !t.q.empty() ) {
    const double * d = t.q.data() + 4 * i;
    exs = ncfma( std::sqrt( e ),
                 ncfma( u, ncfma( u, ncfma( u, d[3], d[2] ), d[1] ), d[0] ),
                 exs );
  }
  if ( t.r.empty() )
    return exs / e;
  const double inv = 1.0 / e;
  return ncfma( t.r[i], inv, exs ) * inv;
}

inline double NCrystal::IsoXSEvaluator::
oneShotEval( const IsoXSTable& tbl, double ekin )
{
  return evalIsoXS( tbl, ekin );
}

inline double NCrystal::IsoXSEvaluator::
evalScaled( double e, double escale, double xsscale ) const ncnoexceptndebug
{
  //Below the table, the 1/sqrt form is evaluated in e (as for the scaled
  //table), since the clamping of e=0 and underflows of escale*e would
  //otherwise differ:
  nc_assert( e >= 0.0 && escale > 0.0 && xsscale >= 0.0 );
  const auto& t = m_table;
  const double x = escale * e;
  if ( x < t.E.front() ) {
    if ( !t.tail[0] )
      return xsscale * t.tail[1];
    const double ee = ncmax( std::numeric_limits<double>::denorm_min(), e );
    return ( t.tail[0] * ( xsscale / std::sqrt( escale ) ) ) / std::sqrt( ee )
      + xsscale * t.tail[1];
  }
  return xsscale * eval( x );
}

#endif
