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

namespace NCRYSTAL_NAMESPACE {

  //Utilities for IsoXSDescr (NCTypes.hh). Its Tabulation t describes a cross
  //section xs(e) (energy e in eV, xs in barn) as:
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
  //
  //IsoXSDescr::breakpts holds the energies of known kinks and discontinuities
  //of the cross section, sorted.

  //Form of a cross section at high energies: xs(E) = B + A/E for E >= emin.
  //Note: emin must be large enough that the form is exact to double precision
  //      (within a few ulp). It becomes the high-energy tail of an
  //      IsoXSDescr::Tabulation (tail[3]=B, tail[2]=A), so only cross sections
  //      with such a form for some finite emin can currently be tabulated
  //      this way.
  struct IsoXSHighEForm final {
    double B, A;
    NeutronEnergy emin;
  };

  //Evaluate the tabulated cross section at e (finding the interval with a
  //binary search, so meant for occasional use):
  double evalIsoXS( const IsoXSDescr::Tabulation&, double e );

  //Merge the descriptions of several components, each with a scale factor, into
  //a description of their summed cross section. The breakpts are the sorted
  //union of all breakpts. The table is the exact sum of the scaled tables
  //(agreeing with the sum of the scaled inputs to within a few ulp), with E the
  //sorted union of all E, and present only if all inputs have a table. Scale
  //factors must be non-negative and finite, and inputs with scale factor 0 are
  //ignored. An empty input gives an empty description:
  IsoXSDescr mergeIsoXSDescr( Span<const std::pair<double,const IsoXSDescr*>> );

  //Explicit copy:
  IsoXSDescr cloneIsoXSDescr( const IsoXSDescr& );

  //The description of xs'(e) = xsscale*xs(escale*e), where xs is the cross
  //section of the given description (escale>0 and xsscale>=0, both finite).
  //For instance, all free-gas cross sections are such scaled versions of a
  //single universal one. It agrees with IsoXSEvaluator::evalScaled to within
  //rounding (at discontinuities, energies within rounding of the jump might
  //evaluate on either side):
  IsoXSDescr scaleIsoXSDescr( const IsoXSDescr&, double escale, double xsscale );

  //Holds an IsoXSDescr (which must have a table) for fast evaluation of its
  //table, giving results bit-identical to evalIsoXS. The interval is found with
  //a LogGuideTable (8 buckets per point, scanning lookup), so evaluation costs
  //one table read, usually at most one comparison step, and a few fused
  //multiply-adds (plus a sqrt in intervals where Q is non-zero).
  class IsoXSEvaluator final : private MoveOnly {
  public:
    explicit IsoXSEvaluator( IsoXSDescr );
    double eval( double e ) const ncnoexceptndebug;//xs(e) in barn, e in eV (e >= 0)
    //The cross section described by scaleIsoXSDescr(descr(),escale,xsscale):
    double evalScaled( double e, double escale, double xsscale ) const ncnoexceptndebug;
    const IsoXSDescr& descr() const noexcept { return m_descr; }
  private:
    IsoXSDescr m_descr;
    VectD m_invw;
    LogGuideTable m_guide;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::IsoXSEvaluator::eval( double e ) const ncnoexceptndebug
{
  nc_assert( e >= 0.0 );
  const auto& t = m_descr.table.value();
  const VectD& E = t.E;
  if ( e < E.front() ) {
    if ( !t.tail[0] )
      return t.tail[1];
    const double ee = std::max( e, std::numeric_limits<double>::denorm_min() );
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

inline double NCrystal::IsoXSEvaluator::evalScaled( double e, double escale,
                                                    double xsscale ) const ncnoexceptndebug
{
  //Below the table, the 1/sqrt form is evaluated in e (as for the scaled
  //description), since the clamping of e=0 and underflows of escale*e would
  //otherwise differ:
  nc_assert( e >= 0.0 && escale > 0.0 && xsscale >= 0.0 );
  const auto& t = m_descr.table.value();
  const double x = escale * e;
  if ( x < t.E.front() ) {
    if ( !t.tail[0] )
      return xsscale * t.tail[1];
    const double ee = std::max( e, std::numeric_limits<double>::denorm_min() );
    return ( t.tail[0] * ( xsscale / std::sqrt( escale ) ) ) / std::sqrt( ee )
      + xsscale * t.tail[1];
  }
  return xsscale * eval( x );
}

#endif
