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

#include "NCrystal/internal/sab/NCSABBetaStripSampler.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCFastSearch.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {

      std::size_t clampIdx( std::size_t i, std::size_t lo, std::size_t hi )
      {
        return i < lo ? lo : ( i > hi ? hi : i );
      }

      //To be able to select the overlay strip efficiently, the alias method by
      //Walker-Vose (cf.  https://en.wikipedia.org/wiki/Alias_method) for
      //picking indices in proportion to weights w with a single random
      //number is used.
      //FIXME: Move to a common utility and use with pickRandIdxByWeight?
      void initAliasTable( const VectD& w, VectD& prob,
                           std::vector<std::uint32_t>& alias )
      {
        const std::size_t n = w.size();
        StableSum sum;
        for ( auto x : w )
          sum.add( x );
        const double scale = n / sum.sum();
        prob.resize( n );
        alias.resize( n );
        std::vector<std::uint32_t> small, large;
        for ( auto i : ncrange( n ) ) {
          prob[i] = w[i] * scale;
          ( prob[i] < 1.0
            ? small
            : large ).push_back(static_cast<std::uint32_t>(i) );
        }
        while ( !small.empty() && !large.empty() ) {
          const auto s = small.back();
          small.pop_back();
          const auto l = large.back();
          alias[s] = l;
          prob[l] -= ( 1.0 - prob[s] );
          if ( prob[l] < 1.0 ) {
            large.pop_back();
            small.push_back( l );
          }
        }
        for ( auto i : large )
          prob[i] = 1.0;
        for ( auto i : small )
          prob[i] = 1.0;//only due to rounding
      }
    }
  }
}

NCS::SABBetaStripSampler::SABBetaStripSampler( shared_obj<const SABData> sd,
                                               double emax_div_kT,
                                               unsigned nStripsPerBetaCell )
  : m_sab( std::move(sd) ),
    m_emax( emax_div_kT )
{
#ifndef NDEBUG
  nc_assert_always( emax_div_kT > 0.0 && std::isfinite(emax_div_kT) );
  nc_assert_always( nStripsPerBetaCell >= 1 && nStripsPerBetaCell <= 1000 );
#endif
  const auto& ag = m_sab->alphaGrid();
  const auto& bg = m_sab->betaGrid();
  m_logS.reserve( m_sab->sab().size() );
  for ( auto v : m_sab->sab() )
    m_logS.push_back( safeLogOrElse( v, 0.0 ) );

  //Strips covering KB(e) for all e<=emax, weighted by
  //Smax*(b1-b0)*sqrt(b1+emax) (the largest alpha range of KB(e) in the strip
  //is 4*sqrt(e)*sqrt(b1+emax)):
  {
    const double blo = ncmax( -m_emax, bg.front() );
    //Beyond this beta, KB(emax) is above the alpha grid:
    const double bcap = getBetaPlus( m_emax, ag.back() );
    const double bhi = ncmin( bg.back(), bcap );
    if ( !( bhi > blo ) )
      return;//no kernel data within KB(emax)
    const std::uint32_t nacells = static_cast<std::uint32_t>( ag.size()-1 );
    const unsigned n = nStripsPerBetaCell;
    VectD weights;
    //Each beta cell (clipped to [blo,bhi]) is split in n equal strips:
    for ( auto ibcell : ncrange( bg.size() - 1 ) ) {
      const auto ib = static_cast<std::uint32_t>( ibcell );
      const double c0 = ncmax( blo, vectAt( bg, ib ) );
      const double c1 = ncmin( bhi, vectAt( bg, ib + 1 ) );
      if ( !( c1 > c0 ) )
        continue;
      for ( auto j : ncrange( n ) ) {
        const double b0 = c0 + ( c1 - c0 ) * j / n;
        const double b1 = ( j + 1 == n ? c1
                            : c0 + ( c1 - c0 ) * ( j + 1 ) / n );
        const auto sc = calcStripBound( b0, b1, ib );
        const double smax = sc.bound.smax;
        if ( !( smax > 0.0 ) )
          continue;
        const double b1pemax = b1 + m_emax;
        //Alpha cells, with one extra cell on each side as rounding margin:
        const auto ia0 = sc.ia0;
        const auto ia1 = sc.ia1;
        m_strips.push_back( { b0, b1, smax, b1pemax, ib,
                              ( ia0 ? ia0 - 1 : ia0 ),
                              ( ia1 + 1 < nacells ? ia1 + 1 : ia1 ) } );
        weights.push_back( smax * std::sqrt( b1pemax ) * ( b1 - b0 ) );
      }
    }
    if ( m_strips.empty() )
      return;
    initAliasTable( weights, m_aliasProb, m_alias );
    StableSum sumw;
    for ( auto w : weights )
      sumw.add( w );
    m_totalWeight = sumw.sum();
  }
}

NCS::SABBetaStripSampler::StripCalc
NCS::SABBetaStripSampler::calcStripBound( double b0, double b1,
                                          std::uint32_t ib ) const
{
  const auto& ag = m_sab->alphaGrid();
  const auto& bg = m_sab->betaGrid();
  const auto& sab = m_sab->sab();
  const std::size_t na = ag.size();
  //Strip: smallest rectangle containing KB(e) for beta in [b0,b1] and all
  //e<=emax (alpha- is smallest at emax and the beta closest to 0, alpha+
  //largest at emax and b1):
  const double a0 = ( b0 > 0.0 ? getAlphaMinus( m_emax, b0 )
                      : ( b1 < 0.0 ? getAlphaMinus( m_emax, b1 ) : 0.0 ) );
  const double a1 = getAlphaPlus( m_emax, b1 );
  //Alpha grid points in (a0,a1) are ia_lo..ia_hi-1, and the cells containing
  //a0 and a1 are next to them:
  const std::size_t ia_lo = fastUpperBoundIdx( ag.data(), na, a0 );
  const std::size_t ia_hi = fastLowerBoundIdx( ag.data(), na, a1 );
  const auto ia0 = static_cast<std::uint32_t>( clampIdx( ia_lo, 1, na-1 ) - 1 );
  const auto ia1 = static_cast<std::uint32_t>( clampIdx( ia_hi, 1, na-1 ) - 1 );
  //S is convex along alpha and linear along beta within each cell, so its
  //maximum in the strip is at a corner of a strip/cell intersection. First the
  //strip corners (S=0 outside the alpha grid):
  auto smaxAtAlpha = [&]( double a, std::size_t ia )
  {
    return ( a >= ag.front() && a <= ag.back()
             ? ncmax( evalSCell( a, b0, ia, ib ), evalSCell( a, b1, ia, ib ) )
             : 0.0 );
  };
  double smax = ncmax( smaxAtAlpha( a0, ia0 ), smaxAtAlpha( a1, ia1 ) );
  //Then the alpha grid points inside, where S is linear in beta between the
  //stored values:
  {
    const double bb0 = vectAt( bg, ib );
    const double inv_db = 1.0 / ( vectAt( bg, ib + 1 ) - bb0 );
    const double t0 = ( b0 - bb0 ) * inv_db;
    const double t1 = ( b1 - bb0 ) * inv_db;
    for ( std::size_t ia = ia_lo; ia < ia_hi; ++ia ) {
      const double s_lo = vectAt( sab, ib * na + ia );
      const double s_hi = vectAt( sab, ( ib + 1 ) * na + ia );
      smax = ncmax( smax, ncmax( nclerp( s_lo, s_hi, t0 ),
                                 nclerp( s_lo, s_hi, t1 ) ) );
    }
  }
  return { { b0, b1, a0, a1, smax }, ia0, ia1 };
}

std::vector<NCS::SABBetaStripSampler::StripBound>
NCS::SABBetaStripSampler::stripBounds() const
{
  std::vector<StripBound> res;
  res.reserve( m_strips.size() );
  for ( auto& s : m_strips )
    res.push_back( calcStripBound( s.b0, s.b1, s.ib ).bound );
  return res;
}

double NCS::SABBetaStripSampler::evalSCell( double alpha, double beta,
                                            std::size_t ia,
                                            std::size_t ib ) const
{
  const auto& ag = m_sab->alphaGrid();
  const auto& bg = m_sab->betaGrid();
  const std::size_t na = ag.size();
  nc_assert( vectAt( ag, ia ) <= alpha && alpha <= vectAt( ag, ia + 1 ) );
  nc_assert( vectAt( bg, ib ) <= beta && beta <= vectAt( bg, ib + 1 ) );
  const auto& sab = m_sab->sab();
  const double a1 = vectAt( ag, ia );
  const double a2 = vectAt( ag, ia + 1 );
  //Along alpha: log-linear, or linear if a value is zero:
  auto edge = [&]( std::size_t i ) {
    const double s1 = vectAt( sab, i );
    const double s2 = vectAt( sab, i + 1 );
    if ( s1 > 0.0 && s2 > 0.0 )
      return interpolate_loglin_fast2_NEW( a1, s1, a2, s2, alpha,
                                           vectAt( m_logS, i ),
                                           vectAt( m_logS, i + 1 ) ).first;
    return interpolate_linlin_NEW( a1, s1, a2, s2, alpha );
  };
  //Along beta: linear.
  return interpolate_linlin_NEW( vectAt( bg, ib ), edge( ib * na + ia ),
                                 vectAt( bg, ib + 1 ),
                                 edge( ( ib + 1 ) * na + ia ), beta );
}

NCS::SABBetaStripSampler::Result
NCS::SABBetaStripSampler::sample( RNG& rng, double e ) const
{
#ifndef NDEBUG
  nc_assert_always( e > 0.0 && e <= m_emax );
  nc_assert_always( hasStrips() );
#endif
  Result res;
  res.ntries = 0;
  while ( true ) {
    ++res.ntries;
    //Guard against endless loops in pathological cases: warn at 1e4
    //tries, abort at 1e5:
    if ( res.ntries == 10000 ) {
      static Msg::WarnFirstTime wf;
      wf("S(alpha,beta) beta-strip sampling needed 10000 tries.");
    }
    nc_assert_always( res.ntries < 100000 );
    //Pick strip (alias method):
    const double u = rng.generate() * m_alias.size();
    std::size_t istrip = static_cast<std::size_t>( u );
    if ( istrip >= m_alias.size() )
      istrip = m_alias.size() - 1;
    if ( u - istrip >= vectAt( m_aliasProb, istrip ) )
      istrip = vectAt( m_alias, istrip );
    const auto& strip = vectAt( m_strips, istrip );
    const double beta = randInterval( rng, strip.b0, strip.b1 );
    const double bpe = beta + e;
    if ( !( bpe > 0.0 ) )
      continue;
    const auto alim = getAlphaLimits( e, beta );
    const double alpha = randInterval( rng, alim.first, alim.second );
    const auto& ag = m_sab->alphaGrid();
    if ( !( alpha >= ag.front() && alpha <= ag.back() ) )
      continue;//S=0
    //Find alpha cell within ia0..ia1:
    const std::size_t na_strip = strip.ia1 - strip.ia0 + 2;
    const std::size_t ia = std::size_t( strip.ia0 ) - 1
      + clampIdx( fastUpperBoundIdx( ag.data() + strip.ia0, na_strip, alpha ),
                  1, na_strip - 1 );
    const double s = evalSCell( alpha, beta, ia, strip.ib );
    //Accept with probability (S/Smax)*sqrt((beta+e)/(b1+emax)), compared
    //squared. With r=S/Smax<=1 nothing can overflow, and r^2 only underflows
    //for r<~1e-154, where acceptance is anyway negligible:
    const double r = s / strip.smax;
    if ( ncsquare( rng.generate() ) * strip.b1plusemax > ncsquare( r ) * bpe )
      continue;
    res.alpha = alpha;
    res.beta = beta;
    return res;
  }
}
