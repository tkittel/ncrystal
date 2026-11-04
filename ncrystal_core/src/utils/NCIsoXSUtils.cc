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

#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCrystal/internal/utils/NCFastSearch.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

#ifndef NDEBUG
namespace NCRYSTAL_NAMESPACE {
  namespace {
    void validateTabulation( const IsoXSDescr::Tabulation& t )
    {
      const VectD& E = t.E;
      nc_assert_always( E.size() >= 2 );
      nc_assert_always( t.p.empty() || t.p.size() == 4 * ( E.size() - 1 ) );
      nc_assert_always( t.q.empty() || t.q.size() == 4 * ( E.size() - 1 ) );
      nc_assert_always( t.r.empty() || t.r.size() + 1 == E.size() );
      nc_assert_always( E.front() > 0.0 && std::isfinite( E.back() ) );
      for ( auto i : ncrange( std::size_t(1), E.size() ) )
        nc_assert_always( vectAt( E, i - 1 ) < vectAt( E, i ) );
    }
  }
}
#endif

double NC::evalIsoXS( const IsoXSDescr::Tabulation& t, double e )
{
  const VectD& E = t.E;
  //Cheap checks only (validateTabulation is O(n), done once in IsoXSEvaluator):
  nc_assert( E.size() >= 2 );
  nc_assert( t.p.empty() || t.p.size() == 4 * ( E.size() - 1 ) );
  nc_assert( t.q.empty() || t.q.size() == 4 * ( E.size() - 1 ) );
  nc_assert( t.r.empty() || t.r.size() + 1 == E.size() );
  nc_assert( E.front() > 0.0 && std::isfinite( E.back() ) );
  nc_assert( e >= 0.0 );
  if ( e < E.front() ) {
    if ( !t.tail[0] )
      return t.tail[1];
    const double ee = ncmax( std::numeric_limits<double>::denorm_min(), e );
    return t.tail[0] / std::sqrt( ee ) + t.tail[1];
  }
  if ( !( e < E.back() ) )
    return t.tail[3] + ( t.tail[2] + t.tail[4] / e ) / e;
  const std::size_t i = fastUpperBoundIdx( E.data(), E.size(), e ) - 1;
  const double u = ( e - vectAt( E, i ) ) * ( 1.0 / ( vectAt( E, i + 1 ) - vectAt( E, i ) ) );
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
  return ncfma( vectAt( t.r, i ), inv, exs ) * inv;
}

double NC::evalIsoXSDerivative( const IsoXSDescr::Tabulation& t, double e )
{
  const VectD& E = t.E;
  nc_assert( E.size() >= 2 );
  nc_assert( e > 0.0 );
  if ( e < E.front() )
    return t.tail[0] ? -0.5 * t.tail[0] / ( e * std::sqrt( e ) ) : 0.0;
  if ( !( e < E.back() ) )
    return -( t.tail[2] + 2.0 * t.tail[4] / e ) / ( e * e );
  //With F(e) = e*xs = P(u) + sqrt(e)*Q(u) + r/e: xs' = ( F' - xs ) / e
  const std::size_t i = fastUpperBoundIdx( E.data(), E.size(), e ) - 1;
  const double invw = 1.0 / ( vectAt( E, i + 1 ) - vectAt( E, i ) );
  const double u = ( e - vectAt( E, i ) ) * invw;
  auto poly = [u]( const double * c )
  {
    return ncfma( u, ncfma( u, ncfma( u, c[3], c[2] ), c[1] ), c[0] );
  };
  auto dpoly = [u,invw]( const double * c )
  {
    return ncfma( u, ncfma( u, 3.0 * c[3], 2.0 * c[2] ), c[1] ) * invw;
  };
  double F = 0.0, dF = 0.0;
  if ( !t.p.empty() ) {
    const double * c = t.p.data() + 4 * i;
    F = poly( c );
    dF = dpoly( c );
  }
  if ( !t.q.empty() ) {
    const double * d = t.q.data() + 4 * i;
    const double se = std::sqrt( e );
    const double Q = poly( d );
    F = ncfma( se, Q, F );
    dF += 0.5 * Q / se + se * dpoly( d );
  }
  if ( !t.r.empty() ) {
    const double ri = vectAt( t.r, i ) / e;
    F += ri;
    dF -= ri / e;
  }
  return ( dF - F / e ) / e;
}

NC::IsoXSDescr NC::scaleIsoXSDescr( const IsoXSDescr& d, double escale,
                                    double xsscale )
{
  //With x=escale*e: e*xs'(e) = (xsscale/escale)*x*xs(x), so in each interval
  //P, Q and r scale with xsscale/escale, xsscale/sqrt(escale) and
  //xsscale/escale^2, and similarly the tails:
  nc_assert_always( escale > 0.0 && std::isfinite( escale ) );
  nc_assert_always( xsscale >= 0.0 && std::isfinite( xsscale ) );
  const double kp = xsscale / escale;
  const double kq = xsscale / std::sqrt( escale );
  const double kr = kp / escale;
  IsoXSDescr res;
  res.breakpts.reserve( d.breakpts.size() );
  for ( auto b : d.breakpts )
    res.breakpts.push_back( b / escale );
  if ( !d.table.has_value() )
    return res;
  const auto& t = d.table.value();
  res.table.emplace();
  auto& c = res.table.value();
  c.E.reserve( t.E.size() );
  for ( auto e : t.E )
    c.E.push_back( e / escale );
  for ( auto i : ncrange( std::size_t(1), c.E.size() ) )
    if ( !( vectAt( c.E, i - 1 ) < vectAt( c.E, i ) ) )
      NCRYSTAL_THROW(CalcError,"scaleIsoXSDescr: energy scale collapses grid points");
  auto scaled = []( const VectD& v, double k )
  {
    VectD r;
    r.reserve( v.size() );
    for ( auto x : v )
      r.push_back( x * k );
    return r;
  };
  c.p = scaled( t.p, kp );
  c.q = scaled( t.q, kq );
  c.r = scaled( t.r, kr );
  c.tail[0] = t.tail[0] * kq;
  c.tail[1] = t.tail[1] * xsscale;
  c.tail[2] = t.tail[2] * kp;
  c.tail[3] = t.tail[3] * xsscale;
  c.tail[4] = t.tail[4] * kr;
  return res;
}

NC::IsoXSDescr NC::cloneIsoXSDescr( const IsoXSDescr& d )
{
  IsoXSDescr res;
  res.breakpts = d.breakpts;
  if ( d.table.has_value() ) {
    const auto& t = d.table.value();
    res.table.emplace();
    auto& c = res.table.value();
    c.E = t.E;
    c.p = t.p;
    c.q = t.q;
    c.r = t.r;
    for ( auto k : ncrange( 5 ) )
      c.tail[k] = t.tail[k];
  }
  return res;
}

NC::IsoXSEvaluator::IsoXSEvaluator( IsoXSDescr d )
  : m_descr( std::move( d ) )
{
  if ( !m_descr.table.has_value() )
    NCRYSTAL_THROW(BadInput,"IsoXSEvaluator needs an IsoXSDescr with a table");
  const VectD& E = m_descr.table.value().E;
#ifndef NDEBUG
  validateTabulation( m_descr.table.value() );
#endif
  m_invw.reserve( E.size() - 1 );
  for ( auto i : ncrange( std::size_t(1), E.size() ) )
    m_invw.push_back( 1.0 / ( vectAt( E, i ) - vectAt( E, i - 1 ) ) );
  m_guide = LogGuideTable( E.data(), E.size(), 8 );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Re-express c0+c1*u+c2*u^2+c3*u^3 with u=a+b*v as a cubic in v, and add s
    //times it to out:
    void addShiftedCubic( const double* c, double a, double b, double s, double* out )
    {
      const double t3 = c[3];
      const double t2 = ncfma( 3.0 * a, t3, c[2] );//c2+3a*c3
      const double t1 = ncfma( a, ncfma( 3.0 * a, t3, 2.0 * c[2] ), c[1] );//c1+a(2c2+3a*c3)
      const double t0 = ncfma( a, ncfma( a, ncfma( a, t3, c[2] ), c[1] ), c[0] );
      out[0] = ncfma( s, t0, out[0] );
      out[1] = ncfma( s * b, t1, out[1] );
      out[2] = ncfma( s * b * b, t2, out[2] );
      out[3] = ncfma( s * b * b * b, t3, out[3] );
    }
    bool allZero( const VectD& v )
    {
      for ( auto x : v )
        if ( x != 0.0 )
          return false;
      return true;
    }
  }
}

NC::IsoXSDescr NC::mergeIsoXSDescr( Span<const std::pair<double,const IsoXSDescr*>> inputs )
{
  IsoXSDescr res;
  bool allHaveTable = true;
  std::size_t nused = 0;
  for ( auto& in : inputs ) {
    nc_assert_always( in.second != nullptr );
    nc_assert_always( in.first >= 0.0 && std::isfinite( in.first ) );
    if ( in.first == 0.0 )
      continue;
    ++nused;
    res.breakpts.insert( res.breakpts.end(), in.second->breakpts.begin(), in.second->breakpts.end() );
    if ( !in.second->table.has_value() )
      allHaveTable = false;
  }
  std::sort( res.breakpts.begin(), res.breakpts.end() );
  res.breakpts.erase( std::unique( res.breakpts.begin(), res.breakpts.end() ), res.breakpts.end() );
  if ( !nused || !allHaveTable )
    return res;

  IsoXSDescr::Tabulation t;
  for ( auto& in : inputs ) {
    if ( in.first == 0.0 )
      continue;
    const auto& ti = in.second->table.value();
#ifndef NDEBUG
    validateTabulation( ti );
#endif
    t.E.insert( t.E.end(), ti.E.begin(), ti.E.end() );
    for ( auto k : ncrange( 5 ) )
      t.tail[k] = ncfma( in.first, ti.tail[k], t.tail[k] );
  }
  std::sort( t.E.begin(), t.E.end() );
  t.E.erase( std::unique( t.E.begin(), t.E.end() ), t.E.end() );
  const std::size_t nint = t.E.size() - 1;
  t.p.assign( 4 * nint, 0.0 );
  t.q.assign( 4 * nint, 0.0 );
  t.r.assign( nint, 0.0 );
  for ( auto& in : inputs ) {
    const double s = in.first;
    if ( s == 0.0 )
      continue;
    const auto& ti = in.second->table.value();
    for ( auto k : ncrange( nint ) ) {
      //Each merged interval lies within one region of each input, since all
      //input breakpoints are merged breakpoints:
      const double lo = vectAt( t.E, k ), w = vectAt( t.E, k + 1 ) - lo;
      double * P = &t.p[4*k];
      double * Q = &t.q[4*k];
      if ( lo < ti.E.front() ) {
        //e*xs = tail0*sqrt(e) + tail1*e, with e = lo + w*v:
        Q[0] = ncfma( s, ti.tail[0], Q[0] );
        P[0] = ncfma( s * ti.tail[1], lo, P[0] );
        P[1] = ncfma( s * ti.tail[1], w, P[1] );
      } else if ( !( lo < ti.E.back() ) ) {
        //e*xs = tail3*e + tail2 + tail4/e:
        P[0] = ncfma( s, ncfma( ti.tail[3], lo, ti.tail[2] ), P[0] );
        P[1] = ncfma( s * ti.tail[3], w, P[1] );
        t.r[k] = ncfma( s, ti.tail[4], t.r[k] );
      } else {
        const std::size_t i = fastUpperBoundIdx( ti.E.data(), ti.E.size(), lo ) - 1;
        const double Ei = vectAt( ti.E, i ), wi = vectAt( ti.E, i + 1 ) - Ei;
        const double a = ( lo - Ei ) / wi, b = w / wi;
        if ( !ti.p.empty() )
          addShiftedCubic( &ti.p[4*i], a, b, s, P );
        if ( !ti.q.empty() )
          addShiftedCubic( &ti.q[4*i], a, b, s, Q );
        if ( !ti.r.empty() )
          t.r[k] = ncfma( s, vectAt( ti.r, i ), t.r[k] );
      }
    }
  }
  if ( allZero( t.p ) )
    t.p.clear();
  if ( allZero( t.q ) )
    t.q.clear();
  if ( allZero( t.r ) )
    t.r.clear();
  res.table = std::move(t);
  return res;
}
