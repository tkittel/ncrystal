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
  if ( e < E.front() )
    return ( t.tail[0] ? t.tail[0] / std::sqrt( e ) : 0.0 ) + t.tail[1];
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
