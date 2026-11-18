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

#include "NCrystal/internal/dyninfoutils/NCVDOSDebyeGroup.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/phys_utils/NCDebyeMSD.hh"
#include "NCrystal/internal/vdos/NCVDOSDebyeXS.hh"
#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;
namespace NCDI = NCrystal::DynInfoUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace DynInfoUtils {
    namespace {
      //(1-3/neff)^2 (0 for neff<=3):
      double neffFactor( double neff )
      {
        nc_assert_always( neff >= 1.0 - 1e-9 );
        return ncsquare( ncmax( 0.0, 1.0 - 3.0 / neff ) );
      }
    }
  }
}

double NCDI::VDOSDebyeGroupData::tolerance( VDOS::VDOSLux vdoslux,
                                            double neff )
{
  static constexpr double tol[6] = { 0.6, 0.4, 0.25, 0.16, 0.1, 0.06 };
  return tol[std::min<unsigned>( 5u, vdoslux.lvl() )] * neffFactor( neff );
}

bool NCDI::VDOSDebyeGroupData::lowersLevels( double neff )
{
  return neff > 3.0;
}

double NCDI::VDOSDebyeGroupData::effectiveNumberOfEntries( const VecEntries& v )
{
  StableSum s, s2;
  for ( auto& e : v ) {
    s.add( e.scale );
    s2.add( e.scale * e.scale );
  }
  nc_assert_always( s2.sum() > 0.0 );
  return ncsquare( s.sum() ) / s2.sum();
}

namespace NCRYSTAL_NAMESPACE {
  namespace DynInfoUtils {
    namespace {

      //An entry with its scale w and x=log(msd):
      struct Pt { std::size_t idx; double w, x; };

      void validateEntries( const std::vector<DebyeTempAndScale>& entries,
                            Temperature temp, AtomMass mass )
      {
        temp.validate();
        mass.validate();
        if ( entries.empty() )
          NCRYSTAL_THROW( BadInput, "VDOSDebyeGroupData: no entries" );
        for ( auto& e : entries ) {
          e.debyeTemperature.validate();
          if ( !( e.scale > 0.0 && std::isfinite( e.scale ) ) )
            NCRYSTAL_THROW2( BadInput, "VDOSDebyeGroupData: invalid scale (must"
                             " be positive and finite): " << e.scale );
        }
      }

      //Cache key: vdoslux, knllux, T, M, and (TDebye,normalised scale) of
      //each entry, sorted:
      using VDGKey = std::tuple<int,int,double,double,
                                std::vector<std::pair<double,double>>>;

      shared_obj<const VDOSDebyeGroupData> createNoCache( const VDGKey& key )
      {
        const VDOS::VDOSLux vdoslux( std::get<0>( key ) );
        const int knllux = std::get<1>( key );
        const Temperature temp{ std::get<2>( key ) };
        const AtomMass mass{ std::get<3>( key ) };
        std::vector<DebyeTempAndScale> entries;
        for ( auto& e : std::get<4>( key ) )
          entries.push_back( { e.second, DebyeTemperature{ e.first } } );
        using GD = VDOSDebyeGroupData;
        const double tol
          = GD::tolerance( vdoslux, GD::effectiveNumberOfEntries( entries ) );
        const auto clusters = GD::clusterEntries( entries, temp, mass, tol );
        const auto nc = clusters.size();

        //Debye temperature of each cluster, from its mean msd (weighted by
        //scale), or its own for single entries:
        std::vector<DebyeTemperature> td;
        VectD wc;
        for ( auto& c : clusters ) {
          StableSum w, wmsd;
          for ( auto i : c ) {
            w.add( entries[i].scale );
            wmsd.add( entries[i].scale
                      * debyeIsotropicMSD( entries[i].debyeTemperature,
                                           temp, mass ) );
          }
          wc.push_back( w.sum() );
          td.push_back( c.size() == 1
                        ? entries[c.front()].debyeTemperature
                        : debyeTempFromIsotropicMSD( wmsd.sum() / w.sum(),
                                                     temp, mass ) );
        }

        //Expensive work as jobs: cross section tables of all entries (thinned,
        //so merging is cheap), and the sampling models (expansions of the
        //Debye models of the clusters). Thinning twice with reltol 2.5e-5
        //keeps cluster tables within 5e-5 of the exact sums:
        const double thin_reltol = 2.5e-5;
        std::vector<IsoXSTable> tables( entries.size() );
        std::vector<std::shared_ptr<const SAB::SABExtended>> sabexts( nc );
        {
          FactoryJobs jobs;
          for ( auto i : ncrange( entries.size() ) ) {
            const auto tdi = entries[i].debyeTemperature;
            IsoXSTable * out = &tables[i];
            jobs.queue( [out,tdi,temp,mass,thin_reltol]()
            {
              *out = thinIsoXSTable( VDOS::VDOSDebyeXS( tdi, temp, mass )
                                     .createXSTabulation(), thin_reltol );
            } );
          }
          //Levels used for all clusters (see lowersLevels):
          const bool lw = GD::lowersLevels( GD::effectiveNumberOfEntries( entries ) );
          auto lower = []( int l ) { return l > 1 ? ncmax( 1, l - 2 ) : l; };
          const int lv = static_cast<int>( vdoslux.lvl() );
          const VDOS::VDOSLux vdoslux_c
            = lw ? vdoslux.reduceLvl( static_cast<unsigned>( lv - lower( lv ) ) )
                 : vdoslux;
          const int knllux_c = lw ? lower( knllux ) : knllux;
          for ( auto ic : ncrange( nc ) ) {
            const auto tdc = td[ic];
            std::shared_ptr<const SAB::SABExtended> * out = &sabexts[ic];
            jobs.queue( [out,tdc,temp,mass,vdoslux_c,knllux_c]()
            {
              auto sabdata
                = extractSABDataFromVDOSDebyeModel( tdc, temp, SigmaBound{ 1.0 },
                                                    mass, vdoslux_c.expansionLux() );
              *out = SAB::createSABExtendedWithCache
                ( knllux_c, std::move( sabdata ),
                  debyeEffectiveTemperature( tdc, temp ), nullptr,
                  SAB::SABExtended::SupportXS::NO );
            } );
          }
          jobs.waitAll();
        }

        //Cluster tables (summed and thinned again), and their sum. With a
        //single cluster, its table is the total:
        std::vector<VDOSDebyeGroupData::Cluster> cls;
        cls.reserve( nc );
        std::vector<IsoXSTable> ctables;
        ctables.reserve( nc );
        for ( auto ic : ncrange( nc ) ) {
          std::vector<std::pair<double,const IsoXSTable*>> parts;
          for ( auto i : clusters[ic] )
            parts.emplace_back( entries[i].scale, &tables[i] );
          ctables.push_back( thinIsoXSTable( mergeIsoXSTables( parts ),
                                             thin_reltol ) );
        }
        tables.clear();
        IsoXSTable total;
        if ( nc == 1 ) {
          total = std::move( ctables.front() );
        } else {
          std::vector<std::pair<double,const IsoXSTable*>> parts;
          for ( auto& t : ctables )
            parts.emplace_back( 1.0, &t );
          total = mergeIsoXSTables( parts );
        }
        for ( auto ic : ncrange( nc ) ) {
          Optional<IsoXSEvaluator> xs;
          if ( nc > 1 )
            xs.emplace( std::move( ctables[ic] ) );
          cls.push_back( { clusters[ic].size(), wc[ic], td[ic], std::move( xs ),
                           std::move( sabexts[ic] ) } );
        }
        PairDD tdrange{ kInfinity, -kInfinity };
        for ( auto& e : entries ) {
          tdrange.first = ncmin( tdrange.first, e.debyeTemperature.dbl() );
          tdrange.second = ncmax( tdrange.second, e.debyeTemperature.dbl() );
        }
        return makeSO<VDOSDebyeGroupData>( temp, mass, entries.size(), tdrange,
                                           std::move( total ),
                                           std::move( cls ) );
      }

      class VDGFactory final
        : public CachedFactoryBase<VDGKey,VDOSDebyeGroupData,10> {
      public:
        const char* factoryName() const final { return "VDOSDebyeGroupFactory"; }
        std::string keyToString( const VDGKey& key ) const final
        {
          std::ostringstream ss;
          ss << "(vdoslux=" << std::get<0>( key )
             << ";knllux=" << std::get<1>( key )
             << ";T=" << std::get<2>( key ) << "K"
             << ";M=" << std::get<3>( key ) << "u"
             << ";nentries=" << std::get<4>( key ).size() << ")";
          return ss.str();
        }
      protected:
        ShPtr actualCreate( const VDGKey& key ) const final
        {
          return createNoCache( key );
        }
      };
    }
  }
}

std::vector<std::vector<std::size_t>>
NCDI::VDOSDebyeGroupData::clusterEntries( const VecEntries& entries,
                                          Temperature temp, AtomMass mass,
                                          double tolerance )
{
  validateEntries( entries, temp, mass );
  if ( !( tolerance >= 0.0 && std::isfinite( tolerance ) ) )
    NCRYSTAL_THROW2( BadInput, "VDOSDebyeGroupData::clusterEntries: invalid"
                     " tolerance: "
                     << tolerance );
  //Points sorted by (x, input order):
  std::vector<Pt> pts;
  pts.reserve( entries.size() );
  for ( auto i : ncrange( entries.size() ) )
    pts.push_back( { i, entries[i].scale,
                     std::log( debyeIsotropicMSD( entries[i].debyeTemperature,
                                                  temp, mass ) ) } );
  std::stable_sort( pts.begin(), pts.end(),
                    []( const Pt& a, const Pt& b ) { return a.x < b.x; } );
  //Cumulative sums of w, w*x, w*x^2 (x relative to the first, for
  //precision), giving the sse of any range [b,e) in O(1):
  const std::size_t n = pts.size();
  VectD cW( n + 1, 0.0 ), cWx( n + 1, 0.0 ), cWxx( n + 1, 0.0 );
  {
    StableSum sW, sWx, sWxx;
    for ( auto i : ncrange( n ) ) {
      const double w = pts[i].w, x = pts[i].x - pts.front().x;
      sW.add( w );
      sWx.add( w * x );
      sWxx.add( w * x * x );
      cW[i+1] = sW.sum();
      cWx[i+1] = sWx.sum();
      cWxx[i+1] = sWxx.sum();
    }
  }
  auto sse = [&cW,&cWx,&cWxx]( std::size_t b, std::size_t e )
  {
    const double wx = cWx[e] - cWx[b];
    return ncmax( 0.0, ( cWxx[e] - cWxx[b] ) - wx * wx / ( cW[e] - cW[b] ) );
  };
  //Possible range ends (no cuts between identical x values):
  std::vector<std::size_t> ends;
  for ( auto i : ncrange( std::size_t(1), n + 1 ) )
    if ( i == n || pts[i].x > pts[i-1].x )
      ends.push_back( i );
  //Dynamic programming over k=1,2,..: cost[j] is the lowest sse of [0,ends[j])
  //split into k ranges, and from[k-1][j] the start of the last of them (ties:
  //first found). Stop when the full range is within the limit, or no longer
  //improves:
  const double limit = tolerance * tolerance * cW[n];
  const std::size_t m = ends.size();
  VectD cost( m );
  for ( auto j : ncrange( m ) )
    cost[j] = sse( 0, ends[j] );
  std::vector<std::vector<std::size_t>> from;
  from.emplace_back( m, 0 );
  while ( cost.back() > limit && from.size() < m ) {
    const std::size_t k = from.size();//ranges so far
    VectD cnew( m, kInfinity );
    std::vector<std::size_t> fnew( m, 0 );
    for ( auto j : ncrange( k, m ) ) {
      for ( auto i : ncrange( k - 1, j ) ) {
        const double c = cost[i] + sse( ends[i], ends[j] );
        if ( c < cnew[j] ) {
          cnew[j] = c;
          fnew[j] = ends[i];
        }
      }
    }
    if ( !( cnew.back() < cost.back() ) )
      break;
    cost.swap( cnew );
    from.push_back( std::move( fnew ) );
  }
  //Ranges, from the back:
  std::vector<std::pair<std::size_t,std::size_t>> rngs;
  {
    std::size_t e = n;
    for ( auto k = from.size(); k > 0; --k ) {
      const std::size_t je = std::lower_bound( ends.begin(), ends.end(), e )
        - ends.begin();
      nc_assert( je < m && ends[je] == e );
      const std::size_t b = ( k == 1 ? 0 : from[k-1][je] );
      rngs.emplace_back( b, e );
      e = b;
    }
    nc_assert( e == 0 );
  }
  //Decreasing weight (summed directly, so equal weights compare equal; ties:
  //lowest x first):
  std::vector<std::pair<double,std::pair<std::size_t,std::size_t>>> wr;
  for ( auto& r : rngs ) {
    StableSum w;
    for ( auto i : ncrange( r.first, r.second ) )
      w.add( pts[i].w );
    wr.emplace_back( w.sum(), r );
  }
  std::sort( wr.begin(), wr.end(),
             []( const std::pair<double,std::pair<std::size_t,std::size_t>>& a,
                 const std::pair<double,std::pair<std::size_t,std::size_t>>& b )
             {
               return a.first > b.first
                 || ( a.first == b.first && a.second.first < b.second.first );
             } );
  for ( auto i : ncrange( wr.size() ) )
    rngs[i] = wr[i].second;
  std::vector<std::vector<std::size_t>> res;
  res.reserve( rngs.size() );
  for ( auto& r : rngs ) {
    res.emplace_back();
    for ( auto i : ncrange( r.first, r.second ) )
      res.back().push_back( pts[i].idx );
  }
  return res;
}

NCDI::VDOSDebyeGroupData::VDOSDebyeGroupData( Temperature temp, AtomMass mass,
                                              std::size_t nentries,
                                              PairDD tdrange, IsoXSTable xs,
                                              std::vector<Cluster> clusters )
  : m_temp( temp ), m_mass( mass ), m_n( nentries ), m_tdrange( tdrange ),
    m_xs( std::move( xs ) ), m_clusters( std::move( clusters ) )
{
  nc_assert_always( !m_clusters.empty() );
#ifndef NDEBUG
  for ( auto& c : m_clusters )
    nc_assert_always( c.xs.has_value() == ( m_clusters.size() > 1 ) );
#endif
}

NC::ScatterOutcomeIsotropic
NCDI::VDOSDebyeGroupData::sampleScatter( RNG& rng, NeutronEnergy ekin ) const
{
  if ( m_clusters.size() == 1 )
    return m_clusters.front().sabext->sampleScatter( rng, ekin );
  //Walk through the clusters (largest first), subtracting their cross
  //sections from a random fraction of the total:
  const double e = ekin.dbl();
  double r = rng.generate() * m_xs.eval( e );
  const auto itE = std::prev( m_clusters.end() );
  auto it = m_clusters.begin();
  for ( ; it != itE; ++it ) {
    r -= it->xs.value().eval( e );
    if ( r < 0.0 )
      break;
  }
  return it->sabext->sampleScatter( rng, ekin );
}

NC::shared_obj<const NCDI::VDOSDebyeGroupData>
NCDI::createVDOSDebyeGroupData( std::vector<DebyeTempAndScale> entries,
                                Temperature temp, AtomMass mass,
                                VDOS::VDOSLux vdoslux, int knllux )
{
  validateEntries( entries, temp, mass );
  if ( !vdoslux.allowGrouping() || knllux == -1 )
    NCRYSTAL_THROW2( BadInput, "VDOSDebyeGroupData: not supported for vdoslux="
                     << vdoslux.raw() << " and knllux=" << knllux );
  //Canonical key: sorted entries, scales normalised:
  std::vector<std::pair<double,double>> v;
  v.reserve( entries.size() );
  for ( auto& e : entries )
    v.emplace_back( e.debyeTemperature.dbl(), e.scale );
  std::sort( v.begin(), v.end() );
  StableSum stot;
  for ( auto& e : v )
    stot.add( e.second );
  const double inv = 1.0 / stot.sum();
  for ( auto& e : v )
    e.second *= inv;
  static VDGFactory s_fact;
  return s_fact.create( VDGKey{ vdoslux.raw(), knllux, temp.dbl(), mass.dbl(),
                                std::move( v ) } );
}
