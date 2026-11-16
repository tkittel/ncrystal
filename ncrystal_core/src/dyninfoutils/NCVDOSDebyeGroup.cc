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

double NCDI::VDOSDebyeGroupData::tolerance( VDOS::VDOSLux vdoslux )
{
  static constexpr double tol[6] = { 0.6, 0.5, 0.4, 0.3, 0.2, 0.1 };
  return tol[std::min<unsigned>( 5u, vdoslux.lvl() )];
}

namespace NCRYSTAL_NAMESPACE {
  namespace DynInfoUtils {
    namespace {

      //An entry with its scale w, x=log(msd), and msd:
      struct Pt { std::size_t idx; double w, x, msd; };

      //Contiguous range [b,e) of points sorted by x (a candidate cluster), with
      //its weight, sum of w*(x-mean)^2, and best split (at index split,
      //reducing sse by red):
      struct PtRange {
        std::size_t b, e;
        double W, sse, red;
        std::size_t split;
      };

      PtRange makePtRange( const std::vector<Pt>& p, std::size_t b,
                           std::size_t e )
      {
        nc_assert( b < e && e <= p.size() );
        PtRange c{ b, e, 0.0, 0.0, 0.0, b };
        StableSum W, Wx;
        for ( auto i = b; i < e; ++i ) {
          W.add( p[i].w );
          Wx.add( p[i].w * p[i].x );
        }
        c.W = W.sum();
        const double mu = Wx.sum() / c.W;
        StableSum v;
        for ( auto i = b; i < e; ++i )
          v.add( p[i].w * ( p[i].x - mu ) * ( p[i].x - mu ) );
        c.sse = v.sum();
        //Split into [b,i) and [i,e) reduces sse by WL*WR/W*(muL-muR)^2. Means
        //from both sides, and no splits between identical x values:
        const std::size_t n = e - b;
        VectD rW( n + 1, 0.0 ), rWx( n + 1, 0.0 );
        for ( auto j = n; j > 0; --j ) {
          rW[j-1] = rW[j] + p[b+j-1].w;
          rWx[j-1] = rWx[j] + p[b+j-1].w * p[b+j-1].x;
        }
        double lW = 0.0, lWx = 0.0;
        for ( std::size_t j = 1; j < n; ++j ) {
          lW += p[b+j-1].w;
          lWx += p[b+j-1].w * p[b+j-1].x;
          if ( !( p[b+j].x > p[b+j-1].x ) )
            continue;
          const double d = lWx / lW - rWx[j] / rW[j];
          const double red = lW * rW[j] / ( lW + rW[j] ) * d * d;
          if ( red > c.red ) {
            c.red = red;
            c.split = b + j;
          }
        }
        return c;
      }

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
        const auto clusters = GD::clusterEntries( entries, temp, mass,
                                                  GD::tolerance( vdoslux ) );
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

        //Expensive work as jobs: cross section tables of all entries, and the
        //sampling models (expansions of the Debye models of the clusters):
        std::vector<IsoXSTable> tables( entries.size() );
        std::vector<std::shared_ptr<const SAB::SABExtended>> sabexts( nc );
        {
          FactoryJobs jobs;
          for ( auto i : ncrange( entries.size() ) ) {
            const auto tdi = entries[i].debyeTemperature;
            IsoXSTable * out = &tables[i];
            jobs.queue( [out,tdi,temp,mass]()
            {
              *out = VDOS::VDOSDebyeXS( tdi, temp, mass ).createXSTabulation();
            } );
          }
          for ( auto ic : ncrange( nc ) ) {
            const auto tdc = td[ic];
            std::shared_ptr<const SAB::SABExtended> * out = &sabexts[ic];
            jobs.queue( [out,tdc,temp,mass,vdoslux,knllux]()
            {
              auto sabdata
                = extractSABDataFromVDOSDebyeModel( tdc, temp, SigmaBound{ 1.0 },
                                                    mass, vdoslux.expansionLux() );
              *out = SAB::createSABExtendedWithCache
                ( knllux, std::move( sabdata ),
                  debyeEffectiveTemperature( tdc, temp ), nullptr,
                  SAB::SABExtended::SupportXS::NO );
            } );
          }
          jobs.waitAll();
        }

        //Cluster tables (summed and thinned), and their sum. With a single
        //cluster, its table is the total:
        constexpr double thin_reltol = 5e-5;
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
  for ( auto i : ncrange( entries.size() ) ) {
    const double msd = debyeIsotropicMSD( entries[i].debyeTemperature,
                                          temp, mass );
    pts.push_back( { i, entries[i].scale, std::log( msd ), msd } );
  }
  std::stable_sort( pts.begin(), pts.end(),
                    []( const Pt& a, const Pt& b ) { return a.x < b.x; } );
  std::vector<PtRange> cls;
  cls.push_back( makePtRange( pts, 0, pts.size() ) );
  StableSum sse;
  sse.add( cls.front().sse );
  //Split the cluster with the largest reduction of sse, until sse is at
  //most tolerance^2 * W (ties: first found):
  const double limit = tolerance * tolerance * cls.front().W;
  while ( sse.sum() > limit ) {
    PtRange* best = nullptr;
    for ( auto& c : cls )
      if ( c.red > 0.0 && ( !best || c.red > best->red ) )
        best = &c;
    if ( !best )
      break;
    const PtRange c = *best;
    *best = makePtRange( pts, c.b, c.split );
    const double sse_left = best->sse;
    cls.push_back( makePtRange( pts, c.split, c.e ) );
    sse.add( -c.sse );
    sse.add( sse_left );
    sse.add( cls.back().sse );
  }
  //Decreasing weight (ties: lowest x first):
  std::stable_sort( cls.begin(), cls.end(),
                    []( const PtRange& a, const PtRange& b )
                    { return a.W > b.W || ( a.W == b.W && a.b < b.b ); } );
  std::vector<std::vector<std::size_t>> res;
  res.reserve( cls.size() );
  for ( auto& c : cls ) {
    res.emplace_back();
    for ( auto i = c.b; i < c.e; ++i )
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
