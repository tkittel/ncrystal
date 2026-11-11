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

#include "NCrystal/internal/sablegacy/NCSABLegacyFactory.hh"
#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/sablegacy/NCSABIntegrator.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SAB {
    namespace {

      //Cache key is (sabdata uid, egrid uid, legacy opts). The thick key also
      //carries a reference to the SAB data.

      using EGridMargin = SABSampler::EGridMargin;

      struct LegacyOptsDecoded final {
        EGridMargin egridMargin = { EGridMargin::default_value };
        LegacyOptsDecoded( LegacySABAlgOpts opts )
        {
          if ( opts == LegacySABAlgOpts::OVERSAMPLE10 )
            egridMargin = EGridMargin{ 10.0 };
        }
      };

      using ScatHelperCacheKey_Thin = std::tuple<UniqueIDValue,
                                                 UniqueIDValue,
                                                 LegacySABAlgOpts>;

      struct ScatHelperCacheKey {
        ScatHelperCacheKey_Thin thin_key;
        std::shared_ptr<const SABData> sabdata_ptr;
      };

      struct ScatHelperCache_KeyThinner {
        using key_type = ScatHelperCacheKey;
        using thinned_key_type = ScatHelperCacheKey_Thin;
        template <class TMap>
        static typename TMap::mapped_type&
        cacheMapLookup( TMap& map, const key_type& key,
                        Optional<thinned_key_type>& tkey )
        {
          if ( !tkey.has_value() )
            tkey = key.thin_key;
          return map[tkey.value()];
        }
      };

      constexpr auto scathelperfact_nstrongrefskept = 20;

      class ScatterHelperFactory
        : public CachedFactoryBase<ScatHelperCacheKey,
                                   SABScatterHelper,
                                   scathelperfact_nstrongrefskept,
                                   ScatHelperCache_KeyThinner> {
      public:
        const char* factoryName() const final { return "ScatterHelperFactory"; }
        std::string keyToString( const ScatHelperCacheKey& key ) const final
        {
          std::ostringstream ss;
          ss<<"(SABData id="<<std::get<0>(key.thin_key).value
            <<";egrid id="<<std::get<1>(key.thin_key).value;
          if ( std::get<2>(key.thin_key) != LegacySABAlgOpts::DEFAULT ) {
            LegacyOptsDecoded o(std::get<2>(key.thin_key));
            LegacyOptsDecoded odef(LegacySABAlgOpts::DEFAULT);
            if ( o.egridMargin.value != odef.egridMargin.value )
              ss << ";egridmargin=" << o.egridMargin.value;
          }
          ss<<")";
          return ss.str();
        }
      protected:
        virtual ShPtr actualCreate( const ScatHelperCacheKey& key ) const final
        {
          nc_assert( key.sabdata_ptr != nullptr );
          nc_assert( key.sabdata_ptr->getUniqueID() == std::get<0>(key.thin_key) );
          return createScatterHelper( key.sabdata_ptr,
                                      egridFromUniqueID(std::get<1>(key.thin_key)),
                                      std::get<2>(key.thin_key) );
        }
      };

      ScatterHelperFactory& getScatterHelperFactory()
      {
        static ScatterHelperFactory s_scathelperfact;
        return s_scathelperfact;
      }

    }
  }
}

std::unique_ptr<const NC::SAB::SABScatterHelper>
NC::SAB::createScatterHelper( shared_obj<const SABData> data,
                              std::shared_ptr<const VectD> energyGrid,
                              LegacySABAlgOpts opts_raw)
{
  nc_assert(!!data);
  LegacyOptsDecoded opts(opts_raw);
  SABIntegrator si( data,
                    energyGrid.get(),
                    nullptr/*default extender*/,
                    opts.egridMargin );
  auto sh = si.createScatterHelper();
  return ncmake_unique<SABScatterHelper>(std::move(sh));
}

NC::shared_obj<const NC::SAB::SABScatterHelper>
NC::SAB::createScatterHelperWithCache( shared_obj<const SABData> sabdataptr,
                                       std::shared_ptr<const VectD> egrid,
                                       LegacySABAlgOpts opts )
{
  ScatHelperCacheKey key;
  std::get<0>(key.thin_key) = sabdataptr->getUniqueID();
  std::get<1>(key.thin_key) = egridToUniqueID(egrid);
  std::get<2>(key.thin_key) = opts;
  key.sabdata_ptr = std::move(sabdataptr);
  return getScatterHelperFactory().create(key);
}
