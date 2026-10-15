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

#include "NCInfoBldOverrides.hh"
#include "NCInfoBldInternalState.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    bool oDataFieldsEqual( const Info::OverrideableDataFields& a,
                           const Info::OverrideableDataFields& b )
    {
      nc_assert(!std::isnan(a.density.dbl()));
      nc_assert(!std::isnan(b.density.dbl()));
      nc_assert(!std::isnan(a.numberdensity.dbl()));
      nc_assert(!std::isnan(b.numberdensity.dbl()));
      if ( a.density != b.density )
        return false;
      if ( a.numberdensity != b.numberdensity )
        return false;
      if (!Cfg::CfgManip::equal(a.cfgData,b.cfgData))
        return false;
      if ( (a.phases == nullptr) != (b.phases == nullptr) )
        return false;
      if ( a.phases != nullptr ) {
        auto& a_pl = *a.phases;
        auto& b_pl = *b.phases;
        if ( a_pl.size() != b_pl.size() )
          return false;
        for ( auto i : ncrange(a_pl.size()) ) {
          auto& pa = a_pl.at(i);
          auto& pb = b_pl.at(i);
          nc_assert(!std::isnan(pa.first));
          nc_assert(!std::isnan(pb.first));
          if ( pa.first != pb.first )
            return false;
          if ( pa.second->getUniqueID() != pb.second->getUniqueID() )
            return false;
        }
      }
      return true;
    }
  }
}

NC::InfoPtr NC::InfoBuilderDetail::overrideInfoFieldsWithCache( InfoPtr orig,
                                                         const Info::OverrideableDataFields& fields )
{
  if ( oDataFieldsEqual( orig->ovrData(), fields ) )
    return orig;
  auto state =  orig->copyInternalState();
  if ( oDataFieldsEqual( state.data->oData.fields, fields ) ) {
    //compatible with underlying Info object:
    state.oData.reset();
    return makeSO<const Info>( std::move(state) );
  }
  auto& cache = state.data->oDataCache;
  NCRYSTAL_LOCK_GUARD(cache.mtx);
  auto& v = cache.entries;
  for ( auto& e : v ) {
    if ( oDataFieldsEqual( e->fields, fields ) ) {
      //already in cache
      state.oData = e;
      return makeSO<const Info>( std::move(state) );
    }
  }
  auto newOData = makeSO<Info::OverrideableData>();
  newOData->fields = fields;
  state.oData = std::move(newOData);
  static constexpr std::size_t nlim = 1000;
  static constexpr std::size_t nlimpreservefrontback = 250;
  static_assert(nlimpreservefrontback*2<=nlim,"");
  if ( v.size() == nlim ) {
    //TODO: unit test this (was tested manually by lowering numbers above from 1000/250 to 3/1)
    //Protect against out-of-control caches while retaining some
    //repeatability. Preserve oldest + youngest and hope middle ones are most
    //unlikely to be requested again:
    std::decay<decltype(v)>::type v2;
    v2.reserve(nlimpreservefrontback*2+1);//+1 since we will do a final push_back below
    for ( auto i : ncrange(nlimpreservefrontback) )
      v2.push_back( std::move(v.at(i)) );
    for ( auto i : ncrange(nlimpreservefrontback) )
      v2.push_back( std::move(v.at((nlim-nlimpreservefrontback)+i)) );
    v = std::move(v2);
  }
  v.push_back( state.oData );
  return makeSO<const Info>( std::move(state) );
}

bool NC::InfoBuilderDetail::phaseListsIdentical( const Info::PhaseList& a,
                                                const Info::PhaseList& b )
{
  if ( a.size() != b.size() )
    return false;
  for ( auto i : ncrange(a.size()) ) {
    auto pa = a.at(i);
    auto pb = b.at(i);
    if ( ! (pa.first == pb.first ) )
      return false;
    if ( pa.second->getUniqueID() != pb.second->getUniqueID() )
      return false;
  }
  return true;
}
