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

#include "NCFilterQuery.hh"
#include "NCrystal/internal/filter/NCFilterTable.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;
namespace NCF = NCrystal::Filter;

void NCF::JSONQuery( std::ostream& os, const Query& query )
{
  auto invalid = [&query](const std::string& reason){
    std::ostringstream ss;
    ss << "Invalid filtertable query: ";
    streamJSON( ss, query );
    ss << " (" << reason << ")";
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  if ( query.size() < 2 || query.front() != "filtertable" )
    invalid("usage: [\"filtertable\",CFGSTR,OPTIONS...] with options like"
            " \"tol=1e-3\", \"wlmax=500\", \"ndense=0\" (auto),"
            " \"discontinuities=1\" and \"algo=greedy\" (or \"algo=dp\")");
  const std::string cfgstr = query.at(1).to_string();

  TableParams params;
  for ( std::size_t i = 2; i < query.size(); ++i ) {
    auto opt = query.at(i).trimmed();
    auto parts = opt.split<2>('=');
    if ( parts.size() != 2 )
      invalid("options must have the form KEY=VALUE");
    auto key = parts.at(0).trimmed();
    auto val = parts.at(1).trimmed();
    auto getdbl = [&]() {
      auto v = val.toDbl();
      if ( !v.has_value() )
        invalid("invalid value for option "+key.to_string());
      return v.value();
    };
    auto getuint = [&]() {
      auto v = val.toInt();
      if ( !v.has_value() || v.value() < 0 || v.value() > 1000000000 )
        invalid("invalid value for option "+key.to_string());
      return static_cast<unsigned>( v.value() );
    };
    if ( key == "tol" ) {
      params.tol = getdbl();
    } else if ( key == "wlmax" ) {
      params.wlmax = getdbl();
    } else if ( key == "ndense" ) {
      params.ndense = getuint();
    } else if ( key == "discontinuities" ) {
      params.discontinuities = ( getuint() != 0 );
    } else if ( key == "algo" ) {
      if ( val == "greedy" )
        params.algo = ReductionAlgo::Greedy;
      else if ( val == "dp" )
        params.algo = ReductionAlgo::DouglasPeucker;
      else
        invalid("algo must be \"greedy\" or \"dp\"");
    } else {
      invalid("unknown option "+key.to_string());
    }
  }

  auto table = createTable( MatCfg( cfgstr ), params );

  streamJSONDictEntry( os, "cfgstr", cfgstr, JSONDictPos::FIRST );
  streamJSONDictEntry( os, "tol", params.tol );
  streamJSONDictEntry( os, "wlmax", params.wlmax );
  streamJSONDictEntry( os, "ndense", ( params.ndense
                                       ? params.ndense
                                       : autoNDense( params.tol, params.wlmax ) ) );
  streamJSONDictEntry( os, "discontinuities", params.discontinuities );
  streamJSONDictEntry( os, "algo", ( params.algo == ReductionAlgo::Greedy
                                     ? "greedy" : "dp" ) );
  streamJSONDictEntry( os, "processes", table.processes );
  streamJSONDictEntry( os, "numberdensity", table.numberDensity );
  streamJSONDictEntry( os, "ndiscontinuities", table.nDiscontinuities );
  streamJSONDictEntry( os, "ndense_total", table.nDense );
  streamJSONDictEntry( os, "nrefine", table.nRefine );
  streamJSONDictEntry( os, "npts", table.wl.size() );
  streamJSONDictEntry( os, "wl", table.wl );
  streamJSONDictEntry( os, "xs", table.xs, JSONDictPos::LAST );
}
