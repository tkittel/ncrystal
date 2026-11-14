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

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCMatCfgQuery.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace MatCfgQuery {
    namespace {
      //Decodes cfg of raw data (named, unless name is empty), as
      //MatCfg::createFromRawData but keeping the name. Returns JSON:
      std::string decodeRawData( std::string&& data, std::string&& pars,
                                 std::string&& name, std::string dtype )
      {
        RawStrData rawdata( std::move(data) );
        if ( dtype.empty() )
          dtype = FactImpl::guessDataType( rawdata, name );
        if ( dtype.empty() )
          NCRYSTAL_THROW(BadInput,"Can not determine format of raw data"
                         " (must be specified explicitly in this case)");
        TextData::DataType dt{ std::move(dtype) };
        auto td = ( name.empty()
                    ? makeSO<const TextData>( std::move(rawdata), std::move(dt) )
                    : makeSO<const TextData>( std::move(rawdata), std::move(dt),
                                              DataSourceName( std::move(name) ) ) );
        return MatCfg( std::move(td), std::move(pars) ).toJSONCfg();
      }
    }
  }
}

void NC::MatCfgQuery::JSONQuery( std::ostream& os, const Query& query )
{
  auto invalid = [&query](const char * reason){
    std::ostringstream ss;
    ss << "Invalid matcfg JSON query: ";
    streamJSON( ss, query );//trick: using streamJSON for easy format.
    if ( reason )
      ss<< " ("<<reason<<')';
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  constexpr auto sv_matcfg = StrView::make("matcfg");
  if ( query.size() < 2 || query.front() != sv_matcfg ) {
    invalid(nullptr);
    return;
  }
  const auto& key = query.at(1).trimmed();
  auto argstr = [&query]( std::size_t i ) { return query.at(i+2).to_string(); };
  const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_decode = StrView::make("decode");
  constexpr auto sv_decodedata = StrView::make("decodedata");
  constexpr auto sv_normalise = StrView::make("normalise");
  constexpr auto sv_doc = StrView::make("doc");

  if ( key == sv_decode ) {
    if ( nargs != 1 )
      invalid("correct usage: [\"matcfg\",\"decode\",CFGSTR]");
    os << MatCfg( argstr(0) ).toJSONCfg();
  } else if ( key == sv_decodedata ) {
    constexpr auto sv_rawdata = StrView::make("RAWDATA");
    if ( nargs != 6 || query.at(4) != sv_rawdata || query.at(6) != sv_rawdata )
      invalid("correct usage: [\"matcfg\",\"decodedata\",DATATYPE,"
              "CFGPARAMS,\"RAWDATA\",NAME,\"RAWDATA\",DATA] (DATATYPE,"
              " CFGPARAMS and NAME can be empty strings, and NAME and DATA"
              " can not contain the character 0x07)");
    os << decodeRawData( argstr(5), argstr(1), argstr(3), argstr(0) );
  } else if ( key == sv_normalise ) {
    if ( nargs != 1 )
      invalid("correct usage: [\"matcfg\",\"normalise\",CFGSTR]");
    streamJSON( os, MatCfg( argstr(0) ).toStrCfg() );
  } else if ( key == sv_doc ) {
    constexpr auto sv_full = StrView::make("full");
    constexpr auto sv_short = StrView::make("short");
    constexpr auto sv_json = StrView::make("json");
    const auto usage
      = "correct usage: [\"matcfg\",\"doc\",MODE] with MODE"
        " \"full\", \"short\" or \"json\"";
    if ( nargs != 1 )
      invalid(usage);
    const auto mode = query.at(2).trimmed();
    if ( mode == sv_json ) {
      MatCfg::genDoc( os, MatCfg::GenDocMode::JSON );
    } else if ( mode == sv_full || mode == sv_short ) {
      std::ostringstream ss;
      MatCfg::genDoc( ss, ( mode == sv_full
                            ? MatCfg::GenDocMode::TXT_FULL
                            : MatCfg::GenDocMode::TXT_SHORT ) );
      streamJSON( os, ss.str() );
    } else {
      invalid(usage);
    }
  } else if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"matcfg\",\"list\"]");
    streamJSON( os, std::array<StrView,4>{ sv_decode, sv_decodedata,
                                           sv_normalise, sv_doc } );
  } else {
    invalid(nullptr);
  }
}
