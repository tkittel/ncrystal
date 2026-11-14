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

#include "NCCIFQuery.hh"
#include "NCQueryImpl.hh"
#include "NCrystal/internal/cif/NCCIFUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

void NC::CIFQuery::JSONQuery( std::ostream& os, const Query& rawquery )
{
  const Query query = QueryImpl::expandRawData( rawquery );
  auto invalid = [&query](const char * reason){
    std::ostringstream ss;
    ss << "Invalid cif query: ";
    streamJSON( ss, query );//trick: use streamJSON for easy format.
    if ( reason )
      ss<< " ("<<reason<<')';
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  constexpr auto sv_cif = StrView::make("cif");
  if ( query.size() < 2 || query.front() != sv_cif ) {
    invalid(nullptr);
    return;
  }
  const auto& key = query.at(1).trimmed();
  auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
  const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_numbers = StrView::make("numbers");
  constexpr auto sv_elements = StrView::make("elements");
  if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"cif\",\"list\"]");
    os << "[\"numbers\", \"elements\"]";
  } else if ( key == sv_numbers ) {
    //["cif","numbers",STR,...]: printed CIF numbers (e.g. "0.0283(4)") as
    //[[value,rounding error],...], with null for strings which are not
    //numbers (cf. parseCIFNumber):
    os << '[';
    for ( auto i : ncrange( nargs ) ) {
      if ( i )
        os << ',';
      const auto x = parseCIFNumber( arg(i) );
      if ( !x.has_value() ) {
        os << "null";
        continue;
      }
      os << '[';
      streamJSON( os, x.value().value );
      os << ',';
      streamJSON( os, x.value().roundingError );
      os << ']';
    }
    os << ']';
  } else if ( key == sv_elements ) {
    //["cif","elements","typesymbol",STR,...] and ["cif","elements","label",
    //FORMULA,STR,...]: elements (or isotopes) of CIF atom type symbols, or of
    //site labels checked against a chemical formula, as [[Z,A],...] (A=0 for
    //elements), with null where not recognised (cf. elementFromCIFTypeSymbol
    //and elementFromCIFLabel):
    const bool fromlabel = ( nargs >= 2 && arg(0) == "label" );
    if ( !fromlabel && !( nargs >= 1 && arg(0) == "typesymbol" ) )
      invalid("correct usage: [\"cif\",\"elements\",\"typesymbol\",STR,...]"
              " or [\"cif\",\"elements\",\"label\",FORMULA,STR,...]");
    os << '[';
    for ( auto i : ncrange( std::size_t( fromlabel ? 2 : 1 ), nargs ) ) {
      if ( i > ( fromlabel ? 2u : 1u ) )
        os << ',';
      const auto e = ( fromlabel ? elementFromCIFLabel( arg(i), arg(1) )
                       : elementFromCIFTypeSymbol( arg(i) ) );
      if ( e.has_value() )
        os << '[' << e.value().Z << ',' << e.value().A << ']';
      else
        os << "null";
    }
    os << ']';
  } else {
    invalid("unknown key");
  }
}
