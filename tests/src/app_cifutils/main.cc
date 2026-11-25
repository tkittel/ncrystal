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

// Tests of NCCIFUtils: numbers with rounding errors, and elements from atom
// type symbols and site labels.

#include "NCrystal/internal/cif/NCCIFUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <iostream>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace {
  void show( const std::string& s )
  {
    auto r = NC::parseCIFNumber( s );
    std::cout << '"' << s << "\" -> ";
    if ( !r.has_value() ) {
      std::cout << "not a number" << std::endl;
      return;
    }
    std::cout << "value " << NC::dbl2shortstr( r.value().value )
              << ", rounding error " << NC::dbl2shortstr( r.value().roundingError )
              << std::endl;
  }
}

void showElement( const char * what, const std::string& s,
                  const NC::Optional<NC::CIFElement>& e )
{
  std::cout << what << " \"" << s << "\" -> ";
  if ( !e.has_value() )
    std::cout << "none" << std::endl;
  else
    std::cout << "Z=" << e.value().Z << " A=" << e.value().A << std::endl;
}

int main()
{
  for ( const char * s : {
      //Decimals (trailing zeros are printed digits):
      "0.0283", "-0.0032", "+0.5", ".34", "-.5", "1.", "12.75", "0.343400",
      " 0.012 ", "007.50",
      //Zeros without uncertainty are exact:
      "0", "-0", "+0", "00", "0e-3", "0E5", "0.000", "-0.000", "0.", ".0",
      "0.0e3",
      //Uncertainties in parentheses are ignored (but zeros with them are not
      //exact):
      "0.345(34)", "0.0283(4)", "1.23e-4(2)", "12(3)", "0(1)", "0.000(2)",
      //Integers and exponents:
      "2344", "-7", "100", "1234.2e-6", "0.0001E2", "3400e-6", "1e-3",
      "5.e-3", "2.83E+2", "-6.02e23", ".5e1", "1e+000000000005",
      "1e-00000000000000001",
      //Extremes (subnormal, underflow, beyond the table, overflow):
      "1.5e-320", "1e-500", "1e-1000000000", "1.2e-50", "0.0e-50(1)",
      "1e39", "1e40", "1e41", "1e300", "1e308", "1e309", "0.0e9999",
      "1e99999999999999999999",
      //Floor of 1e-15 times the value:
      "0.12345678901234567890", "-123456.789012345678", "0.1234567890123",
      //Not numbers:
      "?", ".", "", " ", "-", "+", "e3", ".e3", "1,5", "inf", "-inf", "nan",
      "0x10", "1.2.3", "1e", "1e+", "1.2(3)e-4", "(3)", "1.2(3", "1.2(a)",
      "1 2", "abc", "'0.1'", "1.0d3", "1.2()" } )
    show( s );

  for ( const char * s : { "O", "O2-", "Fe3+", "Fe0", "Mn3.614+", "T", "T+", "D",
                           "D0", "h", "c", "cl", "zn", "MO", "AL", "CO", "Co", "co",
                           "Cl1-", " Na+ ", "OW", "Ow", "X", "Xx", "C_ar", "O2-x",
                           "Fe 3+", "", "?", ".", "1", "2+" } )
    showElement( "type symbol", s, NC::elementFromCIFTypeSymbol( s ) );
  for ( const char * s : { "CaM1", "OA1", "Fe12", "Fe1", "F1", "C1", "CO1", "Co1",
                           "D1", "T2", "H1A", "Ow1", "OW1", "X1", "C", "1", "Na 1",
                           "NaA B2", "Fe1a2", "" } )
    showElement( "label", s, NC::elementFromCIFLabel( s ) );

  //Rounding errors are exact multiples of powers of ten where possible:
  REQUIRE( NC::parseCIFNumber( "2344" ).value().roundingError == 0.5 );
  REQUIRE( NC::parseCIFNumber( "0.000" ).value().roundingError == 0.0 );
  REQUIRE( !NC::parseCIFNumber( "1e41" ).has_value() );
  return 0;
}
