#ifndef NCrystal_NumFmt_hh
#define NCrystal_NumFmt_hh

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

#include "NCrystal/internal/utils/NCSpan.hh"

namespace NCRYSTAL_NAMESPACE {

  //Formatting of numbers for data files (e.g. NCMAT data written by the
  //NCMATComposer), with a given relative tolerance instead of the full
  //precision, so numerical noise from calculations does not end up in the
  //text. The results are idempotent: formatting the value read back from a
  //text gives the same text.

  //Number of significant digits for a relative tolerance: the smallest n with
  //10^(1-n) <= rtol, or 17 (meaning exact) if that is more than 15, the
  //largest number of decimal digits which always survives a round trip
  //through a double. Rounding to n digits then changes a value by at most
  //rtol/2, and a fraction with the same text (cf. fmtSigDigits) differs from
  //it by at most rtol:
  unsigned numFmtDigits( double rtol );

  //Text for x with ndigits significant digits (17: the shortest text reading
  //back exactly), like printf's "%.<ndigits>g" but without trailing zeros or
  //a trailing ".0", in fixed notation where that is shorter (e.g. "300"
  //instead of "3e+02"), as "0" for zeros, and with stripLeadingZero without
  //the zero of e.g. "0.25" (".25"). For 0<x<1 and maxDenominator>0, a
  //fraction "p/q" (q<=maxDenominator, smallest q first) is written instead
  //if it gives the same text, decided from that text alone (and for exact
  //texts, if the double of p/q is exactly that of the text):
  std::string fmtSigDigits( double x, unsigned ndigits, bool stripLeadingZero,
                            unsigned maxDenominator = 0 );

  //Same, with the digits for a relative tolerance (cf. numFmtDigits):
  std::string fmtRelTol( double x, double rtol, bool stripLeadingZero,
                         unsigned maxDenominator = 0 );

  //Texts for all values as from fmtRelTol. For a grid (isGrid and rtol>0),
  //more digits are used where needed (up to exact texts), so the differences
  //between neighbouring values read back are also within rtol of the
  //original differences:
  VectS fmtRelTolVector( Span<const double>, double rtol, bool isGrid,
                         bool stripLeadingZero, unsigned maxDenominator = 0 );

}

#endif
