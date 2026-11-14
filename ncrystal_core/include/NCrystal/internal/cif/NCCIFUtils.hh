#ifndef NCrystal_CIFUtils_hh
#define NCrystal_CIFUtils_hh

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

#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  //Helpers for data in the CIF format.

  //A number in CIF data (e.g. "0.0283", "2344", "1.234e-2(6)") and its
  //rounding error: an upper bound on the difference between the value meant
  //and the double holding it, i.e. half the step of the last printed digit
  //(e.g. 5e-5 for "0.0283", 0.5 for "2344", 5e-7 for "3400e-6") as the
  //smallest double not below that, and at least 1e-15 times the value. A zero
  //without uncertainty in parentheses ("0", "0.", "0.000", "0e-3") is exact
  //(rounding error 0, as it is in practice a constraint or 0 printed with
  //e.g. %.3f). An uncertainty in parentheses after the number (e.g.
  //"0.0283(4)") is otherwise ignored. Leading and trailing whitespace is
  //ignored. Returns NullOpt if the string is not such a number (e.g. "?",
  //".", "1,5", "inf", "1.2(3)e-4"), if the value is not a finite double, or
  //if the last printed digit is a multiple of 1e41 or more (e.g. "1e41"). If
  //it is below 1e-40 (e.g. "1.2e-50"), the rounding error is 5e-41 (a valid
  //but loose bound):
  struct CIFNumber {
    double value;
    double roundingError;
  };
  Optional<CIFNumber> parseCIFNumber( StrView );

  //Element (A=0) or isotope (D and T as Z=1 and A=2 or 3) of a CIF atom type
  //symbol: one or two letters (ignoring case, two if they form an element,
  //e.g. "CO" is cobalt) optionally followed by a charge of digits, decimal
  //points and a sign before or after them (e.g. "O2-", "O-2", "Fe3+", "Fe0",
  //"Mn3.614+", "T+"). Returns NullOpt for anything else (incl. "?" and "."):
  struct CIFElement {
    unsigned Z;
    unsigned A;
  };
  Optional<CIFElement> elementFromCIFTypeSymbol( StrView );

  //Element or isotope deduced from a CIF site label, for data without type
  //symbols: the element (or D or T, as above) named by the first two letters
  //of the label (in their usual case, e.g. "Cl" in "Cl(1)"), or else by its
  //first letter, which also appears in the chemical formula (e.g. the CIF
  //_chemical_formula_sum, cf. tryDecodeChemicalFormula, isotopes counting as
  //their element). Returns NullOpt if neither does, if the formula can not
  //be decoded, or if the label is ambiguous because the first two letters
  //form an element of the formula when ignoring case (e.g. "CL1" for a formula
  //with both C and Cl):
  Optional<CIFElement> elementFromCIFLabel( StrView label, StrView formula );

}

#endif
