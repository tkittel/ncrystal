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

#include "NCrystal/internal/cif/NCCIFUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {

    //Rounding errors 10^k/2 for k=-40..40, as the smallest doubles not below
    //them (generated with exact fractions in Python, and written with
    //repr(..), so the compiler reproduces exactly these doubles):
    constexpr int s_kmin = -40;
    constexpr int s_kmax = 40;
    constexpr double s_rounding[] = {
      5.000000000000001e-41, 5.0000000000000005e-40, 5.0000000000000005e-39,
      5e-38, 5.0000000000000005e-37, 5e-36, 5.000000000000001e-35, 5e-34,
      5e-33, 5e-32, 5e-31, 5.0000000000000004e-30, 5.0000000000000004e-29,
      5e-28, 5e-27, 5e-26, 5.0000000000000005e-25, 5.0000000000000005e-24,
      5e-23, 5.0000000000000005e-22, 5.0000000000000005e-21,
      5.0000000000000005e-20, 5e-19, 5e-18, 5.0000000000000005e-17, 5e-16,
      5.000000000000001e-15, 5e-14, 5.000000000000001e-13,
      5.0000000000000005e-12, 5e-11, 5e-10, 5e-09, 5.0000000000000004e-08,
      5.000000000000001e-07, 5e-06, 5e-05, 0.0005, 0.005, 0.05, 0.5, 5.0,
      50.0, 500.0, 5000.0, 50000.0, 500000.0, 5000000.0, 50000000.0,
      500000000.0, 5000000000.0, 50000000000.0, 500000000000.0,
      5000000000000.0, 50000000000000.0, 500000000000000.0,
      5000000000000000.0, 5e+16, 5e+17, 5e+18, 5e+19, 5e+20, 5e+21,
      5.0000000000000004e+22, 5.0000000000000006e+23, 5e+24, 5e+25, 5e+26,
      5.000000000000001e+27, 5.0000000000000005e+28, 5e+29,
      5.0000000000000004e+30, 5e+31, 5.0000000000000004e+32,
      5.0000000000000003e+33, 5.000000000000001e+34, 5e+35,
      5.0000000000000004e+36, 5.000000000000001e+37, 5.0000000000000005e+38,
      5e+39
    };
    static_assert( sizeof(s_rounding)/sizeof(double) == s_kmax - s_kmin + 1, "" );

    bool isDigit( char c ) { return c >= '0' && c <= '9'; }

    StrView takeDigits( StrView s, std::size_t& i )
    {
      const auto i0 = i;
      while ( i < s.size() && isDigit( s[i] ) )
        ++i;
      return s.substr( i0, i - i0 );
    }

    bool isUpper( char c ) { return c >= 'A' && c <= 'Z'; }
    bool isLower( char c ) { return c >= 'a' && c <= 'z'; }
    char toUpper( char c ) { return isLower( c ) ? char( c - 'a' + 'A' ) : c; }
    char toLower( char c ) { return isUpper( c ) ? char( c - 'A' + 'a' ) : c; }

    //Element or isotope of a name in its usual case (e.g. "Fe", "D"):
    Optional<CIFElement> elementFromName( const std::string& name )
    {
      if ( name == "D" )
        return CIFElement{ 1, 2 };
      if ( name == "T" )
        return CIFElement{ 1, 3 };
      const unsigned z = elementNameToZ( name );
      if ( z )
        return CIFElement{ z, 0 };
      return NullOpt;
    }

    StrView withoutLeadingZeros( StrView s )
    {
      std::size_t i = 0;
      while ( i < s.size() && s[i] == '0' )
        ++i;
      return s.substr( i );
    }
  }
}

NC::Optional<NC::CIFNumber> NC::parseCIFNumber( StrView sv )
{
  //Syntax: [+-]digits[.digits][(e|E)[+-]digits][(digits)], with at least one
  //digit before the exponent:
  const StrView s = sv.trimmed();
  std::size_t i = 0;
  const bool negative = ( i < s.size() && s[i] == '-' );
  if ( i < s.size() && ( s[i] == '-' || s[i] == '+' ) )
    ++i;
  const StrView ipart = takeDigits( s, i );
  StrView fpart;
  if ( i < s.size() && s[i] == '.' ) {
    ++i;
    fpart = takeDigits( s, i );
  }
  if ( ipart.empty() && fpart.empty() )
    return NullOpt;
  bool exp_negative = false;
  StrView exp_digits;
  if ( i < s.size() && ( s[i] == 'e' || s[i] == 'E' ) ) {
    ++i;
    exp_negative = ( i < s.size() && s[i] == '-' );
    if ( i < s.size() && ( s[i] == '-' || s[i] == '+' ) )
      ++i;
    exp_digits = takeDigits( s, i );
    if ( exp_digits.empty() )
      return NullOpt;
  }
  bool has_uncertainty = false;
  if ( i < s.size() && s[i] == '(' ) {
    ++i;
    if ( takeDigits( s, i ).empty() || i >= s.size() || s[i] != ')' )
      return NullOpt;
    ++i;
    has_uncertainty = true;
  }
  if ( i != s.size() )
    return NullOpt;

  //Decimal exponent (huge ones clamped, they can only give 0 or overflow):
  const StrView exp_sig = withoutLeadingZeros( exp_digits );
  constexpr std::int64_t exp_clamp = 1000000000000;
  std::int64_t e = ( exp_sig.size() > 12 ? exp_clamp
                     : ( exp_sig.empty() ? 0 : str2int64( exp_sig ) ) );
  if ( exp_negative )
    e = -e;

  //Value. The decimal exponent of the leading nonzero digit decides about
  //underflow (to zero, as in Python) and overflow (NullOpt):
  const StrView isig = withoutLeadingZeros( ipart );
  const StrView fsig = withoutLeadingZeros( fpart );
  double value = 0.0;
  if ( !isig.empty() || !fsig.empty() ) {
    const std::int64_t lead
      = e + ( !isig.empty() ? static_cast<std::int64_t>( isig.size() ) - 1
              : -static_cast<std::int64_t>( fpart.size() - fsig.size() ) - 1 );
    if ( lead > 308 )
      return NullOpt;
    if ( lead >= -400 ) {
      std::string txt = ( ipart.empty() ? std::string("0") : ipart.to_string() );
      if ( !fpart.empty() )
        txt += '.' + fpart.to_string();
      txt += 'e' + std::to_string( e );
      if ( !safe_str2dbl( txt, value ) ) {
        if ( lead > 0 )
          return NullOpt;//overflow
        value = 0.0;//underflow
      }
      if ( !std::isfinite( value ) )
        return NullOpt;
    }
    if ( negative )
      value = -value;
  }
  if ( value == 0.0 ) {
    value = 0.0;//(no negative zero)
    if ( !has_uncertainty )
      return CIFNumber{ 0.0, 0.0 };
  }

  //Rounding error, half the step of the last printed digit, 10^k/2:
  const std::int64_t k = e - static_cast<std::int64_t>( fpart.size() );
  if ( k > s_kmax )
    return NullOpt;
  //(below the table, its first entry is still an upper bound):
  const double r = s_rounding[ k < s_kmin ? 0 : k - s_kmin ];
  return CIFNumber{ value, ncmax( r, ncabs( value ) * 1e-15 ) };
}

NC::Optional<NC::CIFElement> NC::elementFromCIFTypeSymbol( StrView sv )
{
  const StrView s = sv.trimmed();
  std::size_t nletters = 0;
  while ( nletters < s.size() && ( isUpper( s[nletters] ) || isLower( s[nletters] ) ) )
    ++nletters;
  if ( nletters < 1 || nletters > 2 )
    return NullOpt;
  Optional<CIFElement> res;
  if ( nletters == 2 ) {
    res = elementFromName( std::string{ toUpper( s[0] ), toLower( s[1] ) } );
    if ( !res.has_value() )
      return NullOpt;//(e.g. "OW")
  } else {
    res = elementFromName( std::string( 1, toUpper( s[0] ) ) );
    if ( !res.has_value() )
      return NullOpt;
  }
  //Optional charge, like "2-", "-2", "3.614+", "0" or "+":
  auto isSign = []( char c ) { return c == '+' || c == '-'; };
  std::size_t i = nletters;
  const bool signFirst = ( i < s.size() && isSign( s[i] ) );
  if ( signFirst )
    ++i;
  while ( i < s.size() && ( isDigit( s[i] ) || s[i] == '.' ) )
    ++i;
  if ( !signFirst && i < s.size() && isSign( s[i] ) )
    ++i;
  if ( i != s.size() )
    return NullOpt;
  return res;
}

NC::Optional<NC::CIFElement> NC::elementFromCIFLabel( StrView label,
                                                      StrView formula )
{
  auto terms = tryDecodeChemicalFormula( formula );
  if ( !terms.has_value() )
    return NullOpt;
  auto inFormula = [&terms]( unsigned z ) {
    for ( auto& e : terms.value() )
      if ( e.second.Z() == z )
        return true;
    return false;
  };
  const StrView s = label.trimmed();
  if ( s.empty() )
    return NullOpt;
  if ( s.size() >= 2 ) {
    auto res = elementFromName( s.substr( 0, 2 ).to_string() );
    if ( res.has_value() && inFormula( res.value().Z ) )
      return res;
    //Ambiguous if the two letters only form an element of the formula when
    //ignoring case (e.g. "CL1" with both C and Cl in the formula):
    if ( isUpper( s[1] ) ) {
      auto ci = elementFromName( std::string{ s[0], toLower( s[1] ) } );
      if ( ci.has_value() && inFormula( ci.value().Z ) )
        return NullOpt;
    }
  }
  auto res = elementFromName( s.substr( 0, 1 ).to_string() );
  if ( res.has_value() && inFormula( res.value().Z ) )
    return res;
  return NullOpt;
}
