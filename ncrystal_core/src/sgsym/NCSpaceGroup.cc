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

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCSGSymTable.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    const SGSymTableEntry& entry( std::uint16_t hn )
    {
      return sgsymTableEntry( SpaceGroupSettingIdx{ hn } );
    }

    //Setting index of the first (i.e. default) setting of a given number, by
    //binary search (the table is sorted by number):
    std::uint16_t firstSettingIdx( unsigned number )
    {
      nc_assert( number >= 1 && number <= 230 );
      std::uint16_t lo = SpaceGroupSettingIdx::min_value;
      std::uint16_t hi = SpaceGroupSettingIdx::max_value;
      while ( lo < hi ) {
        const std::uint16_t mid = static_cast<std::uint16_t>( ( lo + hi ) / 2 );
        if ( entry( mid ).number < number )
          lo = static_cast<std::uint16_t>( mid + 1 );
        else
          hi = mid;
      }
      nc_assert( entry( lo ).number == number );
      return lo;
    }

    [[noreturn]] void throwBadSpec( StrView s, const std::string& extra )
    {
      NCRYSTAL_THROW2(BadInput,"Invalid space group specification: \""<<s
                      <<"\" ("<<extra<<")");
    }
  }
}

std::vector<NC::SpaceGroup> NC::SpaceGroup::settingsOfNumber( unsigned number )
{
  nc_assert_always( number >= 1 && number <= 230 );
  std::vector<SpaceGroup> res;
  for ( auto hn = firstSettingIdx( number );
        hn <= SpaceGroupSettingIdx::max_value && entry( hn ).number == number;
        ++hn )
    res.emplace_back( SpaceGroupSettingIdx{ hn } );
  return res;
}

std::string NC::SpaceGroup::listSettings( unsigned number )
{
  std::string res;
  for ( auto& sg : settingsOfNumber( number ) ) {
    if ( !res.empty() )
      res += ", ";
    res += sg.toString();
    if ( sg.isDefaultSetting() && sg.hasMultipleSettings() )
      res += " (default)";
  }
  return res;
}

NC::SpaceGroup::SpaceGroup( SpaceGroupSettingIdx hn )
  : m_idx( DoValidate, hn )
{
}

NC::SpaceGroup::SpaceGroup( StrView s, SpaceGroupParseMode mode )
{
  //Parse "<number>" or "<number>:<choice>":
  const auto parts = s.split<2>( ':' );
  const char * expected = ( "expected \"<number>\" or \"<number>:<choice>\""
                            " with a number in 1..230" );
  if ( parts.size() > 2 || parts.front().empty()
       || parts.front().size() > 3 )
    throwBadSpec( s, expected );
  unsigned number = 0;
  for ( char c : parts.front() ) {
    if ( !( c >= '0' && c <= '9' ) )
      throwBadSpec( s, expected );
    number = 10 * number + static_cast<unsigned>( c - '0' );
  }
  if ( !( number >= 1 && number <= 230 ) )
    throwBadSpec( s, expected );
  const std::uint16_t hn0 = firstSettingIdx( number );
  const bool multi = SpaceGroup( SpaceGroupSettingIdx{ hn0 } ).hasMultipleSettings();
  if ( parts.size() == 1 ) {
    //Bare number: the default setting, which for numbers with several
    //settings is only allowed in relaxed mode:
    if ( multi && mode == SpaceGroupParseMode::Strict )
      throwBadSpec( s, "space group " + std::to_string( number )
                    + " has several settings, so one of them must be"
                    " specified: " + listSettings( number ) );
    m_idx = SpaceGroupSettingIdx{ hn0 };
    return;
  }
  const StrView choice = parts.back();
  if ( !multi )
    throwBadSpec( s, "space group " + std::to_string( number )
                  + " has only one setting, which is written as just \""
                  + std::to_string( number ) + "\"" );
  auto matches = [choice,number]( std::uint16_t h, bool ignore_case )
  {
    const StrView c( entry( h ).choice, std::strlen( entry( h ).choice ) );
    if ( entry( h ).number != number || c.size() != choice.size() )
      return false;
    for ( std::size_t i = 0; i < c.size(); ++i ) {
      char a = choice[i];
      if ( ignore_case && a >= 'A' && a <= 'Z' )
        a = static_cast<char>( a - 'A' + 'a' );
      if ( a != c[i] )
        return false;
    }
    return true;
  };
  for ( auto h = hn0; h <= SpaceGroupSettingIdx::max_value
          && entry( h ).number == number; ++h ) {
    if ( matches( h, false ) ) {
      m_idx = SpaceGroupSettingIdx{ h };
      return;
    }
  }
  //Codes are lower case as in the table, so point out a case-only mismatch:
  for ( auto h = hn0; h <= SpaceGroupSettingIdx::max_value
          && entry( h ).number == number; ++h ) {
    if ( matches( h, true ) )
      throwBadSpec( s, "setting choice codes are case-sensitive and lower"
                    " case, change it to \""
                    + SpaceGroup( SpaceGroupSettingIdx{ h } ).toString()
                    + "\" if that is the intention" );
  }
  throwBadSpec( s, "not a setting of space group " + std::to_string( number )
                + ", the settings are: " + listSettings( number ) );
}

unsigned NC::SpaceGroup::number() const ncnoexceptndebug
{
  return entry( m_idx.get() ).number;
}

const char * NC::SpaceGroup::choice() const ncnoexceptndebug
{
  return entry( m_idx.get() ).choice;
}

const char * NC::SpaceGroup::hallSymbol() const ncnoexceptndebug
{
  return entry( m_idx.get() ).hall;
}

const char * NC::SpaceGroup::hmSymbol() const ncnoexceptndebug
{
  return entry( m_idx.get() ).hm;
}

bool NC::SpaceGroup::isDefaultSetting() const ncnoexceptndebug
{
  const auto hn = m_idx.get();
  return hn == 1 || entry( hn - 1 ).number != entry( hn ).number;
}

bool NC::SpaceGroup::hasMultipleSettings() const ncnoexceptndebug
{
  const auto hn = m_idx.get();
  const auto n = entry( hn ).number;
  return ( ( hn > 1 && entry( hn - 1 ).number == n )
           || ( hn < SpaceGroupSettingIdx::max_value
                && entry( hn + 1 ).number == n ) );
}

std::string NC::SpaceGroup::toString() const
{
  //NB: Only the numbers with a single setting have empty choice codes:
  std::string res = std::to_string( number() );
  if ( *choice() ) {//i.e. choice code is not the empty string
    res += ':';
    res += choice();
  }
  return res;
}

std::ostream& NC::operator<<( std::ostream& os, const SpaceGroup& sg )
{
  return os << sg.toString();
}

NC::SGCrystalSystem NC::SpaceGroup::crystalSystem() const ncnoexceptndebug
{
  const unsigned n = number();
  if ( n <= 2 )
    return SGCrystalSystem::Triclinic;
  if ( n <= 15 )
    return SGCrystalSystem::Monoclinic;
  if ( n <= 74 )
    return SGCrystalSystem::Orthorhombic;
  if ( n <= 142 )
    return SGCrystalSystem::Tetragonal;
  if ( n <= 167 )
    return SGCrystalSystem::Trigonal;
  if ( n <= 194 )
    return SGCrystalSystem::Hexagonal;
  return SGCrystalSystem::Cubic;
}

NC::SGSettingInfo NC::SpaceGroup::settingInfo() const ncnoexceptndebug
{
  SGSettingInfo res{ SGHexFamilyAxes::NotApplicable,
                     SGUniqueAxis::NotApplicable,
                     SGCellChoice::NotApplicable,
                     SGOriginChoice::NotApplicable,
                     SGAxisPermutation::NotApplicable };
  //NB: The choice codes come from the hardwired table (fully verified by
  //tests), so anything unexpected would be a bug:
  const std::string c = choice();
  auto originFromDigit = []( char ch )
  {
    nc_assert( ch == '1' || ch == '2' );
    return ch == '1' ? SGOriginChoice::Choice1 : SGOriginChoice::Choice2;
  };
  switch ( crystalSystem() ) {
  case SGCrystalSystem::Triclinic:
    nc_assert( c.empty() );
    break;
  case SGCrystalSystem::Monoclinic:
    {
      //Codes like "b", "b1", "-c2":
      std::size_t i = 0;
      const bool minus = ( !c.empty() && c[0] == '-' );
      if ( minus )
        ++i;
      nc_assert( i < c.size() );
      const char ax = c[i++];
      nc_assert( ax == 'a' || ax == 'b' || ax == 'c' );
      res.uniqueAxis = ( ax == 'a'
                         ? ( minus ? SGUniqueAxis::minus_a : SGUniqueAxis::a )
                         : ax == 'b'
                         ? ( minus ? SGUniqueAxis::minus_b : SGUniqueAxis::b )
                         : ( minus ? SGUniqueAxis::minus_c : SGUniqueAxis::c ) );
      if ( i < c.size() ) {
        nc_assert( i + 1 == c.size() );
        const char cc = c[i];
        nc_assert( cc >= '1' && cc <= '3' );
        res.cellChoice = ( cc == '1' ? SGCellChoice::Choice1
                           : cc == '2' ? SGCellChoice::Choice2
                           : SGCellChoice::Choice3 );
      }
    }
    break;
  case SGCrystalSystem::Orthorhombic:
    {
      //Codes like "abc", "cab", "1", "2ba-c":
      std::string perm = c;
      if ( !perm.empty() && ( perm[0] == '1' || perm[0] == '2' ) ) {
        res.originChoice = originFromDigit( perm[0] );
        perm = perm.substr( 1 );
      }
      if ( perm.empty() || perm == "abc" )
        res.axisPermutation = SGAxisPermutation::abc;
      else if ( perm == "ba-c" )
        res.axisPermutation = SGAxisPermutation::ba_mc;
      else if ( perm == "cab" )
        res.axisPermutation = SGAxisPermutation::cab;
      else if ( perm == "-cba" )
        res.axisPermutation = SGAxisPermutation::mcba;
      else if ( perm == "bca" )
        res.axisPermutation = SGAxisPermutation::bca;
      else if ( perm == "a-cb" )
        res.axisPermutation = SGAxisPermutation::a_mcb;
      else
        nc_assert( false );
    }
    break;
  case SGCrystalSystem::Tetragonal:
  case SGCrystalSystem::Cubic:
    if ( !c.empty() ) {
      nc_assert( c.size() == 1 );
      res.originChoice = originFromDigit( c[0] );
    }
    break;
  case SGCrystalSystem::Trigonal:
  case SGCrystalSystem::Hexagonal:
    nc_assert( c.empty() || c == "h" || c == "r" );
    res.hexFamilyAxes = ( c == "r" ? SGHexFamilyAxes::Rhombohedral
                          : SGHexFamilyAxes::Hexagonal );
    break;
  }
  return res;
}

std::ostream& NC::operator<<( std::ostream& os, SGCrystalSystem cs )
{
  switch ( cs ) {
  case SGCrystalSystem::Triclinic: return os << "triclinic";
  case SGCrystalSystem::Monoclinic: return os << "monoclinic";
  case SGCrystalSystem::Orthorhombic: return os << "orthorhombic";
  case SGCrystalSystem::Tetragonal: return os << "tetragonal";
  case SGCrystalSystem::Trigonal: return os << "trigonal";
  case SGCrystalSystem::Hexagonal: return os << "hexagonal";
  case SGCrystalSystem::Cubic: return os << "cubic";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGHexFamilyAxes v )
{
  switch ( v ) {
  case SGHexFamilyAxes::NotApplicable: return os << "n/a";
  case SGHexFamilyAxes::Hexagonal: return os << "hexagonal";
  case SGHexFamilyAxes::Rhombohedral: return os << "rhombohedral";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGUniqueAxis v )
{
  switch ( v ) {
  case SGUniqueAxis::NotApplicable: return os << "n/a";
  case SGUniqueAxis::a: return os << "a";
  case SGUniqueAxis::b: return os << "b";
  case SGUniqueAxis::c: return os << "c";
  case SGUniqueAxis::minus_a: return os << "-a";
  case SGUniqueAxis::minus_b: return os << "-b";
  case SGUniqueAxis::minus_c: return os << "-c";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGCellChoice v )
{
  switch ( v ) {
  case SGCellChoice::NotApplicable: return os << "n/a";
  case SGCellChoice::Choice1: return os << "1";
  case SGCellChoice::Choice2: return os << "2";
  case SGCellChoice::Choice3: return os << "3";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGOriginChoice v )
{
  switch ( v ) {
  case SGOriginChoice::NotApplicable: return os << "n/a";
  case SGOriginChoice::Choice1: return os << "1";
  case SGOriginChoice::Choice2: return os << "2";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGAxisPermutation v )
{
  switch ( v ) {
  case SGAxisPermutation::NotApplicable: return os << "n/a";
  case SGAxisPermutation::abc: return os << "abc";
  case SGAxisPermutation::ba_mc: return os << "ba-c";
  case SGAxisPermutation::cab: return os << "cab";
  case SGAxisPermutation::mcba: return os << "-cba";
  case SGAxisPermutation::bca: return os << "bca";
  case SGAxisPermutation::a_mcb: return os << "a-cb";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, const SGSettingInfo& si )
{
  return os << "{hexFamilyAxes=" << si.hexFamilyAxes
            << ", uniqueAxis=" << si.uniqueAxis
            << ", cellChoice=" << si.cellChoice
            << ", originChoice=" << si.originChoice
            << ", axisPermutation=" << si.axisPermutation << "}";
}
