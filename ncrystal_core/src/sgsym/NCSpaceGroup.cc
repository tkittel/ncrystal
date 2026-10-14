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
