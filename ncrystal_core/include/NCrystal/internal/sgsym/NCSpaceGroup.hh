#ifndef NCrystal_SpaceGroup_hh
#define NCrystal_SpaceGroup_hh

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

#include "NCrystal/internal/sgsym/NCSpaceGroupSettingIdx.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  //Immutable class representing one of the 530 space group settings of ITVB
  //2010 Table A1.4.2.7 (see NCSGSymTable.hh), i.e. a space group number
  //(1..230) along with a setting (e.g. an origin choice, a monoclinic unique
  //axis and cell choice, an orthorhombic axis permutation, or hexagonal/
  //rhombohedral axes).
  //
  //String representations are "<number>:<choice>", using the choice codes of
  //the table (e.g. "227:2", "166:r", "14:c1", "62:cab", "62:abc"), or just
  //"<number>" for the 140 numbers with a single setting (whose choice code is
  //empty). This is the only place where such strings are parsed, in one of
  //two modes: Strict accepts exactly these representations (so the setting
  //of a number with several settings must be spelled out), while Relaxed
  //also accepts a bare number for such numbers, selecting their default
  //(i.e. first listed) setting. Choice codes are matched exactly (lower case,
  //as in the table) in both modes.
  //
  //Technically the SpaceGroup class is an internal class wrapping the
  //SpaceGroupSettingIdx from the public API, and which exposes the SGSymTable
  //information through direct methods.

  enum class SpaceGroupParseMode { Relaxed, Strict };

  class SpaceGroup final {
  public:
    explicit SpaceGroup( SpaceGroupSettingIdx );
    SpaceGroup( StrView, SpaceGroupParseMode );
    SpaceGroup( const std::string& s, SpaceGroupParseMode m )
      : SpaceGroup( StrView( s ), m ) {}
    SpaceGroup( const char * s, SpaceGroupParseMode m )
      : SpaceGroup( StrView( s ), m ) {}

    //All settings of a number (default first), and a listing of them for
    //messages, like "227:1 (default), 227:2":
    static std::vector<SpaceGroup> settingsOfNumber( unsigned number );
    static std::string listSettings( unsigned number );

    SpaceGroupSettingIdx settingIdx() const noexcept { return m_idx; }
    unsigned number() const ncnoexceptndebug;//1..230
    const char * choice() const ncnoexceptndebug;//e.g. "2", "r", "c1", "abc",
                                                 //or ""
    const char * hmSymbol() const ncnoexceptndebug;//Hermann-Mauguin symbol in
                                                   //ITVB's computer-entry
                                                   //format, e.g. "F d -3 m:2",
                                                   //"P 1 21/c 1"
    const char * hallSymbol() const ncnoexceptndebug;//e.g. "-f 4vw 2vw 3" (as
                                                     //in ITVB, lattice letter
                                                     //in lower case)
    bool isDefaultSetting() const ncnoexceptndebug;//first listed setting of
                                                   //number?
    bool hasMultipleSettings() const ncnoexceptndebug;//number has other settings?
    std::string toString() const;//e.g. "227:2", "62:cab", "62:abc", or "225"

    //Note: Symmetry operations etc. are available via the SGSymmetry class.

    bool operator==( const SpaceGroup& o ) const noexcept;
    bool operator!=( const SpaceGroup& o ) const noexcept;
    bool operator<( const SpaceGroup& o ) const noexcept;
  private:
    SpaceGroupSettingIdx m_idx;
  };

  //streams as .toString():
  std::ostream& operator<<( std::ostream&, const SpaceGroup& );

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline bool SpaceGroup::operator==( const SpaceGroup& o ) const noexcept
  {
    return m_idx == o.m_idx;
  }
  inline bool SpaceGroup::operator!=( const SpaceGroup& o ) const noexcept
  {
    return m_idx != o.m_idx;
  }
  inline bool SpaceGroup::operator<( const SpaceGroup& o ) const noexcept
  {
    return m_idx < o.m_idx;
  }
}

#endif
