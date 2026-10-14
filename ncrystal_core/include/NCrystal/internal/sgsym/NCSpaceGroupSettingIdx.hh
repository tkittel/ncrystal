#ifndef NCrystal_SpaceGroupSettingIdx_hh
#define NCrystal_SpaceGroupSettingIdx_hh

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

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //Identifies one of the 530 space group settings tabulated in:
  //
  //    Shmueli, U., Hall, S. R. & Grosse-Kunstleve, R. W. (2010).
  //    "Space-group symbols for numeric and symbolic computations."
  //    International Tables for Crystallography, Volume B: Reciprocal
  //    Space, 3rd ed., Appendix A1.4.2, pp. 122-134 (Table A1.4.2.7).
  //
  //That table has 530 rows, and the setting is thus efficiently encoded as the
  //row number, i.e. an integral value from 1 to 530 (both included). We shall
  //refer to this value as the "setting index" (notice that this is instead
  //called the "Hall number" in spglib).

  class SpaceGroupSettingIdx final
    : public EncapsulatedValue<SpaceGroupSettingIdx,std::uint16_t> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
    static constexpr std::uint16_t min_value = 1;
    static constexpr std::uint16_t max_value = 530;
    void validate() const;//throws BadInput if not in 1..530
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline void SpaceGroupSettingIdx::validate() const
  {
    if ( ! ( m_value >= min_value && m_value <= max_value ) )
      NCRYSTAL_THROW2(BadInput,"Invalid space group setting index: "<<m_value
                      <<" (must be in range "<<min_value<<".."<<max_value<<")");
  }
}

#endif
