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

#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include <iostream>
#include <map>

//Basic checks of EqRefl for all space group numbers and both axes choices of
//the rhombohedral groups (the legacy tables of NCrystal 4.4 will be moved here
//next, as a reference to compare against).

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

int main()
{
  std::map<unsigned,unsigned> nrot_count;//nRotations -> number of settings
  for ( int sg = 1; sg <= 230; ++sg ) {
    for ( int irh = 0; irh < 2; ++irh ) {
      const bool rhombo = ( irh == 1 );
      if ( rhombo && !NC::isRhombohedralSpaceGroup( sg ) )
        continue;
      const NC::EqRefl er( sg, rhombo );
      REQUIRE( er.nRotations() >= 1 && er.nRotations() <= 24 );
      nrot_count[ er.nRotations() ] += 1;
      //(h,k,l) is always in its own family, as is (-h,-k,-l):
      REQUIRE( er.getEquivalentReflections( 1, 2, 3 ).contains( NC::HKL( -1, -2, -3 ) ) );
    }
  }
  for ( auto& e : nrot_count )
    std::cout << "  " << e.second << " settings with " << e.first
              << " rotation(s) modulo sign" << std::endl;
  for ( auto bad : { std::make_pair( 0, false ), std::make_pair( 231, false ),
                     std::make_pair( 225, true ) } ) {
    bool caught = false;
    try {
      NC::EqRefl e( bad.first, bad.second );
    } catch ( NC::Error::BadInput& ) {
      caught = true;
    }
    REQUIRE( caught );
  }
  std::cout << "Invalid arguments rejected" << std::endl;
  return 0;
}
