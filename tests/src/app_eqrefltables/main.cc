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

#include "NCEqReflTables.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include <iostream>

//The symmetry-equivalent reflections of EqRefl (derived from the sgsym
//symmetry operations) must be identical to those of the hand-written tables
//of NCrystal 4.4 and earlier (EqReflTables, moved here from the library), for
//all space group numbers and both axes choices of the rhombohedral groups.

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace {
  template<class T1, class T2>
  void expectBad( int sg, bool rhombo )
  {
    for ( int which = 0; which < 2; ++which ) {
      try {
        if ( which == 0 ) {
          T1 e( sg, rhombo );
        } else {
          T2 e( sg, rhombo );
        }
      } catch ( NC::Error::BadInput& ) {
        continue;
      }
      REQUIRE( false );
    }
  }
}

int main()
{
  constexpr int N = 6;
  unsigned long npoints = 0, nsettings = 0;
  for ( int sg = 1; sg <= 230; ++sg ) {
    for ( int irh = 0; irh < 2; ++irh ) {
      const bool rhombo = ( irh == 1 );
      if ( rhombo && !NC::isRhombohedralSpaceGroup( sg ) )
        continue;
      const NC::EqRefl er( sg, rhombo );
      const NC::EqReflTables et( sg, rhombo );
      ++nsettings;
      for ( int h = -N; h <= N; ++h ) {
        for ( int k = -N; k <= N; ++k ) {
          for ( int l = -N; l <= N; ++l ) {
            if ( !h && !k && !l )
              continue;
            const auto a = er.getEquivalentReflections( h, k, l );
            const auto b = et.getEquivalentReflections( h, k, l );
            REQUIRE( a.size() == b.size() );
            auto itb = b.begin();
            for ( auto& e : a )
              REQUIRE( e == *itb++ );
            REQUIRE( er.getEquivalentReflectionsRepresentativeValue( h, k, l )
                     == et.getEquivalentReflectionsRepresentativeValue( h, k, l ) );
            REQUIRE( a.contains( NC::HKL( h, k, l ) )
                     && b.contains( NC::HKL( -h, -k, -l ) ) );
            ++npoints;
          }
        }
      }
    }
  }
  std::cout << "EqRefl and EqReflTables agree for " << nsettings
            << " space group settings and " << npoints << " hkl points"
            << " (|h|,|k|,|l| <= " << N << ")" << std::endl;
  expectBad<NC::EqRefl,NC::EqReflTables>( 0, false );
  expectBad<NC::EqRefl,NC::EqReflTables>( 231, false );
  expectBad<NC::EqRefl,NC::EqReflTables>( 225, true );//not rhombohedral
  std::cout << "Invalid arguments rejected by both" << std::endl;
  return 0;
}
