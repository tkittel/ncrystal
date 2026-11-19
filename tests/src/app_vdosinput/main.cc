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

#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)
#define REQUIREFLTEQ(x, y) nc_assert_always(::NC::floateq((x), (y)))

namespace {
  void prtgrid(NC::Span<const double> g)
  {
    std::cout<<"[";
    bool first = true;
    for (auto e : g) {
      if (first)
        first = false;
      else
        std::cout<<", ";
      std::cout<<NC::fmt(e);
    }
    std::cout<<"]";
  }

  void testCheckAndStandardiseInput()
  {
    using NC::VectD;
    using FZ = NC::VDOS::VDOSEval::FinalZero;
    auto t = []( VectD eg, VectD d, FZ fz )
    {
      std::cout<<"Test checkAndStandardiseInput(";
      prtgrid(eg);
      std::cout<<", ";
      prtgrid(d);
      std::cout<<", "<<(fz==FZ::Require?"Require":"Append")<<")"<<std::endl;
      std::cout<<"   ->   result: ";
      try {
        NC::VDOS::VDOSEval::checkAndStandardiseInput( eg, d, fz );
      } catch ( NC::Error::BadInput& e ) {
        std::cout<<"BadInput: "<<e.what()<<std::endl;
        return;
      }
      prtgrid(eg);
      std::cout<<", ";
      prtgrid(d);
      std::cout<<std::endl;
      REQUIRE( NC::VDOS::VDOSEval::isStandardInput( eg, d ) );
    };
    for ( auto fz : { FZ::Require, FZ::Append } ) {
      t( {0.01,0.05}, {1,2,3,2,1,0}, fz );
      t( {0.01,0.07}, {1,2,3,2,1,0,0,0}, fz );
      t( {0.01,0.02,0.03,0.04,0.05,0.06,0.07},
         {1,2,3,2,1,0,0}, fz );
      t( {0.01,0.05}, {1,2,3,2,1}, fz );
      t( {0.01,0.02,0.04,0.08,0.16}, {1,2,3,2,1}, fz );
      t( {0.01,0.05}, {1,2,3,0,0}, fz );
      t( {0.01,0.05}, {1,2,0,3,0}, fz );
      t( {0.01,0.05}, {0,0,0,0,0}, fz );
      t( {0.01,0.02,0.03}, {1,2,3,2,0}, fz );
    }
    //Dropping excess zeros from a two-value grid keeps the bin width:
    VectD eg{0.001,0.009}, d{1,2,3,4,5,6,0,0,0};
    NC::VDOS::VDOSEval::checkAndStandardiseInput( eg, d, FZ::Require );
    REQUIRE( d.size() == 7 );
    REQUIREFLTEQ( eg.back(), 0.007 );
    REQUIRE( !NC::VDOS::VDOSEval::isStandardInput( {0.01,0.05},
                                                   {1,2,3,0,0,0} ) );
    REQUIRE( !NC::VDOS::VDOSEval::isStandardInput( {0.01,0.05},
                                                   {1,2,3,2,1} ) );
  }
}

int main()
{
  testCheckAndStandardiseInput();
  return 0;
}
