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
#include "NCrystal/internal/sgsym/NCSymOp.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include <iostream>

//Tests of the non-throwing variants of sgsym functions.

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace test_tryfcts {

//The try.. variants (SpaceGroup::tryParse, SymOp::tryParse and
//SGCellConstraints::tryComplete) give the results and messages of the
//throwing versions:
template<class TFct>
std::string thrownMsg( TFct f )
{
  try {
    f();
  } catch ( NC::Error::BadInput& e ) {
    return e.what();
  }
  return {};
}

void run()
{
  for ( auto s : { "70:1", "225", "70", "70:3", "70:A", "bla", "231" } ) {
    const auto r = NC::SpaceGroup::tryParse( s,
                                             NC::SpaceGroupParseMode::Strict );
    const auto t = thrownMsg( [s](){
      NC::SpaceGroup( s, NC::SpaceGroupParseMode::Strict ); } );
    REQUIRE( r.setting.has_value() == r.errorMsg.empty() );
    if ( r.setting.has_value() ) {
      REQUIRE( t.empty() );
      std::cout << "SpaceGroup::tryParse(\"" << s << "\"): "
                << r.setting.value() << std::endl;
    } else {
      REQUIRE( t == r.errorMsg );
      std::cout << "SpaceGroup::tryParse(\"" << s << "\"): "
                << r.errorMsg << std::endl;
    }
  }
  for ( auto s : { "-x,y+1/2,-z", "x,y", "2x,y,z", "x+y,y,z", "x,y,z+1/3" } ) {
    const auto r = NC::SymOp::tryParse( s );
    const auto t = thrownMsg( [s](){ NC::SymOp{ s }; } );
    REQUIRE( r.op.has_value() == r.errorMsg.empty() );
    if ( r.op.has_value() ) {
      REQUIRE( t.empty() );
      std::cout << "SymOp::tryParse(\"" << s << "\"): " << r.op.value()
                << std::endl;
    } else {
      REQUIRE( t == r.errorMsg );
      std::cout << "SymOp::tryParse(\"" << s << "\"): " << r.errorMsg
                << std::endl;
    }
  }
  const auto& cc = NC::SGSymmetry::get( NC::SpaceGroup( "225",
                                         NC::SpaceGroupParseMode::Strict ) )
    .cellConstraints();
  for ( auto c : { NC::CellParameters{ 4.0, 0, 0, 0, 0, 0 },
                   NC::CellParameters{ 4.0, 4.1, 0, 0, 0, 0 } } ) {
    const auto r = cc.tryComplete( c );
    auto c2 = c;
    const auto t = thrownMsg( [&cc,&c2](){ cc.complete( c2 ); } );
    REQUIRE( r.cell.has_value() != r.errorMsg.has_value() );
    if ( r.cell.has_value() ) {
      REQUIRE( t.empty() && r.cell.value().b == c2.b );
      std::cout << "tryComplete: b=" << r.cell.value().b << std::endl;
    } else {
      REQUIRE( t == r.errorMsg.value() );
      std::cout << "tryComplete: " << r.errorMsg.value() << std::endl;
    }
  }
}

}

int main()
{
  test_tryfcts::run();
  return 0;
}
