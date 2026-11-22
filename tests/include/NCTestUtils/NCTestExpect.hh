#ifndef NCrystal_TestExpect_hh
#define NCrystal_TestExpect_hh

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

#include "NCrystal/core/NCException.hh"
#include <iostream>
#include <cstring>

namespace nctest {

  //Run f, which must throw an exception of type TExc (and, if msg is given,
  //with exactly that message). Prints "Caught expected <type>(<message>)",
  //like ensure_error in the Python NCTestUtils/common.py, and throws
  //std::runtime_error if the expected exception was not thrown:
  template<class TExc, class TFct>
  void ensureError( TFct f, const char* msg = nullptr )
  {
    try {
      f();
    } catch ( TExc& e ) {
      if ( msg && std::strcmp( e.what(), msg ) != 0 )
        throw std::runtime_error( std::string( "Expected exception had"
                                               " message \"" ) + e.what()
                                  + "\" instead of \"" + msg + "\"" );
      std::cout << "Caught expected " << e.getTypeName() << "(" << e.what()
                << ")" << std::endl;
      return;
    }
    throw std::runtime_error( "Did not throw exception as required" );
  }

}

#endif
