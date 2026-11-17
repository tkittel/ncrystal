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



#include "NCrystal/interfaces/NCSABData.hh"
#include <iostream>

namespace NC = NCrystal;

namespace {

  //Tries to construct a VDOSData, returns false if rejected with BadInput:
  bool accepted( NC::VectD dens, double emax )
  {
    try {
      NC::VDOSData( NC::PairDD( 0.1, emax ), std::move(dens),
                    NC::Temperature{ 300.0 }, NC::SigmaBound{ 1.0 },
                    NC::AtomMass{ 12.0 } );
    } catch ( NC::Error::BadInput& ) {
      return false;
    }
    return true;
  }

}

int main()
{
  //A VDOS density must end with exactly one zero:
  std::cout << "VDOSData accepts [1 2 0]: " << accepted({1.,2.,0.},0.3)
            << std::endl;
  std::cout << "VDOSData accepts [1 2]: " << accepted({1.,2.},0.2)
            << std::endl;
  std::cout << "VDOSData accepts [1 2 0 0]: " << accepted({1.,2.,0.,0.},0.4)
            << std::endl;
  return 0;
}
