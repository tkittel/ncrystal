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

#include "NCLegacyLatticeChecks.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include <iostream>

//Check that the legacy lattices completion and checks from NCrystal <= 4.4 and
//the next-gen completeAndCheckStructureInfoCell agree.

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace {

  struct Cell { double a, b, c, alpha, beta, gamma; };

  bool operator==( const Cell& x, const Cell& y )
  {
    return x.a == y.a && x.b == y.b && x.c == y.c
      && x.alpha == y.alpha && x.beta == y.beta && x.gamma == y.gamma;
  }

  //Returns true if accepted, in which case cl is the completed cell:
  bool legacy( unsigned sg, Cell& cl )
  {
    try {
      NC::LegacyLattice::checkAndCompleteLattice( sg, cl.a, cl.b, cl.c );
      NC::LegacyLattice::checkAndCompleteLatticeAngles( sg, cl.alpha,
                                                        cl.beta, cl.gamma );
      //(the a=b=c check for rhombohedral axes was done by the callers):
      if ( NC::LegacyLattice::usesRhombohedralAxes( sg, cl.alpha )
           && !( cl.a == cl.b && cl.a == cl.c ) )
        return false;
    } catch ( NC::Error::BadInput& ) {
      return false;
    }
    return true;
  }

  bool modern( unsigned sg, Cell& cl )
  {
    NC::StructureInfo si;
    si.spacegroup = sg;
    si.lattice_a = cl.a;
    si.lattice_b = cl.b;
    si.lattice_c = cl.c;
    si.alpha = cl.alpha;
    si.beta = cl.beta;
    si.gamma = cl.gamma;
    //(the legacy checks assumed the default setting of the number):
    si.spacegroup_setting
      = NC::SpaceGroup( std::to_string( sg )
                        + ( NC::usesRhombohedralAxes( static_cast<int>( sg ),
                                                      cl.alpha ) ? ":r" : "" ),
                        NC::SpaceGroupParseMode::Relaxed ).settingIdx();
    try {
      NC::InfoUtils::completeAndCheckStructureInfoCell( si );
    } catch ( NC::Error::BadInput& ) {
      return false;
    }
    cl = { si.lattice_a, si.lattice_b, si.lattice_c, si.alpha, si.beta,
           si.gamma };
    return true;
  }
}

int main()
{
  const double a = 4.1;
  const double bvals[] = { 4.1, 5.2, 0.0, -1.0 };
  const double cvals[] = { 4.1, 6.3, 0.0 };
  const double alphavals[] = { 90.0, 80.0, 57.0, 125.0, 0.0 };
  const double betavals[] = { 90.0, 95.0, 57.0, 0.0, 185.0 };
  const double gammavals[] = { 90.0, 120.0, 100.0, 57.0, 0.0 };
  unsigned long ncases = 0, naccepted = 0, nskipped = 0;
  for ( unsigned sg = 1; sg <= 230; ++sg ) {
    for ( double b : bvals ) {
      for ( double c : cvals ) {
        for ( double alpha : alphavals ) {
          for ( double beta : betavals ) {
            for ( double gamma : gammavals ) {
              const Cell cell{ a, b, c, alpha, beta, gamma };
              if ( c == 0.0 && NC::usesRhombohedralAxes( static_cast<int>( sg ),
                                                          alpha ) ) {
                //Known difference: the legacy code never completed c for
                //rhombohedral axes (it required all three lengths), while
                //the constraints imply c=a. Unobservable in practice (NCMAT
                //always has three lengths, .laz/.lau forbid these axes):
                ++nskipped;
                continue;
              }
              Cell cl = cell, cm = cell;
              const bool okl = legacy( sg, cl );
              const bool okm = modern( sg, cm );
              ++ncases;
              if ( okl != okm || ( okl && !( cl == cm ) ) ) {
                std::cout << "DISAGREEMENT for space group " << sg << ": cell "
                          << a << ' ' << b << ' ' << c << " / " << alpha << ' '
                          << beta << ' ' << gamma << " legacy "
                          << ( okl ? "accepts" : "rejects" )
                          << ", completeAndCheckStructureInfoCell "
                          << ( okm ? "accepts" : "rejects" ) << std::endl;
                REQUIRE( false );
              }
              if ( okl )
                ++naccepted;
            }
          }
        }
      }
    }
  }
  std::cout << "Legacy lattice checks and completeAndCheckStructureInfoCell"
            << " agree on "
            << ncases << " cases for all 230 space group numbers ("
            << naccepted << " accepted, " << nskipped
            << " rhombohedral-axes cases with c=0 skipped)" << std::endl;
  return 0;
}
