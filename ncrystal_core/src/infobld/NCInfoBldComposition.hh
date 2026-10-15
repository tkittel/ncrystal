#ifndef NCrystal_InfoBldComposition_hh
#define NCrystal_InfoBldComposition_hh

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

#include "NCrystal/internal/infobld/NCInfoBuilder.hh"
#include "NCrystal/internal/utils/NCMath.hh"

//Validation and completion of the composition, temperatures, densities and
//state of matter of the input to InfoBuilder::buildInfo (used by
//NCInfoBuilder.cc).

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilderDetail {

    //The validation steps (throwing BadInput for invalid input), in the order
    //they are applied by InfoBuilder:
    void validateAtomIndexes( const Info::Composition& composition );
    void validateAndCompleteComposition( Optional<Info::Composition>& composition,
                                         const Optional<InfoBuilder::UnitCell>& unitcell,
                                         Optional<DynamicInfoList>& dynamics );
    void validateTemperatures( const Optional<Temperature>& temperature,
                               const Optional<DynamicInfoList>& dynamics );
    AtomMass calculateAverageAtomMass(const Info::Composition& c);
    void validateDensities( Density d, NumberDensity nd );
    void validateAndCompleteDensities( AtomMass averageAtomMass,
                                       const Optional<InfoBuilder::UnitCell>& unitcell,
                                       Optional<Density>& density,
                                       Optional<NumberDensity>& numberDensity );
    void validateAndCompleteStateOfMatter( bool isCrystalline,
                                           const Optional<DynamicInfoList>& dynamics,
                                           Info::StateOfMatter& stateOfMatter );

    //Validate a list of (fraction,object) entries and snap the fractions to
    //sum exactly to unity:
    template<class TList, class TGetFractionFct, class TSetFractionFct, class TGetDescrObjFct>
    void validateFractionListAndSnapToUnity( const char* fieldname, TList& alist,
                                             const TGetFractionFct& getfracfct,
                                             const TSetFractionFct& setfracfct,
                                             const TGetDescrObjFct& getdescrobjfct )
    {
      if ( alist.empty() )
        NCRYSTAL_THROW2(BadInput,"invalid "<<fieldname<<" : no entries!");

      StableSum ftot;
      for (const auto& e : alist ) {
        double fraction = getfracfct(e);
        if ( fraction<=0 || fraction > 1.0 ) {
          NCRYSTAL_THROW2(BadInput,"invalid "<<fieldname<<" fraction for "
                          <<getdescrobjfct(e)<<" : "<<fraction);
        }
        ftot.add( fraction );
      }
      const double ftotval = ftot.sum();
      if ( ftotval >= 1.000000001 || ftotval < 0.999999999 )
        NCRYSTAL_THROW2(BadInput,"invalid "<<fieldname<<" : fractions do not sum to unity");
      const double snapfact = 1.0 / ftotval;
      for (auto& e : alist )
        setfracfct(e,getfracfct(e)*snapfact);
    }

  }
}

#endif
