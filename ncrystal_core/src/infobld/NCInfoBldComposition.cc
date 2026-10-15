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

#include "NCInfoBldComposition.hh"
#include "NCInfoBldStructure.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilderDetail {

      void validateAtomIndexes( const Info::Composition& composition )
      {
        std::set<AtomIndex> seen;
        for ( auto& c : composition ) {
          if ( !seen.insert(c.atom.index).second )
            NCRYSTAL_THROW2(BadInput,"Invalid AtomIndex setup (repeated indices found in composition list)");
          if ( c.atom.index.get() >= composition.size() )
            NCRYSTAL_THROW2(BadInput,"Invalid AtomIndex setup (must be one of 0,...,ncomponents-1)");
        }
      }
      void validateAndCompleteComposition( Optional<Info::Composition>& composition,
                                           const Optional<InfoBuilder::UnitCell>& unitcell,
                                           Optional<DynamicInfoList>& dynamics )
      {
        const bool has_uc_atomlist = unitcell.has_value() && unitcell.value().atomlist.has_value();
        if ( composition.has_value() ) {
          //Composition is set, so just verify that dynamics/unitcell are not set, validate the provided fractions, and return.
          if ( has_uc_atomlist || dynamics.has_value() )
            NCRYSTAL_THROW(BadInput,"Do not set explicit composition on SinglePhaseBuilder when providing unitcell.atomlist or dynamics.");
          validateFractionListAndSnapToUnity( "composition list", composition.value(),
                                              [](const Info::CompositionEntry& e){ return e.fraction; },
                                              [](Info::CompositionEntry& e, double f){ e.fraction = f; },
                                              [](const Info::CompositionEntry& e){ return e.atom; } );
          return;
        }

        //Composition is absent, must be able to get from either unit cell or dynamics.
        if ( !has_uc_atomlist && !dynamics.has_value() ) {
          NCRYSTAL_THROW(BadInput,"SinglePhaseBuilder must have at least one of the following"
                         " pieces of information: composition, atomlist in unit cell, or dynamics.");
        }

        //Create composition based on either unit cell or dynamics (their
        //fractions have already been aligned but we anyway check the more
        //precise unitcell source first):

        Info::Composition new_composition;
        if ( has_uc_atomlist ) {
          //From unit cell:
          auto& atomlist = unitcell.value().atomlist.value();
          const auto ntot = totalNumberOfAtomsInUnitCell( atomlist );
          nc_assert_always( ntot > 0 );
          new_composition.reserve( atomlist.size() );
          for ( const auto& ai : atomlist ) {
            nc_assert( ai.numberPerUnitCell() > 0 );
            new_composition.emplace_back( ai.numberPerUnitCell() / double(ntot), ai.atom() );
          }
        } else {
          //From dynamics
          new_composition.reserve( dynamics.value().size() );
          for ( auto& di : dynamics.value() )
            new_composition.emplace_back( di->fraction(), di->atom() );
        }
        composition = std::move(new_composition);
      }
      void validateTemperatures( const Optional<Temperature>& temperature,
                                 const Optional<DynamicInfoList>& dynamics )
      {
        Temperature tempval{-1.0};
        if ( temperature.has_value() ) {
          tempval = temperature.value();
          if ( ! ( tempval.get() >= Temperature::allowed_range.first
                   && tempval.get() <= Temperature::allowed_range.second ) )
            NCRYSTAL_THROW2(BadInput,"Invalid or out-of-range temperature value provided: "<<tempval);
        }
        if ( !dynamics.has_value() )
          return;
        if ( !temperature.has_value() )
          NCRYSTAL_THROW(BadInput,"Temperature is required whenever dynamic info is supplied");
        nc_assert(tempval.get()>0.0);
        tempval.validate();
        //
        //Check consistency of temperature fields (checking consistency of
        //DI_VDOSDebye Debye temperature fields elsewhere):
        for (auto& di : dynamics.value()) {

          if ( di->temperature() != tempval )
            NCRYSTAL_THROW(BadInput,"temperature info on DynamicInfo object"
                           " is different than temperature value of containing phase!");

          auto di_vdos = dynamic_cast<const DI_VDOS*>(di.get());
          if ( di_vdos && di_vdos->vdosData().temperature() != di_vdos->temperature() )
            NCRYSTAL_THROW(BadInput,"temperature info on VDOSData object provided by DI_VDOS object"
                           " is different than temperature on DI_VDOS object itself!");

          //NB: for DI_ScatKnlDirect we can't validate the temperature field on
          //the SABData object here, since it is built on-demand. However, we do
          //it upon construction in the ensureBuildThenReturnSAB() method.
        }
      }
      AtomMass calculateAverageAtomMass(const Info::Composition& c)
      {
        StableSum sum;
        for ( auto& di : c )
          sum.add( di.fraction * di.atom.data().averageMassAMU().dbl() );
        return AtomMass{ sum.sum() };
      }
      void validateDensities( Density d, NumberDensity nd )
      {
        if ( ! (d.dbl() >= 0.0 ) || std::isinf(d.dbl()) )
          NCRYSTAL_THROW2(BadInput,"Invalid density value: "<<d);
        if ( ! (nd.dbl() >= 0.0 ) || std::isinf(nd.dbl()) )
          NCRYSTAL_THROW2(BadInput,"Invalid number density value: "<<nd);
        if ( d.dbl() == 0.0 || nd.dbl() == 0.0 )
          NCRYSTAL_THROW(BadInput,"Materials with vanishing densities are not presently supported.");
        //General type validation:
        d.validate();
        nd.validate();
      }
      void validateAndCompleteDensities( AtomMass averageAtomMass,
                                         const Optional<InfoBuilder::UnitCell>& unitcell,
                                         Optional<Density>& density,
                                         Optional<NumberDensity>& numberDensity )
      {
        averageAtomMass.validate();
        nc_assert_always( averageAtomMass.get() > 0.0 );

        if ( density.has_value() && numberDensity.has_value() )
          NCRYSTAL_THROW(BadInput,"Do not supply both Density and NumberDensity on"
                         " SinglePhaseBuilder (supply at most one and the other will be calculated).");

        if ( unitcell.has_value() ) {
          //Always take final density/numberdensity values from unit cell when
          //possible - but check compatibility with any number already provided.
          const auto& si = unitcell.value().structinfo;
          nc_assert( si.volume > 0.0 && si.n_atoms > 0 );
          const auto numberDensityUC = NumberDensity{ si.n_atoms / si.volume };
          const auto densityUC = Density( numberDensityUC, averageAtomMass );
          if ( numberDensity.has_value() && !floateq( numberDensityUC.dbl(), numberDensity.value().dbl(), 1e-2 ) )
            NCRYSTAL_THROW2(BadInput,"Provided (" << numberDensity.value() << ") versus calculated-from-unit-cell"
                            " (" << numberDensityUC << ") number density values are incompatible!");
          if ( density.has_value() && !floateq( densityUC.dbl(), density.value().dbl(), 1e-2 ) )
            NCRYSTAL_THROW2(BadInput,"Provided (" << density.value() << ") versus calculated-from-unit-cell"
                            " (" << densityUC << ") density values are incompatible!");
          numberDensity = numberDensityUC;
          density = densityUC;
        } else {
          //If one of density/numberdensity is set, calculate the other:
          if ( density.has_value() && !numberDensity.has_value() )
            numberDensity = NumberDensity( density.value(), averageAtomMass );

          if ( numberDensity.has_value() && !density.has_value() )
            density = Density( numberDensity.value(), averageAtomMass );

          //Final requirement is always that a density must be available:
          if ( !numberDensity.has_value() || !density.has_value() )
            NCRYSTAL_THROW(BadInput,"Density/NumberDensity values must always be supplied directly or it"
                           " must be possible to deduce them (from each other or unit cell information).");
        }

        //Final sanity check:
        nc_assert_always( density.has_value() && numberDensity.has_value() );
        validateDensities( density.value(), numberDensity.value() );
      }
      void validateAndCompleteStateOfMatter( bool isCrystalline,
                                             const Optional<DynamicInfoList>& dynamics,
                                             Info::StateOfMatter& stateOfMatter )
      {
        //Verify and possibly fill out StateOfMatter:
        //
        //  Crystalline fields or presence of VDOS/VDOSDebye dynamics implies
        //  solid. In that case, verify that current state is either Unknown or
        //  Solid (and replace Unknown with Solid).
        //
        bool should_be_solid(isCrystalline);
        if ( dynamics.has_value() ) {
          for ( auto& di : dynamics.value() ) {
            if (should_be_solid)
              break;
            if (dynamic_cast<const DI_VDOSDebye*>(di.get())||dynamic_cast<const DI_VDOS*>(di.get()))
              should_be_solid = true;
          }
        }
        if ( should_be_solid ) {
          if ( !isOneOf( stateOfMatter, Info::StateOfMatter::Unknown, Info::StateOfMatter::Solid ) )
            NCRYSTAL_THROW2(BadInput,"Info objects that are crystalline or have at least one VDOS"
                            " (or VDOSDebye) can not be designated as \""<<Info::toString(stateOfMatter)<<"\"");
          stateOfMatter = Info::StateOfMatter::Solid;
        }
      }

  }
}
