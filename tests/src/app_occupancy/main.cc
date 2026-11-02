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


#include "NCrystal/NCrystal.hh"
#include "NCrystal/internal/infobld/NCInfoBuilder.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)
#define REQUIREFLTEQ(x,y) nc_assert_always(NC::floateq(x,y))

namespace {

  //Rock salt rebuilt from the stdlib NaCl file via InfoBuilder, with the given
  //occupancies of the Na and Cl roles (and optionally explicit dynamics with
  //the given fractions, or a uij value on the first site of the first role):
  NC::InfoPtr buildNaCl( double occ_na, double occ_cl,
                         NC::Optional<NC::PairDD> dynfracs = NC::NullOpt,
                         bool set_uij = false )
  {
    auto ref = NC::createInfo( "stdlib::NaCl_sg225_SodiumChloride.ncmat"
                               ";inelas=0" );
    NC::InfoBuilder::SinglePhaseBuilder bldr;
    bldr.dataSourceName = NC::DataSourceName( "nacl_occu" );
    bldr.temperature = ref->getTemperature();
    NC::InfoBuilder::UnitCell uc;
    uc.structinfo = ref->getStructureInfo();
    uc.atomlist.emplace();
    for ( auto& ai : ref->getAtomInfos() ) {
      const bool is_na = ai.atomData().isElement()
        && ai.atomData().elementName() == "Na";
      NC::AtomInfo::AtomPositions pos = ai.unitCellPositions();
      uc.atomlist.value().emplace_back( ai.indexedAtomData(), std::move(pos),
                                        ai.debyeTemp(), ai.msd(),
                                        is_na ? occ_na : occ_cl );
      if ( set_uij && uc.atomlist.value().size() == 1 ) {
        NC::AtomInfo::Site site;
        site.position = ai.unitCellPositions().front();
        site.multiplicity = ai.numberPerUnitCell();
        site.uij.emplace();
        site.uij.value().fill( 0.0 );
        uc.atomlist.value().back().detail_accessSites().push_back(site);
      }
    }
    if ( dynfracs.has_value() ) {
      bldr.dynamics.emplace();
      for ( auto& ai : ref->getAtomInfos() ) {
        const bool is_na = ai.atomData().elementName() == "Na";
        bldr.dynamics.value().emplace_back(
          new NC::DI_Sterile( is_na ? dynfracs.value().first
                              : dynfracs.value().second,
                              ai.indexedAtomData(), ref->getTemperature() ) );
      }
    }
    bldr.unitcell = std::move( uc );
    bldr.hklPlanes.emplace();
    bldr.hklPlanes.value().dspacingRange = NC::PairDD( 1.0, 2.0 );
    bldr.hklPlanes.value().source = NC::HKLList();
    return NC::InfoBuilder::buildInfoPtr( std::move( bldr ) );
  }

  void dump( const char * descr, const NC::Info& info )
  {
    std::cout << descr << ": n_atoms=" << info.getStructureInfo().n_atoms
              << " averageatomspercell="
              << NC::fmt( info.getAverageAtomsPerCell() )
              << " numberdensity=" << NC::fmt( info.getNumberDensity().dbl() )
              << " density=" << NC::fmt( info.getDensity().dbl() ) << '\n';
    for ( auto& ai : info.getAtomInfos() )
      std::cout << "  " << ai.atomData().elementName()
                << ": numberPerUnitCell=" << ai.numberPerUnitCell()
                << " occupancy=" << NC::fmt( ai.occupancy() ) << '\n';
    for ( auto& c : info.getComposition() )
      std::cout << "  composition " << c.atom.data().elementName()
                << ": " << NC::fmt( c.fraction ) << '\n';
  }

  template<class TFct>
  void expectBadInput( const char * descr, TFct fct )
  {
    bool ok = false;
    try {
      fct();
    } catch ( NC::Error::BadInput& e ) {
      std::cout << descr << ": BadInput: " << e.what() << '\n';
      ok = true;
    }
    REQUIRE( ok );
  }
}

int main()
{
  auto ref = NC::createInfo( "stdlib::NaCl_sg225_SodiumChloride.ncmat"
                             ";inelas=0" );
  auto full = buildNaCl( 1.0, 1.0 );
  dump( "full occupancy", *full );
  //Occupancies of 1 must give bit-identical results:
  REQUIRE( full->getNumberDensity().dbl() == ref->getNumberDensity().dbl() );
  REQUIRE( full->getDensity().dbl() == ref->getDensity().dbl() );
  REQUIRE( full->getComposition().size() == ref->getComposition().size() );
  for ( auto i : NC::ncrange( ref->getComposition().size() ) )
    REQUIRE( full->getComposition().at(i).fraction
             == ref->getComposition().at(i).fraction );

  auto partial = buildNaCl( 0.9, 1.0 );
  dump( "Na occupancy 0.9", *partial );
  //Counts of positions are unchanged, averages are weighted:
  REQUIRE( partial->getStructureInfo().n_atoms == 8 );
  const double vol = partial->getStructureInfo().volume;
  REQUIRE( NC::InfoUtils::occupancyWeightedAtomCount( partial->getAtomInfos() )
           == 7.6 );
  REQUIRE( partial->getAverageAtomsPerCell() == 7.6 );
  REQUIRE( full->getAverageAtomsPerCell() == 8.0 );
  REQUIRE( ref->getAverageAtomsPerCell() == 8.0 );
  REQUIREFLTEQ( partial->getNumberDensity().dbl(), 7.6 / vol );
  for ( auto& c : partial->getComposition() )
    REQUIREFLTEQ( c.fraction, ( c.atom.data().elementName() == "Na"
                                ? 3.6 : 4.0 ) / 7.6 );

  //Dynamics fractions must be the weighted ones:
  buildNaCl( 0.9, 1.0, NC::PairDD( 3.6/7.6, 4.0/7.6 ) );
  expectBadInput( "unweighted dyninfo fractions", []()
  { buildNaCl( 0.9, 1.0, NC::PairDD( 0.5, 0.5 ) ); } );

  //Invalid occupancies:
  for ( double o : { 0.0, -0.5, 1.0000001, 2.0,
                     std::numeric_limits<double>::quiet_NaN() } ) {
    std::cout << "occupancy " << NC::fmt( o );
    expectBadInput( "", [o]() { buildNaCl( o, 1.0 ); } );
  }

  //Anisotropic displacements are reserved:
  expectBadInput( "uij", []() { buildNaCl( 1.0, 1.0, NC::NullOpt, true ); } );

  std::cout << "All OK" << std::endl;
  return 0;
}
