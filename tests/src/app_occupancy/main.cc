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
#include "NCrystal/internal/extd_utils/NCFillHKL.hh"
#include "NCrystal/internal/elincscatter/NCElIncScatter.hh"
#include "NCrystal/internal/powderbragg/NCPowderBragg.hh"
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
    uc.structinfo.n_atoms_average = 0.0;//(recalculated from the atoms)
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
        site.multiplicity = ai.numberOfPositions();
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
    bldr.hklPlanes.value().dspacingRange = NC::PairDD( 0.6, 1e9 );
    bldr.hklPlanes.value().source = NC::InfoBuilder::SinglePhaseBuilder
      ::HKLPlanes::HKLListGenFct( []( const NC::StructureInfo* si,
                                      const NC::AtomInfoList* al,
                                      NC::PairDD )
      {
        nc_assert_always( si && al );
        NC::FillHKLCfg cfg;
        cfg.dcutoff = 0.6;
        return NC::calculateHKLPlanes( *si, *al, cfg );
      } );
    return NC::InfoBuilder::buildInfoPtr( std::move( bldr ) );
  }

  void dump( const char * descr, const NC::Info& info )
  {
    const auto& si = info.getStructureInfo();
    std::cout << descr << ": n_positions=" << si.n_positions
              << " n_atoms_average=" << NC::fmt( si.n_atoms_average )
              << " numberdensity=" << NC::fmt( info.getNumberDensity().dbl() )
              << " density=" << NC::fmt( info.getDensity().dbl() ) << '\n';
    for ( auto& ai : info.getAtomInfos() )
      std::cout << "  " << ai.atomData().elementName()
                << ": numberOfPositions=" << ai.numberOfPositions()
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
  REQUIRE( partial->getStructureInfo().n_positions == 8 );
  const double vol = partial->getStructureInfo().volume;
  REQUIRE( partial->getStructureInfo().n_atoms_average == 7.6 );
  REQUIRE( full->getStructureInfo().n_atoms_average == 8.0 );
  REQUIRE( ref->getStructureInfo().n_atoms_average == 8.0 );
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

  //Structure factors: |F| of every rock salt reflection is linear in the Na
  //occupancy (|F| = |o*A_Na +- A_Cl|, without sign change since b_Cl > b_Na),
  //so |F| at o=0.8 is the average of those at o=1 and o=0.6:
  {
    NC::FillHKLCfg cfg;
    cfg.dcutoff = 0.6;
    auto planes = [&cfg]( double o )
    {
      auto info = buildNaCl( o, 1.0 );
      return NC::calculateHKLPlanes( info->getStructureInfo(),
                                     info->getAtomInfos(), cfg );
    };
    const auto h10 = planes( 1.0 );
    const auto h08 = planes( 0.8 );
    const auto h06 = planes( 0.6 );
    REQUIRE( h10.size() == h08.size() && h10.size() == h06.size() );
    REQUIRE( h10.size() > 10 );
    for ( auto i : NC::ncrange( h10.size() ) ) {
      REQUIRE( h10[i].dspacing == h08[i].dspacing );
      REQUIRE( h10[i].dspacing == h06[i].dspacing );
      REQUIRE( h10[i].multiplicity == h08[i].multiplicity );
      REQUIRE( NC::floateq( 2.0 * std::sqrt( h08[i].fsquared ),
                            std::sqrt( h10[i].fsquared )
                            + std::sqrt( h06[i].fsquared ), 1e-12, 0.0 ) );
    }
    std::cout << "F^2 of " << h10.size()
              << " families linear in |F| versus occupancy: OK\n";
  }

  //Incoherent elastic at very low energies (Debye-Waller factors 1): per unit
  //cell the sum of npos*(o*sigma_inc+o*(1-o)*sigma_coh) over the roles:
  for ( double o : { 1.0, 0.9, 0.5, 0.1 } ) {
    auto info = buildNaCl( o, 1.0 );
    NC::ElIncScatter sc( *info );
    NC::CachePtr cp;
    const double xs_cell
      = sc.crossSectionIsotropic( cp, NC::NeutronEnergy( 1e-12 ) ).dbl()
      * info->getStructureInfo().n_atoms_average;
    double expected = 0.0;
    for ( auto& ai : info->getAtomInfos() ) {
      const double oo = ai.occupancy();
      expected += ai.numberOfPositions()
        * ( oo * ai.atomData().incoherentXS().dbl()
            + oo * ( 1.0 - oo ) * ai.atomData().coherentXS().dbl() );
    }
    std::cout << "Na occupancy " << NC::fmt( o )
              << ": incoherent elastic xs per cell at E->0 = "
              << NC::fmt( xs_cell, "%.6g" ) << " barn\n";
    REQUIRE( NC::floateq( xs_cell, expected, 1e-9, 0.0 ) );
  }

  //Powder Bragg cross sections are per actual atom: times the average number
  //of atoms they must equal lambda^2/(2V)*sum(d*mult*F^2) over planes with
  //2d>lambda:
  for ( double o : { 1.0, 0.7 } ) {
    auto info = buildNaCl( o, 1.0 );
    NC::PowderBragg pb( *info );
    const double natoms = info->getStructureInfo().n_atoms_average;
    const double cellvol = info->getStructureInfo().volume;
    for ( double wl : { 1.0, 2.5, 4.0 } ) {
      double sum = 0.0;
      for ( auto& h : info->hklList() )
        if ( 2.0 * h.dspacing > wl )
          sum += h.dspacing * h.multiplicity * h.fsquared;
      const double expected = wl * wl / ( 2.0 * cellvol ) * sum;
      NC::CachePtr cp;
      const double xs_cell = natoms * pb.crossSectionIsotropic(
        cp, NC::NeutronWavelength( wl ) ).dbl();
      REQUIRE( NC::floateq( xs_cell, expected, 1e-9, 0.0 ) );
    }
    std::cout << "Na occupancy " << NC::fmt( o )
              << ": powder Bragg normalisation OK\n";
  }

  //Anisotropic displacements are reserved:
  expectBadInput( "uij", []() { buildNaCl( 1.0, 1.0, NC::NullOpt, true ); } );

  std::cout << "All OK" << std::endl;
  return 0;
}
