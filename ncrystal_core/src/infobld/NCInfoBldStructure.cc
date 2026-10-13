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

#include "NCInfoBldStructure.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/extd_utils/NCStructureSetting.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilderDetail {

    //Verify that the atom positions actually have the symmetry of the
    //declared space group (in some of its settings), rather than just
    //assuming it. Slight deviations (rounded coordinates in data files)
    //give a warning and the symmetric structure is assumed, while larger
    //ones are an error, as are settings whose symmetry-equivalent
    //reflections differ from those NCrystal assumes (cf. EqRefl):
    SpaceGroupSettingIdx checkStructureSymmetry( const StructureInfo& si,
                                                 const AtomInfoList& atomlist )
    {
      nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
      std::vector<std::vector<Vector>> positions;
      positions.reserve( atomlist.size() );
      for ( auto& ai : atomlist ) {
        positions.emplace_back();
        for ( auto& p : ai.unitCellPositions() )
          positions.back().push_back( p.as<Vector>() );
      }
      const int sg = static_cast<int>( si.spacegroup );
      const bool rhombo = usesRhombohedralAxes( sg, si.alpha );
      //A setting declared by the data source is checked as such, otherwise
      //the settings of the number are searched:
      const Optional<SpaceGroup> declared = declaredStructureSetting( si );
      const auto res = ( declared.has_value()
                         ? NC::checkStructureSymmetry( declared.value(),
                                                       positions )
                         : NC::checkStructureSymmetry( si.spacegroup, rhombo,
                                                       positions ) );
      //Tolerances: rounding noise in data files gives deviations well
      //below tol_warn (all stdlib files are below 3e-16, and positions
      //given with 4 decimals stay below 1e-4), sites rounded to 3 decimals
      //exceed it slightly, while structures beyond tol_error (comparable to
      //the distinct-site threshold of expandSiteToOrbit) simply do not have
      //the declared symmetry:
      constexpr double tol_warn = 5e-4;
      constexpr double tol_error = 1e-2;
      if ( !( res.maxDeviation <= tol_error ) ) {
        if ( declared.has_value() )
          NCRYSTAL_THROW2(BadInput,"Atom positions are not consistent with"
                          " the space group setting "<<declared.value()
                          <<" (largest deviation of a symmetry image from an"
                          " atom of the same type: "
                          <<fmt(res.maxDeviation,"%.2g")<<" in fractional"
                          " coordinates). Either the setting or the"
                          " positions are wrong.");
        NCRYSTAL_THROW2(BadInput,"Atom positions are not consistent with"
                        " space group "<<si.spacegroup<<" in any of its"
                        " settings (largest deviation of a symmetry image"
                        " from an atom of the same type: "
                        <<fmt(res.maxDeviation,"%.2g")<<" in fractional"
                        " coordinates). Either the space group or the"
                        " positions are wrong.");
      }
      if ( !( res.maxDeviation <= tol_warn ) )
        NCRYSTAL_WARN("Atom positions deviate slightly from the symmetry"
                      " of space group "<<si.spacegroup<<" (largest"
                      " deviation of a symmetry image from an atom of the"
                      " same type: "<<fmt(res.maxDeviation,"%.2g")
                      <<" in fractional coordinates). Assuming the"
                      " symmetric structure (specify positions with more"
                      " precision to avoid this warning).");
      const SpaceGroup setting( SpaceGroupSettingIdx{ res.settingIdx } );
      if ( declared.has_value() )
        return setting.settingIdx();
      std::string spec = std::to_string( si.spacegroup );
      if ( rhombo )
        spec += ":r";
      const SpaceGroup assumed( spec, SpaceGroupParseMode::Relaxed );
      if ( setting != assumed
           && !haveSameRotationParts( SGSymmetry::get( setting ),
                                      SGSymmetry::get( assumed ) ) )
        NCRYSTAL_THROW2(BadInput,"Atom positions are only consistent with"
                        " space group "<<si.spacegroup<<" in the setting "
                        <<setting<<", whose symmetry-equivalent reflections"
                        " differ from those of the setting "<<assumed
                        <<" assumed by NCrystal. Transform the structure to"
                        " that setting (e.g. with"
                        " NCMATComposer.refine_crystal_structure()).");
      return setting.settingIdx();
    }

    SpaceGroupSettingIdx settingWithoutPositions( const StructureInfo& si )
    {
      nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
      const int sg = static_cast<int>( si.spacegroup );
      const Optional<SpaceGroup> declared = declaredStructureSetting( si );
      if ( declared.has_value() )
        return declared.value().settingIdx();
      std::string spec = std::to_string( sg );
      if ( usesRhombohedralAxes( sg, si.alpha ) )
        spec += ":r";
      const SpaceGroup setting( spec, SpaceGroupParseMode::Relaxed );
      //(the two settings of an R number are told apart by the cell angles):
      if ( setting.hasMultipleSettings() && !isRhombohedralSpaceGroup( sg ) )
        return SpaceGroupSettingIdx{ 0 };
      return setting.settingIdx();
    }

  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilderDetail {
    namespace {
      namespace details {
        bool atominfo_pos_compare( const NC::AtomInfo::Pos& rh, const NC::AtomInfo::Pos& lh )
        {
          if (rh[0]!=lh[0]) return rh[0] < lh[0];
          if (rh[1]!=lh[1]) return rh[1] < lh[1];
          return rh[2] < lh[2];
        }

        void atominfo_pos_remap( double& x ) {
          const double xorig = x;
          if ( x<0.0 )
            x += 1.0;
          else if ( x>=1.0 )
            x -= 1.0;
          if ( ! (x>=0 && x<1.0) )//must be in [0,1) and not NaN
            NCRYSTAL_THROW2(BadInput,"Invalid coordinate of atom position encountered (out of range or NaN): "<<xorig);
          if (x==0.0)
            x=0.0;//remaps -0 to 0
        }

        void detect_duplicate_positions( AtomInfo::AtomPositions& plist ) {
          //Throws if two positions (all in [0,1)) are closer than a tolerance
          //in all three coordinates, taking periodicity into account (so
          //e.g. x=0.99999 and x=0.0 are close). Sorts plist by x, and then
          //sweeps over it, comparing each position to those with nearby x
          //(also across the x=1 boundary).
          const double tol = 0.0001;//NB: Best if value matches the one in NCrystal/cifutils.py!
          std::sort( plist.begin(), plist.end(), atominfo_pos_compare );
          auto check = [tol]( const AtomInfo::Pos& p1, const AtomInfo::Pos& p2 )
          {
            for ( auto k : ncrange( 3 ) ) {
              const double d = ncabs( p1[k] - p2[k] );
              if ( ncmin( d, 1.0 - d ) >= tol )
                return;
            }
            NCRYSTAL_THROW2(BadInput,"The same atom position used more than"
                            " once: ("<<p1[0]<<", "<<p1[1]<<", "<<p1[2]<<")"
                            " and ("<<p2[0]<<", "<<p2[1]<<", "<<p2[2]<<")");
          };
          const std::size_t n = plist.size();
          for ( auto i : ncrange( n ) ) {
            const auto& p = vectAt( plist, i );
            for ( auto j : ncrange( i+1, n ) ) {
              if ( vectAt( plist, j )[0] - p[0] >= tol )
                break;
              check( p, vectAt( plist, j ) );
            }
            //Across the x=1 boundary (x close to 1 vs. x close to 0):
            for ( auto j : ncrange( i ) ) {
              if ( vectAt( plist, j )[0] + 1.0 - p[0] >= tol )
                break;
              check( p, vectAt( plist, j ) );
            }
          }
        }
      }
    }

      std::size_t totalNumberOfAtomsInUnitCell( const AtomInfoList& atomlist )
      {
        //Derives number of atoms in a unit cell from atom list. Also
        //validatesthat atomlist nonempty, that all ai.numberPerUnitCell() > 0,
        //and returns the sum of those.
        std::size_t ntot(0);
        if ( atomlist.empty() )
          NCRYSTAL_THROW2(BadInput,"AtomInfoList must be non-empty if provided");
        for ( const auto& ai : atomlist ) {
          if ( ai.numberPerUnitCell() == 0 )
            NCRYSTAL_THROW(BadInput,"AtomInfo object should not have numberPerUnitCell()==0");
          nc_assert( ai.numberPerUnitCell() > 0 );
          ntot += ai.numberPerUnitCell();
        }
        nc_assert_always( ntot > 0 );
        return ntot;
      }

      //Fill in the unit cell positions of AtomInfo objects given as sites
      //(one position per orbit of the space group setting):
      void expandSitesToPositions( const StructureInfo& si, AtomInfoList& atomlist )
      {
        if ( si.spacegroup == 0 )
          NCRYSTAL_THROW(BadInput,"Atom positions given as sites require a"
                         " space group");
        const SGSymmetry& sym = SGSymmetry::get( structureSetting( si ) );
        for ( auto& ai : atomlist ) {
          if ( !ai.unitCellPositions().empty() )
            continue;
          auto& sites = ai.detail_accessSites();
          auto& positions = ai.detail_accessPos();
          nc_assert( !sites.empty() );
          for ( auto& site : sites ) {
            const auto orbit = expandSiteToOrbit( sym, site.position.as<Vector>() );
            site.position = orbit.positions.front().as<AtomInfo::Pos>();//(symmetrised)
            site.multiplicity = static_cast<unsigned>( orbit.positions.size() );
            for ( auto& p : orbit.positions )
              positions.push_back( p.as<AtomInfo::Pos>() );
          }
        }
      }

      //Fill in the sites of AtomInfo objects given as positions, by grouping
      //the (sorted) positions into orbits; the representative of an orbit is
      //its first position. Left empty if the grouping is not possible (e.g.
      //positions accepted with a warning only by the structure check):
      void groupPositionsIntoSites( const StructureInfo& si, AtomInfoList& atomlist )
      {
        if ( !si.spacegroup_setting.get() )
          return;
        const SGSymmetry& sym = SGSymmetry::get( SpaceGroup( si.spacegroup_setting ) );
        for ( auto& ai : atomlist ) {
          if ( !ai.sites().empty() )
            continue;
          std::vector<Vector> positions;
          positions.reserve( ai.numberPerUnitCell() );
          for ( auto& p : ai.unitCellPositions() )
            positions.push_back( p.as<Vector>() );
          std::vector<std::vector<unsigned>> orbits;
          try {
            orbits = groupPositionsIntoOrbits( sym, positions );
          } catch ( Error::BadInput& ) {
            continue;
          }
          auto& sites = ai.detail_accessSites();
          sites.reserve( orbits.size() );
          for ( auto& orbit : orbits ) {
            AtomInfo::Site site;
            site.position = vectAt( ai.unitCellPositions(), orbit.front() );
            site.multiplicity = static_cast<unsigned>( orbit.size() );
            sites.push_back( std::move( site ) );
          }
        }
      }

      void validateAndCompleteUnitCell( InfoBuilder::UnitCell& uc )
      {
        auto& si = uc.structinfo;

        if ( si.spacegroup > 230 )
          NCRYSTAL_THROW2(BadInput,"invalid spacegroup number ("<<si.spacegroup<<")");
        if ( si.spacegroup > 0 )
          completeAndCheckCell( si );//(cell constraints of the setting)
        if ( !( si.lattice_a > 0.0 && si.lattice_b > 0.0 && si.lattice_c > 0.0 ) )
          NCRYSTAL_THROW(BadInput,"lattice parameters must be positive numbers");

        if ( ! ( si.alpha > 0 && si.alpha < 180 &&
                 si.beta  > 0 && si.beta  < 180 &&
                 si.gamma > 0 && si.gamma < 180 ) ) {
          NCRYSTAL_THROW(BadInput,"Lattice angles must all be >0 and <180 degrees.");
        }

        if ( ! ( si.alpha > kPi || si.beta > kPi || si.gamma > kPi ) )
          NCRYSTAL_THROW(BadInput,"Lattice angles specified in radians where degrees were expected.");

        if ( uc.atomlist.has_value() ) {
          auto& uc_atomlist = uc.atomlist.value();
          bool anysites = false;
          for ( auto& ai : uc_atomlist )
            if ( ai.unitCellPositions().empty() )
              anysites = true;
          if ( anysites ) {
            expandSitesToPositions( si, uc_atomlist );
            if ( si.n_atoms == 0 )
              si.n_atoms = static_cast<unsigned>( totalNumberOfAtomsInUnitCell( uc_atomlist ) );
          }
          //First validate total number of atoms is >0 and consistently specified:
          auto natompos = totalNumberOfAtomsInUnitCell( uc_atomlist );//throws if empty or have entries with 0 atoms
          nc_assert(natompos>0);
          if ( si.n_atoms != natompos )
            NCRYSTAL_THROW2(BadInput,"Inconsistent total number of atoms deduced from StructureInfo and AtomInfo ("
                            <<si.n_atoms<<" vs. "<<natompos<<")");

          //Sort for reproducibility:
          std::stable_sort(uc_atomlist.begin(),uc_atomlist.end(),
                           []( const NC::AtomInfo& rh, const NC::AtomInfo& lh )
                           { return rh.atomData() < lh.atomData(); } );

          //Check AtomInfo objects for consistency and standardise them a bit:
          AtomInfo::AtomPositions all_positions;
          all_positions.reserve(natompos);

          for ( auto& ai : uc_atomlist ) {
            //Map all atom positions to interval [0,1) (e.g. 1.0 becomes 0.0, -0.3
            //becomes 0.7, etc.). However at most move 1.0:
            for ( auto& pos : ai.detail_accessPos() ) {
              for ( auto& coord : pos )
                details::atominfo_pos_remap( coord );
              all_positions.push_back( pos );
            }

            //Sort positions, ensuring a well defined ordering:
            std::stable_sort( ai.detail_accessPos().begin(),
                              ai.detail_accessPos().end(),
                              details::atominfo_pos_compare );

            nc_assert( ai.numberPerUnitCell() > 0 );//was validated in totalNumberOfAtomsInUnitCell(..)

            //Debye temps and msd's should be specified for all or none of the AtomInfo objects:
            if ( ai.debyeTemp().has_value() != uc_atomlist.front().debyeTemp().has_value() )
              NCRYSTAL_THROW(BadInput,"Inconsistency: Debye temperatures specified for some but not all AtomInfo objects.");

            if ( ai.debyeTemp().has_value() )
              ai.debyeTemp().value().validate();

            if ( ai.msd().has_value() != uc_atomlist.front().msd().has_value() )
              NCRYSTAL_THROW(BadInput,"Inconsistency: MSD values specified for some but not all AtomInfo objects.");

          }

          //Ensure only one atom exists at a given position, within a tolerance:
          details::detect_duplicate_positions(all_positions);
        }

        if ( ! ( si.n_atoms > 0 ) )
          NCRYSTAL_THROW2(BadInput,"StructureInfo n_atoms not a positive number: "<<si.n_atoms);

        //Calculate unit cell volume:
        auto latRot = getLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                     si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
        double expected_volume = latRot.colX().cross( latRot.colY() ).dot( latRot.colZ() );
        if ( si.volume > 0.0 && !floateq(si.volume,expected_volume,1e-3) )
          NCRYSTAL_THROW2(BadInput,"Provided ("<<si.volume<<"Aa3) versus calculated ("
                          <<expected_volume<<"Aa3) unit cell volume are incompatible!");
        si.volume = expected_volume;
        nc_assert( si.volume > 0.0 );

        if ( si.spacegroup > 0 )
          si.spacegroup_setting
            = ( uc.atomlist.has_value()
                ? InfoBuilderDetail::checkStructureSymmetry( si, uc.atomlist.value() )
                : InfoBuilderDetail::settingWithoutPositions( si ) );
        if ( uc.atomlist.has_value() )
          groupPositionsIntoSites( si, uc.atomlist.value() );
      }

  }
}
