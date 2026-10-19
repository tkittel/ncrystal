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
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace InfoBuilderDetail {

    namespace {
      //The setting declared by the data source (spacegroup_setting already
      //set), validated against the number and the axes. NullOpt if none:
      Optional<SpaceGroup> declaredSetting( const StructureInfo& si )
      {
        if ( !si.spacegroup_setting.has_value() )
          return NullOpt;
        si.spacegroup_setting.value().validate();
        const SpaceGroup setting( si.spacegroup_setting.value() );
        if ( setting.number() != si.spacegroup )
          NCRYSTAL_THROW2(BadInput,"Declared space group setting "<<setting
                          <<" is not a setting of space group number "
                          <<si.spacegroup);
        const bool rhombo
          = usesRhombohedralAxes( static_cast<int>( si.spacegroup ),
                                  si.alpha );
        if ( ( StrView( setting.choice() ) == "r" ) != rhombo )
          NCRYSTAL_THROW2(BadInput,"Declared space group setting "<<setting
                          <<" does not match the unit cell axes (hexagonal"
                          " versus rhombohedral)");
        return setting;
      }
    }

    namespace {
      //Tolerances: rounding noise in data files gives deviations well below
      //tol_warn (all stdlib files are O(1e-16), and positions given with 4
      //decimals stay below 1e-4), sites rounded to 3 decimals exceed it
      //slightly, while structures beyond tol_error (comparable to the
      //distinct-site threshold of expandSiteToOrbit) simply do not have the
      //symmetry:
      constexpr double tol_warn = 5e-4;
      constexpr double tol_error = 1e-2;

      std::vector<std::vector<Vector>>
      positionsBySpecies( const AtomInfoList& atomlist )
      {
        std::vector<std::vector<Vector>> positions;
        positions.reserve( atomlist.size() );
        for ( auto& ai : atomlist ) {
          positions.emplace_back();
          for ( auto& p : ai.unitCellPositions() )
            positions.back().push_back( p.as<Vector>() );
        }
        return positions;
      }
    }

    Optional<SpaceGroup>
    deduceSettingFromPositions( const StructureInfo& si,
                                const AtomInfoList& atomlist )
    {
      nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
      nc_assert( !si.spacegroup_setting.has_value() );
      const auto positions = positionsBySpecies( atomlist );
      const bool rhombo
        = usesRhombohedralAxes( static_cast<int>( si.spacegroup ), si.alpha );
      Optional<SpaceGroup> best;
      double bestdev = kInfinity;
      for ( auto& sg : SpaceGroup::settingsOfNumber( si.spacegroup ) ) {
        if ( ( StrView( sg.choice() ) == "r" ) != rhombo )
          continue;
        CellParameters cp{ si.lattice_a, si.lattice_b, si.lattice_c,
                           si.alpha, si.beta, si.gamma };
        try {
          SGSymmetry::get( sg ).cellConstraints().complete( cp );
        } catch ( Error::BadInput& ) {
          continue;
        }
        const double dev
          = NC::checkStructureSymmetry( sg, positions ).maxDeviation;
        if ( dev <= tol_error && dev < bestdev ) {//(first one on ties)
          best.emplace( sg );
          bestdev = dev;
        }
      }
      return best;
    }

    //Verify that the atom positions actually have the symmetry of the
    //declared space group (in the stored setting, else in some of its
    //settings), rather than just assuming it. Slight deviations
    //(rounded coordinates in data files) give a warning and the symmetric
    //structure is assumed, while larger ones are an error:
    SpaceGroupSettingIdx checkStructureSymmetry( const StructureInfo& si,
                                                 const AtomInfoList& atomlist )
    {
      nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
      const auto positions = positionsBySpecies( atomlist );
      const int sg = static_cast<int>( si.spacegroup );
      const bool rhombo = usesRhombohedralAxes( sg, si.alpha );
      //A setting declared by the data source is checked as such, otherwise
      //the settings of the number are searched:
      const Optional<SpaceGroup> declared = declaredSetting( si );
      const auto res = ( declared.has_value()
                         ? NC::checkStructureSymmetry( declared.value(),
                                                       positions )
                         : NC::checkStructureSymmetry( si.spacegroup, rhombo,
                                                       positions ) );
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
      //(the setting found is used as is: all hkl-related code follows the
      //setting stored on the StructureInfo)
      return res.settingIdx.value();
    }

    Optional<SpaceGroupSettingIdx>
    settingWithoutPositions( const StructureInfo& si )
    {
      nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
      const int sg = static_cast<int>( si.spacegroup );
      const Optional<SpaceGroup> declared = declaredSetting( si );
      if ( declared.has_value() )
        return declared.value().settingIdx();
      std::string spec = std::to_string( sg );
      if ( usesRhombohedralAxes( sg, si.alpha ) )
        spec += ":r";
      const SpaceGroup setting( spec, SpaceGroupParseMode::Relaxed );
      //(the two settings of an R number are told apart by the cell angles):
      if ( setting.hasMultipleSettings() && !isRhombohedralSpaceGroup( sg ) )
        return NullOpt;
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

      void validateAndCompleteUnitCell( InfoBuilder::UnitCell& uc )
      {
        auto& si = uc.structinfo;
        if ( si.spacegroup_setting.has_value() && si.spacegroup == 0 ) {
          //The number is implied by a declared setting:
          si.spacegroup_setting.value().validate();
          si.spacegroup = SpaceGroup( si.spacegroup_setting.value() ).number();
        }

        if ( uc.atomlist.has_value() ) {
          auto& uc_atomlist = uc.atomlist.value();
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

        if ( si.spacegroup > 230 )
          NCRYSTAL_THROW2(BadInput,"invalid spacegroup number ("<<si.spacegroup<<")");
        //Without a declared setting, the positions may single out a setting
        //whose cell constraints differ from those of the default one:
        if ( si.spacegroup > 0 ) {
          //(a declared setting is validated before checking the cell with it)
          declaredSetting( si );
          if ( !si.spacegroup_setting.has_value() && uc.atomlist.has_value() ) {
            const auto deduced
              = deduceSettingFromPositions( si, uc.atomlist.value() );
            if ( deduced.has_value() )
              si.spacegroup_setting = deduced.value().settingIdx();
          }
          InfoUtils::completeAndCheckStructureInfoCell( si );
        }
        if ( !( si.lattice_a > 0.0 && si.lattice_b > 0.0 && si.lattice_c > 0.0 ) )
          NCRYSTAL_THROW(BadInput,"lattice parameters must be positive numbers");


        if ( ! ( si.alpha > 0 && si.alpha < 180 &&
                 si.beta  > 0 && si.beta  < 180 &&
                 si.gamma > 0 && si.gamma < 180 ) ) {
          NCRYSTAL_THROW(BadInput,"Lattice angles must all be >0 and <180 degrees.");
        }

        if ( ! ( si.alpha > kPi || si.beta > kPi || si.gamma > kPi ) )
          NCRYSTAL_THROW(BadInput,"Lattice angles specified in radians where degrees were expected.");

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
        if ( si.spacegroup_setting.has_value() )
          si.spacegroup_setting.value().validate();
      }

  }
}
