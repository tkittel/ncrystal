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

    SpaceGroup deduceSettingFromPositions( const StructureInfo& si,
                                           const AtomInfoList& atomlist )
    {
      nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
      nc_assert( !si.spacegroup_setting.has_value() );
      const auto positions = positionsBySpecies( atomlist );
      const CellParameters cell{ si.lattice_a, si.lattice_b, si.lattice_c,
                                 si.alpha, si.beta, si.gamma };
      const auto cands = settingsCompatibleWithCell( si.spacegroup, cell );
      if ( cands.empty() ) {
        //No setting accepts the cell (this throws with the details):
        InfoUtils::getSpaceGroupCandidateForHKLFromStructureInfo( si );
        nc_assert_always( false );
      }
      Optional<SpaceGroup> best;
      double bestdev = kInfinity;
      for ( auto& sg : cands ) {
        const double dev
          = NC::checkStructureSymmetry( sg, positions ).maxDeviation;
        if ( dev <= tol_error && dev < bestdev ) {//(first one on ties)
          best.emplace( sg );
          bestdev = dev;
        }
      }
      if ( best.has_value() )
        return best.value();
      //No setting fits both the cell and the positions. Tell which one of
      //them is at fault:
      const bool rhombo
        = usesRhombohedralAxes( static_cast<int>( si.spacegroup ), si.alpha );
      const auto res = NC::checkStructureSymmetry( si.spacegroup, rhombo,
                                                   positions );
      if ( !( res.maxDeviation <= tol_error ) )
        NCRYSTAL_THROW2(BadInput,"Atom positions are not consistent with"
                        " space group "<<si.spacegroup<<" in any of its"
                        " settings (largest deviation of a symmetry image"
                        " from an atom of the same type: "
                        <<fmt(res.maxDeviation,"%.2g")<<" in fractional"
                        " coordinates). Either the space group or the"
                        " positions are wrong.");
      const SpaceGroup sg( res.settingIdx.value() );
      const auto& cc = SGSymmetry::get( sg ).cellConstraints();
      CellParameters cp = cell;
      std::string reason;
      try {
        cc.complete( cp );
      } catch ( Error::BadInput& e ) {
        reason = e.what();
      }
      NCRYSTAL_THROW2(BadInput,"Atom positions have the symmetry of space"
                      " group setting "<<sg<<", but the unit cell is not"
                      " consistent with it ("<<cc<<"): "<<reason);
    }

    void checkStructureSymmetry( const StructureInfo& si,
                                 const AtomInfoList& atomlist )
    {
      nc_assert_always( si.spacegroup_setting.has_value() );
      const SpaceGroup setting( si.spacegroup_setting.value() );
      const auto res = NC::checkStructureSymmetry( setting,
                                                   positionsBySpecies( atomlist ) );
      if ( !( res.maxDeviation <= tol_error ) )
        NCRYSTAL_THROW2(BadInput,"Atom positions are not consistent with"
                        " the space group setting "<<setting
                        <<" (largest deviation of a symmetry image from an"
                        " atom of the same type: "
                        <<fmt(res.maxDeviation,"%.2g")<<" in fractional"
                        " coordinates). Either the setting or the"
                        " positions are wrong.");
      if ( !( res.maxDeviation <= tol_warn ) )
        NCRYSTAL_WARN("Atom positions deviate slightly from the symmetry"
                      " of space group "<<si.spacegroup<<" (largest"
                      " deviation of a symmetry image from an atom of the"
                      " same type: "<<fmt(res.maxDeviation,"%.2g")
                      <<" in fractional coordinates). Assuming the"
                      " symmetric structure (specify positions with more"
                      " precision to avoid this warning).");
    }

    Optional<SpaceGroupSettingIdx>
    settingWithoutPositions( const StructureInfo& si )
    {
      nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
      const Optional<SpaceGroup> declared = declaredSetting( si );
      if ( declared.has_value() )
        return declared.value().settingIdx();
      const CellParameters cell{ si.lattice_a, si.lattice_b, si.lattice_c,
                                 si.alpha, si.beta, si.gamma };
      const auto cands = settingsCompatibleWithCell( si.spacegroup, cell );
      if ( cands.size() != 1 )
        return NullOpt;
      return cands.front().settingIdx();
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

        //Tolerance for two positions to be the same:
        constexpr double pos_tol = 0.0001;//NB: Best if value matches the
                                          //    one in NCrystal/cifutils.py!

        bool samePosition( const AtomInfo::Pos& p1, const AtomInfo::Pos& p2 )
        {
          //Periodic: e.g. x=0.99999 and x=0.0 are the same.
          for ( auto k : ncrange( 3 ) ) {
            const double d = ncabs( p1[k] - p2[k] );
            if ( ncmin( d, 1.0 - d ) >= pos_tol )
              return false;
          }
          return true;
        }

        //An atom position with the element and coordinates as given:
        struct InputPos {
          AtomInfo::Pos pos;//mapped to [0,1)
          AtomInfo::Pos input;
          const AtomInfo* ai;
        };

        std::string posDescr( const AtomInfo& ai, const AtomInfo::Pos& p )
        {
          std::ostringstream ss;
          ss << ai.atomData().description(false) << " (" << fmt(p[0]) << ", "
             << fmt(p[1]) << ", " << fmt(p[2]) << ")";
          return ss.str();
        }

        void detect_duplicate_positions( std::vector<InputPos>& plist ) {
          //Throws if two positions (all in [0,1)) are the same. Sorts plist
          //by x, and then sweeps over it, comparing each position to those
          //with nearby x (also across the x=1 boundary).
          std::sort( plist.begin(), plist.end(),
                     []( const InputPos& a, const InputPos& b )
                     { return atominfo_pos_compare( a.pos, b.pos ); } );
          auto check = []( const InputPos& p1, const InputPos& p2 )
          {
            if ( samePosition( p1.pos, p2.pos ) )
              NCRYSTAL_THROW2(BadInput,"The same atom position is used twice: "
                              <<posDescr( *p1.ai, p1.input )<<" and "
                              <<posDescr( *p2.ai, p2.input ));
          };
          const std::size_t n = plist.size();
          for ( auto i : ncrange( n ) ) {
            const auto& p = vectAt( plist, i ).pos;
            for ( auto j : ncrange( i+1, n ) ) {
              if ( vectAt( plist, j ).pos[0] - p[0] >= pos_tol )
                break;
              check( vectAt( plist, i ), vectAt( plist, j ) );
            }
            //Across the x=1 boundary (x close to 1 vs. x close to 0):
            for ( auto j : ncrange( i ) ) {
              if ( vectAt( plist, j ).pos[0] + 1.0 - p[0] >= pos_tol )
                break;
              check( vectAt( plist, i ), vectAt( plist, j ) );
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
          std::vector<details::InputPos> all_positions;
          all_positions.reserve(natompos);

          for ( auto& ai : uc_atomlist ) {
            //Map all atom positions to interval [0,1) (e.g. 1.0 becomes 0.0, -0.3
            //becomes 0.7, etc.). However at most move 1.0:
            for ( auto& pos : ai.detail_accessPos() ) {
              const auto input = pos;
              for ( auto& coord : pos )
                details::atominfo_pos_remap( coord );
              all_positions.push_back( { pos, input, &ai } );
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
        if ( si.spacegroup > 0 ) {
          //(a declared setting is validated before checking the cell with it)
          declaredSetting( si );
          //Without a declared setting, atom positions determine it:
          if ( !si.spacegroup_setting.has_value() && uc.atomlist.has_value() )
            si.spacegroup_setting
              = deduceSettingFromPositions( si, uc.atomlist.value() ).settingIdx();
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

        if ( si.spacegroup > 0 ) {
          if ( uc.atomlist.has_value() )
            InfoBuilderDetail::checkStructureSymmetry( si, uc.atomlist.value() );
          else
            si.spacegroup_setting = InfoBuilderDetail::settingWithoutPositions( si );
        }
        if ( si.spacegroup_setting.has_value() )
          si.spacegroup_setting.value().validate();
      }

  }
}
