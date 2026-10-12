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
#include "NCrystal/internal/utils/NCLatticeUtils.hh"

namespace NC = NCrystal;

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

        checkAndCompleteLattice( si.spacegroup, si.lattice_a,
                                 si.lattice_b,  si.lattice_c );

        checkAndCompleteLatticeAngles( si.spacegroup, si.alpha, si.beta, si.gamma );
        if ( usesRhombohedralAxes( static_cast<int>(si.spacegroup), si.alpha )
             && !( si.lattice_a == si.lattice_b
                   && si.lattice_a == si.lattice_c ) )
          NCRYSTAL_THROW2(BadInput,"Spacegroup ("<<si.spacegroup<<") in"
                          " rhombohedral axes requires a=b=c");


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
      }

  }
}
