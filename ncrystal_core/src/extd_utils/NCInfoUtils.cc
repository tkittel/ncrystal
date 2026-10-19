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

#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NC = NCrystal;

const NC::DynamicInfo*
NC::InfoUtils::findDynInfo( const Info& info,
                            Optional<std::string> atomDisplayLabel )
{
  const DynamicInfo* res = nullptr;
  for ( auto& di : info.getDynamicInfoList() ) {
    if ( atomDisplayLabel.has_value()
         && ( info.displayLabel(di->atom().index)
              != atomDisplayLabel.value() ) ) {
      continue;
    }
    if (res) {
      NCRYSTAL_THROW(BadInput,"Could not uniquely determine DynInfo object."
                     " Perhaps providing an explicit atom label will help.");
    } else {
      nc_assert_always( di.get() != nullptr );
      res = di.get();
    }
  }
  if (!res)
    NCRYSTAL_THROW(BadInput,"Could not find requested DynInfo object.");
  return res;
}

NC::SpaceGroup
NC::InfoUtils::getSpaceGroupFromStructureInfo( const StructureInfo& si )
{
  nc_assert_always( si.spacegroup >= 1 && si.spacegroup <= 230 );
  if ( si.spacegroup_setting.has_value() )
    return SpaceGroup( si.spacegroup_setting.value() );
  const CellParameters cell{ si.lattice_a, si.lattice_b, si.lattice_c,
                             si.alpha, si.beta, si.gamma };
  const auto cands = settingsCompatibleWithCell( si.spacegroup, cell );
  if ( cands.size() == 1 )
    return cands.front();
  if ( cands.empty() )
    getSpaceGroupCandidateForHKLFromStructureInfo( si );//throws (details)
  std::ostringstream msg;
  msg << "The space group setting is unknown (the unit cell can be given in"
         " several settings of space group " << si.spacegroup << ": ";
  for ( auto& c : cands )
    msg << ( &c == &cands.front() ? "" : ", " ) << c;
  msg << ")";
  NCRYSTAL_THROW(BadInput,msg.str());
}

NC::SpaceGroup
NC::InfoUtils::getSpaceGroupCandidateForHKLFromStructureInfo( const StructureInfo& si )
{
  nc_assert_always( si.spacegroup >= 1 && si.spacegroup <= 230 );
  if ( si.spacegroup_setting.has_value() )
    return SpaceGroup( si.spacegroup_setting.value() );
  const CellParameters cell{ si.lattice_a, si.lattice_b, si.lattice_c,
                             si.alpha, si.beta, si.gamma };
  const auto cands = settingsCompatibleWithCell( si.spacegroup, cell );
  if ( cands.empty() ) {
    //Report the problem with the first setting with the axes of the cell:
    const bool rhombo
      = usesRhombohedralAxes( static_cast<int>( si.spacegroup ), si.alpha );
    unsigned nsettings = 0;
    std::string reason;
    for ( auto& sg : SpaceGroup::settingsOfNumber( si.spacegroup ) ) {
      if ( ( StrView( sg.choice() ) == "r" ) != rhombo )
        continue;
      if ( nsettings++ )
        continue;
      const auto& cc = SGSymmetry::get( sg ).cellConstraints();
      CellParameters cp = cell;
      try {
        cc.complete( cp );
      } catch ( Error::BadInput& e ) {
        std::ostringstream ss;
        ss << sg << " (" << cc << "): " << e.what();
        reason = ss.str();
      }
    }
    if ( nsettings == 1 )
      NCRYSTAL_THROW2(BadInput,"Unit cell not consistent with space group "
                      << reason);
    NCRYSTAL_THROW2(BadInput,"Unit cell not consistent with any setting of"
                    " space group " << si.spacegroup << " (e.g. " << reason
                    << ")");
  }
  const SGSymmetry& sym0 = SGSymmetry::get( cands.front() );
  for ( auto& s : cands ) {
    if ( haveSameRotationParts( SGSymmetry::get( s ), sym0 ) )
      continue;
    std::ostringstream msg;
    msg << "The space group setting is unknown, and the unit cell can be"
           " given in settings of space group " << si.spacegroup
        << " with different symmetry-equivalent reflections (";
    for ( auto& c : cands )
      msg << ( &c == &cands.front() ? "" : ", " ) << c;
    msg << ")";
    NCRYSTAL_THROW(BadInput,msg.str());
  }
  return cands.front();
}

void NC::InfoUtils::completeAndCheckStructureInfoCell( StructureInfo& si )
{
  nc_assert( si.spacegroup >= 1 && si.spacegroup <= 230 );
  const SpaceGroup setting
    = getSpaceGroupCandidateForHKLFromStructureInfo( si );
  const auto& cc = SGSymmetry::get( setting ).cellConstraints();
  CellParameters cp{ si.lattice_a, si.lattice_b, si.lattice_c,
                     si.alpha, si.beta, si.gamma };
  try {
    cc.complete( cp );
  } catch ( Error::BadInput& e ) {
    NCRYSTAL_THROW2(BadInput,"Unit cell not consistent with space group "
                    << setting << " (" << cc << "): " << e.what());
  }
  si.lattice_a = cp.a;
  si.lattice_b = cp.b;
  si.lattice_c = cp.c;
  si.alpha = cp.alpha;
  si.beta = cp.beta;
  si.gamma = cp.gamma;
}
