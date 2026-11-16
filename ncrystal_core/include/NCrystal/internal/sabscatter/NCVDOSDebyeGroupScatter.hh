#ifndef NCrystal_VDOSDebyeGroupScatter_hh
#define NCrystal_VDOSDebyeGroupScatter_hh

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

#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/dyninfoutils/NCVDOSDebyeGroup.hh"

namespace NCRYSTAL_NAMESPACE {

  class VDOSDebyeGroupScatter final : public ProcImpl::ScatterIsotropicMat {
  public:

    //Scattering from a group of VDOSDebye entries with the same mass and
    //temperature (see NCVDOSDebyeGroup.hh). The SigmaBound is the bound cross
    //section per atom of the group (i.e. the summed fraction*boundXS of its
    //entries divided by their summed fraction, which should then be the scale
    //of this process in a ProcComposition):

    const char * name() const noexcept override
    {
      return "VDOSDebyeGroupScatter";
    }

    VDOSDebyeGroupScatter( shared_obj<const DynInfoUtils::VDOSDebyeGroupData>,
                           SigmaBound );

    CrossSect crossSectionIsotropic( CachePtr&, NeutronEnergy ) const override;

    ScatterOutcomeIsotropic
    sampleScatterIsotropic( CachePtr&, RNG&, NeutronEnergy ) const override;

    Optional<IsoXSTable> createIsoXSTable() const override;

  protected:
    Optional<std::string> specificJSONDescription() const override;

  private:
    shared_obj<const DynInfoUtils::VDOSDebyeGroupData> m_data;
    SigmaBound m_sb;
  };

}

#endif
