#ifndef NCrystal_RotatedScatter_hh
#define NCrystal_RotatedScatter_hh

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
#include "NCrystal/internal/utils/NCRotMatrix.hh"

namespace NCRYSTAL_NAMESPACE {

  class RotatedScatter final : public ProcImpl::ScatterAnisotropicMat {
  public:
    //A scatter process seen from a rotated frame. This allows e.g. a single
    //crystal model set up in the frame of the crystal to be shared between
    //all orientations of that crystal. The rotation must take directions in
    //the frame of the wrapped process to the frame of this one.
    RotatedScatter( ProcImpl::ProcPtr, const RotMatrix& );
    const char * name() const noexcept override { return "RotatedScatter"; }
    const ProcImpl::ProcPtr& process() const noexcept { return m_proc; }
    EnergyDomain domain() const noexcept override;
    CrossSect crossSection( CachePtr&, NeutronEnergy,
                            const NeutronDirection& ) const override;
    ScatterOutcome sampleScatter( CachePtr&, RNG&, NeutronEnergy,
                                  const NeutronDirection& ) const override;
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    bool isPureElasticScatter() const override;
#endif
  protected:
    Optional<std::string> specificJSONDescription() const override;
  private:
    ProcImpl::ProcPtr m_proc;
    RotMatrix m_rot, m_invrot;
  };

}

#endif
