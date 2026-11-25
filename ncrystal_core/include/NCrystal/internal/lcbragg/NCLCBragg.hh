#ifndef NCrystal_LCBragg_hh
#define NCrystal_LCBragg_hh

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
#include "NCrystal/interfaces/NCSCOrientation.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"

namespace NCRYSTAL_NAMESPACE {

  class PlaneProvider;
  class Info;

  class LCBragg final : public ProcImpl::ScatterAnisotropicMat {
  public:

    //This class models the Bragg diffraction in a layered crystal such as
    //pyrolytic graphite. The lcaxis is the axis (in the crystal frame) around
    //which the crystallites are randomly rotated. The mode parameter can be
    //used to modify which model is used internally to implement the layered
    //crystal:
    //
    //     mode=0: LCHelper
    //     mode>0: LCBraggRef(nsample=mode)
    //     mode<0: LCBraggRndmRef(nsample=-mode)
    //
    //For a description of the prec and ntrunc parameters, see NCGaussMos.hh.
    LCBragg( const Info&,
             const SCOrientation&,
             MosaicityFWHM,
             const LCAxis& lcaxis,
             int mode = 0,
             double delta_d = 0,
             PlaneProvider * plane_provider = 0,
             double prec=1e-3,
             double ntrunc=0.0 );

    //The model itself is set up in the frame of the crystal, and can be
    //shared between LCBragg instances for crystals differing only in
    //orientation. Parameters are as for the constructor above:
    class Model;
    using ModelPtr = shared_obj<const Model>;
    static ModelPtr createModel( const Info&,
                                 MosaicityFWHM,
                                 const LCAxis& lcaxis,
                                 int mode = 0,
                                 double delta_d = 0,
                                 PlaneProvider * plane_provider = 0,
                                 double prec=1e-3,
                                 double ntrunc=0.0 );

    //Construct from a model and the rotation from the crystal to the lab frame
    //(cf. getCrystal2LabRot):
    LCBragg( ModelPtr, const RotMatrix& crystal2lab );
    const ModelPtr& model() const { return m_model; }
    const RotMatrix& crystal2lab() const { return m_rot; }

    const char * name() const noexcept override { return "LCBragg"; }

    virtual ~LCBragg();

    //There is a maximum wavelength at which Bragg diffraction is possible, so
    //lower bound will reflect this (upper bound is infinity):
    EnergyDomain domain() const noexcept override;

    CrossSect crossSection(CachePtr&, NeutronEnergy, const NeutronDirection& ) const override;
    ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    bool isPureElasticScatter() const override { return true; }
#endif

  protected:
    Optional<std::string> specificJSONDescription() const override;
  private:
    ModelPtr m_model;
    RotMatrix m_rot, m_invrot;//crystal to lab and vice versa
  };
}

#endif
