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

#include "NCrystal/internal/lcbragg/NCLCBragg.hh"
#include "NCrystal/internal/extd_utils/NCLCUtils.hh"
#include "NCrystal/internal/extd_utils/NCLCRefModels.hh"
#include "NCrystal/internal/scbragg/NCSCBragg.hh"
#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

  class LCBragg::Model final : private ::NC::MoveOnly {
  public:

    //Everything here is in the frame of the crystal.

    Model( LCAxis lcaxis, int mode, const Info& cinfo,
           PlaneProvider * plane_provider, MosaicityFWHM mosaicity,
           double delta_d, double prec, double ntrunc )
      : m_ekin_low(-1),
        m_mos(mosaicity),
        m_lcaxis(lcaxis.as<Vector>().unit()),
        m_mode(mode)
    {
      if (!cinfo.hasStructureInfo())
        NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks structure information.");
      const StructureInfo& si = cinfo.getStructureInfo();
      const LCAxis lcaxis_unit = lcaxis.as<Vector>().unit().as<LCAxis>();

      if (mode==0) {
        nc_assert_always(delta_d==0);//mode=0 does not currently support delta_d!=0

        std::unique_ptr<PlaneProvider> stdpp;
        if (!plane_provider) {
          stdpp = createStdPlaneProvider(&cinfo);
          plane_provider = stdpp.get();
        }

        const double v0natoms = si.volume * si.n_atoms_average;
        m_lchelper = NC::ncmake_unique<LCHelper>( lcaxis_unit,
                                                  lcaxis_unit,
                                                  mosaicity,
                                                  v0natoms,
                                                  plane_provider,
                                                  prec, ntrunc );

        m_ekin_low = wl2ekin( m_lchelper->braggThreshold() );

      } else {
        static const double identity[9] = { 1., 0., 0., 0., 1., 0., 0., 0., 1. };
        auto scbragg = makeSO<SCBragg>( SCBragg::createModel( cinfo, mosaicity,
                                                              delta_d,
                                                              plane_provider,
                                                              prec, ntrunc ),
                                        RotMatrix( identity ) );
        if (mode>0) {
          m_scmodel = std::make_shared<LCBraggRef>(scbragg, lcaxis_unit, mode);
        } else {
          int nsample = -mode;
          nc_assert(nsample>0);
          m_scmodel = std::make_shared<LCBraggRndmRot>(scbragg, lcaxis_unit, nsample);
        }
        m_ekin_low = m_scmodel->domain().elow.get();
        nc_assert(ncisinf(m_scmodel->domain().ehigh.get()));
      }

      nc_assert(m_ekin_low>0);
      nc_assert_always(bool(m_lchelper)!=bool(m_scmodel!=nullptr));
    }

    double m_ekin_low;
    std::unique_ptr<LCHelper> m_lchelper;
    ProcImpl::OptionalProcPtr m_scmodel;
    //For descriptions:
    MosaicityFWHM m_mos;
    Vector m_lcaxis;
    int m_mode;
    UniqueID m_uid;
  };

}

NC::LCBragg::ModelPtr NC::LCBragg::createModel( const Info& ci,
                                                MosaicityFWHM mosaicity,
                                                const LCAxis& lcaxis,
                                                int mode, double delta_d,
                                                PlaneProvider * plane_provider,
                                                double prec, double ntrunc )
{
  return makeSO<const Model>( lcaxis, mode, ci, plane_provider, mosaicity,
                              delta_d, prec, ntrunc );
}

NC::LCBragg::LCBragg( ModelPtr model, const RotMatrix& crystal2lab )
  : m_model( std::move(model) ),
    m_rot( Matrix( MatrixAllowCopy, crystal2lab ) ),
    m_invrot( ~crystal2lab )
{
  nc_assert_always( ncabs( crystal2lab.determinant() - 1.0 ) < 1e-6 );
}

NC::LCBragg::LCBragg( const Info& ci, const SCOrientation& sco, MosaicityFWHM mosaicity,
                      const LCAxis& lcaxis, int mode, double delta_d, PlaneProvider * plane_provider,
                      double prec, double ntrunc)
  : LCBragg( createModel( ci, mosaicity, lcaxis, mode, delta_d, plane_provider,
                          prec, ntrunc ),
             getCrystal2LabRot( sco, getReciprocalLatticeRot(
                                  ci.getStructureInfo() ) ) )
{
}

NC::LCBragg::~LCBragg() = default;

NC::EnergyDomain NC::LCBragg::domain() const noexcept
{
  return { NeutronEnergy{m_model->m_ekin_low}, NeutronEnergy{kInfinity} };
}

NC::CrossSect NC::LCBragg::crossSection(NC::CachePtr& cp, NC::NeutronEnergy ekin, const NC::NeutronDirection& indir_lab ) const
{
  if ( ekin.get() < m_model->m_ekin_low )
    return CrossSect{ 0.0 };

  const Vector indir = m_invrot * indir_lab.as<Vector>();
  if (! m_model->m_scmodel ) {
    NeutronWavelength wl{ekin};
    if (!(wl.get()>0.0))
      return CrossSect{ 0.0 };
    return CrossSect{ m_model->m_lchelper->crossSection( accessCache<LCHelper::Cache>(cp), wl.get(), indir.unit() ) };
  } else {
    return CrossSect{ m_model->m_scmodel->crossSection( cp, ekin, indir.as<NeutronDirection>() ) };
  }
}

NC::ScatterOutcome NC::LCBragg::sampleScatter(NC::CachePtr& cp, NC::RNG& rng, NC::NeutronEnergy ekin, const NC::NeutronDirection& indir_lab ) const
{
  if ( ekin.get() < m_model->m_ekin_low )
    return { ekin, indir_lab };

  const Vector indir = m_invrot * indir_lab.as<Vector>();
  if (! m_model->m_scmodel ) {
    NeutronWavelength wl{ekin};
    if (!(wl.get()>0.0))
      return { ekin, indir_lab };
    const Vector indirv = indir.unit();
    Vector outdir;
    m_model->m_lchelper->genScatter( accessCache<LCHelper::Cache>(cp), rng, wl.get(), indirv, outdir );
    if ( outdir == indirv )
      return { ekin, indir_lab };//unscattered, preserve exactly
    return { ekin, ( m_rot * outdir ).as<NeutronDirection>() };
  } else {
    const auto d = indir.as<NeutronDirection>();
    auto res = m_model->m_scmodel->sampleScatter( cp, rng, ekin, d );
    if ( res.direction == d )
      res.direction = indir_lab;//unscattered, preserve exactly
    else
      res.direction = ( m_rot * res.direction.as<Vector>() ).as<NeutronDirection>();
    return res;
  }
}

NC::Optional<std::string> NC::LCBragg::specificJSONDescription() const
{
  const Model& m = m_model;
  //Plane sets (only for the standard mode, lcmode=0):
  Optional<std::size_t> nplanesets;
  double dmin(-1),dmax(-1);
  if ( m.m_lchelper ) {
    const auto& ps = m.m_lchelper->planeSets();
    nplanesets = ps.size();
    if ( !ps.empty() ) {
      dmin = 0.5*ps.back().twodsp;
      dmax = 0.5*ps.front().twodsp;
    }
  }
  const Vector& a = m.m_lcaxis;
  std::ostringstream ss;
  {
    std::ostringstream tmp;
    if ( nplanesets.has_value() )
      tmp << "nplanesets="<<nplanesets.value()
          <<";dmin="<<dmin
          <<"Aa;dmax="<<dmax<<"Aa;";
    tmp << "mos="<<m.m_mos
        << ";lcaxis=("<<a.x()<<','<<a.y()<<','<<a.z()<<')';
    if ( m.m_mode )
      tmp << ";lcmode="<<m.m_mode;
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "nplanesets", nplanesets );
  streamJSONDictEntry( ss, "dmin", dmin );
  streamJSONDictEntry( ss, "dmax", dmax );
  streamJSONDictEntry( ss, "mos", m.m_mos.dbl() );
  streamJSONDictEntry( ss, "lcaxis", VectD{ a.x(), a.y(), a.z() } );
  streamJSONDictEntry( ss, "lcmode", m.m_mode );
  streamJSONDictEntry( ss, "model_uid", m.m_uid.getUniqueID().value );
  streamJSONDictEntry( ss, "crystal2lab", rotMatrixRows( m_rot ),
                       JSONDictPos::LAST );
  return ss.str();
}
