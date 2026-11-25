
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

#include "NCrystal/internal/scbragg/NCSCBragg.hh"
#include "NCrystal/internal/phys_utils/NCGaussMos.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/interfaces/NCSCOrientation.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"
#include <functional>//std::greater
namespace NC=NCrystal;

class NC::SCBragg::Model final : private ::NC::MoveOnly {
public:

  //Everything here is in the frame of the crystal.

  typedef std::map<std::pair<uint64_t,uint64_t>,std::vector<Vector>,
                   std::greater<std::pair<uint64_t,uint64_t> > > SCBraggSortMap;

  Model( const NC::Info&, MosaicityFWHM, double dd,
         PlaneProvider * plane_provider, double prec, double ntrunc );

  double setupFamilies( const Info& cinfo,
                        PlaneProvider * plane_provider,
                        double V0numAtom );

  class Cache : public CacheBase {
  public:
    void invalidateCache() override { ekin = -1.0; }
    //cache signature:
    double ekin = -1.0;//Start with invalid cache
    Vector dir;
    Vector dir_raw;//before normalisation
    //cache contents:
    double wl;
    VectD xs_commul;
    std::vector<GaussMos::ScatCache> scatcache;
  };

  void genScat( Cache&, RNG&, Vector& outdir ) const;
  void updateCache( Cache&, NeutronEnergy, const Vector& ) const;

  double m_threshold_ekin;
  GaussMos::PlaneList m_planes;//in families sharing d-spacing and fsquared.
  GaussMos m_gm;
  UniqueID m_uid;//for descriptions
};

NC::SCBragg::Model::Model(const NC::Info& cinfo, MosaicityFWHM mosaicity,
                          double dd, PlaneProvider * plane_provider,
                          double prec, double ntrunc)
  : m_threshold_ekin(kInfinity),
    m_gm(mosaicity,prec,ntrunc)
{
  m_gm.setDSpacingSpread(dd);

  //Always needs structure info:
  if (!cinfo.hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Structure information.");

  //Setup based on structure info:
  double V0numAtom = ( cinfo.getStructureInfo().n_atoms_average
                       * cinfo.getStructureInfo().volume );

  double maxdsp = setupFamilies( cinfo, plane_provider, V0numAtom );

  m_threshold_ekin = wl2ekin(maxdsp * 2.0);

}

NC::SCBragg::ModelPtr NC::SCBragg::createModel( const NC::Info& cinfo,
                                                MosaicityFWHM mosaicity,
                                                double dd,
                                                PlaneProvider * plane_provider,
                                                double prec, double ntrunc )
{
  return makeSO<const Model>( cinfo, mosaicity, dd, plane_provider,
                              prec, ntrunc );
}

NC::SCBragg::SCBragg( ModelPtr model, const RotMatrix& crystal2lab )
  : m_model( std::move(model) ),
    m_rot( Matrix( MatrixAllowCopy, crystal2lab ) ),
    m_invrot( ~crystal2lab )
{
  nc_assert_always( ncabs( crystal2lab.determinant() - 1.0 ) < 1e-6 );
}

NC::SCBragg::SCBragg( const NC::Info& cinfo,
                      const SCOrientation& sco,
                      MosaicityFWHM mosaicity,
                      double dd,
                      PlaneProvider * plane_provider,
                      double prec, double ntrunc)
  : SCBragg( createModel( cinfo, mosaicity, dd, plane_provider, prec, ntrunc ),
             getCrystal2LabRot( sco, getReciprocalLatticeRot(
                                  cinfo.getStructureInfo() ) ) )
{
}

NC::SCBragg::~SCBragg() = default;

double NC::SCBragg::Model::setupFamilies( const NC::Info& cinfo,
                                          NC::PlaneProvider * plane_provider,
                                          double V0numAtom )
{
  //expand crystal info
  nc_assert_always(cinfo.hasHKLInfo());
  nc_assert_always(cinfo.hasStructureInfo());
  nc_assert(m_planes.nFamilies()==0);

  //collect all planes, sorted by (dsp,fsq). To avoid issues connected to
  //floating point number keys, we store dspacing/fsquared as integers, keeping
  //precision down to 1e-10 angstrom and 1e-10 barn respectively.
  SCBraggSortMap planes;
  //but also try to avoid rounding issues when floating point values are not misbehaving:
  std::map<uint64_t,double> origvals_dsp;
  std::map<uint64_t,double> origvals_fsq;

  const double two30 = 1073741824.0;//2^30 ~= 1.07e9

  std::unique_ptr<PlaneProvider> ppguard;
  if (!plane_provider) {
    //fall back to standard plane provider
    ppguard = createStdPlaneProvider(&cinfo);//NB: cinfo must outlive ppguard!
    plane_provider = ppguard.get();
  } else {
    //use supplied plane provider - we must reset looping since we don't know
    //its current state:
    plane_provider->prepareLoop();
  }

  double maxdspacing(0);

  Optional<PlaneProvider::Plane> opt_plane;
  while ( ( opt_plane = plane_provider->getNextPlane() ).has_value() ) {
    auto& pl = opt_plane.value();
    if (pl.dspacing>maxdspacing)
      maxdspacing = pl.dspacing;

    nc_assert(pl.dspacing>0.0&&pl.fsq>0.0&&pl.dspacing<1e7&&pl.fsq<1e7);
    uint64_t ui_dsp = (uint64_t)(pl.dspacing*two30+0.5);
    uint64_t ui_fsq = (uint64_t)(pl.fsq*two30+0.5);

    //a bit messy, but nice to preserve values when possible:
    std::map<uint64_t,double>::iterator itOrig = origvals_dsp.find(ui_dsp);
    if (itOrig==origvals_dsp.end()) {
      origvals_dsp[ui_dsp] = pl.dspacing;
    } else if (ncabs(pl.dspacing-itOrig->second)>1e-12) {
      itOrig->second = -1;//multiple values observed ...!
    }
    itOrig = origvals_fsq.find(ui_fsq);
    if (itOrig==origvals_fsq.end()) {
      origvals_fsq[ui_fsq] = pl.fsq;
    } else if (ncabs(pl.fsq-itOrig->second)>1e-12) {
      itOrig->second = -1;//multiple values observed ...!
    }

    static_assert( std::is_same<SCBraggSortMap::key_type,
                                std::pair<uint64_t,uint64_t> >::value, "" );
    SCBraggSortMap::key_type key(ui_dsp,ui_fsq);

    SCBraggSortMap::iterator it = planes.find(key);
    if ( it != planes.end() ) {
      it->second.push_back(pl.demi_normal);
    } else {
      std::pair<SCBraggSortMap::key_type,std::vector<Vector> > newentry;
      newentry.first = key;
      newentry.second.push_back(pl.demi_normal);
      planes.insert(it,newentry);
    }
  }

  //Transfer to m_planes, in order of decreasing d-spacing:
  for ( auto& e : planes ) {
    auto itOrig = origvals_dsp.find(e.first.first);
    nc_assert(itOrig!=origvals_dsp.end());
    const double dsp = (itOrig->second > 0 ? itOrig->second : e.first.first / two30);

    itOrig = origvals_fsq.find(e.first.second);
    nc_assert(itOrig!=origvals_fsq.end());
    const double fsq = (itOrig->second > 0 ? itOrig->second : e.first.second / two30);

    nc_assert(fsq>0&&dsp>0);
    m_planes.addFamily( 0.5/dsp, fsq/V0numAtom, e.second );
  }

  return maxdspacing;
}


namespace NCRYSTAL_NAMESPACE {
  inline double SCBragg_cacheRound(double x) {
    //Cut off input at 15 decimals, which should be a negligible effect on any
    //realistic value of ekin in eV, but ensures that we don't
    //get call-order irreproducibilities.
    return std::floor(x*1e15+0.5)*1e-15;
  }
}

void NC::SCBragg::Model::updateCache( Cache& cache, NeutronEnergy ekin_raw, const NC::Vector& dir ) const
{
  //We check the cache validity on the rounded ekin value, but for simplicity we
  //keep the direction as it is. We could consider rounding the direction as
  //well...
  //NB: We used to check co-alignment of angles via a dot-product, but that is
  //actually numerically imprecise for small angles, leading to occurances of
  //cache validity where it should have been invalid.
  double ekin = SCBragg_cacheRound(ekin_raw.get());
  if ( cache.ekin==ekin ) {
    //Cache is valid if the angle between dir and cache.dir is less than
    //1e-12. First try to decide without the slow angle_highres:
    if ( dir == cache.dir_raw )
      return;//valid, identical direction
    const double dp = dir.dot(cache.dir);
    const bool angle_above_1em10 = ( !( dp > 0.0 )
                                     || dir.cross(cache.dir).mag2() > 1e-20*dp*dp );
    if ( !angle_above_1em10 && dir.angle_highres(cache.dir)<1.0e-12 )
      return;//valid
  }

  //Cache not valid!
  cache.ekin = -1.0;//in case of exceptions below
  cache.dir_raw = dir;
  cache.dir = dir;
  cache.dir.normalise();

  //Energy or direction is new, we must recalculate.

  cache.ekin = ekin;
  cache.wl = ekin2wl(ekin);
  nc_assert(cache.wl>=0);
  cache.scatcache.clear();
  cache.xs_commul.clear();
  if (cache.wl==0)
    return;//done, all cross-sections will be zero

  m_gm.calcCrossSections( m_planes, cache.wl, cache.dir,
                          cache.scatcache, cache.xs_commul );

  nc_assert(cache.xs_commul.empty()||cache.xs_commul.back()>0.0);
}

void NC::SCBragg::Model::genScat( Cache& cache, RNG& rng, NC::Vector& outdir ) const
{
  nc_assert(!cache.xs_commul.empty());
  nc_assert(cache.xs_commul.back()>0.0);
  nc_assert(cache.xs_commul.size()==cache.scatcache.size());

  std::size_t idx = RandIdxByWeightPicker::oneShotPick(rng,cache.xs_commul);
  nc_assert(idx<cache.scatcache.size());
  GaussMos::ScatCache& chosen_scatcache = cache.scatcache[idx];

  m_gm.genScat( rng, chosen_scatcache, cache.wl, cache.dir, outdir );
}

NC::EnergyDomain NC::SCBragg::domain() const noexcept
{
  return { NeutronEnergy{m_model->m_threshold_ekin}, NeutronEnergy{kInfinity} };
}

NC::CrossSect NC::SCBragg::crossSection(CachePtr& cp, NeutronEnergy ekin, const NeutronDirection& dir ) const
{
  if ( ekin.get() <= m_model->m_threshold_ekin )
    return CrossSect{ 0.0 };
  auto& cache = accessCache<Model::Cache>(cp);
  m_model->updateCache( cache, ekin, m_invrot * dir.as<Vector>() );
  return CrossSect{ cache.xs_commul.empty() ? 0.0 : cache.xs_commul.back() };
}

NC::ScatterOutcome NC::SCBragg::sampleScatter( CachePtr& cp, RNG& rng, NeutronEnergy ekin, const NeutronDirection& indir ) const
{
  if ( ekin.get() <= m_model->m_threshold_ekin ) {
    //Scatterings not actually possible at this configuration, so don't change
    //state:
    return { ekin, indir };
  }

  auto& cache = accessCache<Model::Cache>(cp);
  m_model->updateCache( cache, ekin, m_invrot * indir.as<Vector>() );

  if ( cache.xs_commul.empty() || cache.xs_commul.back()<=0.0 ) {
    //Again, scatterings are not actually possible here:
    return { ekin, indir };
  }

  Vector outdir;
  m_model->genScat( cache, rng, outdir );
  return { ekin, ( m_rot * outdir ).as<NeutronDirection>() };
}

NC::Optional<std::string> NC::SCBragg::specificJSONDescription() const
{
  auto nfam = m_model->m_planes.nFamilies();
  auto mos = m_model->m_gm.mosaicityFWHM();

  double dmin(-1),dmax(-1);
  if ( nfam ) {
    dmin = 0.5/m_model->m_planes.maxInv2dsp();
    dmax = 0.5/m_model->m_planes.minInv2dsp();
  }

  std::ostringstream ss;
  {
    std::ostringstream tmp;
    tmp << "nfamilies="<<nfam
        <<";dmin="<<dmin
        <<"Aa;dmax="<<dmax
        <<"Aa;mos="<<mos;
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "nfamilies", nfam );
  streamJSONDictEntry( ss, "dmin", dmin );
  streamJSONDictEntry( ss, "dmax", dmax );
  streamJSONDictEntry( ss, "mos", mos.dbl() );
  streamJSONDictEntry( ss, "model_uid", m_model->m_uid.getUniqueID().value );
  streamJSONDictEntry( ss, "crystal2lab", rotMatrixRows( m_rot ),
                       JSONDictPos::LAST );
  return ss.str();
}
