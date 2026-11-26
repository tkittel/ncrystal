#ifndef NCrystal_GaussMos_hh
#define NCrystal_GaussMos_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/phys_utils/NCGaussOnSphere.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCSpan.hh"

namespace NCRYSTAL_NAMESPACE {

  class GaussMos {
  public:

    //Helper class which implements a mosaic model, assuming a (truncated)
    //Gaussian distribution of crystallite orientations. The provided mosaicity
    //should be in radians, and will be interpreted to be either the sigma or
    //FWHM of the (untruncated) Gaussian mosaicity distribution, depending on
    //the value of the mosaicity_is_fhwm parameter. The precision ("prec")
    //parameter is passed on to the internal GaussOnSphere object (see
    //NCGaussOnSphere.hh for details). If ntrunc is 0,
    //GaussOnSphere::estimateNTruncFromPrec will be used to estimate the
    //truncation range from prec, otherwise the provided value will be used
    //directly.
    GaussMos( MosaicityFWHM, double prec = 1e-3, double ntrunc=0 );
    GaussMos( MosaicitySigma, double prec = 1e-3, double ntrunc=0 );
    ~GaussMos();

    //Mosaicity can also be changed after construction:
    void setMosaicity( MosaicityFWHM );
    void setMosaicity( MosaicitySigma );

    //Modify certain aspects:
    void setTruncationN(double);
    void setPrecision(double);

    void setDSpacingSpread(double);//Enable dspacing deviation in non-ideal crystal [default value of 0.0 means no deviation]

    //Access parameters of Gaussian mosaicity distribution:
    MosaicityFWHM mosaicityFWHM() const;
    MosaicitySigma mosaicityGaussSigma() const;
    double mosaicityGaussNormFact() const;//gauss(x)=normfact*exp(-k*x*x)
    double mosaicityTruncationN() const;
    double mosaicityTruncationAngle() const;
    double mosaicityCosTruncationAngle() const;
    double mosaicitySinTruncationAngle() const;
    double precision() const;
    const GaussOnSphere& gos() const { return m_gos; }

    //Wavelengths and 1/(2*dspacing) values are rounded to 15 decimals before
    //use, which has a negligible effect on realistic values, but ensures that
    //results do not depend on the order of calls. Results are >0.0, so 1e-200
    //does not become 0.0:
    static double cacheRound( double x )
    {
      return std::floor( ncmax( x, 1e-15 ) * 1e15 + 0.5 ) * 1e-15;
    }

    //Before calculating cross-sections, the relevant interaction parameters for
    //the neutron and plane family must be set in an InteractionPars object
    //(reusing it might save calculations under certain conditions). The xsfact
    //parameter should typically be set to fsquared/(v0*natoms), and it will
    //simply scale the resulting cross-sections. It must hold that
    //neutron_wavelength<=2*dspacing (otherwise cross-sections are trivially
    //zero and GaussMos should not be invoked in any case).
    class InteractionPars;

    //Calculate cross-section, given (the cosine of) an angle between the
    //neutron direction and a given normal. This function does NOT apply any
    //initial truncation cutoff!!  It is assumed that this function will only be
    //called after the caller has verified that the truncation cutoff is
    //satisfied (but note that the truncationN value still plays a role for the
    //normalisation of the mosaicity distribution used):
    double calcRawCrossSectionValue( InteractionPars& ip,
                                     double cos_angle_indir_normal ) const;

    //Cross-sections for a large number of demi-normals can be found in one go,
    //by first collecting them in a PlaneList. This method DOES implement the
    //Gauss truncation internally (and exactly), and ignores planes with
    //neutron_wavelength>=2*dspacing. Intermediate results for individual
    //deminormals with non-zero contributions will be appended to the passed-in
    //cache vector, while the commulative values of the corresponding
    //individual cross-sections will be appended to the passed in xs_commul
    //vector. This setup allows the caller to subsequently select a plane to
    //scatter on with a binary search in the xs_commul vector, and subsequently
    //use the ScatCache object at the same index in the cache vector, to
    //generate scatterings.
    //
    //The method is fast even for many planes, since exact calculations are
    //only carried out for the few planes surviving a single-precision
    //preselection (which is conservative, so results are unaffected by it).
    class ScatCache;
    class PlaneList;
    void calcCrossSections( const PlaneList&,
                            double neutron_wavelength,
                            const Vector& neutron_indir,
                            std::vector<ScatCache>& cache,
                            VectD& xs_commul ) const;

    //Scatterings can only be generated once appropriate info has been found via
    //previous calls to cross-section methods, and with relevant info embedded
    //into ScatCache objects (of course, they will only be relevant for the
    //neutron state for which they were created):
    class ScatCache {
    public:
      ScatCache() = default;//invalid cache object
      ScatCache(const Vector& pn, double i2d);
      void set(const Vector& pn, double i2d);
      void clear();//invalidate
      bool isValid() const;
      const Vector& plane_normal() const;
      double plane_inv2d() const;
    private:
      Vector m_plane_normal;//signed normal
      double m_plane_inv2d = 0.0;
    };

    //Generate scatterings. Needs a valid ScatCache object, a source of random
    //numbers, and the neutron state parameters (must be consistent with the
    //ones used when the ScatCache object was set up):
    void genScat( RNG&, const ScatCache&,
                  double neutron_wavelength, const Vector& neutron_indir,
                  Vector& neutron_outdir ) const;

  private:
    GaussOnSphere m_gos;
    MosaicityFWHM m_mos_fwhm = MosaicityFWHM{-99.0};
    double m_mos_truncN = 0.0;
    MosaicitySigma m_mos_sigma = MosaicitySigma{-99};
    double m_prec = 0.0;
    double m_delta_d = 0.0;
    float m_presel_eps = 0.0f;//tolerance of PlaneList preselection
    void updateDerivedValues();
    double calcRawCrossSectionValueInit( InteractionPars&, double ) const;
  };

  class GaussMos::PlaneList final : private MoveOnly {
  public:
    //List of (demi) normals of crystal planes, grouped in families sharing
    //d-spacing and xsfact (cf. InteractionPars). Families must be added in
    //order of increasing inv2dsp=1/(2*dspacing).
    void addFamily( double inv2dsp, double xsfact,
                    Span<const Vector> deminormals );
    std::size_t nFamilies() const { return m_families.size(); }
    std::size_t nPlanes() const { return m_planes.size(); }
    double minInv2dsp() const;//of first family (0 if empty)
    double maxInv2dsp() const;//of last family (0 if empty)
  private:
    friend class GaussMos;
    struct Family {
      double inv2dsp = 0.0;
      double xsfact = 0.0;
      double inv2dsp_rounded = 0.0;//as needed by InteractionPars::setRounded
      std::size_t nbegin = 0;//index of first normal
    };
    struct Plane { Vector normal; std::uint32_t ifam = 0; };
    std::vector<Family> m_families;
    std::vector<Plane> m_planes;//of all families, one after the other
    //Single-precision copies for preselection, padded to whole chunks:
    std::vector<float> m_nx, m_ny, m_nz, m_inv2dsp;
  };

  class GaussMos::InteractionPars {
  public:
    InteractionPars(double neutron_wavelength, double inv2dsp, double xsfact);
    InteractionPars() = default;
    bool isValid() const;
    void set(double neutron_wavelength, double inv2dsp, double xsfact);
  private:
    friend class GaussMos;
    //As set(), but with wavelength and inv2dsp already rounded as set() would,
    //and without trying to reuse values from the previous call (same results,
    //but faster when most calls are for different values):
    void setRoundedNoReuse( double wavelength, double inv2dsp, double xsfact );
    //Frequently accessed members first:
    double m_Q = 0.0;//Once initialised, m_Q = m_Qprime * m_xsfact
    double m_sin_perfect_theta = 0.0;
    double m_cos_perfect_theta = 0.0;
    //Less often accessed
    double m_wl = -1.0;
    double m_wl3 = 0.0;
    double m_inv2dsp = -1.0;
    double m_cos_perfect_theta_sq = 0.0;
    double m_alpha = 0.0;
    double m_Qprime = 0.0;
    double m_xsfact = 0.0;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline MosaicityFWHM GaussMos::mosaicityFWHM() const { return m_mos_fwhm; }
  inline MosaicitySigma GaussMos::mosaicityGaussSigma() const { nc_assert(m_mos_sigma.dbl()==m_gos.getSigma()); return MosaicitySigma{ m_gos.getSigma() }; }
  inline double GaussMos::mosaicityGaussNormFact() const { return m_gos.getNormFactor(); }
  inline double GaussMos::mosaicityTruncationN() const { return m_mos_truncN; }
  inline double GaussMos::mosaicityTruncationAngle() const { return m_gos.getTruncangle(); }
  inline double GaussMos::mosaicityCosTruncationAngle() const { return m_gos.getCosTruncangle(); }
  inline double GaussMos::mosaicitySinTruncationAngle() const { return m_gos.getSinTruncangle(); }
  inline double GaussMos::precision() const { return m_gos.getPrecisionParameter(); }

  inline GaussMos::InteractionPars::InteractionPars(double wl, double inv2dsp, double xsfact)
  {
    set(wl, inv2dsp, xsfact);
  }

  inline void GaussMos::InteractionPars::setRoundedNoReuse( double wl,
                                                            double inv2dsp,
                                                            double xsfact )
  {
    nc_assert(wl>0);
    nc_assert(inv2dsp>0);
    nc_assert(xsfact>0.);
    m_xsfact = xsfact * 0.5;//absorb factor 1/2 here, to avoid a multiplication later
    m_wl = wl;
    m_wl3 = m_wl*m_wl*m_wl;
    //New thetabragg and invalid Q/Qprime/alpha/cos_perfect_theta values:
    m_inv2dsp = inv2dsp;
    //Clamped, since rounding can give values just above 1 at the threshold:
    m_sin_perfect_theta = ncmin( wl * inv2dsp, 1.0 );
    nc_assert(valueInInterval(0.,1.,m_sin_perfect_theta));
    m_cos_perfect_theta_sq = 1 - m_sin_perfect_theta*m_sin_perfect_theta;
    m_Q = m_Qprime = m_cos_perfect_theta = -1;//invalidate
    m_alpha = -99;//invalidate
  }

  inline bool GaussMos::InteractionPars::isValid() const
  {
    return m_wl>0;
  }

  inline double GaussMos::PlaneList::minInv2dsp() const
  {
    return m_families.empty() ? 0.0 : m_families.front().inv2dsp;
  }

  inline double GaussMos::PlaneList::maxInv2dsp() const
  {
    return m_families.empty() ? 0.0 : m_families.back().inv2dsp;
  }

  inline GaussMos::ScatCache::ScatCache(const Vector& pn, double i2d) : m_plane_normal(pn), m_plane_inv2d(i2d) { nc_assert(i2d>0&&pn.isUnitVector()); }
  inline void GaussMos::ScatCache::set(const Vector& pn, double i2d) { nc_assert(i2d>0&&pn.isUnitVector()); m_plane_normal = pn; m_plane_inv2d = i2d; }
  inline void GaussMos::ScatCache::clear() { m_plane_inv2d=0; }
  inline bool GaussMos::ScatCache::isValid() const { return m_plane_inv2d>0; }
  inline const Vector& GaussMos::ScatCache::plane_normal() const { nc_assert(isValid()); return m_plane_normal; }
  inline double GaussMos::ScatCache::plane_inv2d() const { nc_assert(isValid()); return m_plane_inv2d; }

  inline double GaussMos::calcRawCrossSectionValue( InteractionPars& ip, double cos_angle_indir_normal ) const
  {
    nc_assert(ip.isValid());
    nc_assert(ncabs(cos_angle_indir_normal)<=1.+1e-10);
    cos_angle_indir_normal = ncclamp(cos_angle_indir_normal,-1.0,1.0);
    if (ip.m_Q>0.) {
      double sin_angle_indir_normal = std::sqrt(1.0-cos_angle_indir_normal*cos_angle_indir_normal);//>0 since angle is in 0..pi.
      return ip.m_Q * m_gos.circleIntegral( cos_angle_indir_normal, sin_angle_indir_normal, ip.m_sin_perfect_theta, ip.m_cos_perfect_theta );
    }
    return calcRawCrossSectionValueInit(ip,cos_angle_indir_normal);
  }

  inline double GaussMos::calcRawCrossSectionValueInit(InteractionPars& ip, double cos_angle_indir_normal ) const
  {
    nc_assert(ip.m_Q<=0.);
    if ( ip.m_Qprime == -1 ) {
      //First call after getting new wavelength or d-spacing.
      ip.m_cos_perfect_theta = std::sqrt(ip.m_cos_perfect_theta_sq);// = sinalpha
      double tmp2 = ip.m_cos_perfect_theta*ip.m_sin_perfect_theta;
      if (tmp2>0) {
        //usual case:
        ip.m_Qprime = ip.m_wl3 / tmp2;
      } else {
        //Avoid a zero division when wavelength is 0 or 2*dspacing. The actual
        //cross-section will be 0 or infinity at the two limits (special value -2.0
        //designates infinity).
        //Values are for wl~=2dsp and wl~=0 respectively:
        ip.m_Qprime = (ip.m_sin_perfect_theta>0.5&&ip.m_xsfact) ? -2.0 : 0.0;
      }
    }
    if (ip.m_Qprime>0.) {
      ip.m_Q = ip.m_Qprime * ip.m_xsfact;
      nc_assert(ip.m_Q>0.);
      return calcRawCrossSectionValue(ip,cos_angle_indir_normal);
    }

    if (ip.m_Qprime) {
      nc_assert(ip.m_Qprime==-2);
      return kInfinity;//Q=inf, W factor irrelevant (assuming called within truncation radius)
    } else {
      return 0.0;//Q=0, W factor irrelevant
    }
  }
}

#endif
