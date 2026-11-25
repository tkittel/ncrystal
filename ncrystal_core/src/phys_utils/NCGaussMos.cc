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

#include "NCrystal/internal/phys_utils/NCGaussMos.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {
  inline double GaussMos_cacheRound(double x) {
    //Cut off input at 15 decimals, which should be a negligible effect on any
    //realistic value of wavelength or 1/2dspacing, but ensures that we don't
    //get call-order irreproducibilities. The ncmax fct ensures that results
    //will be >0.0, so the code won't have to worry about FPE's when e.g. 1e-200
    //gets rounded to 0.0.
    return std::floor(ncmax(x,1e-15)*1e15+0.5)*1e-15;
  }
}

NC::GaussMos::GaussMos( MosaicityFWHM mosaicity, double prec, double ntrunc )
  : m_mos_truncN( ntrunc==0 ? GaussOnSphere::estimateNTruncFromPrec(prec) : ntrunc),
    m_prec(prec)
{
  double override_ntrunc = ncgetenv_dbl("GAUSSMOS_OVERRIDE_NTRUNC");
  if (override_ntrunc)
    m_mos_truncN = override_ntrunc;
  nc_assert(prec>=0);
  nc_assert(m_mos_truncN>=0);
  //Set mosaicity and trigger one call to updateDerivedValues:
  setMosaicity(mosaicity);
  nc_assert( m_mos_fwhm.dbl() != -99 );
  nc_assert( m_mos_sigma.dbl() != -99 );
}

NC::GaussMos::GaussMos( MosaicitySigma mosaicity, double prec, double ntrunc )
  : GaussMos( mosaicity.fwhm(), prec, ntrunc )
{
  //Same, but keeping exact MosaicitySigma value:
  m_mos_sigma = mosaicity;
  m_mos_sigma.validate();
}

NC::GaussMos::~GaussMos() = default;

void NC::GaussMos::updateDerivedValues()
{
  double truncangle = m_mos_truncN*m_mos_sigma.dbl();
  if ( ! (truncangle < kPiHalf) )
    NCRYSTAL_THROW(BadInput,"Mosaicity too large, truncation angle (sigma*Ntrunc) must be less than pi/2");
  m_gos.set(m_mos_sigma.dbl(), truncangle, m_prec );
  //Bound on ||n.k|-wl/2d| for contributing normals n, plus safety margin for
  //single-precision evaluations:
  m_presel_eps = static_cast<float>( 2.0*std::sin(0.5*truncangle) + 1e-5 );
}

void NC::GaussMos::setMosaicity( MosaicityFWHM mosaicity )
{
  mosaicity.validate();
  nc_assert_always(mosaicity.get()>0);
  m_mos_fwhm = mosaicity;
  m_mos_sigma = mosaicity.sigma();
  updateDerivedValues();
}

void NC::GaussMos::setMosaicity( MosaicitySigma mosaicity )
{
  //Same, but keeping exact MosaicitySigma value:
  auto mf = mosaicity.fwhm();
  setMosaicity(mf);
  m_mos_sigma = mosaicity;
}

void NC::GaussMos::setTruncationN(double N)
{
  nc_assert(N>0);
  if (m_mos_truncN != N) {
    m_mos_truncN = N;
    updateDerivedValues();
  }
}

void NC::GaussMos::setPrecision(double p)
{
  nc_assert(p>=0);
  if (m_prec != p) {
    m_prec = p;
    updateDerivedValues();
  }
}

void NC::GaussMos::setDSpacingSpread(double dd)
{
  if (dd==m_delta_d)
    return;
  NCRYSTAL_THROW(LogicError,"GaussMos::setDSpacingSpread not actually implemented and debugged fully yet");//TODO: implement+test this?
  nc_assert(dd>=0&&dd<0.99);
  m_delta_d = dd;
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Planes are preselected in chunks of this size, or in small chunks when
    //few planes remain (so few planes are always cheap):
    constexpr std::size_t gmpl_chunk = 64;
    constexpr std::size_t gmpl_chunksmall = 16;

    //Neutron parameters for the preselection (eps is the tolerance).
    struct GMPLNeutron { float kx, ky, kz, wl, eps; };

    //Preselection of the planes in a chunk which might contribute. A plane
    //with normal n and s=wl/(2*dspacing) only contributes to scattering of a
    //neutron with direction k, if ||n.k|-s|<=2*sin(truncangle/2). The
    //single-precision evaluation is safe when eps includes a safety
    //margin. Indices (offset by idx0) of selected planes are written to
    //idxout (must have room for a full chunk), and their number is returned.
    //
    //Performance depends crucially on the first loop being auto-vectorised,
    //which is why it has a fixed length and no conditionals (needed for
    //MSVC). Selected planes are those with a negative value, i.e. a sign bit.
    template<std::size_t NCHUNK>
    inline unsigned gmpl_preselect( const float* nx, const float* ny,
                                    const float* nz, const float* inv2dsp,
                                    const GMPLNeutron& neutron,
                                    std::size_t idx0, std::uint32_t* idxout )
    {
      static_assert( NCHUNK%8 == 0 && sizeof(float) == 4
                     && std::numeric_limits<float>::is_iec559, "" );
      const float kx(neutron.kx), ky(neutron.ky), kz(neutron.kz);
      const float wl(neutron.wl), eps(neutron.eps);
      float q[NCHUNK];
      for ( std::size_t j = 0; j < NCHUNK; ++j ) {
        const float dot = nx[j]*kx + ny[j]*ky + nz[j]*kz;
        q[j] = std::fabs( std::fabs(dot) - wl*inv2dsp[j] ) - eps;
      }
      //Check sign bits of 8 values at a time:
      std::uint64_t w[NCHUNK/2];
      std::memcpy( w, q, sizeof(w) );
      constexpr std::uint64_t signbits = 0x8000000080000000ull;
      std::uint64_t grp[NCHUNK/8];
      std::uint64_t any = 0;
      for ( std::size_t g = 0; g < NCHUNK/8; ++g ) {
        grp[g] = ( w[4*g] | w[4*g+1] | w[4*g+2] | w[4*g+3] ) & signbits;
        any |= grp[g];
      }
      if ( !any )
        return 0;
      unsigned n = 0;
      for ( std::size_t g = 0; g < NCHUNK/8; ++g ) {
        if ( !grp[g] )
          continue;
        for ( std::size_t j = 8*g; j < 8*g+8; ++j ) {
          std::uint32_t b;
          std::memcpy( &b, &q[j], sizeof(b) );
          idxout[n] = static_cast<std::uint32_t>(idx0+j);
          n += ( b >> 31 );//branch-free
        }
      }
      return n;
    }
  }
}

void NC::GaussMos::PlaneList::addFamily( double inv2dsp, double xsfact,
                                         Span<const Vector> deminormals )
{
  nc_assert_always( inv2dsp > 0.0 && inv2dsp < 1e30 && xsfact > 0.0 );
  nc_assert_always( m_families.empty()
                    || inv2dsp >= m_families.back().inv2dsp );
  const std::size_t n0 = m_planes.size();
  const std::size_t n1 = n0 + deminormals.size();
  nc_assert_always( n1 < std::numeric_limits<std::uint32_t>::max()-gmpl_chunk );
  const auto ifam = static_cast<std::uint32_t>( m_families.size() );
  m_families.push_back( Family{ inv2dsp, xsfact,
                                GaussMos_cacheRound(inv2dsp), n0 } );
  //Remove padding, add planes, and pad again (with planes never selected):
  for ( auto v : { &m_nx, &m_ny, &m_nz, &m_inv2dsp } )
    v->resize( n0 );
  for ( auto& n : deminormals ) {
    nc_assert( n.isUnitVector() );
    m_planes.push_back( Plane{ n, ifam } );
    m_nx.push_back( static_cast<float>( n.x() ) );
    m_ny.push_back( static_cast<float>( n.y() ) );
    m_nz.push_back( static_cast<float>( n.z() ) );
    m_inv2dsp.push_back( static_cast<float>( inv2dsp ) );
  }
  const std::size_t npad = ( (n1+gmpl_chunk-1) / gmpl_chunk ) * gmpl_chunk;
  for ( auto v : { &m_nx, &m_ny, &m_nz } )
    v->resize( npad, 0.0f );
  m_inv2dsp.resize( npad, 1e15f );
}

void NC::GaussMos::calcCrossSections( const PlaneList& pl,
                                      double wl,
                                      const NC::Vector& indir,
                                      std::vector<ScatCache>& cache,
                                      VectD& xs_commul ) const
{
  nc_assert(wl>0);
  nc_assert(indir.isUnitVector());

  //Only the first nplanes planes fulfill wl<2d:
  const double inv2dcutoff = (1.0-2*std::numeric_limits<double>::epsilon())/wl;
  const auto& fams = pl.m_families;
  auto itfamE = std::lower_bound( fams.begin(), fams.end(), inv2dcutoff,
                                  []( const PlaneList::Family& f, double v )
                                  { return f.inv2dsp < v; } );
  const std::size_t nplanes = ( itfamE == fams.end()
                                ? pl.m_planes.size() : itfamE->nbegin );

  const double cta = m_gos.getCosTruncangle();
  const GMPLNeutron neutron{ static_cast<float>(indir.x()),
                             static_cast<float>(indir.y()),
                             static_cast<float>(indir.z()),
                             static_cast<float>(wl), m_presel_eps };
  const float * const nx( pl.m_nx.data() );
  const float * const ny( pl.m_ny.data() );
  const float * const nz( pl.m_nz.data() );
  const float * const inv2dsp( pl.m_inv2dsp.data() );
  const double wl_rounded = GaussMos_cacheRound(wl);
  InteractionPars ip;
  constexpr auto nofam = std::numeric_limits<std::uint32_t>::max();
  std::uint32_t ifam = nofam;
  //Cross sections are summed up separately in each family:
  double xscommul = xs_commul.empty() ? 0.0 : xs_commul.back();
  double xsoffset(0.0), xssum(0.0);
  //Planes are preselected in batches, since it is faster to not interleave
  //preselection and exact calculations too finely:
  constexpr unsigned nbatch = 4*gmpl_chunk;
  std::uint32_t idx[nbatch];
  for ( std::size_t i0 = 0; i0 < nplanes; ) {
    unsigned nsel = 0;
    while ( i0 < nplanes && nsel <= nbatch - gmpl_chunk ) {
      if ( nplanes - i0 > gmpl_chunk - gmpl_chunksmall ) {
        nsel += gmpl_preselect<gmpl_chunk>( nx+i0, ny+i0, nz+i0, inv2dsp+i0,
                                            neutron, i0, idx + nsel );
        i0 += gmpl_chunk;
      } else {
        nsel += gmpl_preselect<gmpl_chunksmall>( nx+i0, ny+i0, nz+i0,
                                                 inv2dsp+i0, neutron, i0,
                                                 idx + nsel );
        i0 += gmpl_chunksmall;
      }
    }
    for ( unsigned isel = 0; isel < nsel; ++isel ) {
      const std::size_t i = idx[isel];
      if ( i >= nplanes )
        break;
      const auto& plane = pl.m_planes[i];
      const auto& fam = fams[plane.ifam];
      //New family => new sum (selecting via arrays to avoid branches which
      //can not be predicted):
      const std::size_t newfam = ( plane.ifam != ifam ? 1 : 0 );
      ifam = plane.ifam;
      const double xsoffset_select[2] = { xsoffset, xscommul };
      const double xssum_select[2] = { xssum, 0.0 };
      xsoffset = xsoffset_select[newfam];
      xssum = xssum_select[newfam];
      ip.setRoundedNoReuse( wl_rounded, fam.inv2dsp_rounded, fam.xsfact );
      const Vector& normal = plane.normal;
      const double dot = normal.dot(indir);
      const double sdotcptsq = (1.0 - dot * dot)*ip.m_cos_perfect_theta_sq;
      const double ds = dot * ip.m_sin_perfect_theta;

      //The normal (anti-normal) contributes when sqrt(sdotcptsq) exceeds
      //cta+ds (cta-ds), i.e. when it is within the truncated Gauss. This is
      //always the case for negative values, which includes sdotcptsq=0, i.e.
      //a neutron exactly along a normal near back-scattering. First check
      //the one of them most likely to contribute (the anti-normal if ds>0):
      const double absds = ncabs(ds);
      const double A0 = cta - absds;
      if ( A0 >= 0.0 && sdotcptsq <= A0*A0 )
        continue;//neither contributes
      const double A1 = cta + absds;
      if ( !( sdotcptsq > A1*A1 ) ) {
        //Usual case, just one of them contributes.
        const double xs = calcRawCrossSectionValue( ip, ncabs(dot) );
        if (xs) {
          xs_commul.push_back( xscommul = xsoffset + (xssum += xs) );
          cache.emplace_back( normal * ( ds > 0.0 ? -1.0 : 1.0 ),
                              ip.m_inv2dsp );
        }
      } else {
        //Both contribute, first the anti-normal, then the normal:
        for ( double sign : { -1.0, 1.0 } ) {
          const double xs = calcRawCrossSectionValue( ip, -sign*dot );
          if (xs) {
            xs_commul.push_back( xscommul = xsoffset + (xssum += xs) );
            cache.emplace_back( normal * sign, ip.m_inv2dsp );
          }
        }
      }
    }
  }
}

void NC::GaussMos::calcRawCrossSectionValues( InteractionPars& ip,
                                              std::size_t n,
                                              const double* cosvals,
                                              double* results ) const
{
  if ( !n )
    return;
  //The first value separately, which also initialises ip.m_Q if needed:
  results[0] = calcRawCrossSectionValue( ip, cosvals[0] );
  if ( !( ip.m_Q > 0.0 ) ) {
    //Special case (vanishing or infinite cross sections).
    for ( std::size_t i = 1; i < n; ++i )
      results[i] = calcRawCrossSectionValue( ip, cosvals[i] );
    return;
  }
  constexpr std::size_t nchunk = 32;
  double cg[nchunk], ca[nchunk], sa[nchunk];
  for ( auto& e : ca )
    e = ip.m_sin_perfect_theta;
  for ( auto& e : sa )
    e = ip.m_cos_perfect_theta;
  const double Q = ip.m_Q;
  for ( std::size_t i0 = 1; i0 < n; i0 += nchunk ) {
    const std::size_t m = ncmin( n - i0, nchunk );
    for ( std::size_t i = 0; i < m; ++i ) {
      nc_assert(ncabs(cosvals[i0+i])<=1.+1e-10);
      cg[i] = ncclamp( cosvals[i0+i], -1.0, 1.0 );
    }
    m_gos.circleIntegrals( m, cg, ca, sa, results + i0 );
    for ( std::size_t i = 0; i < m; ++i )
      results[i0+i] *= Q;
  }
}

void NC::GaussMos::genScat( RNG& rng, const ScatCache& cache, double wl_raw, const NC::Vector& indir, NC::Vector& outdir) const
{
  nc_assert(wl_raw>0.);
  nc_assert(cache.plane_inv2d()>0.);
  double wl = GaussMos_cacheRound(wl_raw);
  double inv2d = GaussMos_cacheRound(cache.plane_inv2d());
  nc_assert(wl>0);
  nc_assert(inv2d>0);

  //NB: indir has opposite sign convention from that used in derivation of equations. Hence a few extra minus signs in this method.
  nc_assert(cache.isValid());
  nc_assert(indir.isUnitVector());
  double sinthetabragg = wl * inv2d;
  nc_assert(ncabs(sinthetabragg)<=1.0);

  if ( sinthetabragg==0. ) {
    //completely forward scattering
    outdir = indir;
    return;
  }

  //NB: In principle we could cache sa and sg on the cache object, for increased
  //efficiency when scattering from the same neutron state constantly. However,
  //this is a small improvement, coming with mem/cpu costs, and would not help
  //the usual stepping-MC case.
  double ca = sinthetabragg;
  double sa = std::sqrt(1.0-ca*ca);
  double cg = ncclamp(-(indir.dot(cache.plane_normal())),-1.0,1.0);
  double sg = std::sqrt(1.0-cg*cg);
  double ct,st;
  if (!m_gos.genPointOnCircle( rng, cg,sg,ca,sa,ct,st) ) {
    //something went wrong, or numerical imprecision caused us to be called at
    //vanishing cross-section. Assume the latter.
    outdir = indir;
    return;
  }

  //genPointOnCircle selected an actual position of the normal for us (e.g. it
  //selected an actual crystallite). It did not do that in the lab frame, but in
  //a frame where the neutron sits at (0,0,1) rather than -indir, while the
  //nominal normal sits at (sg,0,cg) rather than cache.plane_normal, and the
  //selected normal sits at (sa*ct,sa*st,ca). Perform the reflection of the
  //neutron in that frame as well (using double-angle formulas to get
  //sin(2*alpha) and cos(2*alpha)):
  double s2a = 2*sa*ca;
  double c2a = ca*ca-sa*sa;
  outdir.set(s2a*ct,s2a*st,c2a);

  //Now rotate back to the lab frame:
  rotateToFrame( sg, cg, cache.plane_normal(), -indir, outdir, &rng );

  //For numerical safety:
  outdir.normalise();
  return;
}

void NC::GaussMos::InteractionPars::set(double wl_raw, double inv2dsp_raw, double xsfact) {
  nc_assert(wl_raw>0);
  nc_assert(inv2dsp_raw>0);
  nc_assert(xsfact>0.);
  double wl = GaussMos_cacheRound(wl_raw);
  double inv2dsp = GaussMos_cacheRound(inv2dsp_raw);
  if ( wl==m_wl && inv2dsp==m_inv2dsp ) {
    //great, we don't have to invalidate the Qprime/alpha values, but m_Q might have changed:
    m_xsfact = xsfact * 0.5;//absorb factor 1/2 here, to avoid a multiplication later
    m_Q = (m_Qprime > 0.0 ? m_Qprime * m_xsfact : -1/*will trigger recalc later*/);
  } else {
    //wl or inv2dsp changed:
    setRoundedNoReuse( wl, inv2dsp, xsfact );
  }
}
