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

#include "NCrystal/internal/phys_utils/NCGaussOnSphere.hh"
#include "NCrystal/internal/utils/NCRomberg.hh"
#include "NCrystal/internal/utils/NCGaussLegendre.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  class CosExpansionRadiusFct : public Fct1D {
    double m_prec;
  public:
    CosExpansionRadiusFct(double prec) : m_prec(prec) {};
    virtual ~CosExpansionRadiusFct(){}
    virtual double eval(double x) const
    {
      return (1.0-(1.0-0.5*x*x)/cos_mpi2pi2(x)) - m_prec;
    }
  };
  double gos_cosexpansionradius(double target_precision) {
    //Returns largest x for which (cos(x)-(1-0.5*x^2))/cos(x) is smaller than the
    //target_precision (but clips the target precision to at worst 0.999999
    //corresponding to a |x| smaller than roughly 81degrees.
    nc_assert_always(target_precision>0);
    CosExpansionRadiusFct f(ncmin(0.999999,target_precision));
    return findRoot(&f,0.0,0.999999*kPiHalf,1e-13);
  }
}

namespace NCRYSTAL_NAMESPACE {

  class GOSNormInt : public Romberg {
  public:
    GOSNormInt(double mos): Romberg(), m_c(-0.5/(mos*mos)) {}
    virtual ~GOSNormInt(){}
    double evalFunc(double x) const { return ncmax(0.0, std::sin( x ) * std::exp( m_c * x * x)); }
    virtual bool accept(unsigned level, double prev_estimate, double estimate,double,double) const
    {
      return level>12 || ncabs(prev_estimate-estimate)<1e-12*ncabs(estimate);
    }
  private:
    const double m_c;
  };

  //Numbers of Gauss-Legendre points (per segment) available for the
  //integration of circles:
  constexpr unsigned gos_gl_npts[7] = { 8, 10, 12, 14, 16, 20, 24 };

  //Maximum relative errors observed in the integration of circles, for
  //truncation angles of N sigma (rows: N<=3,3.5,..,7,8,10) and the numbers of
  //points in gos_gl_npts (columns). Found by comparison with 400 points, for
  //many circles and mosaicities of 0.1-40 degrees (see the developer app
  //devel/simplebuild/static/NCDev/app_glerrs, which regenerates the values):
  constexpr double gos_gl_maxN[11] = { 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5,
                                       7.0, 8.0, 10.0 };
  constexpr double gos_gl_errs[11][7] = {
    { 7e-07, 8e-09, 4e-11, 5e-13, 2e-15, 9e-16, 1e-15 },
    { 1e-06, 7e-08, 9e-10, 5e-12, 1e-13, 1e-15, 1e-15 },
    { 1e-05, 1e-07, 6e-09, 5e-11, 5e-13, 1e-15, 1e-15 },
    { 5e-05, 7e-07, 1e-08, 4e-10, 3e-12, 2e-15, 2e-15 },
    { 1e-04, 4e-06, 3e-08, 1e-09, 3e-11, 6e-15, 2e-15 },
    { 1e-04, 9e-06, 3e-07, 3e-09, 1e-10, 3e-14, 2e-15 },
    { 1e-04, 1e-05, 9e-07, 3e-08, 5e-10, 6e-13, 2e-15 },
    { 6e-04, 2e-05, 1e-06, 8e-08, 4e-09, 1e-12, 2e-15 },
    { 1e-03, 9e-05, 3e-06, 1e-07, 9e-09, 1e-11, 1e-14 },
    { 3e-03, 3e-04, 2e-05, 2e-06, 7e-08, 8e-11, 3e-13 },
    { 3e-03, 7e-04, 1e-04, 1e-05, 7e-07, 3e-09, 3e-11 } };

  //Coefficients of series for delta^2 in w=1-cos(delta), i.e.
  //2^(n+1)/(n^2*binom(2n,n)) for n=1..11:
  constexpr double gos_deltasq_c[11] = {
    2.0, 1.0/3.0, 4.0/45.0, 1.0/35.0, 16.0/1575.0, 8.0/2079.0, 32.0/21021.0,
    4.0/6435.0, 256.0/984555.0, 128.0/1154725.0, 512.0/10669659.0 };

  class SLTFct_EvalCosX : public Fct1D {
    double m_norm, m_expval;
  public:
    SLTFct_EvalCosX(double norm, double expval) : m_norm(norm), m_expval(expval) {}
    virtual ~SLTFct_EvalCosX(){}
    virtual double eval(double cosx) const
    {
      nc_assert(ncabs(cosx)<1.000001);
      double x = std::acos(ncmin(1.0,ncmax(-1.0,cosx)));
      return m_norm * std::exp(m_expval*x*x);
    }
  };



}

double NC::GaussOnSphere::calcNormFactor( double sigma, double trunc_angle )
{
  //Find by integration, but limit to at most 20*sigma for numerical
  //stability (contributions outside will likely not affect the result
  //anyway):
  GOSNormInt normcalc(sigma);
  double integral = k2Pi*normcalc.integrate(0.0,ncmin(trunc_angle,20.0*sigma));
  nc_assert(integral>0.0);
  return 1.0/integral;
}


NC::GaussOnSphere::GaussOnSphere()
{
  nc_assert(!isValid());
}

NC::GaussOnSphere::GaussOnSphere( double sigma, double trunc_angle, double prec  )
{
  nc_assert(!isValid());
  set(sigma,trunc_angle,prec);
  nc_assert(isValid());
}

NC::GaussOnSphere::~GaussOnSphere()
{
#ifndef NDEBUG
  if (m_stats.genpointworst)
    produceStatReport("destructed");
#endif
}

void NC::GaussOnSphere::produceStatReport(const char * callpt)
{
  markused(callpt);
#ifndef NDEBUG
  uint64_t worst = (m_stats.genpointcalled?(uint64_t)m_stats.genpointworst:0);
  NCRYSTAL_MSG("GaussOnSphere(sigma="<<m_sigma<<", truncangle="
               <<m_truncangle/m_sigma<<"sigma, prec="<<m_prec<<") "
               <<callpt<<". Used "<<m_stats.genpointtries
               <<" tries to generate "<<m_stats.genpointcalled
               <<" pts on circles (acceptance rate: "
               <<(m_stats.genpointtries?m_stats.genpointcalled*100.0/m_stats.genpointtries:0.0)
               <<"%). Worst case used "<<worst<<" tries."
               << " Performed "<<m_stats.circleintnumber
               <<" numerical circle integrations using an average of "
               <<(m_stats.circleintnumber?double(m_stats.circleintevals)/m_stats.circleintnumber:0.0)
               << " function evaluations each time (worst case used "
               <<m_stats.circleintworst<<" evaluations).");
#else
  markused(callpt);
#endif
}

void NC::GaussOnSphere::set(double sigma, double trunc_angle, double prec ) {
  nc_assert_always(sigma>0);
  nc_assert_always(trunc_angle>0);
  nc_assert_always(trunc_angle<kPiHalf);
  if ( ! (valueInInterval(0.9999e-7,0.10000001,prec) || valueInInterval(1.0,10000.0,prec)) )
    NCRYSTAL_THROW(BadInput,"prec must either be in the range [1e-7,1e-1] or in the range [1,10000].");
  if (prec<=1)
    prec = ncclamp(prec,1e-7,0.1);
  else
    prec = ncclamp(prec,1.0,10000.0);

  if (m_truncangle==trunc_angle&&m_sigma==sigma&&m_prec == prec)
    return;

#ifndef NDEBUG
  if (m_stats.genpointworst) {
    produceStatReport("settings changed.");
    m_stats.genpointworst = m_stats.genpointtries = m_stats.genpointcalled = 0;
    m_stats.circleintworst = m_stats.circleintnumber = m_stats.circleintevals = 0;
  }
  if (ncgetenv_bool("DEBUG_GAUSSONSPHERE")) {
    m_stats.genpointworst = 1;//enable stat collection
    m_stats.genpointtries = m_stats.genpointcalled = 0;
    m_stats.circleintworst = m_stats.circleintnumber = m_stats.circleintevals = 0;
  }
#endif
  unsigned nlt = 0;
  m_truncangle=trunc_angle;
  m_sigma=sigma;
  m_prec=prec;
  sincos(trunc_angle,m_cta,m_sta);
  {
    const double st2 = std::sin(0.5*trunc_angle);
    m_wtrunc = 2.0*st2*st2;
  }
  m_expfact = -0.5/(sigma*sigma);
  m_norm = calcNormFactor(sigma,trunc_angle);

  //analyse prec and trunc_angle parameters to establish some overall settings:
  nc_assert(prec>0);
  m_cf_s4 = 0.75 * sigma * sigma * sigma * sigma;
  m_inv_sigma2 = 1.0 / ( sigma * sigma );
  if (prec<1.0) {
    m_cf_lim = 3.0 * prec * sigma * sigma;
    double circleint_approx_maxangle = gos_cosexpansionradius(prec*0.5);//0.5 is safety
    nc_assert_always(circleint_approx_maxangle<kPiHalf);
    //The closed form takes cos(delta_0), whose rounding (~eps) gives relative
    //errors of ~eps/sigma^2 (relevant for sub-arcsec mosaicities):
    const bool cf_rounding_ok
      = ( 4.0 * std::numeric_limits<double>::epsilon() < 0.5 * prec * sigma * sigma );
    if (m_truncangle>circleint_approx_maxangle||circleint_approx_maxangle<=1e-10
        || !cf_rounding_ok ) {
      //approximation formula never valid at requested precision:
      m_circleint_k1 = kInfinity;
      m_circleint_k2 = 0.0;
    } else {
      m_circleint_k1 = cos_mpi2pi2(circleint_approx_maxangle);
      m_circleint_k2 = m_cta - 1e-5;
      //NB: 1e-5 above is for numerical stability, so we are sure to always fall
      //back to numerical integration when gamma~=0
    }
    m_numint_accuracy = ncclamp(prec * 0.1,1e-7,1e-4);
    nlt = std::max<unsigned>(40,std::min<unsigned>(10000,unsigned(10.0/std::sqrt(prec)+0.5)));
  } else {
    nc_assert_always(prec>=1.0);
    //prec>1: always use approximation formula, use prec pts for lookup tables (but at least 20)
    nlt = std::max<unsigned>(20,unsigned(prec+0.5));
    m_numint_accuracy = 1e-4;//won't actually be used
    m_circleint_k1 = 0.0;
    m_circleint_k2 = kInfinity;
    m_cf_lim = kInfinity;
  }
  nc_assert_always(nlt>=20&&nlt<=10000);
  //Settings for Gauss-Legendre integration of circles. Use segments of at
  //most 10 sigma, and points per segment giving a maximum error at least 10
  //times below the requested accuracy (or else the most points available):
  {
    const double ntrunc = trunc_angle / sigma;
    m_numint_nseg = static_cast<unsigned>( ncmax( 1.0, std::ceil( ntrunc / 10.0 ) ) );
    const double neff = ntrunc / m_numint_nseg;
    unsigned irow = 0;
    while ( irow + 1 < 11 && neff > gos_gl_maxN[irow] )
      ++irow;
    unsigned icol = 0;
    while ( icol + 1 < 7 && gos_gl_errs[irow][icol] > 0.1*m_numint_accuracy )
      ++icol;
    m_numint_npts = gos_gl_npts[icol];
    //Terms needed in series for delta^2 for relative precision 1e-17 (if
    //there are enough):
    m_numint_nseries = 0;
    for ( auto k : ncrange( 1u, 11u ) ) {
      if ( gos_deltasq_c[k] * std::pow( m_wtrunc, k ) < 2e-17 ) {
        m_numint_nseries = k;
        break;
      }
    }
  }

  //Prepare lookup table for N*exp(-k*delta^2) by cos(delta) value, thus
  //avoiding calls to both acos and exp. The lookup table essentially
  //covers the interval from ltmin = cos(trunc_angle*1.01 to 1.0.

  double ltmin(0.), ltinvdelta(0.);
  SplinedLookupTable lt_evalcosx;
  {
    double nltm1 = nlt-1.0;
    ltmin = cos_mpi2pi2(ncmin(kPiHalf,trunc_angle));
    ltinvdelta = nltm1/(1.0-ltmin);

    SLTFct_EvalCosX sltfct(m_norm,m_expfact);
    //determine h, "epsilon" parameter for derivative estimation. It can not be
    //too small, or numerical impression results.
    double h = 0.01/ltinvdelta;
    double minh = ncmin(1e-4,(1.0-ltmin)*0.1);
    h=ncmax(h,minh);
    double fprime_b = estimateSingleSidedDerivative(&sltfct, 1.0,-h);
    double fprime_a = estimateSingleSidedDerivative(&sltfct, ltmin,h);
    lt_evalcosx.set(&sltfct,ltmin,1.0,fprime_a,fprime_b,nlt,
                    "gosdensityfromcosx",
                    "Gaussian density on sphere as a function of cos(x)");
  }

  SplinedLookupTable lt_cf;
  {
    const double cta = m_cta;
    class CFFct final : public Fct1D {
      const GaussOnSphere* m_gos;
    public:
      CFFct( const GaussOnSphere* gos ) : m_gos(gos) {}
      double eval( double x ) const override { return m_gos->closedFormF(x); }
    } fct_cf(this);
    //determine h, "epsilon" parameter for derivative estimation. It can not be
    //too small, or numerical impression results.
    double h = 0.01*(1.0-cta)/nlt;
    double minh = ncmin(1e-4,(1.0-cta)*0.1);
    h=ncmax(h,minh);
    lt_cf.set(&fct_cf,cta,1.0,
              estimateSingleSidedDerivative(&fct_cf, cta, h),
              estimateSingleSidedDerivative(&fct_cf, 1.0, -h),nlt,
              "gosclosedform",
              "Closed form of circle integral (cf. closedFormF)");
  }

  m_lt_cf.swap(lt_cf);
  m_lt_evalcosx.swap(lt_evalcosx);
  setupTable();
}

void NC::GaussOnSphere::setupTable()
{
  //Where the closed-form approximation is not valid, circle integrals would
  //need numerical integration. Instead, a table is set up from such
  //integrals, and used with 4x4-point Lagrange interpolation. It covers
  //alpha in [1.2*truncangle,pi/2] and gamma-alpha in [-truncangle,
  //truncangle], with coordinates which need no transcendental functions and
  //resolve the important regions: z=sin(gamma-alpha)/sin(truncangle)
  //uniform in [-1,1], and r=q/(q+kappa), q=tan(alpha/2), uniform. Tabulated
  //is log(g/(sqrt(1-z^2)*sqrt(sin(alpha)/sin(gamma)))), which is smooth. The
  //table is validated at cell centres and refined until the maximum error is
  //10 times below the requested accuracy (or else not used), ignoring
  //negligible values.
  m_table = Table();
  const double amin = 1.2 * m_truncangle;
  if ( !( amin < kPiHalf - 0.05 ) )
    return;
  {
    //Not needed if closed form valid at the peak for alpha=amin:
    const double sa = std::sin( amin ), ca = std::cos( amin );
    if ( closedFormValid( 1.0, sa * sa, ca * ca ) )
      return;
  }
  const double qmin = std::tan( 0.5 * amin );
  const double target = 0.1 * m_numint_accuracy;
  for ( unsigned na = 64, nz = 256; na <= 256; na *= 2, nz *= 2 ) {
    Table t;
    t.na = na;
    t.nz = nz;
    t.kappa = qmin;
    t.r0 = 0.5;//qmin/(qmin+kappa)
    t.dr = ( 1.0 / ( 1.0 + t.kappa ) - t.r0 ) / ( na - 1 );
    t.dz = 2.0 / ( nz - 1 );
    t.sintrunc = m_sta;
    t.cosamin = std::cos( amin );
    //angles from coordinates:
    const double kappa = t.kappa;
    const double sintrunc = t.sintrunc;
    auto angles = [kappa,sintrunc]( double r, double z, double& cg,
                                    double& sg, double& ca, double& sa )
    {
      const double q = kappa * r / ( 1.0 - r );
      const double q2 = q * q;
      sa = 2.0 * q / ( 1.0 + q2 );
      ca = ( 1.0 - q2 ) / ( 1.0 + q2 );
      const double sd = z * sintrunc;
      const double cd = std::sqrt( 1.0 - sd * sd );
      sg = sa * cd + ca * sd;
      cg = ca * cd - sa * sd;
    };
    t.g.resize( std::size_t(na) * nz );
    for ( unsigned ia = 0; ia < na; ++ia ) {
      for ( unsigned iz = 0; iz < nz; ++iz ) {
        //Nodes at |z|=1 are moved slightly inwards, to avoid cancellations in
        //the numerical integration (its upper limit there is t'~0):
        const double z = ncclamp( -1.0 + iz * t.dz, -1.0 + 1e-8, 1.0 - 1e-8 );
        double cg, sg, ca, sa;
        angles( t.r0 + ia * t.dr, z, cg, sg, ca, sa );
        const double w = circleIntegralNumerical( cg, sg, ca, sa );
        if ( !( w > 0.0 ) )
          return;//density underflows (truncangle above ~38 sigma), no table
        vectAt( t.g, ia*nz+iz ) = std::log( w / ( std::sqrt( 1.0 - z*z )
                                        * std::sqrt( sa / sg ) ) );
      }
    }
    m_table = std::move( t );
    //Validate at the centres of cells of every second row. Errors are
    //relative to max(g,wfloor), with wfloor 0.4% of typical peak values of
    //norm*sqrt(2pi)*sigma (below that, e.g. at the truncation edges, only
    //absolute errors matter):
    const double wfloor = 1e-2 * m_norm * m_sigma;
    double maxerr = 0.0;
    for ( unsigned ia = 0; ia + 1 < na; ia += 2 ) {
      for ( unsigned iz = 0; iz + 1 < nz; ++iz ) {
        double cg, sg, ca, sa;
        angles( m_table.r0 + ( ia + 0.5 ) * m_table.dr,
                -1.0 + ( iz + 0.5 ) * m_table.dz, cg, sg, ca, sa );
        const double ref = circleIntegralNumerical( cg, sg, ca, sa );
        const double val = circleIntegralTable( cg, sg, ca, sa );
        maxerr = ncmax( maxerr, ncabs( val - ref ) / ncmax( ref, wfloor ) );
      }
    }
    if ( maxerr <= target )
      return;//ok
    m_table = Table();
  }
}

std::pair<unsigned,unsigned> NC::GaussOnSphere::tableSize() const
{
  return { m_table.na, m_table.nz };
}

double NC::GaussOnSphere::circleIntegralTable( double cg, double sg,
                                               double ca, double sa ) const
{
  const Table& t = m_table;
  nc_assert( !t.g.empty() );
  const double z = ncclamp( ( sg * ca - cg * sa ) / t.sintrunc, -1.0, 1.0 );
  const double q = sa / ( 1.0 + ca );
  const double fr = ( q / ( q + t.kappa ) - t.r0 ) / t.dr;
  const double fz = ( z + 1.0 ) / t.dz;
  const int ia = ncclamp( static_cast<int>( fr ), 1, int(t.na) - 3 );
  const int iz = ncclamp( static_cast<int>( fz ), 1, int(t.nz) - 3 );
  //Lagrange weights for nodes at -1,0,1,2 relative to ia (iz):
  auto weights = []( double u, double* w )
  {
    const double um1(u-1.0), um2(u-2.0), up1(u+1.0);
    w[0] = -u * um1 * um2 * ( 1.0 / 6.0 );
    w[1] = up1 * um1 * um2 * 0.5;
    w[2] = -up1 * u * um2 * 0.5;
    w[3] = up1 * u * um1 * ( 1.0 / 6.0 );
  };
  double wa[4], wz[4];
  weights( fr - ia, wa );
  weights( fz - iz, wz );
  const double * g = &t.g[ std::size_t(ia-1) * t.nz + (iz-1) ];
  double sum = 0.0;
  for ( unsigned i = 0; i < 4; ++i, g += t.nz )
    sum += wa[i] * ( wz[0]*g[0] + wz[1]*g[1] + wz[2]*g[2] + wz[3]*g[3] );
  return std::sqrt( ncmax( 0.0, 1.0 - z * z ) ) * std::sqrt( sa / sg )
    * std::exp( sum );
}

double NC::GaussOnSphere::circleIntegralSlow( double cg, double sg, double ca, double sa ) const
{
  const double sasg = sa*sg; nc_assert(sasg>=0.0);
  const double cacg = ca*cg;
  const double cd = cacg+sasg;

  if (cd<=m_cta)
    return 0.0;

  if (sasg<1e-14) {
    //degenerate case, either sa~=0 and the result should vanish since the
    //circle has vanishing circumference, or sg~=0 and the density will be
    //roughly identical on all points of the circle and we can perform the
    //integral by evaluating the Gaussian at x=alpha and multiplying by the
    //circle's circumference.
    return k2Pi*sa*evalCosX(ca);
  }

  //Use table if available and applicable, otherwise numerical integration:
  if ( !m_table.g.empty() && ca >= 0.0 && ca <= m_table.cosamin )
    return circleIntegralTable( cg, sg, ca, sa );
  return circleIntegralNumerical( cg, sg, ca, sa );
}

double NC::GaussOnSphere::circleIntegralNumerical( double cg, double sg,
                                                   double ca, double sa ) const
{
  const double sasg = sa*sg;
  const double cacg = ca*cg;
  const double cd = cacg+sasg;
  if ( cd <= m_cta || !( sasg > 0.0 ) )
    return 0.0;
  //Use Gauss-Legendre quadrature of Eq. 24 in the NCrystal elastic paper
  //(CPC 267 (2021) 108082) over the angle t along the circle (in [0,t'],
  //doubled by symmetry), with the exact density. With delta the angle to the
  //centre of the Gaussian, Eq. 23 gives w=1-cos(delta)=w0+2*sasg*sin^2(t/2),
  //where w0=1-cd is found via sin(gamma-alpha) to avoid cancellations for
  //small angles.
  const double sd = sg*ca - cg*sa;
  const double w0 = sd*sd/(1.0+cd);
  const double y = ( m_wtrunc - w0 ) / ( 2.0*sasg );//sin^2(t'/2)
  if ( !( y > 0.0 ) )
    return 0.0;
  const double tp = ( y >= 1.0 ? kPi : 2.0*std::asin( std::sqrt( y ) ) );
  if ( tp<=1e-12 )
    return 0.0;
  //Integrating over u=t/halfdt in [0,nseg], i.e. segments of unit length:
  const double halfdt = 0.5 * tp / m_numint_nseg;
  const double twosasg = 2.0 * sasg;
  auto f = [this,w0,twosasg,halfdt]( double u )
  {
    const double sh = std::sin( halfdt * u );
    return std::exp( m_expfact * numIntDeltaSq( w0 + twosasg*sh*sh ) );
  };
  const double integral = integrateGaussLegendre( f, 0.0, m_numint_nseg,
                                                  m_numint_npts,
                                                  m_numint_nseg );
#ifndef NDEBUG
  if ( m_stats.genpointworst>0 ) {
    const uint64_t nevals = m_numint_npts * m_numint_nseg;
    m_stats.circleintworst = std::max<uint64_t>(m_stats.circleintworst,nevals);
    ++m_stats.circleintnumber;
    m_stats.circleintevals += nevals;
  }
#endif
  //sa is the radius of circle (Jacobian):
  return 4.0 * sa * m_norm * halfdt * integral;
}

double NC::GaussOnSphere::closedFormF( double x ) const
{
  //Closed form of circle integrals, g=sqrt(u*sa/sg)*F(x)*(1+u*E[v^2]/(4sasg)),
  //with x=cos(delta_0) and u=x-cos(truncangle). With s=sqrt(2sasg)*sin(t/2),
  //Eq. 23 of the elastic paper gives w=1-cos(delta)=1-x+s^2 exactly, and
  //dt=sqrt(2/sasg)ds/cos(t/2) with 1/cos(t/2)=1/sqrt(1-s^2/(2sasg)) which is
  //~1+s^2/(4sasg). With s=sqrt(u)*v, F(x)=2sqrt(2)*int_0^1 W(delta(w))dv is
  //smooth at u=0 (unlike S=sqrt(u)F, splined before), and E[v^2] is the mean
  //of v^2 with weight W, approximated in closedFormCorr. Unlike Eq. 35 there,
  //this uses the exact density and includes the correction. Gauss-Legendre
  //with segments of 2 sigma in s (where W falls like exp(-s^2/sigma^2)).
  const double u = ncmax( 0.0, x - m_cta );
  const double w0 = 1.0 - x;
  const unsigned nseg = static_cast<unsigned>( ncmax( 1.0, std::ceil( 0.5 * std::sqrt( u ) / m_sigma ) ) );
  auto f = [this,w0,u]( double v )
  {
    return std::exp( m_expfact * numIntDeltaSq( w0 + u * v * v ) );
  };
  return 2.0 * kSqrt2 * m_norm
    * integrateGaussLegendre( f, 0.0, 1.0, 16, nseg );
}

double NC::GaussOnSphere::numIntDeltaSq( double w ) const
{
  //Squared angle, delta^2, from w=1-cos(delta).
  nc_assert( w >= 0.0 );
  if ( m_numint_nseries ) {
    double r = gos_deltasq_c[m_numint_nseries-1];
    for ( unsigned k = m_numint_nseries-1; k > 0; --k )
      r = r * w + gos_deltasq_c[k-1];
    return r * w;
  }
  const double b = 2.0 * std::asin( std::sqrt( 0.5 * w ) );
  return b * b;
}

void NC::GaussOnSphere::circleIntegrals( std::size_t n, const double* cg,
                                         const double* ca, const double* sa,
                                         double* results ) const
{
  //As circleIntegral, but for speed some square roots and divisions are
  //carried out in a separate loop. The tiny value prevents division by zero,
  //and never affects results (cf. the check on sasg).
  nc_assert(isValid());
  constexpr std::size_t nchunk = 32;
  constexpr double tiny = 1e-300;
  double sg[nchunk], sqrt_sa_over_sg[nchunk];
  while ( n ) {
    const std::size_t m = ncmin( n, nchunk );
    for ( std::size_t i = 0; i < m; ++i ) {
      nc_assert( ncabs(cg[i]) <= 1.0 && sa[i] >= 0.0 );
      nc_assert(ncabs(ca[i]*ca[i]+sa[i]*sa[i]-1.0)<1e-6);
      sg[i] = std::sqrt( 1.0 - cg[i]*cg[i] );
      sqrt_sa_over_sg[i] = std::sqrt( sa[i] / ( sg[i] + tiny ) );
    }
    for ( std::size_t i = 0; i < m; ++i ) {
      const double sasg = sa[i]*sg[i];
      const double cacg = ca[i]*cg[i];
      const double cd = cacg+sasg;
      if ( closedFormValid( cd, sasg, cacg ) ) {
        //closed-form approximation gives accurate result here:
        results[i] = closedFormValue( cd, sasg, sqrt_sa_over_sg[i] );
      } else {
        //special case, needs more careful treatment:
        results[i] = circleIntegralSlow( cg[i], sg[i], ca[i], sa[i] );
      }
    }
    n -= m;
    cg += m;
    ca += m;
    sa += m;
    results += m;
  }
}

bool NC::GaussOnSphere::genPointOnCircle( RNG& rng, double cg, double sg, double ca, double sa,
                                          double& ct, double& st ) const
{
  nc_assert(isValid());
  nc_assert(ncabs(cg*cg+sg*sg-1.0)<1e-6);
  nc_assert(ncabs(ca*ca+sa*sa-1.0)<1e-6);
  const double sasg = sa*sg; nc_assert(sasg>=0.0);
  const double cacg = ca*cg;
  const double cd = cacg+sasg;

  if (cd<=m_cta)
    return false;

  if (sasg<1e-14) {
    //Degenerate case, either sa~=0 and the density integrated over the circle
    //should vanish since the circle has vanishing circumference (so return
    //false), or sg~=0 and the density will be roughly identical on all points
    //of the circle, so pick random t:
    if (sa<1e-7)
      return false;
    std::tie(ct,st) = randPointOnUnitCircle(rng);
    return true;
  }

  //Find valid range of t values:
  nc_assert(sasg>0);
  double cos_tmax =  (m_cta-cacg)/sasg;
  if (cos_tmax>=1.0)//vanishing length of circle inside truncation zone
    return false;
  double tmax = ( cos_tmax<=-1.0 ? kPi : std::acos(cos_tmax) );

  //The highest contribution is at t=0, at which cos(delta) = cd. Generate t via MC-rejection.
  double densitymax = evalCosXInRange(cd)*1.00000001;//1.00000001 is overlay safety
  const int maxtriesplus1(1001);//we should usually use *much* fewer tries than this (averaging around 3-6 depending on parameters).
  int triesleft = maxtriesplus1;
  while (--triesleft) {
    ct = cos_mpipi( rng.generate()*tmax );//generate t uniformly in allowed range
    double cd_at_t = sasg*ct+cacg;
    double density_at_t = evalCosXInRange(cd_at_t);
    if ( density_at_t > densitymax ) {
      static std::atomic<bool> first(true);
      if ( first.exchange(false) ) {
        NCRYSTAL_WARN("Problems sampling with rejection method during GaussOnSphere::genPointOnCircle "
                      "invocation. Overlay value was not larger than actual cross-section value at sampled point "
                      "(overshot by factor of "<<(densitymax?density_at_t/densitymax:kInfinity)<<"). Further warnings"
                      " of this type will not be emitted.");
      }
    }
    if ( density_at_t > densitymax * rng.generate())
      break;
  }
#ifndef NDEBUG
  if (m_stats.genpointworst) {
    ++m_stats.genpointcalled;
    uint64_t triesused = maxtriesplus1-triesleft;
    m_stats.genpointtries += triesused;;
    m_stats.genpointworst = std::max<uint64_t>(m_stats.genpointworst,triesused);
  }
#endif
  if (triesleft<=0) {
    static std::atomic<bool> first(true);
    if ( first.exchange(false) ) {
      NCRYSTAL_WARN("Problems sampling with rejection method during GaussOnSphere::genPointOnCircle "
                    "invocation. Did not accept sampled value after "<<maxtriesplus1-1<<" attempts. Further warnings"
                    " of this type will not be emitted.");
    }
    return false;
  }
  st = std::sqrt(1.0-ct*ct);
  st = (rng.coinflip()?st:-st);//pick t in [-pi,pi], not just in [0,pi]
  return true;
}

double NC::GaussOnSphere::estimateNTruncFromPrec( double prec, double minval, double maxval )
{
  nc_assert(minval>0&&maxval>minval&&prec>=0.0);
  //special cases: reference mode always selects ntrunc=maxval, approx mode always ntrunc=minval
  if (prec==0)
    return maxval;
  if (prec>=1.0)
    return minval;
  nc_assert_always(prec>0.0&&prec<1.0);
#if 0
  //Volume within a radius of N*sigma from the center of a 1D Gaussian in the
  //plane is 1-erf(N/sqrt(2)), which can not be inverted easily, so we use
  //resort to root-finding:
  class EstNTruncFct : public Fct1D {
    double m_prec;
  public:
    EstNTruncFct(double prec) : m_prec(prec) {};
    virtual ~EstNTruncFct(){}
    virtual double eval(double x) const
    {
      return ncerfc(x*kInvSqrt2) - m_prec;
    }
  };
  EstNTruncFct f(prec);
  if (f.eval(maxval)>=0)
    return maxval;
  if (f.eval(minval)<=0)
    return minval;
  double n = findRoot(&f,minval,maxval,1e-13);
#else
  //Volume within a radius of N*sigma from the center of a 2D Gaussian in the
  //plane is 1-exp(-N^2/2), which can be use to find the following exact
  //formula:
  double n = std::sqrt(-2.0*std::log(ncmax(1e-300,prec)));
#endif
  //Returned clamped value and with a safety factor of 1.1:
  return ncclamp(1.1*n,minval,maxval);
}
