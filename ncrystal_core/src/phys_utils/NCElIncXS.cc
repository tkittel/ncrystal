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

#include "NCrystal/internal/phys_utils/NCElIncXS.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

NC::ElIncXS::ElIncXS( const VectD& elm_msd,
                      const VectD& elm_bixs,
                      const VectD& elm_scale )
{
  set( elm_msd, elm_bixs, elm_scale );
}

//FIXME: Revisit the libm usage (std::exp here) for platform reproducibility:
//stable_expm1 etc. are slower (especially in wasm without fast fma), the
//NCMath.hh approximations too. Maybe tabulate (1-exp(-t))/t, e.g. with a
//spline from NCSpline.hh.
double NC::ElIncXS::eval_1mexpmtdivt(double t)
{
  //safe eval of (1-exp(-t))/t for t>=0.0
  nc_assert(t>=0.0);
  if ( t < 0.01 ) {
    //Taylor expansion, sum_k (-t)^k/(k+1)!, for numerical stability. Terms up
    //to t^6 make the truncation error (<t^7/8!=2.5e-19 at t=0.01) negligible,
    //so there is no jump at t=0.01:
    double p = 1.0 / 5040.0;
    p = ncfma( t, p, -1.0 / 720.0 );
    p = ncfma( t, p, 1.0 / 120.0 );
    p = ncfma( t, p, -1.0 / 24.0 );
    p = ncfma( t, p, 1.0 / 6.0 );
    p = ncfma( t, p, -0.5 );
    return ncfma( t, p, 1.0 );
  }
  if ( t > 37.0 ) {
    //Limiting behaviour, exact to double precision since exp(-37)<1e-16:
    return 1.0 / t;
  }
  //evaluate at intermediate t using actual formula (precision is fine without
  //std::expm1 since t>0.01, and std::exp is typically faster):
  return ( 1.0 - std::exp( -t ) ) / t;
}

NC::CrossSect NC::ElIncXS::evaluateMonoAtomic(NeutronEnergy ekin, double meanSqDisp, SigmaBound bound_incoh_xs)
{
  nc_assert(ekin.get()>=0.0&&meanSqDisp>=0.0&&bound_incoh_xs.get()>=0.0);
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  return CrossSect{ bound_incoh_xs.dbl() * eval_1mexpmtdivt(kkk * meanSqDisp * ekin.dbl()) };
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Piecewise cubic Hermite interpolation of f(t)=1-exp(-t), with nodes t
    //from 2e-10 to beyond -log(1e-10), and coefficients p (in u, see
    //NCIsoXSUtils.hh). On [a,a+h], |f-H| <= h^4*max|f^(4)|/384 =
    //h^4*exp(-a)/384, and f >= f(a), so h=(192e-10*f(a)/exp(-a))^(1/4) gives a
    //relative error below 0.5e-10. Below and above the nodes, f(t)=t and
    //f(t)=1 are good to 0.5e-10 and 0.4e-10:
    struct ElIncUnivTable {
      VectD t, p;
    };

    const ElIncUnivTable& elIncUnivTable()
    {
      static const ElIncUnivTable tab = []()
      {
        constexpr double t_low = 1e-10;
        constexpr double t_high = 24.0;//exp(-t_high)<0.4e-10
        ElIncUnivTable res;
        res.t.push_back( t_low );
        while ( res.t.back() < t_high ) {
          const double a = res.t.back();
          const double fa = -stable_expm1( -a );
          const double ea = stable_exp( -a );
          const double h = std::sqrt( std::sqrt( 192e-10 * fa / ea ) );
          const double b = a + h;
          const double fb = -stable_expm1( -b );
          const double d0 = h * ea;
          const double d1 = h * stable_exp( -b );
          res.p.insert( res.p.end(), { fa, d0,
                                       ncfma( 3.0, fb - fa,
                                              -ncfma( 2.0, d0, d1 ) ),
                                       ncfma( 2.0, fa - fb, d0 + d1 ) } );
          res.t.push_back( b );
        }
        return res;
      }();
      return tab;
    }
  }
}

NC::IsoXSTable NC::ElIncXS::tabulateMonoAtomic( double meanSqDisp,
                                                            SigmaBound bound_incoh_xs )
{
  nc_assert( meanSqDisp >= 0.0 && meanSqDisp < 1e6 );
  nc_assert( bound_incoh_xs.dbl() >= 0.0 && bound_incoh_xs.dbl() < 1e12 );
  //e*xs = sigma/c * f(c*e), with f(t)=1-exp(-t) from the universal table:
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  const double sigma = bound_incoh_xs.dbl();
  const double c = kkk * meanSqDisp;
  IsoXSTable res;
  if ( meanSqDisp < 1e-290 ) {
    //Constant cross section (for msd<1e-290 good to 1e-10 below 1e277eV, and
    //keeping E=t/c and the coefficients sigma/c*p, |p|<4, finite):
    res.E = { 1.0, 2.0 };
    res.p = { sigma, sigma, 0.0, 0.0 };
    res.tail[1] = res.tail[3] = sigma;
    return res;
  }
  const auto& ut = elIncUnivTable();
  const double k = sigma / c;
  res.E.reserve( ut.t.size() );
  for ( auto t : ut.t )
    res.E.push_back( t / c );
  res.p.reserve( ut.p.size() );
  for ( auto v : ut.p )
    res.p.push_back( v * k );
  res.tail[1] = sigma;
  res.tail[2] = k;
  return res;
}

NC::IsoXSTable NC::ElIncXS::tabulate() const
{
  //Same construction as the universal table, for the sum over elements: with
  //t=c*e (c for the largest msd), e*xs = K*f(t) where f(t) = sum_j
  //a_j*(1-exp(-g_j*t)) + b*t (sum a_j = 1, 0<g_j<=1, b*t for elements with
  //zero msd). Then |f^(4)| = sum_j a_j*g_j^4*exp(-g_j*t) decreases and f
  //increases, so h=(192e-10*f(a)/|f^(4)(a)|)^(1/4) again gives a relative
  //error below 0.5e-10. Nodes end once exp(-g_min*t)<0.4e-10.
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  constexpr double msd_zero = 1e-290;//as in tabulateMonoAtomic
  if ( m_elm_data.size() == 1 )
    return tabulateMonoAtomic( m_elm_data.front().first,
                               SigmaBound{ m_elm_data.front().second } );
  StableSum sig_all, sig_zero, sum_k;
  double msdmax = 0.0;
  for ( auto& em : m_elm_data ) {
    sig_all.add( em.second );
    if ( em.first < msd_zero )
      sig_zero.add( em.second );
    else if ( em.second > 0.0 )
      msdmax = ncmax( msdmax, em.first );
  }
  if ( !msdmax )
    return tabulateMonoAtomic( 0.0, SigmaBound{ sig_all.sum() } );
  const double c = kkk * msdmax;
  struct Term { double g, a; };
  std::vector<Term> terms;
  for ( auto& em : m_elm_data ) {
    if ( em.first < msd_zero || !( em.second > 0.0 ) )
      continue;
    terms.push_back( { em.first / msdmax, em.second / ( kkk * em.first ) } );
    sum_k.add( terms.back().a );
  }
  const double K = sum_k.sum();
  double gmin = 1.0;
  for ( auto& tm : terms ) {
    tm.a /= K;
    gmin = ncmin( gmin, tm.g );
  }
  const double b = sig_zero.sum() / ( K * c );
  //f, f' and |f^(4)| at t:
  struct FVals { double f, fp, f4; };
  auto evalF = [&terms,b]( double t )
  {
    StableSum f, fp, f4;
    for ( auto& tm : terms ) {
      const double ex = stable_exp( -tm.g * t );
      f.add( -tm.a * stable_expm1( -tm.g * t ) );
      fp.add( tm.a * tm.g * ex );
      const double g2 = tm.g * tm.g;
      f4.add( tm.a * g2 * g2 * ex );
    }
    f.add( b * t );
    return FVals{ f.sum(), fp.sum() + b, f4.sum() };
  };
  const double t_low = 1e-10;
  const double t_high = 24.0 / gmin;
  VectD tt, pp;
  tt.push_back( t_low );
  FVals va = evalF( t_low );
  while ( tt.back() < t_high ) {
    const double ta = tt.back();
    nc_assert( va.f > 0.0 && va.f4 > 0.0 );
    const double h = std::sqrt( std::sqrt( 192e-10 * va.f / va.f4 ) );
    const double tb = ta + h;
    const FVals vb = evalF( tb );
    const double d0 = h * va.fp;
    const double d1 = h * vb.fp;
    pp.insert( pp.end(), { va.f, d0,
                           ncfma( 3.0, vb.f - va.f, -ncfma( 2.0, d0, d1 ) ),
                           ncfma( 2.0, va.f - vb.f, d0 + d1 ) } );
    tt.push_back( tb );
    va = vb;
  }
  IsoXSTable res;
  res.E.reserve( tt.size() );
  for ( auto t : tt )
    res.E.push_back( t / c );
  res.p.reserve( pp.size() );
  for ( auto v : pp )
    res.p.push_back( v * K );
  res.tail[1] = sig_all.sum();
  res.tail[2] = K;
  res.tail[3] = sig_zero.sum();
  return res;
}

NC::ElIncXS::~ElIncXS() = default;

void NC::ElIncXS::set( const VectD& elm_msd,
                       const VectD& elm_bixs,
                       const VectD& elm_scale )
{
  //sanity check element data:
  nc_assert_always(elm_msd.size()==elm_bixs.size());
  nc_assert_always(elm_msd.size()==elm_scale.size());
  for ( auto i : ncrange( elm_msd.size() ) ) {
    nc_assert_always(elm_msd.at(i)>=0.0&&elm_msd.at(i)<1e6);
    nc_assert_always(elm_bixs.at(i)>=0.0&&elm_bixs.at(i)<1e6);
    nc_assert_always(elm_scale.at(i)>=0.0&&elm_scale.at(i)<=1e6);
  }

  m_elm_data.clear();//releases all memory since it is SmallVector
  m_elm_data.reserve_hint(elm_bixs.size());
  for ( auto i : ncrange( elm_msd.size() ) )
    m_elm_data.emplace_back( elm_msd[i], elm_bixs[i]*elm_scale[i] );
}

NC::CrossSect NC::ElIncXS::evaluate(NeutronEnergy ekin) const
{
  //NB: The cross-section code here must be consistent with code in
  //evaluateMany().
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  double e = kkk*ekin.dbl();
  double xs = 0.0;
  for ( auto& elmdata : m_elm_data )
    xs = ncfma( elmdata.second, eval_1mexpmtdivt( elmdata.first * e ), xs );
  return CrossSect{ xs };
}

void NC::ElIncXS::evaluateMany( Span<const double> ekin, Span<double> tgt) const
{
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  nc_assert(ekin.size()==tgt.size());
  const std::size_t n = tgt.size();

  constexpr auto nbuf = 2048;
  double buf[nbuf];
  double buf_e[nbuf];

  const double * ekin_ptr = ekin.data();
  double * tgt_ptr = tgt.data();

  std::fill( tgt_ptr, tgt_ptr + n, 0.0 );

  std::size_t nleft = n;
  while ( nleft ) {
    std::size_t nstep = std::min<std::size_t>(nleft,nbuf);
    for( auto i : ncrange(nstep) )
      buf_e[i] = kkk * ekin_ptr[i];
    for ( auto& elmdata : m_elm_data ) {
      for( auto i : ncrange(nstep) )
        buf[i] = elmdata.first * buf_e[i];
      for( auto i : ncrange(nstep) )
        buf[i] = eval_1mexpmtdivt( buf[i] );//NB: One day we want to vectorize
                                            //this!
      for( auto i : ncrange(nstep) )
        tgt_ptr[i] = ncfma( elmdata.second, buf[i], tgt_ptr[i] );
    }
    tgt_ptr += nstep;
    ekin_ptr += nstep;
    nleft -= nstep;
  }
}

NC::ElIncXS::ElIncXS( const ElIncXS& a, double scale_a, const ElIncXS& b, double scale_b )
{
  //m_elm_data is (msd,boundincohxs*scale)
  nc_assert(scale_a>=0.0);
  using ElmDataVector = decltype(m_elm_data);
  ElmDataVector tmp;
  tmp.reserve_hint( a.m_elm_data.size() + b.m_elm_data.size() );
  auto addData = [&tmp](const ElmDataVector& data, double scale )
  {
    for ( auto& e : data ) {
      double strength = e.second * scale;
      if (strength)
        tmp.emplace_back( e.first, strength );
    }
  };
  addData(a.m_elm_data,scale_a);
  addData(b.m_elm_data,scale_b);
  std::sort(tmp.begin(),tmp.end());
  ElmDataVector result;
  result.reserve_hint(tmp.size());

  for ( const auto& e : tmp ) {
    if ( result.empty() || !floateq(result.back().first,e.first,1e-15) )
      result.emplace_back( e );
    else
      result.back().second += e.second;
  }
  result.shrink_to_fit();
  std::swap(m_elm_data,result);
}

std::size_t NC::ElIncXS::nElements() const
{
  return m_elm_data.size();
}
