#ifndef NCrystal_SABUtils_hh
#define NCrystal_SABUtils_hh

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

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/sab/NCScatKnlData.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    //Transform ScatKnlData to the standard unscaled and asymmetric
    //S(alpha,beta) format. Internally this calls validateScatKnlData(...) on
    //the input and, in debug builds only, also on the output. Thus, creating
    //SABData objects via this function alleviates the need for littering code
    //elsewhere with validate calls. Parts of the grid not filled by the input
    //kernel are trimmed away, with a WARNING message in case anything is
    //trimmed.
    SABData transformKernelToStdFormat(ScatKnlData&&);

    //If beta grid is defined as [0,b1,b2,..,bn] it is assumed that it is a
    //space-saving shorthand for [-bn,...,-b2,-b1,0,b1,b2,..,bn]. Expand, using
    //S(alpha,-beta)=S(alpha,beta):
    void expandBetaAndSABToAllBetas( Span<const double> halfbetagrid,
                                     Span<const double> alphagrid,
                                     Span<const double> sab_for_halfbetagrid,
                                     VectD& complete_betagrid,
                                     VectD& complete_sab );

    //Access S(alpha|beta_i) slice (const and mutable version):
    Span<const double> sliceSABAtBetaIdx( Span<const double> sab, std::size_t nalpha, std::size_t beta_idx);
    Span<double> sliceSABAtBetaIdxMutable( Span<double> sab, std::size_t nalpha, std::size_t beta_idx);

    //Or just calculate the index:
    std::size_t calcSABIdx( std::size_t nalpha, std::size_t alpha_idx, std::size_t beta_idx);

    //Interpolation between (a,fa) and (b,fb) at x in [a,b], exact at a and b.
    //interpLin is linear, interpLogLinWCL linear in log(f) (requires fa,fb>0)
    //and interpLogLinOrLinWCL is log-linear, or linear if fa or fb is 0. WCL
    //(with cached logs) versions take logfa=log(fa), logfb=log(fb) and return
    //{f,log(f)} (log(f)=0 when linear):
    double interpLin(double a, double fa, double b, double fb, double x);
    PairDD interpLogLinWCL(double a, double fa, double b, double fb, double x, double logfa, double logfb);
    PairDD interpLogLinOrLinWCL(double a, double fa, double b, double fb, double x, double logfa, double logfb);

    //S(alpha,beta) at (alpha,beta) inside grid cell (ia,ib), i.e.
    //alphaGrid[ia]<=alpha<=alphaGrid[ia+1] and similarly for beta: log-linear
    //in alpha (linear if a corner is 0), then linear in beta. Exact at the
    //corners, and robust for any finite S values (incl. subnormal ones). The
    //WCL version takes the createLogSTable(sd) table, the other computes the
    //logs of the 4 corners:
    double evalSABCellWCL( const SABData& sd, const VectD& logS,
                           std::size_t ia, std::size_t ib,
                           double alpha, double beta );
    double evalSABCell( const SABData& sd, std::size_t ia, std::size_t ib,
                        double alpha, double beta );

    //Same, given the cell edges and the S and log(S) values s[i], ls[i] at
    //the 4 corners (alpha,beta) = (a0,b0), (a1,b0), (a0,b1), (a1,b1):
    double evalSABCellCorners( double a0, double a1, double b0, double b1,
                               const double* s, const double* ls,
                               double alpha, double beta );

    //S(alpha,beta) of the table at any point, with the interpolation used for
    //kernels (log-linear in alpha, linear if a corner value is 0, then linear
    //in beta, as evalSABCell). Returns 0 outside the grid. Locates the cell by
    //binary search at each call, so not efficient for many evaluations:
    double evalSABOneShot( const SABData&, double alpha, double beta );

    //The next two functions can be used to evaluate a SAB integral over a bin
    //of alpha-values, assuming a linear relationship in log(s).  The analytical
    //solution when s1,s2 are both nonzero is (a2-a1)*(s2-s1)/log(s2/s1), but
    //the evaluation will take performance and numerical precision issues into
    //account. In case either s1 or s2 is zero, the function reverts to a simple
    //trapezoidal integral.

    //Version without pre-calculated logs. Returns full (better than 1e-15) precision results.
    double integrateLogLinOrLin(double a1,double s1, double a2 , double s2 );

    //Faster version with cached logs (WCL). Returns slightly reduced (better than 1e-14) precision results.
    double integrateLogLinOrLinWCL(double a1,double s1, double a2 , double s2 , double logs1, double logs2);

    //log(s) if s>0.0 else valueIfNonPositive. Use instead of the naive
    //ternary: unsafe under clang, see definition below.
    double safeLogOrElse( double s, double valueIfNonPositive );

    //safeLogOrElse(S,0.0) for all values of the S(alpha,beta) table (same
    //indexing as SABData::sab()), to be shared by all users of a kernel:
    shared_obj<const VectD> createLogSTable( const SABData& );


  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::Span<const double> NCrystal::SABUtils::sliceSABAtBetaIdx( NCrystal::Span<const double> sab, std::size_t nalpha, std::size_t beta_idx)
{
  nc_assert( (beta_idx+1) * nalpha <= static_cast<std::size_t>(sab.size()) );
  return {sab.begin() + beta_idx*nalpha, sab.begin() + (beta_idx+1)*nalpha};
}

inline NCrystal::Span<double> NCrystal::SABUtils::sliceSABAtBetaIdxMutable( NCrystal::Span<double> sab, std::size_t nalpha, std::size_t beta_idx)
{
  nc_assert( (beta_idx+1) * nalpha <= static_cast<std::size_t>(sab.size()) );
  return {sab.begin() + beta_idx*nalpha, sab.begin() + (beta_idx+1)*nalpha};
}

inline std::size_t NCrystal::SABUtils::calcSABIdx( std::size_t nalpha, std::size_t alpha_idx, std::size_t beta_idx)
{
  return beta_idx * nalpha + alpha_idx;
}

inline double NCrystal::SABUtils::interpLin(double a, double fa,
                                            double b, double fb,
                                            double x)
{
  nc_assert( fa >= 0.0 && fb >= 0.0 );
  nc_assert( x >= a && x <= b );
  nc_assert( std::isfinite(fa) && std::isfinite(fb) );
  const double bma = b - a;
  nc_assert( bma > 0.0 );
  const double midpoint = 0.5 * ( b + a );
  if ( x < midpoint ) {
    //choose forms most numerically stable for x near a
    const double r = (x-a) / bma;
    return nclerp(fa,fb,r);
  } else {
    //choose forms most numerically stable for x near b
    const double s = (b-x) / bma;
    return nclerp(fb,fa,s);
  }
}

inline NCrystal::PairDD
NCrystal::SABUtils::interpLogLinWCL(double a, double fa,
                                    double b, double fb, double x,
                                    double logfa, double logfb)
{
  nc_assert( fa > 0.0 && fb > 0.0 );
  nc_assert( x >= a && x <= b );
  nc_assert( std::isfinite(logfa) && std::isfinite(logfb) );
  nc_assert( std::isfinite(fa) && std::isfinite(fb) );
  const double bma = b - a;
  nc_assert( bma > 0.0 );
  const double midpoint = 0.5 * ( b + a );
  auto flexexp = [](double xx)
  {
    //exp_smallarg_approx has 7th order Taylor, which might be faster than a
    //library call. FIXME: std::exp might differ slightly across platforms
    //(libm implementations), consider stable_exp (check speed):
    return ncabs(xx)<0.02 ? exp_smallarg_approx(xx) : std::exp(xx);
  };
  //f = f0*exp(dlf), or exp(log f) when exp(dlf) could overflow or become
  //subnormal (only for extreme ratios fb/fa, above ~e^1400):
  auto fval = [&flexexp]( double f0, double dlf, double logf )
  {
    if ( ncabs( dlf ) < 700.0 )
      return f0 * flexexp( dlf );
    return std::exp( logf );
  };
  PairDD res;
  if ( x < midpoint ) {
    //choose forms most numerically stable (and with smallest flexexp arg) for x
    //near a
    const double r = (x-a) / bma;
    const double dl = logfb-logfa;
    const double dlf = dl*r;
    res.second = ncfma( dl, r, logfa );
    res.first = fval( fa, dlf, res.second );
  } else {
    //choose forms most numerically stable (and with smallest flexexp arg) for x
    //near b
    const double s = (b-x) / bma;
    const double dl = logfa-logfb;
    const double dlf = dl*s;
    res.second = ncfma( dl, s, logfb );
    res.first = fval( fb, dlf, res.second );
  }
  return res;
}

inline NCrystal::PairDD
NCrystal::SABUtils::interpLogLinOrLinWCL( double a, double fa, double b,
                                          double fb, double x,
                                          double logfa, double logfb )
{
  nc_assert( fa >= 0.0 && fb >= 0.0 );
  if ( fa > 0.0 && fb > 0.0 )
    return interpLogLinWCL( a, fa, b, fb, x, logfa, logfb );
  return { interpLin( a, fa, b, fb, x ), 0.0 };
}

inline double NCrystal::SABUtils::integrateLogLinOrLin(double a1,double s1, double a2 , double s2 )
{
  nc_assert( s1 >= 0.0 && s2 >= 0.0 );
  nc_assert( a2 > a1 );

  //With ds=s2-s1, ps=s1+s2, da=a2-a1, y=ds/ps, one gets 2*s2=ps+ds,
  //2*s1=ps-ds and thus the analytical result can be expressed as:
  //  da*ps*y/log( (ps+ds) / (ps-ds) ) = da*ps* [ y/log((1+y)/(1-y)) ]
  //This is an even function in y. As s1,s2>0, |y|<=1. When |y| is
  //sufficiently small, the expression can be evaluated via a Taylor expansion:
  //
  // y/log((1+y)/(1-y)) = [1/2]-[1/6]y^2-[2/45]y^4-[22/945]y^6-[214/14175]y^8-[5098/467775]y^10-[5359534/638512875]y^12
  //
  //The Taylor expansion is not only faster, as it avoids a log call, but it
  //is also numerically safer when s1~=s2.
  const double da = a2-a1;
  const double ps = s1+s2;
  const double ds = s2-s1;
  if ( ncabs(ds) < 0.10*ps ) {
    //both s1 and s2 are nonzero and s1 ~= s2 so we can use the Taylor expansion
    nc_assert( ps > 0. );
    const double y = ds / ps;
    const double c1           = 0.166666666666666666666666666666666666666666667   ; // = 1/6
    const double c2           = 0.0444444444444444444444444444444444444444444444  ; // = 2/45
    const double c3           = 0.0232804232804232804232804232804232804232804233  ; // = 22/945
    const double c4           = 0.0150970017636684303350970017636684303350970018  ; // = 214/14175
    const double c5           = 0.0108984020095131206242317353428464539575650687  ; // = 5098/467775
    const double c6           = 0.00839377592816746255899695053134206573360012513 ; // = 5359534/638512875
    const double ysq = y*y;
    //Horner's method with ncfma for extra reproducibility:
    double p = c6;
    p = ncfma( ysq, p, c5 );
    p = ncfma( ysq, p, c4 );
    p = ncfma( ysq, p, c3 );
    p = ncfma( ysq, p, c2 );
    p = ncfma( ysq, p, c1 );
    return da*ps*ncfma( -ysq, p, 0.5 );
  }
  //Evaluate via analytical expression or fall-back to trapezoidal integration:
  return ncmin(s1,s2)<1e-300 ? 0.5*da*ps : da*ds/std::log(s2/s1);
}

inline double NCrystal::SABUtils::safeLogOrElse( double s, double valueIfNonPositive )
{
#if !defined(__clang__)
  return ( s > 0.0 ? std::log(s) : valueIfNonPositive );
#else
  //Workaround for buggy clang behaviour (seen in clang 21.1.8): clang -O3
  //optimises the ternary above into an unconditional log(s), so causing a
  //std::log(0.0) call despite our guard -> FPE. The volatile keeps clang from
  //replacing the log argument with s, and the result is selected on s itself
  //(log(denorm_min) is finite, so it can not be used to detect s<=0):
  volatile double vs = ( s > 0.0 ? s : 1.0 );
  const double lv = std::log( vs );
  return ( s > 0.0 ? lv : valueIfNonPositive );
#endif
}

inline double NCrystal::SABUtils::integrateLogLinOrLinWCL(double a1,double s1, double a2 , double s2 , double logs1, double logs2)
{
  //Similar implementation as integrateLogLinOrLin, but due to the
  //precalculated logs, the Taylor expansion will be used more rarely and
  //with fewer terms (it is not needed for speed reasons, just numerical
  //stability). We also reorder, to put the most likely code first:
  nc_assert(!ncisnan(logs1)&&!ncisnan(logs2));
  nc_assert( s1 >= 0.0 && s2 >= 0.0 );
  nc_assert( a2 > a1 );
  const double da = a2-a1;
  const double ps = s1+s2;
  const double ds = s2-s1;
  if ( ncmin(s1,s2) < 1e-300 )//1e-300, because we want to stay clear of log(s)=-inf issues.
    return 0.5*da*ps;//trapezoidal integration
  if ( ncabs(ds) > 0.006*ps ) {
    nc_assert_always(!ncisinf(logs1)&&!ncisinf(logs2));//should have been caught by checks above
    return da*ds/(logs2-logs1);
  }
  nc_assert( ps > 0. );
  const double y = ds / ps;
  const double c1           = 0.166666666666666666666666666666666666666666667   ; // = 1/6
  const double c2           = 0.0444444444444444444444444444444444444444444444  ; // = 2/45
  const double c3           = 0.0232804232804232804232804232804232804232804233  ; // = 22/945
  const double ysq = y*y;
  //Horner's method with ncfma for extra reproducibility:
  double p = c3;
  p = ncfma( ysq, p, c2 );
  p = ncfma( ysq, p, c1 );
  return da*ps*ncfma( -ysq, p, 0.5 );
}

inline double NCrystal::SABUtils::evalSABCellCorners
( double a0, double a1, double b0, double b1, const double* s,
  const double* ls, double alpha, double beta )
{
  const double e0 = interpLogLinOrLinWCL( a0, s[0], a1, s[1], alpha,
                                          ls[0], ls[1] ).first;
  const double e1 = interpLogLinOrLinWCL( a0, s[2], a1, s[3], alpha,
                                          ls[2], ls[3] ).first;
  return interpLin( b0, e0, b1, e1, beta );
}

inline double NCrystal::SABUtils::evalSABCellWCL( const SABData& sd,
                                                  const VectD& logS,
                                                  std::size_t ia,
                                                  std::size_t ib,
                                                  double alpha, double beta )
{
  const auto& ag = sd.alphaGrid();
  const auto& bg = sd.betaGrid();
  const auto& S = sd.sab();
  const std::size_t na = ag.size();
#ifndef NDEBUG
  nc_assert_always( logS.size() == S.size() );
  nc_assert_always( ia + 1 < na && ib + 1 < bg.size() );
  nc_assert_always( ag[ia] <= alpha && alpha <= ag[ia+1] );
  nc_assert_always( bg[ib] <= beta && beta <= bg[ib+1] );
#endif
  const std::size_t i0 = calcSABIdx( na, ia, ib );
  const std::size_t i1 = calcSABIdx( na, ia, ib + 1 );
  const double s[4] = { vectAt( S, i0 ), vectAt( S, i0 + 1 ),
                        vectAt( S, i1 ), vectAt( S, i1 + 1 ) };
  const double ls[4] = { vectAt( logS, i0 ), vectAt( logS, i0 + 1 ),
                         vectAt( logS, i1 ), vectAt( logS, i1 + 1 ) };
  return evalSABCellCorners( vectAt( ag, ia ), vectAt( ag, ia + 1 ),
                             vectAt( bg, ib ), vectAt( bg, ib + 1 ),
                             s, ls, alpha, beta );
}

#endif
