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

#include "NCrystal/internal/extd_utils/NCSABAnalyser.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/core/NCMem.hh"
#include <mutex>
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABAnalyser {
    namespace {

      //Per-row validity gates. These are NOT tuning parameters: each
      //value follows from an error-propagation or convergence argument
      //(and the floor pair from how kernel files encode placeholder
      //floors), and a sensitivity scan varying each x2 up/down (orders
      //of magnitude for the floors) over 263 ENDF8-, VDOS- and
      //tests/data-derived kernels moved accepted results by <<0.1%:
      constexpr double opt_alpha_min = 0.1;       //1/alpha blowup guard
      constexpr double opt_trunc_relbias = 1e-4;  //max est. truncation bias
      constexpr unsigned opt_resolution_minpts = 8;//points per +-2sigma
      constexpr double opt_floor_relmin = 10.0;   //S<=relmin*min(S>0) or...
      constexpr double opt_floor_relmax = 1e-30;  //...S<=relmax*max(S)
      constexpr double opt_msd_m0_low = 0.05;     //1/N and 1/(1-N) error
      constexpr double opt_msd_m0_high = 0.95;    //  amplification guards
      constexpr double opt_center_tol = 0.05;     //first moment sanity check
                                                  //(good kernels: <0.3%)
      //msd rows: require the trapezoidal m0 to be converged, tested by
      //halving the beta resolution (error ratio 4:1), relative to the
      //deficit 1-m0 which sets the msd error scale. Catches coarse
      //kernels whose sharp small-alpha quasi-elastic peak is
      //under-sampled (bias invisible to the spread metric):
      constexpr double opt_msd_m0_convergence = 0.15;

      //Median and relative IQR spread of a scratch vector (modified):
      double calcMedian( VectD& v )
      {
        nc_assert( !v.empty() );
        std::sort( v.begin(), v.end() );
        const std::size_t n = v.size();
        return ( n % 2 == 1
                 ? vectAt( v, n/2 )
                 : 0.5 * ( vectAt( v, n/2 - 1 ) + vectAt( v, n/2 ) ) );
      }

      double calcRelIQR( VectD& v, double median )
      {
        //v must already be sorted; interpolation-free quartiles are
        //fine for a robustness metric:
        nc_assert( !v.empty() && median != 0.0 );
        const std::size_t n = v.size();
        const double q1 = vectAt( v, n/4 );
        const double q3 = vectAt( v, (3*n)/4 );
        return ( q3 - q1 ) / median;
      }

      //Zeroth beta-moment of a row using only every second retained
      //point, for the msd trapezoid-convergence self-test (an
      //under-resolved quasi-elastic peak converges badly):
      double calcRowM0Stride2( const VectD& beta, const VectD& sab,
                               std::size_t nalpha, std::size_t ia,
                               double floor )
      {
        StableSum m0;
        double bprev(0.0), wprev(0.0);
        unsigned nkept(0);
        bool have_prev(false), take(true);
        const std::size_t nb = beta.size();
        for ( std::size_t j = 0; j < nb; ++j ) {
          const double w = vectAt( sab, j*nalpha + ia );
          if ( !( w > floor ) )
            continue;
          if ( !take ) {
            take = true;
            continue;
          }
          take = false;
          ++nkept;
          const double b = vectAt( beta, j );
          if ( have_prev )
            m0.add( 0.5*( b - bprev )*( w + wprev ) );
          bprev = b; wprev = w; have_prev = true;
        }
        return nkept >= 2 ? m0.sum() : -1.0;
      }

      //Integrals over [b0,b1] of (b-c)^2*S(b) with S(b) linear between w0 and w1
      //(d0=b0-c, d1=b1-c), i.e. the second moment of the kernel as it is
      //interpolated between grid points:
      double segmentM2( double b0, double b1, double w0, double w1,
                        double d0, double d1 )
      {
        const double dd = d0*d1;
        return ( b1 - b0 ) * ( w0*( 3.0*d0*d0 + 2.0*dd + d1*d1 )
                               + w1*( d0*d0 + 2.0*dd + 3.0*d1*d1 ) ) / 12.0;
      }

      //Row moments over the floor-masked polyline (exact integrals of the
      //linear-in-beta interpolation, so independent of how densely a given
      //curve is sampled):
      struct RowMoments {
        double m0 = -1.0, mean = -1.0, var = -1.0;
        unsigned npts = 0;
        unsigned npts_core = 0;//points within +-2 std.dev. of mean
        //Raw (not m0-normalised) first moment, and the variance of the full
        //S(alpha,beta) incl. the elastic delta function of weight 1-m0 at
        //beta=0, from the second moment about the recoil centre mu:
        double raw_m1 = 0.0, full_var = -1.0, mu = 0.0;
        //Estimated contributions to m0 and to the central second moment
        //(unnormalised) of the parts of the row beyond the beta grid, from
        //exponential extrapolation of the outermost two points where the
        //data reaches a grid end (if not decaying there, e.g. noise near the
        //floor, the edge value is assumed to continue for a full grid width):
        double tail_m0 = 0.0, tail_m2 = 0.0;
      };

      //Add the tail estimate beyond grid end point (b0,w0), with (b1,w1)
      //the next point inwards:
      void addTail( RowMoments& rm, double b0, double w0,
                    double b1, double w1, double gridwidth )
      {
        const double lambda = ( w0 < w1
                                ? ncmin( gridwidth, ncabs( b1 - b0 )
                                         / std::log( w1 / w0 ) )
                                : gridwidth );
        const double s = ( b0 < b1 ? -1.0 : 1.0 );//outward direction
        const double d0 = b0 - rm.mu;
        rm.tail_m0 += w0 * lambda;
        rm.tail_m2 += w0 * lambda * ( d0*d0 + 2.0*s*d0*lambda
                                      + 2.0*lambda*lambda );
      }

      RowMoments calcRowMoments( const VectD& beta, const VectD& sab,
                                 std::size_t nalpha, std::size_t ia,
                                 double floor, double mu )
      {
        RowMoments res;
        res.mu = mu;
        //First pass: m0 and m1 over the masked polyline:
        StableSum m0, m1;
        double bprev(0.0), wprev(0.0);
        bool have_prev = false;
        const std::size_t nb = beta.size();
        std::size_t j_first(nb), j_second(nb), j_last(nb), j_prevlast(nb);
        for ( std::size_t j = 0; j < nb; ++j ) {
          const double w = vectAt( sab, j*nalpha + ia );
          if ( !( w > floor ) )
            continue;
          const double b = vectAt( beta, j );
          if ( j_first == nb )
            j_first = j;
          else if ( j_second == nb )
            j_second = j;
          j_prevlast = j_last;
          j_last = j;
          if ( res.npts ) {
            const double db = b - bprev;
            m0.add( 0.5*db*( w + wprev ) );
            m1.add( db*( wprev*( 2.0*bprev + b ) + w*( bprev + 2.0*b ) )
                    / 6.0 );
          }
          bprev = b;
          wprev = w;
          have_prev = true;
          ++res.npts;
        }
        (void)have_prev;
        if ( res.npts < 10 || !( m0.sum() > 0.0 ) )
          return res;
        res.m0 = m0.sum();
        res.raw_m1 = m1.sum();
        res.mean = m1.sum() / res.m0;
        //Second pass: central second moment, second moment about mu:
        StableSum m2, m2mu;
        have_prev = false;
        for ( std::size_t j = 0; j < nb; ++j ) {
          const double w = vectAt( sab, j*nalpha + ia );
          if ( !( w > floor ) )
            continue;
          const double b = vectAt( beta, j );
          const double d = b - res.mean;
          const double dmu = b - mu;
          if ( have_prev ) {
            const double dprev = bprev - res.mean;
            const double dmuprev = bprev - mu;
            m2.add( segmentM2( bprev, b, wprev, w, dprev, d ) );
            m2mu.add( segmentM2( bprev, b, wprev, w, dmuprev, dmu ) );
          }
          bprev = b;
          wprev = w;
          have_prev = true;
        }
        res.var = m2.sum() / res.m0;
        const double dm1 = res.raw_m1 - mu;
        res.full_var = m2mu.sum() - dm1*dm1 + mu*mu*( 1.0 - res.m0 );
        //Core-resolution count (points within +-2 std.dev. of the mean):
        const double sd = std::sqrt( ncmax( 0.0, res.var ) );
        for ( std::size_t j = 0; j < nb; ++j ) {
          const double w = vectAt( sab, j*nalpha + ia );
          if ( w > floor && ncabs( vectAt( beta, j ) - res.mean ) < 2.0*sd )
            ++res.npts_core;
        }
        //Tails beyond the grid ends (only where data reaches them):
        auto sv = [&sab,nalpha,ia](std::size_t j)
        { return vectAt( sab, j*nalpha + ia ); };
        const double gw = beta.back() - beta.front();
        if ( j_first == 0 )
          addTail( res, vectAt( beta, j_first ), sv( j_first ),
                   vectAt( beta, j_second ), sv( j_second ), gw );
        if ( j_last + 1 == nb )
          addTail( res, vectAt( beta, j_last ), sv( j_last ),
                   vectAt( beta, j_prevlast ), sv( j_prevlast ), gw );
        return res;
      }

      //Linear (first) and 4-point Lagrange cubic (second) interpolation of
      //the curve (x,y) at xx, with x strictly ascending, at least 4 points,
      //and x.front()<xx<x.back():
      PairDD interpLinCubic( const VectD& x, const VectD& y, double xx )
      {
#ifndef NDEBUG
        nc_assert_always( x.size() == y.size() && x.size() >= 4 );
        nc_assert_always( xx > x.front() && xx < x.back() );
#endif
        const std::size_t n = x.size();
        const std::size_t i = static_cast<std::size_t>
          ( std::upper_bound( x.begin(), x.end(), xx ) - x.begin() );
        const double x0 = vectAt( x, i-1 );
        const double y0 = vectAt( y, i-1 );
        const double lin = y0 + ( xx - x0 ) / ( vectAt( x, i ) - x0 )
          * ( vectAt( y, i ) - y0 );
        const std::size_t i0 = ( i < 2 ? 0 : ncmin( i - 2, n - 4 ) );
        double cub = 0.0;
        for ( auto k : ncrange( i0, i0 + 4 ) ) {
          double l = 1.0;
          for ( auto m : ncrange( i0, i0 + 4 ) )
            if ( m != k )
              l *= ( xx - vectAt( x, m ) ) / ( vectAt( x, k ) - vectAt( x, m ) );
          cub += l * vectAt( y, k );
        }
        return { lin, cub };
      }

      //Detailed balance deviation of an alpha row. Detailed balance means
      //log(S(alpha,-x)) = x + log(S(alpha,x)). The two curves are compared
      //at the x>0 grid points of each, interpolating the other curve both
      //linearly and with a cubic and keeping the smaller relative
      //deviation, to limit the impact of interpolation errors. Returns the
      //integral of S times that deviation, relative to the integral of S
      //over the whole row (so deviations in far tails count little), or -1
      //if the curves do not overlap. Only points above the floor are used.
      double calcRowDetailedBalanceDev( const VectD& beta, const VectD& sab,
                                        std::size_t nalpha, std::size_t ia,
                                        double floor )
      {
#ifndef NDEBUG
        nc_assert_always( ia < nalpha );
        nc_assert_always( sab.size() == nalpha * beta.size() );
#endif
        VectD xa, ya, xb, yb;//curves (-beta,log(S)) and (beta,beta+log(S))
        for ( auto j : ncrange( beta.size() ) ) {
          const double w = vectAt( sab, j*nalpha + ia );
          if ( !( w > floor ) )
            continue;
          const double b = vectAt( beta, j );
          const double lw = std::log( w );
          xa.push_back( -b );
          ya.push_back( lw );
          xb.push_back( b );
          yb.push_back( b + lw );
        }
        if ( xa.size() < 4 )
          return -1.0;
        std::reverse( xa.begin(), xa.end() );
        std::reverse( ya.begin(), ya.end() );
        const double lo = ncmax( 0.0, ncmax( xa.front(), xb.front() ) );
        const double hi = ncmin( xa.back(), xb.back() );
        if ( !( lo < hi ) )
          return -1.0;
        //Trapezoidal integrals, nodes weighted by half the adjacent spacings:
        auto trapw = []( const VectD& x, std::size_t i )
        {
          const std::size_t n = x.size();
          return 0.5 * ( vectAt( x, i + 1 < n ? i + 1 : i )
                         - vectAt( x, i > 0 ? i - 1 : i ) );
        };
        StableSum total;
        for ( auto i : ncrange( xa.size() ) )
          total.add( trapw( xa, i ) * std::exp( vectAt( ya, i ) ) );
        StableSum num;
        auto compare = [lo,hi,&num,&trapw]( const VectD& xn, const VectD& yn,
                                            const VectD& xo, const VectD& yo )
        {
          for ( auto i : ncrange( xn.size() ) ) {
            const double x = vectAt( xn, i );
            if ( !( x > lo && x < hi ) )
              continue;
            const double y = vectAt( yn, i );
            const PairDD o = interpLinCubic( xo, yo, x );
            const double d = ncmin( ncabs( std::expm1( o.first - y ) ),
                                    ncabs( std::expm1( o.second - y ) ) );
            num.add( trapw( xn, i ) * std::exp( y ) * d );
          }
        };
        compare( xa, ya, xb, yb );
        compare( xb, yb, xa, ya );
        //Each side covers the overlap once, hence the factor 2:
        return num.sum() / ( 2.0 * total.sum() );
      }

    }
  }
}

NC::SABAnalyser::Result
NC::SABAnalyser::analyse( const SABData& sab, const Options& opt )
{
  const auto& agrid = sab.alphaGrid();
  const auto& bgrid = sab.betaGrid();
  const auto& S = sab.sab();
  const std::size_t na = agrid.size();
  nc_assert_always( S.size() == na * bgrid.size() );
  const double T = sab.temperature().dbl();
  const double kT = sab.temperature().kT();
  const double massratio
    = sab.elementMassAMU().dbl() / const_neutron_mass_amu;

  Result res;

  //Global floor-marker level (placeholder values in far wings are not
  //data):
  double smin_pos = kInfinity, smax = 0.0;
  for ( auto v : S ) {
    if ( v > 0.0 )
      smin_pos = ncmin( smin_pos, v );
    smax = ncmax( smax, v );
  }
  if ( !( smax > 0.0 ) )
    return res;//empty/degenerate kernel
  const double floor = ncmax( smin_pos * opt_floor_relmin,
                              smax * opt_floor_relmax );

  std::shared_ptr<Diagnostics> diag;
  if ( opt.collect_diagnostics ) {
    diag = std::make_shared<Diagnostics>();
    diag->alpha = agrid;
    diag->m0.resize( na, -1.0 );
    diag->mean.resize( na, -1.0 );
    diag->variance.resize( na, -1.0 );
    diag->teff_row.resize( na, -1.0 );
    diag->msd_row.resize( na, -1.0 );
    diag->detailed_balance_dev.resize( na, -1.0 );
    diag->teff_row_status.resize( na, RowStatus::NoData );
    diag->msd_row_status.resize( na, RowStatus::NoData );
    diag->floor_value = floor;
  }

  //msd constant: 2W = alpha * alpha2x_factor * kT * msd (exactly as in
  //NCVDOSExpand.cc):
  constexpr double alpha2x_factor
    = 2.0*const_neutron_mass_evc2/(constant_hbar*constant_hbar);

  VectD teff_rows, msd_rows, center_rows;
  teff_rows.reserve( na );
  msd_rows.reserve( na );
  center_rows.reserve( na );

  for ( auto ia : ncrange( na ) ) {
    const double a = vectAt( agrid, ia );
    auto rm = calcRowMoments( bgrid, S, na, ia, floor, -a / massratio );
    if ( diag != nullptr )
      vectAt( diag->detailed_balance_dev, ia )
        = calcRowDetailedBalanceDev( bgrid, S, na, ia, floor );
    auto setstat = [&diag,ia]( RowStatus st, bool for_teff )
    {
      if ( diag != nullptr ) {
        if ( for_teff )
          vectAt( diag->teff_row_status, ia ) = st;
        else
          vectAt( diag->msd_row_status, ia ) = st;
      }
    };
    if ( rm.m0 < 0.0 ) {
      setstat( RowStatus::NoData, true );
      setstat( RowStatus::NoData, false );
      continue;
    }
    res.max_row_m0 = ncmax( res.max_row_m0, rm.m0 );
    if ( diag != nullptr ) {
      vectAt( diag->m0, ia ) = rm.m0;
      vectAt( diag->mean, ia ) = rm.mean;
      vectAt( diag->variance, ia ) = rm.var;
    }
    //msd row:
    if ( !( a > opt_alpha_min ) ) {
      setstat( RowStatus::AlphaTooLow, false );
    } else if ( !( rm.m0 > opt_msd_m0_low
                   && rm.m0 < opt_msd_m0_high ) ) {
      setstat( RowStatus::M0OutOfWindow, false );
    } else if ( !( rm.tail_m0 <= opt_trunc_relbias * ( 1.0 - rm.m0 )
                   * ( -std::log1p( -rm.m0 ) ) ) ) {
      //relative msd bias is tail_m0/((1-m0)*(-log(1-m0))):
      setstat( RowStatus::WingsNotCovered, false );
    } else if ( ncabs( calcRowM0Stride2( bgrid, S, na, ia, floor )
                       - rm.m0 )
                > opt_msd_m0_convergence * ( 1.0 - rm.m0 ) ) {
      setstat( RowStatus::Unresolved, false );
    } else {
      const double c = -std::log1p( -rm.m0 ) / a;
      const double msd_row = c / ( alpha2x_factor * kT );
      msd_rows.push_back( msd_row );
      setstat( RowStatus::Ok, false );
      if ( diag != nullptr )
        vectAt( diag->msd_row, ia ) = msd_row;
    }
    //Teff row:
    if ( !( a > opt_alpha_min ) ) {
      setstat( RowStatus::AlphaTooLow, true );
      continue;
    }
    if ( !( rm.m0 >= opt_msd_m0_high ) ) {
      //elastic part negligible (stored kernels, e.g. ENDF liquids or ice,
      //can violate the sum rules at lower alpha):
      setstat( RowStatus::M0OutOfWindow, true );
      continue;
    }
    if ( !( rm.full_var > 0.0
            && rm.tail_m2 <= opt_trunc_relbias * rm.full_var ) ) {
      setstat( RowStatus::WingsNotCovered, true );
      continue;
    }
    if ( rm.npts_core < opt_resolution_minpts ) {
      setstat( RowStatus::Unresolved, true );
      continue;
    }
    //First moment sum rule (raw first moment is exactly -alpha/A):
    const double c_row = rm.raw_m1 * massratio / ( -a );
    if ( !( ncabs( c_row - 1.0 ) <= opt_center_tol ) ) {
      setstat( RowStatus::CenterOffRecoil, true );
      continue;
    }
    //Second moment sum rule (full variance is 2*alpha*Teff/(A*T)):
    const double teff_row = massratio * T * rm.full_var / ( 2.0 * a );
    teff_rows.push_back( teff_row );
    center_rows.push_back( c_row );
    setstat( RowStatus::Ok, true );
    if ( diag != nullptr )
      vectAt( diag->teff_row, ia ) = teff_row;
  }

  res.teff_nrows = static_cast<unsigned>( teff_rows.size() );
  res.msd_nrows = static_cast<unsigned>( msd_rows.size() );
  if ( !teff_rows.empty() ) {
    const double med = calcMedian( teff_rows );//NB: sorts
    res.teff = Temperature{ med };
    res.teff_relspread = calcRelIQR( teff_rows, med );
    res.recoil_center_ratio = calcMedian( center_rows );
  }
  if ( !msd_rows.empty() ) {
    const double med = calcMedian( msd_rows );//NB: sorts
    res.msd = med;
    res.msd_relspread = calcRelIQR( msd_rows, med );
  }
  res.diagnostics = std::move( diag );
  return res;
}

namespace NCRYSTAL_NAMESPACE {
  namespace SABAnalyser {
    namespace {

      //Acceptance policies for the turn-key estimates, per quantity since the
      //two consumers have opposite failure asymmetries:
      //
      //Teff: refusal is NOT neutral - the consumer then extends with free gas
      //at T, an error of (Teff_true-T), usually far larger than any credible
      //estimate's error (e.g. H in 100K solid ethanol: fallback ~12x off vs a
      //refused 3% estimate). Hence a loose policy: over 501 pooled
      //truth-records (ENDF8- and JENDL5-converted, stdlib and tests/data
      //kernels incl.  temperature extremes), nrows>=10 && relspread<=0.1 never
      //accepted an estimate worse than the fallback, worst error 5.3% (p95 well
      //under 1%); garbage only appears below the row floor (140% error at
      //nrows=2).
      //
      //msd: refusal IS the status quo (no elastic component), while a wrong
      //value adds wrong physics, so strictness is correct. The 0.03 spread cut
      //refuses the self-flagged tail, e.g. U in UO2 at 5K (spread 0.035, error
      //3%), for a worst accepted error of 0.77% (the msd values remain
      //provisional pending ENDF-grade ground truth). Kernels with a row
      //integrating to more than 1+policy_msd_max_row_m0_excess are refused:
      //their beta grids do not resolve all peaks (e.g. crude VDOS expansions,
      //or rotating molecules), which also makes the m0 values used for msd
      //unreliable:
      constexpr unsigned policy_teff_min_nrows = 10;
      constexpr double policy_teff_max_relspread = 0.1;
      constexpr unsigned policy_msd_min_nrows = 20;
      constexpr double policy_msd_max_relspread = 0.03;
      constexpr double policy_msd_sanity_max = 10.0;//[Aa^2]
      constexpr double policy_msd_max_row_m0_excess = 0.01;

      struct TeffMSDCacheDB {
        std::mutex mtx;
        //Both quantities are strictly positive, so -1.0 encodes
        //absence:
        struct Entry { double teff, msd; };
        std::map<std::uint64_t,Entry> map;
        bool cleanup_registered = false;
      };

      TeffMSDCacheDB& getTeffMSDCacheDB()
      {
        static TeffMSDCacheDB db;
        return db;
      }

      TeffMSD decodeEntry( const TeffMSDCacheDB::Entry& e )
      {
        TeffMSD res;
        if ( e.teff >= 0.0 )
          res.effectiveTemperature = Temperature{ e.teff };
        if ( e.msd >= 0.0 )
          res.msd = e.msd;
        return res;
      }

    }
  }
}

NC::SABAnalyser::TeffMSD
NC::SABAnalyser::estimateTeffMSD( const SABData& sab )
{
  auto& db = getTeffMSDCacheDB();
  const std::uint64_t uid = sab.getUniqueID().value;
  {
    NCRYSTAL_LOCK_GUARD(db.mtx);
    if ( !db.cleanup_registered ) {
      db.cleanup_registered = true;
      registerCacheCleanupFunction( []()
      {
        auto& thedb = getTeffMSDCacheDB();
        NCRYSTAL_LOCK_GUARD(thedb.mtx);
        thedb.map.clear();
      });
    }
    auto it = db.map.find( uid );
    if ( it != db.map.end() )
      return decodeEntry( it->second );
  }
  //Not cached; analyse outside the lock (a racing duplicate
  //computation is benign):
  auto r = analyse( sab );
  TeffMSD res;
  const double T = sab.temperature().dbl();
  if ( r.teff.has_value()
       && r.teff_nrows >= policy_teff_min_nrows
       && r.teff_relspread <= policy_teff_max_relspread ) {
    double tv = r.teff.value().dbl();
    //Same semantics as sct_checkedTeff, but discarding rather than
    //throwing:
    if ( std::isfinite( tv ) && tv >= 0.999*T ) {
      tv = ncmax( tv, T );
      if ( tv >= Temperature::allowed_range.first
           && tv <= Temperature::allowed_range.second )
        res.effectiveTemperature = Temperature{ tv };
    }
  }
  if ( r.msd.has_value()
       && r.msd_nrows >= policy_msd_min_nrows
       && r.msd_relspread <= policy_msd_max_relspread
       && r.max_row_m0 <= 1.0 + policy_msd_max_row_m0_excess
       && std::isfinite( r.msd.value() )
       && r.msd.value() > 0.0
       && r.msd.value() < policy_msd_sanity_max )
    res.msd = r.msd.value();
  {
    NCRYSTAL_LOCK_GUARD(db.mtx);
    //Reset if somehow reaching an unreasonable size, on principle:
    if ( db.map.size() >= 5000 )
      db.map.clear();
    db.map[uid]
      = TeffMSDCacheDB::Entry{ ( res.effectiveTemperature.has_value()
                                 ? res.effectiveTemperature.value().dbl()
                                 : -1.0 ),
                               res.msd.has_value() ? res.msd.value()
                                                   : -1.0 };
  }
  return res;
}
