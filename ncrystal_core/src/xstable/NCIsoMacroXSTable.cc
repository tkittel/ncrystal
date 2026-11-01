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

#include "NCrystal/internal/xstable/NCIsoMacroXSTable.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"

//FIXME: Get back to reproducibility for the IsoMacroXSTable tables: the
//selected table points depend on floating point details, which might differ
//between platforms and compilers (e.g. due to FMA contractions or
//auto-vectorisation).

namespace NC = NCrystal;
namespace NCX = NCrystal::XSTable;

namespace NCRYSTAL_NAMESPACE {
  namespace XSTable {
    namespace {

      //Relative distance from a discontinuity, at which the cross sections
      //just below and above it are evaluated. Discontinuities closer than
      //twice this are treated as a single one.
      constexpr double edge_eps = 1e-10;

      //The table is built with a slightly smaller tolerance than requested,
      //so the final verification at random points (which are not part of the
      //dense sample) has a margin:
      constexpr double tol_build_factor = 0.98;

      //Relative jump below which a known discontinuity is treated as a kink:
      constexpr double kink_tol = 1e-6;

      //Neighbouring dense sample points closer than this (relative) are not
      //refined further. If the tolerance is not met between such points, the
      //cross section has a discontinuity (or an extremely sharp feature) which
      //was not known in advance, and an exception is thrown:
      constexpr double min_rel_spacing = 1e-9;

      //The tolerance is relative to max(xs,abs_floor), with xs in barn per
      //atom (an absolute floor is needed for cross sections vanishing at
      //wavelength 0):
      constexpr double xs_abs_floor = 1e-12;

      //Max number of refinement iterations (each iteration at least halves
      //the spacing of the points where the tolerance is not met):
      constexpr unsigned max_refine_iterations = 100;

      //Leaf processes without a table in their isoXSDescr(), whose cross
      //sections are known to have no features which the table creation could
      //miss: they are smooth, apart from the boundaries of their energy domains
      //(which are included explicitly). Null processes and processes with
      //tables (smooth between their breakpts) need no such validation. Any
      //other process (e.g. a process with absorption resonances, or the UCN
      //processes) gives a warning, and must be studied before being added here.
      bool isValidatedUntabulated( const std::string& name )
      {
        static const std::set<std::string> validated = { "SABScatter", "AbsOOV" };
        return validated.count( name ) > 0;
      }

      //Collect the names of the leaf processes (of ProcComposition trees), and
      //of those among them which need a warning (see isValidatedUntabulated):
      void inspectProcess( const ProcImpl::Process& proc,
                           std::set<std::string>& names,
                           std::set<std::string>& unvalidated )
      {
        auto comp = dynamic_cast<const ProcImpl::ProcComposition*>( &proc );
        if ( comp ) {
          for ( auto& c : comp->components() )
            inspectProcess( *c.process, names, unvalidated );
          return;
        }
        names.insert( proc.name() );
        if ( proc.isNull() )
          return;
        auto iso = dynamic_cast<const ProcImpl::ScatterIsotropicMat*>( &proc );
        if ( iso && iso->isoXSDescr().table.has_value() )
          return;
        if ( !isValidatedUntabulated( proc.name() ) )
          unvalidated.insert( proc.name() );
      }

      //The wavelengths of the known discontinuities and kinks of the cross
      //section of a process: the breakpts of its isoXSDescr() (e.g. the Bragg
      //edges), or otherwise the boundaries of its energy domain:
      void collectBreakpoints( const ProcImpl::Process& proc, VectD& wl )
      {
        IsoXSDescr d;
        if ( auto comp = dynamic_cast<const ProcImpl::ProcComposition*>( &proc ) ) {
          d = comp->isoXSDescr();
        } else if ( auto iso = dynamic_cast<const ProcImpl::ScatterIsotropicMat*>( &proc ) ) {
          d = iso->isoXSDescr();
        } else {
          const auto domain = proc.domain();
          if ( !domain.isNull() )
            d.breakpts = { domain.elow.dbl(), domain.ehigh.dbl() };
        }
        wl.reserve( wl.size() + d.breakpts.size() );
        for ( auto e : d.breakpts )
          if ( e > 0.0 && std::isfinite( e ) )
            wl.push_back( ekin2wl( e ) );
      }

      //Sorted discontinuities within (0,wlmax), with clusters of nearby ones
      //represented by their largest member (for Bragg edges, the cross section
      //drops to the value above the whole cluster there).
      VectD mergeDiscontinuities( VectD d, double wlmax )
      {
        std::sort( d.begin(), d.end() );
        VectD res;
        res.reserve( d.size() );
        for ( auto v : d ) {
          if ( !( v > 0.0 && v < wlmax*(1.0-4*edge_eps) ) )
            continue;
          if ( !res.empty() && v <= res.back() * ( 1.0 + 2*edge_eps ) )
            res.back() = v;//merge with previous
          else
            res.push_back( v );
        }
        return res;
      }

      //The cubic Hermite interpolation between (x0,y0) and (x1,y1), with the
      //slopes s0 and s1 there:
      inline double hermite( double x0, double x1, double y0, double y1,
                             double s0, double s1, double x )
      {
        const double h = x1 - x0;
        const double u = ( x - x0 ) / h;
        const double v = 1.0 - u;
        return ( v * v * ( ( 1.0 + 2.0 * u ) * y0 + u * h * s0 )
                 + u * u * ( ( 3.0 - 2.0 * u ) * y1 - v * h * s1 ) );
      }

      //Try to cover the sample points from i to j with one interval, given the
      //slope s0 at i. The slope s1 at j is chosen by a weighted least squares
      //fit to the interior points, and limited (see below) such that the
      //interval and the following one (at least up to the next sample point)
      //can not become negative: h*s0 >= -3*y0 and h*s1 <= 3*y1 for an interval
      //of width h are sufficient for that (Schmidt and Hess 1988). Returns
      //whether all interior points are reproduced within tol:
      bool fitInterval( const VectD& x, const VectD& y, std::size_t i,
                        std::size_t j, double s0, double tol, double abs_floor,
                        double& s1 )
      {
        nc_assert( j > i && x[j] > x[i] );
        const double h = x[j] - x[i];
        //The interpolation is A(u)+s1*B(u), with B(u)=-u^2*(1-u)*h:
        double num = 0.0, den = 0.0;
        for ( std::size_t k = i + 1; k < j; ++k ) {
          const double u = ( x[k] - x[i] ) / h;
          const double v = 1.0 - u;
          const double A = ( v * v * ( ( 1.0 + 2.0 * u ) * y[i] + u * h * s0 )
                             + u * u * ( 3.0 - 2.0 * u ) * y[j] );
          const double B = -u * u * v * h;
          const double w = 1.0 / ncmax( y[k], abs_floor );
          num += w * w * B * ( y[k] - A );
          den += w * w * B * B;
        }
        s1 = ( den > 0.0 ? num / den : ( y[j] - y[i] ) / h );
        //Keep s1 within the slopes of the sample chords next to j (where the
        //slope of a smooth curve is, unless it has an inflection point there),
        //so the slopes of consecutive intervals can not run away:
        if ( x[j-1] < x[j] ) {
          double clo = ( y[j] - y[j-1] ) / ( x[j] - x[j-1] );
          double chi = clo;
          if ( j + 1 < x.size() && x[j+1] > x[j] ) {
            const double c = ( y[j+1] - y[j] ) / ( x[j+1] - x[j] );
            clo = ncmin( clo, c );
            chi = ncmax( chi, c );
          }
          s1 = ncclamp( s1, clo, chi );
        }
        s1 = ncmin( s1, 3.0 * y[j] / h );
        if ( j + 1 < x.size() && x[j+1] > x[j] )
          s1 = ncmax( s1, -3.0 * y[j] / ( x[j+1] - x[j] ) );
        if ( h * s0 < -3.0 * y[i] * ( 1.0 + 1e-12 ) )
          return false;
        for ( std::size_t k = i + 1; k < j; ++k ) {
          const double yk = hermite( x[i], x[j], y[i], y[j], s0, s1, x[k] );
          if ( !( ncabs( yk - y[k] ) <= tol * ncmax( y[k], abs_floor ) ) )
            return false;
        }
        return true;
      }

    }
  }
}

NCX::ReducedHermitePoints NCX::reduceHermitePoints( const VectD& x,
                                                    const VectD& y,
                                                    double tol,
                                                    double abs_floor )
{
  //From the current point (with its slope fixed), find the most distant point
  //which can be reached with an interval satisfying the tolerance. An
  //exponential search followed by a binary search is used, which assumes that
  //feasibility is (approximately) monotonic in the end point. Every selected
  //interval is explicitly verified, so the tolerance is always respected at all
  //the sample points. An interval can span small discontinuities (both values
  //must then be within the tolerance), otherwise a new point is started at
  //the discontinuity, with a fresh slope.
  const std::size_t n = x.size();
  nc_assert_always( n >= 2 );
  auto startSlope = [&x,&y,n]( std::size_t i )
  {
    //Slope towards the next sample point (never violating the constraint,
    //since y is non-negative):
    return ( i + 1 < n && x[i+1] > x[i]
             ? ( y[i+1] - y[i] ) / ( x[i+1] - x[i] )
             : 0.0 );
  };
  ReducedHermitePoints res;
  res.idx.push_back( 0 );
  res.slope.push_back( startSlope( 0 ) );
  std::size_t i = 0;
  while ( i + 1 < n ) {
    if ( !( x[i+1] > x[i] ) ) {
      //Discontinuity:
      ++i;
      res.idx.push_back( i );
      res.slope.push_back( startSlope( i ) );
      continue;
    }
    auto reach = [&]( double s0, double& s_good )
    {
      std::size_t good = i + 1;
      (void)fitInterval( x, y, i, good, s0, tol, abs_floor, s_good );
      std::size_t bad = n;//first known infeasible end point (n: none)
      std::size_t step = 2;
      double s;
      while ( true ) {
        std::size_t j = ncmin( i + step, n - 1 );
        if ( j <= good )
          break;
        if ( x[j] > x[i] && fitInterval( x, y, i, j, s0, tol, abs_floor, s ) ) {
          good = j;
          s_good = s;
          if ( j == n - 1 )
            break;
          step *= 2;
        } else {
          bad = j;
          break;
        }
      }
      if ( bad < n ) {
        while ( bad - good > 1 ) {
          std::size_t mid = good + ( bad - good ) / 2;
          if ( fitInterval( x, y, i, mid, s0, tol, abs_floor, s ) ) {
            good = mid;
            s_good = s;
          } else {
            bad = mid;
          }
        }
      }
      return good;
    };
    double s_good;
    std::size_t good = reach( res.slope.back(), s_good );
    if ( good == i + 1 ) {
      //The slope inherited from the previous interval does not reach beyond
      //the next sample point, and could make the interpolation between them
      //poor. Restart with a fresh slope instead, by repeating the point (with
      //the same value):
      const double s0 = startSlope( i );
      if ( s0 != res.slope.back() ) {
        res.idx.push_back( i );
        res.slope.push_back( s0 );
        good = reach( s0, s_good );
      }
    }
    res.idx.push_back( good );
    res.slope.push_back( s_good );
    i = good;
  }
  return res;
}

unsigned NCX::autoNDense( double tol, double wlmax )
{
  //10000 points per decade for tol=1e-3. The segments of the table get
  //shorter for smaller tolerances (roughly as sqrt(tol) in smooth regions), so
  //the dense sample must get denser:
  nc_assert_always( tol > 0.0 && wlmax > denseWavelengthMin() );
  const double ndecades = std::log10( wlmax / denseWavelengthMin() );
  const double n = ndecades * 10000.0 * std::sqrt( 1e-3 / tol );
  return static_cast<unsigned>( ncclamp( n, 1000.0, 1e7 ) );
}

double NCX::evalIsoMacroXSTable( const double* wl, const double* xs,
                                 const double* slope, std::size_t n,
                                 double wavelength )
{
  //Evaluate the table. This is both used by NCrystal itself, but it is also a
  //reference implementation for people using our tables elsewhere.
  //
  //Note for NCrystal developers: Do not change this reference implementation in
  //                              a way that changes its behaviour. It might
  //                              have been copied elsewhere outside our
  //                              code-base!
  if ( n == 0 )
    return 0.0;
  if ( !( wavelength > wl[0] ) || n == 1 )
    return xs[0];
  if ( wavelength >= wl[n-1] ) {
    //Linear extrapolation of the value and slope of the last point, clamped
    //at 0:
    if ( slope[n-1] == 0.0 )
      return xs[n-1];
    const double v = xs[n-1] + ( wavelength - wl[n-1] ) * slope[n-1];
    return v > 0.0 ? v : 0.0;
  }
  //Last point with wl[i] <= wavelength (the second of a pair at a
  //discontinuity):
  const std::size_t i
    = static_cast<std::size_t>( std::upper_bound( wl, wl + n,
                                                  wavelength ) - wl ) - 1;
  nc_assert( i + 1 < n && wl[i+1] > wl[i] );
  return hermite( wl[i], wl[i+1], xs[i], xs[i+1], slope[i], slope[i+1],
                  wavelength );
}

NCX::IsoMacroXSTable
NCX::createIsoMacroXSTable( const std::function<double(double)>& xs_of_wl,
                            VectD discontinuities_in,
                            const IsoMacroXSTableParams& params )
{
  if ( !( params.wlmax > 10.0 * denseWavelengthMin() && std::isfinite( params.wlmax ) ) )
    NCRYSTAL_THROW2( BadInput, "IsoMacroXSTable: invalid wlmax " << params.wlmax
                     << " (must be larger than " << 10.0 * denseWavelengthMin()
                     << " Aa)" );
  if ( !( params.tol > 0.0 && params.tol < 1.0 ) )
    NCRYSTAL_THROW2( BadInput, "IsoMacroXSTable: invalid tolerance "
                     << params.tol << " (must be in the range (0,1))" );
  if ( params.ndense == 1 )
    NCRYSTAL_THROW( BadInput, "IsoMacroXSTable: ndense must be at least 2" );
  const unsigned ndense = ( params.ndense
                            ? params.ndense
                            : autoNDense( params.tol, params.wlmax ) );
  const double tol = params.tol;
  const double tol_build = tol_build_factor * tol;

  auto evalXS = [&xs_of_wl]( double wl )
  {
    const double v = xs_of_wl( wl );
    if ( !( std::isfinite( v ) && v >= 0.0 ) )
      NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: invalid cross section " << v
                       << " at wavelength " << wl << " Aa" );
    return v;
  };

  const VectD discontinuities = ( params.discontinuities
                                  ? mergeDiscontinuities( std::move( discontinuities_in ),
                                                          params.wlmax )
                                  : VectD() );

  //Dense sample: wavelength 0, log-spaced points, and points just below and
  //above each discontinuity. The x-values are the wavelengths, where a
  //discontinuity contributes two identical x-values (with the cross sections
  //below and above it).
  VectD wl_eval;//wavelengths at which to evaluate
  VectD x;//the corresponding x values of the dense sample
  wl_eval.reserve( 1 + ndense + 2 * discontinuities.size() );
  x.reserve( 1 + ndense + 2 * discontinuities.size() );
  wl_eval.push_back( 0.0 );//NB: the limit is calculated below
  x.push_back( 0.0 );
  {
    const VectD grid = geomspace( denseWavelengthMin(), params.wlmax, ndense );
    auto itD = discontinuities.begin();
    for ( auto wl : grid ) {
      for ( ; itD != discontinuities.end() && *itD <= wl*(1.0+2*edge_eps); ++itD ) {
        const double e = *itD;
        if ( x.size() > 1 && x.back() >= e*(1.0-2*edge_eps) ) {
          //previous grid point too close to the discontinuity, remove it:
          x.pop_back();
          wl_eval.pop_back();
        }
        wl_eval.push_back( e * ( 1.0 - edge_eps ) );
        x.push_back( e );
        wl_eval.push_back( e * ( 1.0 + edge_eps ) );
        x.push_back( e );
      }
      if ( x.size() > 1 && wl <= x.back()*(1.0+2*edge_eps) )
        continue;//grid point too close to the previous discontinuity
      wl_eval.push_back( wl );
      x.push_back( wl );
    }
    nc_assert_always( itD == discontinuities.end() );
  }

  VectD y;
  y.reserve( wl_eval.size() );
  y.push_back( 0.0 );
  for ( std::size_t i = 1; i < wl_eval.size(); ++i )
    y.push_back( evalXS( wl_eval[i] ) );
  wl_eval.clear();
  const double abs_floor = xs_abs_floor;

  //Known discontinuities where the cross section turns out to be continuous
  //(kinks, where the values just below and above agree to 1e-6) are not
  //counted:
  std::size_t ndisc = 0;
  for ( std::size_t i = 0; i + 1 < x.size(); ++i )
    if ( x[i+1] == x[i]
         && ncabs( y[i+1] - y[i] ) > kink_tol * ncmax( ncmax( y[i], y[i+1] ), abs_floor ) )
      ++ndisc;

  //The limit for wavelength -> 0, from linear extrapolation to 0 from very
  //short wavelengths (which is exact for any cross section of the form a+b*wl
  //there), at three scales. When the limit exists, the differences between
  //the estimates at consecutive scales decrease rapidly (e.g. by a factor of
  //100 for a term c*wl^2), while they do not for cross sections without a
  //limit (e.g. with 1/wl, log(wl) or sqrt(wl) behaviour):
  {
    auto extrapTo0 = [&evalXS]( double h ) { return 2.0 * evalXS( h ) - evalXS( 2.0 * h ); };
    const double lim_a = extrapTo0( 1e-7 );
    const double lim_b = extrapTo0( 1e-8 );
    const double lim_c = extrapTo0( 1e-9 );
    const double d_ab = ncabs( lim_a - lim_b );
    const double d_bc = ncabs( lim_b - lim_c );
    const double lim_tol = tol * ncmax( ncabs( lim_c ), abs_floor );
    //Either the estimates agree, or they converge (geometrically, with a
    //remaining error of at most d_bc*0.05/0.95 in lim_c):
    const bool converges = ( d_bc <= 1e-3 * lim_tol
                             || ( d_bc <= 0.05 * d_ab && 0.06 * d_bc <= lim_tol ) );
    if ( !converges )
      NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: the cross section does not"
                       " converge for wavelength -> 0 (estimated limits "
                       << lim_a << ", " << lim_b << " and " << lim_c << " barn)" );
    if ( lim_c < -abs_floor )
      NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: invalid limit " << lim_c
                       << " of the cross section for wavelength -> 0" );
    y.front() = ( lim_c > 0.0 ? lim_c : 0.0 );
  }

  //Reduce the dense sample, and verify the result at the midpoints between
  //all neighbouring sample points. Midpoints where the tolerance is not met
  //are added to the sample, and the procedure is repeated.
  ReducedHermitePoints red;
  std::size_t nrefine = 0;
  while ( true ) {
    red = reduceHermitePoints( x, y, tol_build, abs_floor );
    const auto& idx = red.idx;
    VectD add_x, add_y;
    add_x.reserve( x.size() );//at most one point per interval
    add_y.reserve( x.size() );
    std::size_t seg = 0;//table segment is from x[idx[seg]] to x[idx[seg+1]]
    for ( std::size_t k = 0; k + 1 < x.size(); ++k ) {
      while ( idx[seg+1] <= k )
        ++seg;
      if ( !( x[k+1] > x[k] ) )
        continue;//known discontinuity
      const std::size_t i = idx[seg];
      const std::size_t j = idx[seg+1];
      const double xm = 0.5 * ( x[k] + x[k+1] );
      const double ym = evalXS( xm );
      const double yc = hermite( x[i], x[j], y[i], y[j],
                                 red.slope[seg], red.slope[seg+1], xm );
      if ( ncabs( yc - ym ) <= tol_build * ncmax( ym, abs_floor ) )
        continue;
      if ( !( x[k+1] - x[k] > min_rel_spacing * x[k+1] ) )
        NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: the cross section has an"
                         " unexpected discontinuity (or an extremely sharp"
                         " feature) at wavelength " << xm << " Aa" );
      add_x.push_back( xm );
      add_y.push_back( ym );
    }
    if ( add_x.empty() )
      break;
    if ( ++nrefine > max_refine_iterations )
      NCRYSTAL_THROW( CalcError, "IsoMacroXSTable: the refinement of the table"
                      " did not converge" );
    //Merge the new points into the (sorted) sample:
    VectD new_x, new_y;
    new_x.reserve( x.size() + add_x.size() );
    new_y.reserve( x.size() + add_x.size() );
    std::size_t ia = 0;
    for ( std::size_t k = 0; k < x.size(); ++k ) {
      for ( ; ia < add_x.size() && add_x[ia] < x[k]; ++ia ) {
        new_x.push_back( add_x[ia] );
        new_y.push_back( add_y[ia] );
      }
      new_x.push_back( x[k] );
      new_y.push_back( y[k] );
    }
    nc_assert_always( ia == add_x.size() );
    x.swap( new_x );
    y.swap( new_y );
  }

  IsoMacroXSTable res;
  res.wl.reserve( red.idx.size() );
  res.xs.reserve( red.idx.size() );
  for ( auto i : red.idx ) {
    res.wl.push_back( x[i] );
    res.xs.push_back( y[i] );
  }
  res.slope = std::move( red.slope );
  nc_assert_always( res.wl.size() >= 2 && res.wl.front() == 0.0
                    && res.wl.back() == params.wlmax
                    && res.wl[res.wl.size()-2] < params.wlmax );

  //Final verification at random points in each segment (in addition to the
  //sample points, which are verified by construction):
  {
    RandXRSRImpl rng;//(fixed default seed, for reproducible verification)
    constexpr unsigned npts_per_segment = 4;
    for ( std::size_t i = 0; i + 1 < res.wl.size(); ++i ) {
      const double x0 = res.wl[i];
      const double x1 = res.wl[i+1];
      if ( !( x1 > x0 ) )
        continue;
      for ( unsigned k = 0; k < npts_per_segment; ++k ) {
        const double xv = x0 + ( 0.01 + 0.98 * rng.generate() ) * ( x1 - x0 );
        const double yv = evalXS( xv );
        const double yt = evalIsoMacroXSTable( res.wl.data(), res.xs.data(),
                                               res.slope.data(), res.wl.size(),
                                               xv );
        if ( !( ncabs( yt - yv ) <= tol * ncmax( yv, abs_floor ) ) )
          NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: verification failed at"
                           " wavelength " << xv << " Aa (table: " << yt
                           << " barn, exact: " << yv << " barn)" );
      }
    }
  }

  res.numberDensity = 0.0;
  res.nDiscontinuities = ndisc;
  res.nDense = x.size();
  res.nRefine = nrefine;
  return res;
}

NCX::IsoMacroXSTable
NCX::createIsoMacroXSTable( const MatCfg& cfg,
                            const IsoMacroXSTableParams& params )
{
  auto info = FactImpl::createInfo( cfg );
  auto scat = FactImpl::createScatter( cfg );
  auto absn = FactImpl::createAbsorption( cfg );
  if ( scat->isOriented() || absn->isOriented() )
    NCRYSTAL_THROW( BadInput, "IsoMacroXSTable: only isotropic materials are"
                    " supported" );

  std::set<std::string> names, unvalidated;
  inspectProcess( *scat, names, unvalidated );
  inspectProcess( *absn, names, unvalidated );
  for ( auto& n : unvalidated )
    NCRYSTAL_WARN( "IsoMacroXSTable: the physics process \"" << n << "\" has"
                   " not been validated for IsoMacroXSTable (features of its"
                   " cross section might be missed)" );
  VectD discontinuities;
  collectBreakpoints( *scat, discontinuities );
  collectBreakpoints( *absn, discontinuities );

  CachePtr cp_scat, cp_abs;
  auto xs_of_wl = [&scat,&absn,&cp_scat,&cp_abs]( double wl )
  {
    const NeutronEnergy ekin{ wl2ekin( wl ) };
    return ( scat->crossSectionIsotropic( cp_scat, ekin ).dbl()
             + absn->crossSectionIsotropic( cp_abs, ekin ).dbl() );
  };

  IsoMacroXSTable res = createIsoMacroXSTable( xs_of_wl,
                                               std::move( discontinuities ),
                                               params );
  res.numberDensity = info->getNumberDensity().dbl();
  res.processes.assign( names.begin(), names.end() );
  return res;
}
