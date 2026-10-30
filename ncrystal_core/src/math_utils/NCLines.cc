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

#include "NCrystal/internal/math_utils/NCLines.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCTabInterp.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace MathUtils {
    namespace {

      // log( e^{-y} I_k(y) ) for k = 0..kmax (y > 0; -inf where zero). For
      // y < 1 from the power series e^{-y} (y/2)^k/k! sum_m (y^2/4)^m k!/(m!
      // (k+m)!) in log form (no underflow, as the weights multiply these by
      // e^{k betaE/2}, which can be huge). Otherwise by Miller's backward
      // recurrence I_{k-1} = (2k/y) I_k + I_{k+1} (per-step growth <= 2k, with
      // rescaling), started well above kmax and y, and normalised exactly by
      // I_0 + 2 sum_{k>0} I_k = e^y.
      VectD logScaledBesselI( double y, std::size_t kmax )
      {
#ifndef NDEBUG
        nc_assert_always( y > 0.0 && std::isfinite(y) );
#endif
        VectD res( kmax + 1, -kInfinity );
        if ( y < 1.0 ) {
          const double q = 0.25 * y * y;
          const double logy2 = std::log( 0.5 * y );
          double logkfact = 0.0;//log(k!)
          for ( std::size_t k = 0; k <= kmax; ++k ) {
            if ( k )
              logkfact += std::log( static_cast<double>(k) );
            double term = 1.0, sum = 1.0;
            for ( unsigned m = 1; m < 100; ++m ) {
              term *= q / ( m * ( static_cast<double>(k) + m ) );
              sum += term;
              if ( term < 1e-17 * sum )
                break;
            }
            res[k] = -y + static_cast<double>(k) * logy2 - logkfact
              + std::log( sum );
          }
          return res;
        }
        const double kk = ncmax( static_cast<double>(kmax), y );
        const double nstart_d = kk + 30.0 + 10.0 * std::sqrt( kk );
        nc_assert_always( nstart_d < 1e8 );
        const auto nstart = static_cast<std::size_t>( nstart_d );
        VectD val( kmax + 1, 0.0 );
        //Track a common log-scale so nothing overflows or underflows:
        double logscale = 0.0;
        double ip1 = 0.0, i0 = 1.0, sum = 0.0;//I_{k+1}, I_k, sum_{j>=k} I_j
        VectD scale_at( kmax + 1, 0.0 );
        for ( std::size_t k = nstart; k > 0; --k ) {
          if ( k <= kmax ) {
            val[k] = i0;
            scale_at[k] = logscale;
          }
          sum += i0;
          const double im1 = ( 2.0 * static_cast<double>(k) / y ) * i0 + ip1;
          ip1 = i0;
          i0 = im1;
          if ( i0 > 1e200 ) {
            //exact rescaling by 2^-664 (~1e-200):
            ip1 = std::ldexp( ip1, -664 );
            i0 = std::ldexp( i0, -664 );
            sum = std::ldexp( sum, -664 );
            logscale += 664.0 * std::log( 2.0 );
          }
        }
        //Actual (unnormalised) values are val[k]*exp(scale_at[k]-logscale)
        //relative to the final scale; normalisation 1/(I_0 + 2 sum):
        const double lognorm = -std::log( i0 + 2.0 * sum );
        res[0] = std::log( i0 ) + lognorm;
        for ( std::size_t k = 1; k <= kmax; ++k )
          if ( val[k] > 0.0 )
            res[k] = std::log( val[k] ) + ( scale_at[k] - logscale ) + lognorm;
        return res;
      }
    }
  }
}

NC::VectD NC::MathUtils::shiftAndAdd( Span<const double> xg,
                                      Span<const double> g,
                                      Span<const double> shifts,
                                      Span<const double> weights,
                                      Span<const double> z,
                                      ConvolveInterp interp )
{
#ifndef NDEBUG
  nc_assert_always( xg.size() >= 2 && xg.size() == g.size() );
  nc_assert_always( shifts.size() == weights.size() );
  for ( std::size_t i = 1; i < xg.size(); ++i )
    nc_assert_always( xg[i] > xg[i-1] );
#endif
  const detail::TabInterp gint( xg, g, interp );
  VectD res;
  res.reserve( z.size() );
  for ( auto zz : z ) {
    StableSum sum;
    for ( std::size_t k = 0; k < shifts.size(); ++k ) {
      const double x = zz - shifts[k];
      if ( x >= xg.front() && x <= xg.back() )
        sum.add( weights[k] * gint( x ) );
    }
    res.push_back( sum.sum() );
  }
  return res;
}

NC::MathUtils::Lines NC::MathUtils::oscillatorLines( double A, double betaE,
                                                     double reltol )
{
#ifndef NDEBUG
  nc_assert_always( A >= 0.0 && std::isfinite(A) );
  nc_assert_always( betaE > 0.0 && std::isfinite(betaE) );
  nc_assert_always( reltol >= 0.0 && reltol < 1.0 );
#endif
  Lines res;
  if ( !( A > 0.0 ) ) {
    res.shifts.push_back( 0.0 );
    res.weights.push_back( 1.0 );
    return res;
  }
  const double y = A / std::sinh( 0.5 * betaE );
  //log(w_k) = -A tanh(betaE/4) + k betaE/2 + log(e^{-y} I_k(y)). Grow kmax
  //(energy loss side, k>0, largest weights) until the weights are negligible:
  const double logw0 = -A * std::tanh( 0.25 * betaE );
  const double logtol = std::log( ncmax( reltol, 1e-300 ) );
  std::size_t kmax = static_cast<std::size_t>( y + 10.0 * std::sqrt( y ) )
    + static_cast<std::size_t>( A + 10.0 * std::sqrt( A ) ) + 20;
  VectD lie;
  double lmax;
  while ( true ) {
    lie = logScaledBesselI( y, kmax );
    lmax = -kInfinity;
    for ( std::size_t k = 0; k <= kmax; ++k )
      lmax = ncmax( lmax, logw0 + 0.5 * static_cast<double>(k) * betaE
                    + lie[k] );
    const double lend = logw0 + 0.5 * static_cast<double>(kmax) * betaE
      + lie[kmax];
    if ( lend < lmax + logtol - 5.0 )
      break;
    if ( kmax > 10000000 )
      NCRYSTAL_THROW2(CalcError,"oscillatorLines: too many lines needed for"
                      " A="<<A<<" and betaE="<<betaE);
    kmax *= 2;
  }
  //Gain side k<0 (shift +|k| betaE, weight e^{-|k| betaE} relative), from
  //large positive shifts down to zero, then the loss side:
  auto logw = [&]( long k ) {
    const std::size_t ak = static_cast<std::size_t>( k < 0 ? -k : k );
    return logw0 + 0.5 * static_cast<double>(k) * betaE + lie[ak];
  };
  const long kmaxl = static_cast<long>( kmax );
  for ( long k = kmaxl; k >= -kmaxl; --k ) {
    const double lw = logw( k );
    if ( lw - lmax < logtol )
      continue;
    res.shifts.push_back( -static_cast<double>(k) * betaE );
    res.weights.push_back( std::exp( lw ) );
  }
  return res;
}
