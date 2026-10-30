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

#include "NCrystal/internal/math_utils/NCFourierInt.hh"
#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <complex>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace MathUtils {
    namespace {

      using Cplx = std::complex<double>;
      using VectC = std::vector<Cplx>;

      // g_m(w) = Int_0^h u^m e^{iwu} du for m=0..3 (into g[0..3]). Series for
      // |wh|<2 (the upward recurrence loses precision there), else recurrence.
      void pieceIntegrals( double w, double h, Cplx* g )
      {
        const double x = w * h;
        if ( ncabs(x) < 2.0 ) {
          //Int_0^h u^m e^{iwu} du = h^(m+1) sum_k (ix)^k/(k!(m+k+1))
          for ( unsigned m = 0; m < 4; ++m ) {
            Cplx term(1.0,0.0), acc(1.0/(m+1),0.0);
            for ( unsigned k = 1; k < 60; ++k ) {
              term *= Cplx(0.0,x/k);
              const Cplx c = term / double(m+k+1);
              acc += c;
              if ( std::norm(c) < 1e-40 * std::norm(acc) )
                break;
            }
            g[m] = acc * std::pow(h,m+1);
          }
          return;
        }
        const Cplx e( std::cos(x), std::sin(x) );
        const Cplx inv_iw( 0.0, -1.0/w );
        double hm = 1.0;
        g[0] = ( e - 1.0 ) * inv_iw;
        for ( unsigned m = 1; m < 4; ++m ) {
          hm *= h;
          g[m] = ( hm * e - double(m) * g[m-1] ) * inv_iw;
        }
      }

      // Second derivatives of the cubic spline through y (uniform spacing h),
      // natural at the end, and at the start either natural or clamped with
      // first derivative zero.
      VectD splineSecondDerivs( Span<const double> y, double h,
                                bool clamped_start )
      {
        //Rows: start (2M_0+M_1 = k(y_1-y_0) if clamped, else M_0 = 0),
        //interior M_{i-1}+4M_i+M_{i+1} = k(y_{i+1}-2y_i+y_{i-1}) and
        //M_{n-1} = 0, solved with the Thomas algorithm:
        const std::size_t n = y.size();
        nc_assert( n >= 2 );
        const double k = 6.0 / ( h * h );
        VectD sub(n,1.0), diag(n,4.0), sup(n,1.0), rhs(n,0.0);
        diag[0] = clamped_start ? 2.0 : 1.0;
        sup[0] = clamped_start ? 1.0 : 0.0;
        rhs[0] = clamped_start ? k * ( y[1] - y[0] ) : 0.0;
        for ( std::size_t i = 1; i + 1 < n; ++i )
          rhs[i] = k * ( y[i+1] - 2.0 * y[i] + y[i-1] );
        sub[n-1] = 0.0;
        diag[n-1] = 1.0;
        for ( std::size_t i = 1; i < n; ++i ) {
          const double f = sub[i] / diag[i-1];
          diag[i] -= f * sup[i-1];
          rhs[i] -= f * rhs[i-1];
        }
        VectD res(n);
        res[n-1] = rhs[n-1] / diag[n-1];
        for ( std::size_t i = n - 1; i-- > 0; )
          res[i] = ( rhs[i] - sup[i] * res[i+1] ) / diag[i];
        return res;
      }

      // Polynomial coefficients a[m][n] (m=0..3) of the cubic spline pieces.
      void splineCoeffs( Span<const double> y, double h, bool clamped_start,
                         VectD* a )
      {
        const auto M = splineSecondDerivs( y, h, clamped_start );
        const std::size_t np = y.size() - 1;
        for ( unsigned m = 0; m < 4; ++m )
          a[m].resize(np);
        for ( std::size_t i = 0; i < np; ++i ) {
          a[0][i] = y[i];
          a[1][i] = ( y[i+1] - y[i] ) / h - h * ( 2.0*M[i] + M[i+1] ) / 6.0;
          a[2][i] = 0.5 * M[i];
          a[3][i] = ( M[i+1] - M[i] ) / ( 6.0 * h );
        }
      }

      // Sum_n c_n e^{iwnh} by Horner's scheme (one cos/sin call per w).
      Cplx phaseSumDirect( const VectC& c, double w, double h )
      {
        const Cplx z( std::cos(w*h), std::sin(w*h) );
        Cplx acc(0.0,0.0);
        for ( std::size_t i = c.size(); i-- > 0; )
          acc = acc * z + c[i];
        return acc;
      }

      // Periodic phase sums on the FFT grid w_k = 2pi k/(M h), k=0..M-1.
      VectC phaseSumsFFT( const VectC& c, std::size_t M )
      {
        nc_assert( c.size() <= M );
        VectC buf( c.begin(), c.end() );
        FastConvolve fc;
        fc.fft( buf, true, M );//NB: "forward" uses exp(+i..)
        nc_assert( buf.size() == M );
        return buf;
      }

      Cplx interpPeriodic( const VectC& s, double w, double h, bool cubic )
      {
        const double Md = static_cast<double>( s.size() );
        double p = w * h * Md / k2Pi;
        p -= Md * std::floor( p / Md );//into [0,M)
        double fl = std::floor(p);
        const double u = p - fl;
        auto idx = [&s](double ii)
        {
          auto i = static_cast<std::ptrdiff_t>(ii);
          const auto n = static_cast<std::ptrdiff_t>(s.size());
          i %= n;
          if ( i < 0 )
            i += n;
          return s[static_cast<std::size_t>(i)];
        };
        if ( !cubic )
          return idx(fl) * (1.0-u) + idx(fl+1.0) * u;
        //4-point Lagrange through fl-1..fl+2:
        const double wm1 = -u * ( u - 1.0 ) * ( u - 2.0 ) / 6.0;
        const double w0 = ( u + 1.0 ) * ( u - 1.0 ) * ( u - 2.0 ) / 2.0;
        const double w1 = -( u + 1.0 ) * u * ( u - 2.0 ) / 2.0;
        const double w2 = ( u + 1.0 ) * u * ( u - 1.0 ) / 6.0;
        return idx(fl-1.0)*wm1 + idx(fl)*w0 + idx(fl+1.0)*w1 + idx(fl+2.0)*w2;
      }
    }
  }
}

NC::VectD NC::MathUtils::fourierIntegral( double h, Span<const double> fr,
                                          Span<const double> fi,
                                          Span<const double> w,
                                          FourierScheme scheme,
                                          FourierEval evalmode,
                                          unsigned oversample )
{
#ifndef NDEBUG
  nc_assert_always( h > 0.0 && std::isfinite(h) );
  nc_assert_always( fr.size() >= 2 );
  nc_assert_always( fi.empty() || fi.size() == fr.size() );
  nc_assert_always( oversample >= 1 );
#endif
  const std::size_t N = fr.size();
  const bool has_im = !fi.empty();

  //Terms: F(w) = (1/pi) Re sum_m G_m(w) sum_n c[m]_n e^{iw t_n}, where q =
  //fr - i*fi so Re[q e^{iwt}] = fr cos(wt) + fi sin(wt). G_m(w) is 1 for the
  //trapezoid schemes and the piece integral g_m(w) for the others:
  const bool piecewise = ( scheme == FourierScheme::PiecewiseLinear
                           || scheme == FourierScheme::CubicSpline );
  std::vector<VectC> c;
  if ( !piecewise ) {
    c.resize(1);
    c[0].resize(N);
    for ( std::size_t n = 0; n < N; ++n ) {
      const double wt = ( n == 0 || n + 1 == N ) ? 0.5 * h : h;
      c[0][n] = Cplx( fr[n], has_im ? -fi[n] : 0.0 ) * wt;
    }
  } else if ( scheme == FourierScheme::PiecewiseLinear ) {
    c.assign( 2, VectC(N-1) );
    for ( std::size_t n = 0; n + 1 < N; ++n ) {
      const Cplx q0( fr[n], has_im ? -fi[n] : 0.0 );
      const Cplx q1( fr[n+1], has_im ? -fi[n+1] : 0.0 );
      c[0][n] = q0;
      c[1][n] = ( q1 - q0 ) / h;
    }
  } else {
    nc_assert( scheme == FourierScheme::CubicSpline );
    VectD ar[4], ai[4];
    splineCoeffs( fr, h, true, ar );
    if ( has_im )
      splineCoeffs( fi, h, false, ai );
    c.assign( 4, VectC(N-1) );
    for ( unsigned m = 0; m < 4; ++m )
      for ( std::size_t n = 0; n + 1 < N; ++n )
        c[m][n] = Cplx( ar[m][n], has_im ? -ai[m][n] : 0.0 );
  }

  //Endpoint correction for TrapezoidEC: -h^2/12 [D(T)-D(0)] with
  //D(t) = (q'(t) + i w q(t)) e^{iwt}, q' from 2nd order one-sided differences
  //(first order if N==2):
  Cplx q0, q0p, qT, qTp;
  if ( scheme == FourierScheme::TrapezoidEC ) {
    auto q = [&fr,&fi,has_im](std::size_t i)
    {
      return Cplx( fr[i], has_im ? -fi[i] : 0.0 );
    };
    q0 = q(0);
    qT = q(N-1);
    if ( N >= 3 ) {
      q0p = ( -3.0*q(0) + 4.0*q(1) - q(2) ) / ( 2.0*h );
      qTp = ( 3.0*q(N-1) - 4.0*q(N-2) + q(N-3) ) / ( 2.0*h );
    } else {
      q0p = qTp = ( q(1) - q(0) ) / h;
    }
    //fr is even so fr'(0)=0 exactly (an estimate would leave an error which
    //does not decay with w):
    q0p.real(0.0);
  }
  const double T = h * static_cast<double>(N-1);

  std::vector<VectC> fftsums;
  if ( evalmode != FourierEval::Direct ) {
    std::size_t M = 1;
    while ( M < static_cast<std::size_t>(oversample) * N )
      M *= 2;
    fftsums.reserve( c.size() );
    for ( auto& cm : c )
      fftsums.push_back( phaseSumsFFT( cm, M ) );
  }
  const bool cubic_interp = ( evalmode == FourierEval::FFTCubic );

  VectD res;
  res.reserve( w.size() );
  Cplx g[4];
  for ( auto ww : w ) {
    if ( piecewise )
      pieceIntegrals( ww, h, g );
    Cplx acc(0.0,0.0);
    for ( std::size_t m = 0; m < c.size(); ++m ) {
      const Cplx s = ( fftsums.empty()
                       ? phaseSumDirect( c[m], ww, h )
                       : interpPeriodic( fftsums[m], ww, h, cubic_interp ) );
      acc += ( piecewise ? g[m] * s : s );
    }
    if ( scheme == FourierScheme::TrapezoidEC ) {
      const Cplx iw( 0.0, ww );
      const Cplx eT( std::cos(ww*T), std::sin(ww*T) );
      acc -= ( h * h / 12.0 ) * ( ( qTp + iw * qT ) * eT - ( q0p + iw * q0 ) );
    }
    res.push_back( acc.real() * kInvPi );
  }
  return res;
}
