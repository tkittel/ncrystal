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

#include "NCrystal/internal/math_utils/NCConvolve.hh"
#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace MathUtils {
    namespace {

      // Piecewise linear interpolation, zero outside [x.front(),x.back()].
      double interpLin( Span<const double> x, Span<const double> y,
                        double xv )
      {
        if ( !( xv >= x.front() && xv <= x.back() ) )
          return 0.0;
        auto it = std::upper_bound( x.begin(), x.end(), xv );
        if ( it == x.end() )
          return y.back();
        auto i = static_cast<std::size_t>( std::distance( x.begin(), it ) );
        nc_assert( i >= 1 );
        const double u = ( xv - x[i-1] ) / ( x[i] - x[i-1] );
        return y[i-1] + u * ( y[i] - y[i-1] );
      }

      // Interpolation of tabulated (x,y), zero outside [x.front(),x.back()].
      class TabInterp final {
      public:
        TabInterp( Span<const double> x, Span<const double> y,
                   ConvolveInterp mode )
          : m_x(x), m_y(y), m_mode(mode)
        {
          if ( m_x.size() < 4 && m_mode == ConvolveInterp::Lagrange4 )
            m_mode = ConvolveInterp::Linear;
          if ( m_mode == ConvolveInterp::Spline )
            initSpline();
        }

        double operator()( double xv ) const
        {
          if ( !( xv >= m_x.front() && xv <= m_x.back() ) )
            return 0.0;
          if ( m_mode == ConvolveInterp::Linear )
            return interpLin( m_x, m_y, xv );
          //Segment i with x_i <= xv <= x_{i+1}:
          const std::size_t n = m_x.size();
          auto it = std::upper_bound( m_x.begin(), m_x.end(), xv );
          std::size_t i = static_cast<std::size_t>(
                                       std::distance( m_x.begin(), it ) );
          i = ncmin( ncmax( i, std::size_t(1) ), n - 1 ) - 1;
          if ( m_mode == ConvolveInterp::Spline ) {
            const double h = m_x[i+1] - m_x[i];
            const double a = ( m_x[i+1] - xv ) / h;
            const double b = 1.0 - a;
            return a * m_y[i] + b * m_y[i+1]
              + ( ( a*a*a - a ) * m_M[i] + ( b*b*b - b ) * m_M[i+1] )
              * ( h * h / 6.0 );
          }
          //Lagrange through nodes j0..j0+3, the 4 nodes around segment i:
          const std::size_t j0 = ( i == 0 ? 0 : ncmin( i - 1, n - 4 ) );
          double res = 0.0;
          for ( std::size_t j = j0; j < j0 + 4; ++j ) {
            double l = 1.0;
            for ( std::size_t k = j0; k < j0 + 4; ++k )
              if ( k != j )
                l *= ( xv - m_x[k] ) / ( m_x[j] - m_x[k] );
            res += l * m_y[j];
          }
          return res;
        }

      private:
        Span<const double> m_x, m_y;
        ConvolveInterp m_mode;
        VectD m_M;//spline second derivatives

        void initSpline()
        {
          //Natural spline (M_0 = M_{n-1} = 0), Thomas algorithm:
          const std::size_t n = m_x.size();
          m_M.assign( n, 0.0 );
          if ( n < 3 )
            return;
          VectD diag(n,1.0), sup(n,0.0), sub(n,0.0), rhs(n,0.0);
          for ( std::size_t i = 1; i + 1 < n; ++i ) {
            const double h0 = m_x[i] - m_x[i-1];
            const double h1 = m_x[i+1] - m_x[i];
            sub[i] = h0;
            diag[i] = 2.0 * ( h0 + h1 );
            sup[i] = h1;
            rhs[i] = 6.0 * ( ( m_y[i+1] - m_y[i] ) / h1
                             - ( m_y[i] - m_y[i-1] ) / h0 );
          }
          for ( std::size_t i = 1; i < n; ++i ) {
            const double f = sub[i] / diag[i-1];
            diag[i] -= f * sup[i-1];
            rhs[i] -= f * rhs[i-1];
          }
          m_M[n-1] = rhs[n-1] / diag[n-1];
          for ( std::size_t i = n - 1; i-- > 0; )
            m_M[i] = ( rhs[i] - sup[i] * m_M[i+1] ) / diag[i];
        }
      };

      VectD convolveDirect( Span<const double> xf, Span<const double> f,
                            Span<const double> xg, Span<const double> g,
                            Span<const double> z, ConvolveInterp interp )
      {
        const TabInterp gint( xg, g, interp );
        const std::size_t n = xf.size();
        VectD wf(n);
        for ( std::size_t i = 0; i < n; ++i ) {
          const double lo = xf[ i ? i - 1 : 0 ];
          const double hi = xf[ i + 1 < n ? i + 1 : n - 1 ];
          wf[i] = 0.5 * ( hi - lo ) * f[i];
        }
        VectD res;
        res.reserve( z.size() );
        for ( auto zz : z ) {
          //Only x with z-x inside [xg.front(),xg.back()] contribute:
          auto itB = std::lower_bound( xf.begin(), xf.end(), zz - xg.back() );
          auto itE = std::upper_bound( itB, xf.end(), zz - xg.front() );
          StableSum sum;
          for ( auto it = itB; it != itE; ++it ) {
            auto i = static_cast<std::size_t>( std::distance( xf.begin(), it ) );
            sum.add( wf[i] * gint( zz - *it ) );
          }
          res.push_back( sum.sum() );
        }
        return res;
      }

      // Exact convolution of piecewise linear f and g, with f'' as deltas:
      // h(z) = Int f''(x) G2(z-x) dx, G1/G2 the 1st/2nd antiderivatives of g.
      // Splitting at p = z-yM (yM = xg.back(), beyond which G2 is linear)
      // gives h(z) = G1tot*f(p-) + G2(yM)*f'(p-) + sum over f nodes x_i >= p of
      // kappa_i*G2(z-x_i) (kappa_i: slope changes), plus the step terms
      // f_0*G1(z-x_0) and -f_{N-1}*G1(z-x_{N-1}) when x_0, x_{N-1} >= p.
      VectD convolveExact( Span<const double> xf, Span<const double> f,
                           Span<const double> xg, Span<const double> g,
                           Span<const double> z )
      {
        const std::size_t nf = xf.size();
        const std::size_t ng = xg.size();
        VectD sf(nf-1);//slopes of f
        for ( std::size_t i = 0; i + 1 < nf; ++i )
          sf[i] = ( f[i+1] - f[i] ) / ( xf[i+1] - xf[i] );
        VectD kappa(nf);
        for ( std::size_t i = 0; i < nf; ++i )
          kappa[i] = ( i + 1 < nf ? sf[i] : 0.0 ) - ( i ? sf[i-1] : 0.0 );
        //Antiderivatives of g at its nodes, and slopes:
        VectD sg(ng-1), G1(ng), G2(ng);
        G1[0] = G2[0] = 0.0;
        for ( std::size_t j = 0; j + 1 < ng; ++j ) {
          const double d = xg[j+1] - xg[j];
          sg[j] = ( g[j+1] - g[j] ) / d;
          G1[j+1] = G1[j] + d * ( g[j] + 0.5 * d * sg[j] );
          G2[j+1] = G2[j] + d * ( G1[j] + d * ( 0.5 * g[j] + d * sg[j] / 6.0 ) );
        }
        const double G1tot = G1.back();
        const double G2M = G2.back();
        const double yM = xg.back();
        //G1 and G2 at y <= yM:
        auto evalG = [&]( double y, double& g1, double& g2 )
        {
          if ( !( y > xg.front() ) ) {
            g1 = g2 = 0.0;
            return;
          }
          auto it = std::upper_bound( xg.begin(), xg.end(), y );
          auto j = static_cast<std::size_t>( std::distance( xg.begin(), it ) );
          j = ncmin( j, ng - 1 ) - 1;
          const double u = y - xg[j];
          g1 = G1[j] + u * ( g[j] + 0.5 * u * sg[j] );
          g2 = G2[j] + u * ( G1[j] + u * ( 0.5 * g[j] + u * sg[j] / 6.0 ) );
        };

        VectD res;
        res.reserve( z.size() );
        for ( auto zz : z ) {
          const double p = zz - yM;
          //k: first node with x_k >= p.
          auto itk = std::lower_bound( xf.begin(), xf.end(), p );
          auto k = static_cast<std::size_t>( std::distance( xf.begin(), itk ) );
          StableSum sum;
          if ( k > 0 && k < nf ) {
            //p inside segment k-1 (left limits of f and f'):
            const double u = ( p - xf[k-1] ) / ( xf[k] - xf[k-1] );
            sum.add( G1tot * ( f[k-1] + u * ( f[k] - f[k-1] ) ) );
            sum.add( G2M * sf[k-1] );
          }
          //Local terms, nodes with p <= x_i < z - xg.front():
          const double xi_end = zz - xg.front();
          double g1, g2;
          for ( std::size_t i = k; i < nf && xf[i] < xi_end; ++i ) {
            evalG( zz - xf[i], g1, g2 );
            sum.add( kappa[i] * g2 );
            if ( i == 0 )
              sum.add( f[0] * g1 );
            if ( i + 1 == nf )
              sum.add( -f[nf-1] * g1 );
          }
          res.push_back( sum.sum() );
        }
        return res;
      }

      VectD convolveFFT( Span<const double> xf, Span<const double> f,
                         Span<const double> xg, Span<const double> g,
                         Span<const double> z, double step,
                         ConvolveInterp interp )
      {
        auto sample = [step,interp]( Span<const double> x,
                                     Span<const double> y )
        {
          const TabInterp yint( x, y, interp );
          const double w = x.back() - x.front();
          const double nd = std::ceil( w / step ) + 1.0;
          if ( !( nd < 1e8 ) )
            NCRYSTAL_THROW2(BadInput,"convolve: FFT step "<<step
                            <<" too small for grid width "<<w);
          const auto n = static_cast<std::size_t>( nd );
          VectD res(n);
          for ( std::size_t i = 0; i < n; ++i )
            res[i] = yint( x.front() + step * static_cast<double>(i) );
          return res;
        };
        VectD a = sample( xf, f );
        VectD b = sample( xg, g );
        VectD y;
        FastConvolve().convolve( a, b, y, step );
        //Output on grid x0 + k*step, interpolated at z with the same mode:
        const double x0 = xf.front() + xg.front();
        VectD xy( y.size() );
        for ( std::size_t k = 0; k < y.size(); ++k )
          xy[k] = x0 + step * static_cast<double>(k);
        const TabInterp yint( xy, y, interp );
        VectD res;
        res.reserve( z.size() );
        for ( auto zz : z )
          res.push_back( yint( zz ) );
        return res;
      }
    }
  }
}

NC::VectD NC::MathUtils::convolve( Span<const double> xf, Span<const double> f,
                                   Span<const double> xg, Span<const double> g,
                                   Span<const double> z, ConvolveScheme scheme,
                                   double fft_step, ConvolveInterp interp )
{
#ifndef NDEBUG
  nc_assert_always( xf.size() >= 2 && xf.size() == f.size() );
  nc_assert_always( xg.size() >= 2 && xg.size() == g.size() );
  for ( std::size_t i = 1; i < xf.size(); ++i )
    nc_assert_always( xf[i] > xf[i-1] );
  for ( std::size_t i = 1; i < xg.size(); ++i )
    nc_assert_always( xg[i] > xg[i-1] );
  nc_assert_always( fft_step >= 0.0 && std::isfinite(fft_step) );
#endif
  if ( scheme == ConvolveScheme::DirectTrapezoid )
    return convolveDirect( xf, f, xg, g, z, interp );
  if ( scheme == ConvolveScheme::ExactLinear ) {
    if ( interp != ConvolveInterp::Linear )
      NCRYSTAL_THROW(BadInput,"convolve: the exact scheme requires linear"
                     " interpolation");
    //Use the narrower function as g (fewer local terms per z):
    if ( xg.back() - xg.front() > xf.back() - xf.front() )
      return convolveExact( xg, g, xf, f, z );
    return convolveExact( xf, f, xg, g, z );
  }
  nc_assert( scheme == ConvolveScheme::FFT );
  if ( fft_step == 0.0 ) {
    fft_step = kInfinity;
    for ( std::size_t i = 1; i < xf.size(); ++i )
      fft_step = ncmin( fft_step, xf[i] - xf[i-1] );
    for ( std::size_t i = 1; i < xg.size(); ++i )
      fft_step = ncmin( fft_step, xg[i] - xg[i-1] );
  }
  return convolveFFT( xf, f, xg, g, z, fft_step, interp );
}
