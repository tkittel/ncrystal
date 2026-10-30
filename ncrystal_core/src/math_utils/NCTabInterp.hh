#ifndef NCrystal_TabInterp_hh
#define NCrystal_TabInterp_hh

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
#include "NCrystal/internal/utils/NCSpline.hh"

//Private header of math_utils: interpolation of tabulated functions (shared
//by NCConvolve.cc and NCLines.cc).

namespace NCRYSTAL_NAMESPACE {
  namespace MathUtils {
    namespace detail {

      // Piecewise linear interpolation, zero outside [x.front(),x.back()].
      inline double interpLin( Span<const double> x, Span<const double> y,
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
          if ( m_mode == ConvolveInterp::PCHIP )
            m_pchip.emplace( m_x, m_y );
        }

        double operator()( double xv ) const
        {
          if ( !( xv >= m_x.front() && xv <= m_x.back() ) )
            return 0.0;
          if ( m_mode == ConvolveInterp::Linear )
            return interpLin( m_x, m_y, xv );
          if ( m_mode == ConvolveInterp::PCHIP )
            return m_pchip.value().eval( xv );
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
        Optional<PCHIPInterp> m_pchip;

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


    }
  }
}

#endif
