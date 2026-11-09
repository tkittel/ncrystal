#ifndef NCrystal_Spline_FMA_hh
#define NCrystal_Spline_FMA_hh

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

#include "NCrystal/internal/utils/NCFMADispatch.hh"
#include "NCrystal/internal/utils/NCSpline.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace {

    //Hyman filter of slope m at a node, given the secant slopes dl and dr on
    //its two sides: keeps the sign of the data where it is monotone, limits it
    //to 3 times the smaller secant slope (the Fritsch-Carlson bound), and
    //zeroes it at local extrema of the data, which keeps each interval within
    //its node values:
    ncforceinline double hymanFilter( double m, double dl, double dr )
    {
      if ( !( dl * dr > 0.0 ) )
        return 0.0;
      const double lim = 3.0 * ncmin( ncabs( dl ), ncabs( dr ) );
      if ( m * dl < 0.0 )
        return 0.0;
      return ( ncabs( m ) > lim ? ( dl > 0.0 ? lim : -lim ) : m );
    }

    //HymanInterp's slope at node i of n>=4 values y on a uniform grid with
    //spacing h, with the derivatives of the 4-node cubics in closed form
    //(d0..d3: times 6h, at the 1st..4th node of the cubic through y[j..j+3]):
    ncforceinline double hymanSlopeUniform( const double* y, std::size_t n,
                                            double h, std::size_t i )
    {
      auto d0 = [y]( std::size_t j )
      {
        return ncfma( 2.0, y[j+3], ncfma( 18.0, y[j+1],
                                          -ncfma( 11.0, y[j], 9.0 * y[j+2] ) ) );
      };
      auto d1 = [y]( std::size_t j )
      {
        return ncfma( 3.0, ncfma( 2.0, y[j+2], -y[j+1] ),
                      -ncfma( 2.0, y[j], y[j+3] ) );
      };
      auto d2 = [y]( std::size_t j )
      {
        return ncfma( 3.0, ncfma( -2.0, y[j+1], y[j+2] ),
                      ncfma( 2.0, y[j+3], y[j] ) );
      };
      auto d3 = [y]( std::size_t j )
      {
        return ncfma( 11.0, y[j+3], ncfma( 9.0, y[j+1],
                                           -ncfma( 2.0, y[j], 18.0 * y[j+2] ) ) );
      };
      //Same stencils as HymanInterp (the average of the two cubics with the
      //node in their middle where available):
      double d6h;
      if ( i >= 2 && i + 2 < n )
        d6h = 0.5 * ( d1( i - 1 ) + d2( i - 2 ) );
      else if ( i == 1 )
        d6h = d1( 0 );
      else if ( i + 2 == n )
        d6h = d2( n - 4 );
      else
        d6h = ( i == 0 ? d0( 0 ) : d3( n - 4 ) );
      const double dl = ( i > 0 ? y[i] - y[i-1] : y[1] - y[0] ) / h;
      const double dr = ( i + 1 < n ? ( y[i+1] - y[i] ) / h : dl );
      return hymanFilter( d6h / ( 6.0 * h ), dl, dr );
    }

  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    void, hymanSlopesUniform,
    ( const double* y, std::size_t n, double h, double* slopes )
  )
  NCRYSTAL_FMADISPATCH_DECLARATOR_C(void,detail_hymanSlopesUniform,
                                    hymanSlopesUniform)
  ( const double* y, std::size_t n, double h, double* slopes )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(hymanSlopesUniform,(y,n,h,slopes));
    for ( std::size_t i = 0; i < n; ++i )
      slopes[i] = hymanSlopeUniform( y, n, h, i );
  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    void, hymanCubicsUniform,
    ( const double* y, std::size_t n, double h, double* coeffs )
  )
  NCRYSTAL_FMADISPATCH_DECLARATOR_C(void,detail_hymanCubicsUniform,
                                    hymanCubicsUniform)
  ( const double* y, std::size_t n, double h, double* coeffs )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(hymanCubicsUniform,(y,n,h,coeffs));
    //Cubic Hermite interpolation with the Hyman slopes, as coefficients of
    //c0+u*(c1+u*(c2+u*c3)) in each interval (u in [0,1]):
    double hs1 = h * hymanSlopeUniform( y, n, h, 0 );
    for ( std::size_t i = 0; i + 1 < n; ++i ) {
      const double hs0 = hs1;
      hs1 = h * hymanSlopeUniform( y, n, h, i + 1 );
      double * c = coeffs + 4 * i;
      c[0] = y[i];
      c[1] = hs0;
      c[2] = ncfma( 3.0, y[i+1] - y[i], -ncfma( 2.0, hs0, hs1 ) );
      c[3] = ncfma( 2.0, y[i] - y[i+1], hs0 + hs1 );
    }
  }
}

#endif
