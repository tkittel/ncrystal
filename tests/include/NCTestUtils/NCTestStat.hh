#ifndef NCrystal_TestStat_hh
#define NCrystal_TestStat_hh

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

// Statistics utilities for tests.

namespace nctest {

  using VectD = NCrystal::VectD;

  //Kolmogorov-Smirnov distance: the maximum distance between the empirical
  //CDFs of two non-empty samples.
  inline double ksDistance( VectD x, VectD y )
  {
    nc_assert_always( !x.empty() && !y.empty() );
    std::sort( x.begin(), x.end() );
    std::sort( y.begin(), y.end() );
    const double nx = static_cast<double>( x.size() );
    const double ny = static_cast<double>( y.size() );
    //Walk through all values of both samples in increasing order, comparing
    //the CDFs after each distinct value:
    std::size_t i = 0, j = 0;
    double dmax = 0.0;
    while ( i < x.size() || j < y.size() ) {
      const double v = ( j == y.size()
                         || ( i < x.size() && x[i] <= y[j] ) ) ? x[i] : y[j];
      while ( i < x.size() && x[i] <= v )
        ++i;
      while ( j < y.size() && y[j] <= v )
        ++j;
      dmax = NCrystal::ncmax( dmax, NCrystal::ncabs( i / nx - j / ny ) );
    }
    return dmax;
  }

  //P-value for two non-empty samples coming from the same distribution, from
  //their Kolmogorov-Smirnov distance (asymptotic two-sided Kolmogorov
  //distribution).
  inline double ksPValue( const VectD& x, const VectD& y )
  {
    const double n = static_cast<double>( x.size() );
    const double m = static_cast<double>( y.size() );
    const double lam = ksDistance( x, y ) * std::sqrt( n * m / ( n + m ) );
    if ( !( lam > 0.0 ) )
      return 1.0;
    //p/2 = sum_{k>=1} (-1)^{k-1} exp(-2 k^2 lam^2):
    const double m2lam2 = -2.0 * lam * lam;
    double s = 0.0;
    double sgn = -1.0;
    for ( int k = 1; k <= 200; ++k ) {
      sgn = -sgn;
      const double e = std::exp( m2lam2 * k * k );
      s += sgn * e;
      if ( k > 20 && e < 1e-10 * NCrystal::ncabs( s ) )
        break;
    }
    return NCrystal::ncclamp( 2.0 * s, 0.0, 1.0 );
  }

}

#endif
