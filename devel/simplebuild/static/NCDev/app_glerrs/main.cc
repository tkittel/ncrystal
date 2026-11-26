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

// Regenerates the gos_gl_errs table in NCGaussOnSphere.cc: maximum relative
// errors of Gauss-Legendre integration of circles with the numbers of points
// in gos_gl_npts, for truncation angles of N sigma. Found by comparison with
// 400 points, for many circles and mosaicities of 0.1-40 degrees.

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include <cstdio>

namespace NC = NCrystal;

namespace {

  //n-point rule on [0,1], with Newton iterations for the roots of P_n:
  void glRule( unsigned n, NC::VectD& x, NC::VectD& w )
  {
    x.resize(n);
    w.resize(n);
    for ( auto i : NC::ncrange(n) ) {
      double z = std::cos( NC::kPi*(i+0.75)/(n+0.5) );
      double dp = 1.0;
      for ( int it = 0; it < 100; ++it ) {
        double p1 = 1.0, p2 = 0.0;
        for ( unsigned j = 1; j <= n; ++j ) {
          const double p3 = p2;
          p2 = p1;
          p1 = ( (2.0*j-1.0)*z*p2 - (j-1.0)*p3 ) / j;
        }
        dp = n * ( z*p1 - p2 ) / ( z*z - 1.0 );
        const double dz = p1 / dp;
        z -= dz;
        if ( NC::ncabs( dz ) < 1e-16 )
          break;
      }
      x[i] = 0.5 * ( z + 1.0 );
      w[i] = 1.0 / ( ( 1.0 - z*z ) * dp * dp );
    }
  }

  //Circle integral (up to constant factors) as in GaussOnSphere:
  double integ( double expfact, double wtrunc, double cg, double sg,
                double ca, double sa, const NC::VectD& x, const NC::VectD& w )
  {
    const double sasg = sa*sg, cd = ca*cg + sasg;
    const double sd = sg*ca - cg*sa;
    const double w0 = sd*sd / ( 1.0 + cd );
    const double y = ( wtrunc - w0 ) / ( 2.0*sasg );
    if ( !( y > 0.0 ) )
      return 0.0;
    const double tmax = ( y >= 1.0 ? NC::kPi : 2.0*std::asin( std::sqrt(y) ) );
    NC::StableSum sum;
    for ( auto i : NC::ncrange( x.size() ) ) {
      const double sh = std::sin( 0.5*tmax*x[i] );
      const double b = 2.0*std::asin( std::sqrt( 0.5*( w0 + 2*sasg*sh*sh ) ) );
      sum.add( w[i] * std::exp( expfact*b*b ) );
    }
    return sum.sum() * tmax;
  }
}

int main()
{
  const unsigned npts[] = { 8, 10, 12, 14, 16, 20, 24 };
  std::vector<NC::VectD> xx(7), ww(7);
  for ( auto k : NC::ncrange(7) )
    glRule( npts[k], xx[k], ww[k] );
  NC::VectD xr, wr;
  glRule( 400, xr, wr );
  for ( double N : { 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0 } ) {
    double worst[7] = {};
    for ( double mosdeg : { 0.1, 1.0, 5.0, 20.0, 40.0 } ) {
      const double sigma = mosdeg*NC::kDeg/2.354820045030949;
      const double ta = N*sigma;
      if ( !( ta < NC::kPiHalf*0.99 ) )
        continue;
      const double expfact = -0.5/(sigma*sigma);
      const double st2 = std::sin(0.5*ta), wtrunc = 2*st2*st2;
      auto rng = NC::createBuiltinRNG(1);
      for ( auto ic : NC::ncrange(3000) ) {
        const double alpha = ( ic%10==0 ? rng->generate()*2*ta
                               : (0.001 + 0.999*rng->generate())*NC::kPiHalf );
        const double gamma = alpha + (2*rng->generate()-1.0)*ta;
        if ( gamma <= 0 || gamma >= NC::kPi || alpha <= 0 )
          continue;
        const double cg(std::cos(gamma)), sg(std::sin(gamma));
        const double ca(std::cos(alpha)), sa(std::sin(alpha));
        const double ref = integ(expfact,wtrunc,cg,sg,ca,sa,xr,wr);
        //Ignore negligible values (below 1e-4 of a typical peak value):
        if ( !( ref > 1e-4*std::sqrt(NC::k2Pi)*sigma ) )
          continue;
        for ( auto k : NC::ncrange(7) ) {
          const double v = integ(expfact,wtrunc,cg,sg,ca,sa,xx[k],ww[k]);
          worst[k] = NC::ncmax( worst[k], NC::ncabs(v-ref)/ref );
        }
      }
    }
    std::printf("    {");
    for ( auto k : NC::ncrange(7) )
      std::printf(" %.0e%s", worst[k], k+1<7 ? "," : " },\n");
  }
  return 0;
}
