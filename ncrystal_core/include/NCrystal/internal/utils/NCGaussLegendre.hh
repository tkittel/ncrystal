#ifndef NCrystal_GaussLegendre_hh
#define NCrystal_GaussLegendre_hh

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


#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  //Gauss-Legendre quadrature rules with n points on [-1,1], for n in
  //{6,8,10,12,14,16,20,24}. As n is even, the nodes come in pairs (-x,+x)
  //with the same weight, so only the n/2 positive nodes, x (ascending), and
  //their weights, w, are provided (the weights thus sum to 1). Computed with
  //mpmath at 40 digits.
  struct GaussLegendreRule {
    const double * x;
    const double * w;
    unsigned nhalf;
  };
  inline GaussLegendreRule gaussLegendreRule( unsigned n );

  //Integral of f over [a,b] with npts-point Gauss-Legendre quadrature (npts
  //as for gaussLegendreRule) in each of nseg segments of equal width:
  template <class Func>
  double integrateGaussLegendre( const Func& f, double a, double b,
                                 unsigned npts, unsigned nseg = 1 );

  //8-point Gauss-Legendre quadrature (exact for polynomials up to degree 15):
  template <class Func>
  double integrateGaussLegendre8( const Func& f, double a, double b );

}

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace detail_gl {
    inline constexpr double x6[3] = {
      0.2386191860831969, 0.6612093864662645, 0.932469514203152 };
    inline constexpr double w6[3] = {
      0.46791393457269104, 0.3607615730481386, 0.17132449237917036 };
    inline constexpr double x8[4] = {
      0.1834346424956498, 0.525532409916329, 0.7966664774136267,
      0.9602898564975363 };
    inline constexpr double w8[4] = {
      0.362683783378362, 0.31370664587788727, 0.22238103445337448,
      0.10122853629037626 };
    inline constexpr double x10[5] = {
      0.14887433898163122, 0.4333953941292472, 0.6794095682990244,
      0.8650633666889845, 0.9739065285171717 };
    inline constexpr double w10[5] = {
      0.29552422471475287, 0.26926671930999635, 0.21908636251598204,
      0.1494513491505806, 0.06667134430868814 };
    inline constexpr double x12[6] = {
      0.1252334085114689, 0.3678314989981802, 0.5873179542866175,
      0.7699026741943047, 0.9041172563704749, 0.9815606342467192 };
    inline constexpr double w12[6] = {
      0.24914704581340277, 0.2334925365383548, 0.20316742672306592,
      0.16007832854334622, 0.10693932599531843, 0.04717533638651183 };
    inline constexpr double x14[7] = {
      0.10805494870734367, 0.31911236892788974, 0.5152486363581541,
      0.6872929048116855, 0.827201315069765, 0.9284348836635735,
      0.9862838086968123 };
    inline constexpr double w14[7] = {
      0.2152638534631578, 0.2051984637212956, 0.18553839747793782,
      0.15720316715819355, 0.12151857068790319, 0.08015808715976021,
      0.03511946033175186 };
    inline constexpr double x16[8] = {
      0.09501250983763744, 0.2816035507792589, 0.45801677765722737,
      0.6178762444026438, 0.755404408355003, 0.8656312023878318,
      0.9445750230732326, 0.9894009349916499 };
    inline constexpr double w16[8] = {
      0.1894506104550685, 0.18260341504492358, 0.16915651939500254,
      0.14959598881657674, 0.12462897125553388, 0.09515851168249279,
      0.062253523938647894, 0.027152459411754096 };
    inline constexpr double x20[10] = {
      0.07652652113349734, 0.22778585114164507, 0.37370608871541955,
      0.5108670019508271, 0.636053680726515, 0.7463319064601508,
      0.8391169718222188, 0.912234428251326, 0.9639719272779138,
      0.9931285991850949 };
    inline constexpr double w20[10] = {
      0.15275338713072584, 0.14917298647260374, 0.14209610931838204,
      0.13168863844917664, 0.11819453196151841, 0.10193011981724044,
      0.08327674157670475, 0.06267204833410907, 0.04060142980038694,
      0.017614007139152118 };
    inline constexpr double x24[12] = {
      0.06405689286260563, 0.1911188674736163, 0.3150426796961634,
      0.4337935076260451, 0.5454214713888396, 0.6480936519369755,
      0.7401241915785544, 0.820001985973903, 0.8864155270044011,
      0.9382745520027328, 0.9747285559713095, 0.9951872199970213 };
    inline constexpr double w24[12] = {
      0.12793819534675216, 0.1258374563468283, 0.12167047292780339,
      0.1155056680537256, 0.10744427011596563, 0.09761865210411388,
      0.08619016153195327, 0.0733464814110803, 0.05929858491543678,
      0.04427743881741981, 0.028531388628933663, 0.0123412297999872 };
  }

  inline GaussLegendreRule gaussLegendreRule( unsigned n )
  {
    switch ( n ) {
    case 6: return { detail_gl::x6, detail_gl::w6, 3 };
    case 8: return { detail_gl::x8, detail_gl::w8, 4 };
    case 10: return { detail_gl::x10, detail_gl::w10, 5 };
    case 12: return { detail_gl::x12, detail_gl::w12, 6 };
    case 14: return { detail_gl::x14, detail_gl::w14, 7 };
    case 16: return { detail_gl::x16, detail_gl::w16, 8 };
    case 20: return { detail_gl::x20, detail_gl::w20, 10 };
    case 24: return { detail_gl::x24, detail_gl::w24, 12 };
    default: break;
    };
    NCRYSTAL_THROW2(LogicError,"No Gauss-Legendre rule with "<<n<<" points");
  }

  template <class Func>
  inline double integrateGaussLegendre( const Func& f, double a, double b,
                                        unsigned npts, unsigned nseg )
  {
    nc_assert( nseg >= 1 );
    const auto gl = gaussLegendreRule( npts );
    //Half-width and centres of segments:
    const double h = 0.5 * ( b - a ) / nseg;
    double sum = 0.0;
    for ( auto iseg : ncrange(nseg) ) {
      const double m = ( nseg == 1 ? 0.5 * ( a + b ) : a + ( 2*iseg + 1 ) * h );
      for ( auto i : ncrange(gl.nhalf) )
        sum += gl.w[i] * ( f( m - h * gl.x[i] ) + f( m + h * gl.x[i] ) );
    }
    return sum * h;
  }

  template <class Func>
  inline double integrateGaussLegendre8( const Func& f, double a, double b )
  {
    return integrateGaussLegendre( f, a, b, 8 );
  }

}

#endif
