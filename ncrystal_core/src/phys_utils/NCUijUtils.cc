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

#include "NCrystal/internal/phys_utils/NCUijUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;
namespace NCU = NCrystal::UijUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Index into Uij of component (i,j):
    unsigned uijIdx( unsigned i, unsigned j )
    {
      static constexpr unsigned idx[9] = { 0, 3, 4, 3, 1, 5, 4, 5, 2 };
      return idx[3*i+j];
    }
  }
}

NCU::CellMetric NCU::cellMetric( double a, double b, double c,
                                 double alpha, double beta, double gamma )
{
#ifndef NDEBUG
  nc_assert_always( a > 0.0 && b > 0.0 && c > 0.0 );
  nc_assert_always( alpha > 0.0 && alpha < 180.0 );
  nc_assert_always( beta > 0.0 && beta < 180.0 );
  nc_assert_always( gamma > 0.0 && gamma < 180.0 );
#endif
  //Exact zero for the (common) 90 degree angles:
  auto cosd = []( double x ) { return x == 90.0 ? 0.0 : std::cos( x * kDeg ); };
  const double ca = cosd( alpha ), cb = cosd( beta ), cg = cosd( gamma );
  CellMetric m;
  m.g = { a*a, a*b*cg, a*c*cb,
          a*b*cg, b*b, b*c*ca,
          a*c*cb, b*c*ca, c*c };
  //Reciprocal lengths from the diagonal of G^-1 (cofactors over det):
  const auto& g = m.g;
  const double c00 = g[4]*g[8] - g[5]*g[7];
  const double c11 = g[0]*g[8] - g[2]*g[6];
  const double c22 = g[0]*g[4] - g[1]*g[3];
  const double det = g[0]*c00 - g[1]*( g[3]*g[8] - g[5]*g[6] )
    + g[2]*( g[3]*g[7] - g[4]*g[6] );
  nc_assert_always( det > 0.0 );
  m.recilen = { std::sqrt( c00 / det ), std::sqrt( c11 / det ),
                std::sqrt( c22 / det ) };
  return m;
}

NCU::Uij NCU::cifToFrac( const Uij& u, const CellMetric& m )
{
  const auto& r = m.recilen;
  return { u[0]*r[0]*r[0], u[1]*r[1]*r[1], u[2]*r[2]*r[2],
           u[3]*r[0]*r[1], u[4]*r[0]*r[2], u[5]*r[1]*r[2] };
}

NCU::Uij NCU::fracToCif( const Uij& u, const CellMetric& m )
{
  const auto& r = m.recilen;
  return { u[0]/(r[0]*r[0]), u[1]/(r[1]*r[1]), u[2]/(r[2]*r[2]),
           u[3]/(r[0]*r[1]), u[4]/(r[0]*r[2]), u[5]/(r[1]*r[2]) };
}

NCU::Uij NCU::transformFrac( const SymOp& op, const Uij& u )
{
  static constexpr unsigned ri[6] = { 0, 1, 2, 0, 0, 1 };
  static constexpr unsigned rj[6] = { 0, 1, 2, 1, 2, 2 };
  Uij res;
  for ( unsigned n = 0; n < 6; ++n ) {
    //(R U R^T)_ij = sum_kl R_ik U_kl R_jl, with integer R:
    StableSum s;
    for ( unsigned k = 0; k < 3; ++k ) {
      const int rik = op.rot( ri[n], k );
      if ( !rik )
        continue;
      for ( unsigned l = 0; l < 3; ++l ) {
        const int rjl = op.rot( rj[n], l );
        if ( rjl )
          s.add( rik * rjl * u[uijIdx( k, l )] );
      }
    }
    res[n] = s.sum();
  }
  return res;
}

double NCU::ueqFromFrac( const Uij& u, const CellMetric& m )
{
  StableSum s;
  for ( unsigned i = 0; i < 3; ++i )
    for ( unsigned j = 0; j < 3; ++j )
      s.add( u[uijIdx( i, j )] * m.g[3*i+j] );
  return s.sum() / 3.0;
}

bool NCU::isPositiveDefinite( const Uij& u )
{
  //Sylvester's criterion (leading principal minors positive):
  const double m2 = u[0]*u[1] - u[3]*u[3];
  const double m3 = u[0]*( u[1]*u[2] - u[5]*u[5] )
    - u[3]*( u[3]*u[2] - u[5]*u[4] )
    + u[4]*( u[3]*u[5] - u[1]*u[4] );
  return u[0] > 0.0 && m2 > 0.0 && m3 > 0.0;
}

double NCU::maxAbs( const Uij& u )
{
  double res = 0.0;
  for ( auto e : u )
    res = ncmax( res, ncabs( e ) );
  return res;
}

std::array<double,3> NCU::cartEigenvalues( const Uij& u, const CellMetric& m )
{
  //Eigenvalues of U_frac*G (similar to U_cart), from the characteristic
  //polynomial (all roots real and positive for positive definite U_frac):
  std::array<double,9> M;
  for ( unsigned i = 0; i < 3; ++i )
    for ( unsigned j = 0; j < 3; ++j )
      M[3*i+j] = ( u[uijIdx(i,0)] * m.g[j] + u[uijIdx(i,1)] * m.g[3+j]
                   + u[uijIdx(i,2)] * m.g[6+j] );
  const double tr = M[0] + M[4] + M[8];
  const double c2 = ( M[0]*M[4] - M[1]*M[3] ) + ( M[0]*M[8] - M[2]*M[6] )
    + ( M[4]*M[8] - M[5]*M[7] );
  const double det = M[0]*( M[4]*M[8] - M[5]*M[7] )
    - M[1]*( M[3]*M[8] - M[5]*M[6] ) + M[2]*( M[3]*M[7] - M[4]*M[6] );
  //Depressed cubic t^3+p*t+q=0 for lambda=t+tr/3, solved trigonometrically:
  const double s = tr / 3.0;
  const double p = c2 - tr * s;
  const double q = -( 2.0 * s * s * s - c2 * s + det );
  std::array<double,3> res;
  if ( !( p < -1e-12 * s * s ) ) {
    res = { s, s, s };//(numerically) isotropic
    return res;
  }
  const double r = std::sqrt( -p / 3.0 );
  const double arg = ncclamp( ( 3.0 * q ) / ( 2.0 * p * r ), -1.0, 1.0 );
  const double phi = std::acos( arg ) / 3.0;
  for ( unsigned k = 0; k < 3; ++k )
    res[k] = s + 2.0 * r * std::cos( phi - k * ( 2.0 * kPi / 3.0 ) );
  std::sort( res.begin(), res.end() );
  return res;
}

std::vector<NC::PairDD> NCU::directionalMSDs( const Uij& u,
                                              const CellMetric& m )
{
  //6 point Gauss-Legendre on [-1,1]:
  static constexpr double x[6] = { -0.9324695142031521, -0.6612093864662645,
                                   -0.2386191860831969, 0.2386191860831969,
                                   0.6612093864662645, 0.9324695142031521 };
  static constexpr double w[6] = { 0.1713244923791704, 0.3607615730481386,
                                   0.4679139345726910, 0.4679139345726910,
                                   0.3607615730481386, 0.1713244923791704 };
  const auto ev = cartEigenvalues( u, m );
  std::vector<PairDD> res;
  res.reserve( 36 );
  for ( unsigned i = 0; i < 6; ++i ) {
    const double ct = 0.5 * ( x[i] + 1.0 );//cos(theta) in [0,1]
    const double st2 = 1.0 - ct * ct;
    for ( unsigned j = 0; j < 6; ++j ) {
      const double ph = 0.25 * kPi * ( x[j] + 1.0 );//phi in [0,pi/2]
      const double cp = std::cos( ph );
      const double cp2 = cp * cp;
      const double msd = ( ev[0] * st2 * cp2 + ev[1] * st2 * ( 1.0 - cp2 )
                           + ev[2] * ct * ct );
      res.emplace_back( msd, 0.25 * w[i] * w[j] );
    }
  }
  return res;
}
