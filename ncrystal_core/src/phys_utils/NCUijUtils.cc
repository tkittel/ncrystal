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
