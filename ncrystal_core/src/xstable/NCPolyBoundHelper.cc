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

#include "NCPolyBoundHelper.hh"
#include "NCrystal/internal/utils/NCSpan.hh"

namespace NC = NCrystal;
namespace NCX = NCrystal::XSTable;

double NCX::Poly9::eval( double z ) const
{
  double r = m_c[m_deg];
  for ( int i = m_deg - 1; i >= 0; --i )
    r = r * z + m_c[i];
  return r;
}

NCX::Poly9& NCX::Poly9::operator+=( const Poly9& o )
{
  for ( int i = 0; i <= o.m_deg; ++i )
    m_c[i] += o.m_c[i];
  m_deg = std::max( m_deg, o.m_deg );
  return *this;
}

NCX::Poly9& NCX::Poly9::operator-=( const Poly9& o )
{
  for ( int i = 0; i <= o.m_deg; ++i )
    m_c[i] -= o.m_c[i];
  m_deg = std::max( m_deg, o.m_deg );
  return *this;
}

NCX::Poly9& NCX::Poly9::operator*=( double f )
{
  for ( int i = 0; i <= m_deg; ++i )
    m_c[i] *= f;
  return *this;
}

NCX::Poly9 NCX::operator+( Poly9 a, const Poly9& b )
{
  a += b;
  return a;
}

NCX::Poly9 NCX::operator-( Poly9 a, const Poly9& b )
{
  a -= b;
  return a;
}

NCX::Poly9 NCX::operator*( Poly9 a, double f )
{
  a *= f;
  return a;
}

NCX::Poly9 NCX::operator*( const Poly9& a, const Poly9& b )
{
  Poly9 r;
  r.m_deg = a.m_deg + b.m_deg;
  nc_assert_always( r.m_deg <= Poly9::maxDegree );
  for ( int i = 0; i <= a.m_deg; ++i )
    for ( int j = 0; j <= b.m_deg; ++j )
      r.m_c[i+j] += a.m_c[i] * b.m_c[j];
  return r;
}

NCX::Poly9 NCX::composeCubic( const double * c, const Poly9& u )
{
  Poly9 r( c[3] );
  for ( int k = 2; k >= 0; --k )
    r = r * u + Poly9( c[k] );
  return r;
}

void NCX::bernsteinCoefficients( const Poly9& p, double * b )
{
  //b_i = sum_{j<=i} binom(i,j)/binom(n,j) * a_j:
  constexpr int n = Poly9::maxDegree;
  static const double binom[n+1][n+1] = {
    {1},{1,1},{1,2,1},{1,3,3,1},{1,4,6,4,1},{1,5,10,10,5,1},
    {1,6,15,20,15,6,1},{1,7,21,35,35,21,7,1},{1,8,28,56,70,56,28,8,1},
    {1,9,36,84,126,126,84,36,9,1} };
  for ( int i = 0; i <= n; ++i ) {
    double s = 0.0;
    for ( int j = 0; j <= i && j <= p.degree(); ++j )
      s += binom[i][j] / binom[n][j] * p.coeff( j );
    b[i] = s;
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace XSTable {
    namespace {
      bool bernsteinNonNegativeImpl( const double * b, int depth )
      {
        constexpr int n = Poly9::maxDegree;
        if ( b[0] < 0.0 || b[n] < 0.0 )
          return false;
        bool allNonNeg = true;
        for ( int i = 1; i < n; ++i ) {
          if ( b[i] < 0.0 ) {
            allNonNeg = false;
            break;
          }
        }
        if ( allNonNeg )
          return true;
        if ( depth >= 30 )
          return false;
        //Subdivide at the midpoint (de Casteljau):
        double left[n+1], right[n+1], tmp[n+1];
        std::copy( b, b + n + 1, tmp );
        for ( int r = 0; r <= n; ++r ) {
          left[r] = tmp[0];
          right[n-r] = tmp[n-r];
          for ( int i = 0; i < n - r; ++i )
            tmp[i] = 0.5 * ( tmp[i] + tmp[i+1] );
        }
        return ( bernsteinNonNegativeImpl( left, depth + 1 )
                 && bernsteinNonNegativeImpl( right, depth + 1 ) );
      }
    }
  }
}

bool NCX::bernsteinNonNegative( const double * b )
{
  return bernsteinNonNegativeImpl( b, 0 );
}

namespace NCRYSTAL_NAMESPACE {
  namespace XSTable {
    namespace detail {
      //For testing (declared in the tests):
      VectD polyMultiply( Span<const double> a, Span<const double> b );
      VectD polyComposeCubic( Span<const double> cubic, Span<const double> u );
      VectD polyBernsteinCoefficients( Span<const double> p );
      bool polyBernsteinNonNegative( Span<const double> bernstein );
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace XSTable {
    namespace {
      Poly9 polyFromCoeffs( Span<const double> c )
      {
        nc_assert_always( !c.empty()
                          && c.size() <= static_cast<std::size_t>( Poly9::maxDegree + 1 ) );
        Poly9 res( c[0] );
        Poly9 zk( 1.0 );
        const Poly9 z( 0.0, 1.0 );
        for ( std::size_t i = 1; i < c.size(); ++i ) {
          zk = zk * z;
          res += zk * c[i];
        }
        return res;
      }
      VectD coeffsOf( const Poly9& p )
      {
        VectD res;
        for ( int i = 0; i <= p.degree(); ++i )
          res.push_back( p.coeff( i ) );
        return res;
      }
    }
  }
}

NC::VectD NCX::detail::polyMultiply( Span<const double> a, Span<const double> b )
{
  return coeffsOf( polyFromCoeffs( a ) * polyFromCoeffs( b ) );
}

NC::VectD NCX::detail::polyComposeCubic( Span<const double> cubic,
                                         Span<const double> u )
{
  nc_assert_always( cubic.size() == 4 );
  return coeffsOf( composeCubic( cubic.data(), polyFromCoeffs( u ) ) );
}

NC::VectD NCX::detail::polyBernsteinCoefficients( Span<const double> p )
{
  VectD res( Poly9::maxDegree + 1 );
  bernsteinCoefficients( polyFromCoeffs( p ), res.data() );
  return res;
}

bool NCX::detail::polyBernsteinNonNegative( Span<const double> bernstein )
{
  nc_assert_always( bernstein.size()
                    == static_cast<std::size_t>( Poly9::maxDegree + 1 ) );
  return bernsteinNonNegative( bernstein.data() );
}
