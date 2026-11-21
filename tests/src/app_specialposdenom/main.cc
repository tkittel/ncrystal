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

#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <iostream>
#include <map>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

//Strict upper bound on the denominators of the coordinates of special
//positions, in all settings. Fractional coordinates are written as fractions
//p/q with q up to SgFmt::maxFractionalCoordinateDenominator when they are
//such fractions within the tolerance (cf. SgFmt::fractionalCoordinate), so
//this bound must stay within it for all (snapped) special positions to be
//written as fractions. A point fixed by an operation (R,t) satisfies
//(R-I)x=-(t+n) with n integer. On the fixed set of a site symmetry, x_i is
//fixed iff e_i is a combination of rows a_k of the stacked R-I, and then
//x_i = sum_k c_k b_k (b_k=-(t_k+n_k), any independent rows with
//e_i = sum_k c_k a_k), so its denominator divides lcm_k(den(c_k) den_t(a_k)),
//den_t(a) being the lcm of the denominators of the translations of the
//operations with that row:

namespace {
  struct Frac {
    std::int64_t n, d;
    Frac( std::int64_t nn = 0, std::int64_t dd = 1 ) : n( nn ), d( dd )
    {
      nc_assert_always( d != 0 );
      if ( d < 0 ) {
        n = -n;
        d = -d;
      }
      const std::int64_t g = gcd( n < 0 ? -n : n, d );
      if ( g > 1 ) {
        n /= g;
        d /= g;
      }
    }
    static std::int64_t gcd( std::int64_t a, std::int64_t b )
    {
      while ( b ) {
        const std::int64_t t = a % b;
        a = b;
        b = t;
      }
      return a;
    }
  };
  Frac operator-( const Frac& a, const Frac& b )
  {
    return Frac( a.n * b.d - b.n * a.d, a.d * b.d );
  }
  Frac operator*( const Frac& a, const Frac& b )
  {
    return Frac( a.n * b.n, a.d * b.d );
  }
  Frac operator/( const Frac& a, const Frac& b )
  {
    nc_assert_always( b.n != 0 );
    return Frac( a.n * b.d, a.d * b.n );
  }

  using Row = std::array<int,3>;

  //Solve sum_k c_k rows[k] = target (3 equations, k unknowns) exactly.
  //Returns false if inconsistent or if the rows are dependent:
  bool solve( const std::vector<Row>& rows, const Row& target,
              std::vector<Frac>& c )
  {
    const std::size_t k = rows.size();
    //Augmented 3 x (k+1) matrix, equation m: sum_j rows[j][m] c_j = t[m]:
    std::vector<std::vector<Frac>> M( 3, std::vector<Frac>( k + 1 ) );
    for ( auto m : NC::ncrange( 3u ) ) {
      for ( auto j : NC::ncrange( k ) )
        M[m][j] = Frac( rows[j][m] );
      M[m][k] = Frac( target[m] );
    }
    std::size_t r = 0;
    std::vector<std::size_t> pivcol;
    for ( std::size_t col = 0; col < k; ++col ) {
      std::size_t piv = r;
      while ( piv < 3 && M[piv][col].n == 0 )
        ++piv;
      if ( piv == 3 )
        return false;//dependent rows
      std::swap( M[r], M[piv] );
      for ( auto m : NC::ncrange( 3u ) ) {
        if ( m == r || M[m][col].n == 0 )
          continue;
        const Frac f = M[m][col] / M[r][col];
        for ( auto j : NC::ncrange( k + 1 ) )
          M[m][j] = M[m][j] - f * M[r][j];
      }
      pivcol.push_back( col );
      ++r;
    }
    for ( std::size_t m = r; m < 3; ++m )
      if ( M[m][k].n != 0 )
        return false;//inconsistent
    c.assign( k, Frac() );
    for ( auto i : NC::ncrange( r ) )
      c[pivcol[i]] = M[i][k] / M[i][pivcol[i]];
    return true;
  }
}

void run()
{
  std::int64_t overall = 1;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const auto& sym = NC::SGSymmetry::get(
                             NC::SpaceGroup{ NC::SpaceGroupSettingIdx{ hn } } );
    std::map<Row,std::int64_t> dent;
    for ( auto& op : sym.operations() ) {
      for ( auto i : NC::ncrange( 3u ) ) {
        Row a;
        bool nonzero = false;
        for ( auto j : NC::ncrange( 3u ) ) {
          a[j] = op.rot( i, j ) - ( i == j ? 1 : 0 );
          if ( a[j] )
            nonzero = true;
        }
        if ( !nonzero )
          continue;
        const std::int64_t td = Frac( op.trans( i ), NC::SymOp::tdenom ).d;
        auto& e = dent[a];
        e = ( e ? e / Frac::gcd( e, td ) * td : td );
      }
    }
    std::vector<Row> rows;
    for ( auto& e : dent )
      rows.push_back( e.first );
    const std::size_t nr = rows.size();
    std::vector<Frac> c;
    auto consider = [&]( const std::vector<Row>& S )
    {
      for ( auto i : NC::ncrange( 3u ) ) {
        Row ei = { 0, 0, 0 };
        ei[i] = 1;
        if ( !solve( S, ei, c ) )
          continue;
        std::int64_t d = 1;
        for ( auto j : NC::ncrange( S.size() ) ) {
          const std::int64_t f = c[j].d * dent[S[j]];
          d = d / Frac::gcd( d, f ) * f;
        }
        overall = NC::ncmax( overall, d );
      }
    };
    for ( std::size_t a = 0; a < nr; ++a ) {
      consider( { rows[a] } );
      for ( std::size_t b = a + 1; b < nr; ++b ) {
        consider( { rows[a], rows[b] } );
        for ( std::size_t cc = b + 1; cc < nr; ++cc )
          consider( { rows[a], rows[b], rows[cc] } );
      }
    }
  }
  std::cout << "Largest possible denominator of special position"
            << " coordinates in all settings: " << overall << std::endl;
  REQUIRE( overall <= NC::SgFmt::maxFractionalCoordinateDenominator );
}

//Tests of SgFmt::fractionalCoordinate: within the tolerance (periodically,
//as "0" may stand for a value just below 1), idempotent, exact for atol=0, all
//fractions p/q up to the maximum denominator written as such, and rounding
//noise in small coordinates removed:
namespace test_fmt {

namespace {
  double parse( const std::string& t )
  {
    const NC::StrView sv( t );
    const auto pos = t.find( '/' );
    if ( pos == std::string::npos )
      return NC::str2dbl( sv );
    return ( NC::str2dbl( sv.substr( 0, pos ) )
             / NC::str2dbl( sv.substr( pos + 1 ) ) );
  }
  std::uint64_t s_lcg = 7;
  double rnd()
  {
    s_lcg = s_lcg * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>( s_lcg >> 11 ) / 9007199254740992.0;
  }
}

void run()
{
  using NC::SgFmt::fractionalCoordinate;
  constexpr unsigned maxden = NC::SgFmt::maxFractionalCoordinateDenominator;
  std::size_t n = 0;
  for ( double atol : { 1e-13, 1e-6, 1e-3, 0.0 } ) {
    for ( unsigned i = 0; i < 100000; ++i ) {
      double x;
      switch ( i % 4 ) {
      case 0: x = rnd(); break;
      case 1: x = std::floor( rnd() * 1e4 ) / 1e4; break;
      case 2: {
        const double q = std::floor( rnd() * 48 ) + 1;
        x = std::floor( rnd() * q ) / q;
        x += ( rnd() - 0.5 ) * 4e-16;
        break;
      }
      default: x = 1.0 - rnd() * 1e-12; break;
      }
      x = x - std::floor( x );
      if ( !( x < 1.0 ) )
        x = 0.0;
      const std::string t = fractionalCoordinate( x, atol );
      const double y = parse( t );
      REQUIRE( y >= 0.0 && y < 1.0 );
      REQUIRE( fractionalCoordinate( y, atol ) == t );
      const double d = NC::ncabs( y - x );
      const double pd = NC::ncmin( d, 1.0 - d );
      if ( atol == 0.0 )
        REQUIRE( y == x || ( t == "0" && x == 0.0 ) );
      else
        REQUIRE( pd <= atol * ( 1.0 + 1e-9 ) );
      ++n;
    }
  }
  for ( unsigned q = 2; q <= maxden; ++q ) {
    for ( unsigned p = 1; p < q; ++p ) {
      unsigned a = p, b = q;
      while ( b ) {
        const unsigned tmp = a % b;
        a = b;
        b = tmp;
      }
      if ( a != 1 )
        continue;
      std::ostringstream ss;
      ss << p << '/' << q;
      for ( double atol : { 1e-13, 1e-6, 0.0 } )
        REQUIRE( fractionalCoordinate( static_cast<double>( p ) / q, atol )
                 == ss.str() );
    }
  }
  std::cout << "fractionalCoordinate: " << n << " random coordinates OK,"
            << " fractions OK" << std::endl;
  for ( double x : { 0.0, 0.25, 1.0 / 3, 0.008389999999999898, 0.0458,
                     std::nextafter( 1.0, 0.0 ), 1.0 - 1e-12,
                     0.123456789012345 } )
    std::cout << "  " << NC::fmt( x ) << " -> "
              << fractionalCoordinate( x, 1e-13 ) << std::endl;
  REQUIRE( fractionalCoordinate( 0.008389999999999898, 1e-13 ) == ".00839" );
}

}

int main()
{
  run();
  test_fmt::run();
  return 0;
}
