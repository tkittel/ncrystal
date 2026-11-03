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

#include "NCSGExactOrbit.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {

    constexpr std::int64_t kMax = std::int64_t(1) << 53;

    //Whether the correctly rounded double value of p/q (|p|,q <= 2^53, q>0)
    //is below (-1), equal to (0) or above (+1) x:
    int roundedCmp( std::int64_t p, std::int64_t q, double x )
      ncnoexceptndebug
    {
      nc_assert( q > 0 && q <= kMax && p <= kMax && p >= -kMax );
      const double v = static_cast<double>( p ) / static_cast<double>( q );
      return v < x ? -1 : ( v > x ? 1 : 0 );
    }

    //The fraction p/q=(a+k*c)/(b+k*d) for k>=1 (k=1 is the mediant of a/b
    //and c/d, see rationalFromDouble), or false if p or q would exceed the
    //limits:
    bool stepValue( std::int64_t a, std::int64_t b,
                    std::int64_t c, std::int64_t d,
                    std::int64_t k, std::int64_t maxden,
                    std::int64_t& p, std::int64_t& q ) ncnoexceptndebug
    {
      nc_assert( k >= 1 && b >= 1 && d >= 1 && maxden >= b );
      if ( k > ( maxden - b ) / d )
        return false;
      q = b + k * d;//<= maxden
      const std::int64_t ac = ( c < 0 ? -c : c );
      if ( ac != 0 && k > ( kMax - ( a < 0 ? -a : a ) ) / ac )
        return false;
      p = a + k * c;
      return true;
    }
  }
}

NC::Optional<std::pair<std::int64_t,std::int64_t>>
NC::detail::rationalFromDouble( double x, std::int64_t maxDenominator )
{
  using result_t = std::pair<std::int64_t,std::int64_t>;
  if ( !std::isfinite( x ) || !( ncabs( x ) <= double( kMax ) )
       || maxDenominator < 1 )
    return NullOpt;
  maxDenominator = ncmin( maxDenominator, kMax );
  const double fl = std::floor( x );
  if ( fl == x )
    return result_t( static_cast<std::int64_t>( fl ), 1 );

  //Search of the Stern-Brocot tree (Graham, Knuth, Patashnik, "Concrete
  //Mathematics", 2nd ed., sec. 4.5, or
  //https://en.wikipedia.org/wiki/Stern%E2%80%93Brocot_tree): the fractions
  //which round to x form an interval, initially bracketed by the integers
  //a/b=floor(x)/1 and c/d=(floor(x)+1)/1. The "mediant" (a+c)/(b+d) always
  //lies between a/b and c/d, and has the smallest denominator of all
  //fractions between them. If it rounds to x it is the answer, otherwise it
  //replaces the bound on its side, and so on. Runs of replacements on the
  //same side are done at once, by finding the largest k for which
  //(a+k*c)/(b+k*d) (or (c+k*a)/(d+k*b)) stays on that side (exponential and
  //then binary search), so the number of iterations is logarithmic:
  std::int64_t a = static_cast<std::int64_t>( fl ), b = 1;
  std::int64_t c = a + 1, d = 1;
  while ( true ) {
    std::int64_t p, q;
    if ( !stepValue( a, b, c, d, 1, maxDenominator, p, q ) )
      return NullOpt;
    const int s = roundedCmp( p, q, x );
    if ( s == 0 )
      return result_t( p, q );//(mediants are always in lowest terms)
    //Moving the lower bound up (s<0), or the upper bound down (s>0):
    const std::int64_t lo_n = ( s < 0 ? a : c ), lo_d = ( s < 0 ? b : d );
    const std::int64_t hi_n = ( s < 0 ? c : a ), hi_d = ( s < 0 ? d : b );
    auto sameSide = [&]( std::int64_t k, std::int64_t& pk,
                         std::int64_t& qk )
    {
      return ( stepValue( lo_n, lo_d, hi_n, hi_d, k, maxDenominator, pk, qk )
               && roundedCmp( pk, qk, x ) == s );
    };
    std::int64_t kgood = 1;//(k=1 is the mediant just tested)
    std::int64_t kbad = 2;
    std::int64_t pk, qk;
    while ( sameSide( kbad, pk, qk ) ) {
      kgood = kbad;
      if ( kbad > kMax / 2 )
        return NullOpt;
      kbad *= 2;
    }
    while ( kbad - kgood > 1 ) {
      const std::int64_t kmid = kgood + ( kbad - kgood ) / 2;
      if ( sameSide( kmid, pk, qk ) )
        kgood = kmid;
      else
        kbad = kmid;
    }
    nc_assert_always( stepValue( lo_n, lo_d, hi_n, hi_d, kgood,
                                 maxDenominator, pk, qk ) );
    if ( s < 0 ) {
      a = pk;
      b = qk;
    } else {
      c = pk;
      d = qk;
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    //Exact counterparts of the above for rational coordinates, represented
    //by integer numerators over a common denominator:
    using IPos = std::array<std::int64_t,3>;

    //Largest denominator of rational coordinates recovered from doubles,
    //enough for all decimals with up to 8 digits after the decimal point:
    constexpr std::int64_t site_max_denominator = 100000000;

    //Limit on the common denominator, so all numerators in
    //expandSiteToOrbitExact stay below 2^53 (in magnitude):
    constexpr std::int64_t site_max_num = std::int64_t(1) << 53;
    constexpr std::int64_t site_max_common_denominator = site_max_num / 1024;

#ifndef NDEBUG
    //For debug checks of the bounds:
    bool numOK( std::int64_t n ) noexcept
    {
      return n < site_max_num && n > -site_max_num;
    }
#endif

    std::int64_t wrapNum( std::int64_t n, std::int64_t den ) ncnoexceptndebug
    {
      nc_assert( den > 0 && numOK( n ) );
      n %= den;
      return n < 0 ? n + den : n;
    }

    //Periodic difference of numerators, in [-den/2,den/2]:
    std::int64_t periodicDiffNum( std::int64_t a, std::int64_t b,
                                  std::int64_t den ) ncnoexceptndebug
    {
      const std::int64_t r = wrapNum( a - b, den );
      return 2 * r > den ? r - den : r;
    }

    double periodicDist( const IPos& a, const IPos& b, std::int64_t den )
      ncnoexceptndebug
    {
      nc_assert( den > 0 );
      std::int64_t m = 0;
      for ( unsigned i = 0; i < 3; ++i ) {
        const std::int64_t d = periodicDiffNum( a[i], b[i], den );
        m = ncmax( m, d < 0 ? -d : d );
      }
      return double( m ) / double( den );
    }

    //Correctly rounded (since |numerators| and den are at most 2^53):
    Vector toVector( const IPos& p, std::int64_t den ) ncnoexceptndebug
    {
      nc_assert( den > 0 && den <= site_max_num && numOK( p[0] )
                 && numOK( p[1] ) && numOK( p[2] ) );
      const double d = static_cast<double>( den );
      return Vector( static_cast<double>( p[0] ) / d,
                     static_cast<double>( p[1] ) / d,
                     static_cast<double>( p[2] ) / d );
    }

    //Apply op, with translations in units of 1/24 scaled to the denominator
    //(tscale=den/24):
    IPos applyOp( const SymOp& op, const IPos& p, std::int64_t tscale )
      ncnoexceptndebug
    {
      nc_assert( tscale >= 1 && tscale <= site_max_num / SymOp::tdenom );
      IPos r;
      for ( unsigned i = 0; i < 3; ++i ) {
        std::int64_t v = op.trans( i ) * tscale;
        for ( unsigned j = 0; j < 3; ++j ) {
          nc_assert( op.rot( i, j ) >= -1 && op.rot( i, j ) <= 1 );
          nc_assert( numOK( p[j] ) );
          v += op.rot( i, j ) * p[j];
        }
        nc_assert( numOK( v ) );
        r[i] = v;
      }
      return r;
    }
  }
}

//Same as expandSiteToOrbit (and with the same tolerances), but in exact
//integer arithmetic, so the resulting positions are the correctly
//rounded values of the exact positions (e.g. 1/10 and not 1-0.9). The
//site coordinates are recovered as rationals from the doubles (see
//rationalFromDouble), and represented with a common denominator D (a
//multiple of 24, D*nstab after symmetrisation). Bounds, given D <=
//2^53/1024: symmetrised sums are below 192*1.5*D, and images of numerators
//below D*nstab<=192*D are below 4*192*D. Returns false if the coordinates
//are not recognised as rationals, or if D would exceed the limit:
bool NC::detail::expandSiteToOrbitExact( const SGSymmetry& sym,
                                         const Vector& input_site,
                                         SGSiteOrbit& res )
{
  static_assert( site_max_common_denominator * 1024 <= site_max_num, "" );
  const double coords[3] = { input_site.x(), input_site.y(),
                             input_site.z() };
  std::pair<std::int64_t,std::int64_t> rc[3];//numerators and denominators
  std::int64_t den = SymOp::tdenom;
  for ( unsigned i = 0; i < 3; ++i ) {
    auto r = detail::rationalFromDouble( coords[i], site_max_denominator );
    if ( !r.has_value() )
      return false;
    rc[i] = r.value();
    nc_assert( rc[i].second >= 1 );
    //den = lcm(den,rc[i].second), if within the limit:
    std::int64_t a = den, b = rc[i].second;
    while ( b ) {
      const std::int64_t tmp = a % b;
      a = b;
      b = tmp;
    }
    nc_assert( a >= 1 );
    const std::int64_t f = rc[i].second / a;
    if ( den > site_max_common_denominator / f )
      return false;
    den *= f;
  }
  nc_assert( den % SymOp::tdenom == 0 );
  IPos site;
  for ( unsigned i = 0; i < 3; ++i ) {
    //(|numerator| <= 2^53 and the factor is at most den, so we reduce
    //modulo the denominator before scaling, to not overflow):
    const std::int64_t d = rc[i].second;
    site[i] = wrapNum( rc[i].first, d ) * ( den / d );
  }
  const Vector site_dbl = toVector( site, den );
  const auto& ops = sym.operations();

  //Symmetrise:
  IPos sum = { 0, 0, 0 };
  std::int64_t nstab = 0;
  for ( auto& op : ops ) {
    const IPos img = applyOp( op, site, den / SymOp::tdenom );
    const double dist = periodicDist( img, site, den );
    if ( dist < site_merge_tol ) {
      for ( unsigned i = 0; i < 3; ++i ) {
        sum[i] += site[i] + periodicDiffNum( img[i], site[i], den );
        nc_assert( numOK( sum[i] ) );
      }
      ++nstab;
    } else if ( dist < site_distinct_tol ) {
      throwNearSpecial( site_dbl, dist );
    }
  }
  nc_assert_always( nstab >= 1 && nstab <= 192 );
  const std::int64_t den2 = den * nstab;
  IPos symsite;
  for ( unsigned i = 0; i < 3; ++i )
    symsite[i] = wrapNum( sum[i], den2 );

  //Expand:
  std::vector<IPos> positions;
  for ( auto& op : ops ) {
    IPos img = applyOp( op, symsite, den2 / SymOp::tdenom );
    for ( auto& e : img )
      e = wrapNum( e, den2 );
    bool found = false;
    for ( auto& p : positions ) {
      const double dist = periodicDist( img, p, den2 );
      if ( dist < site_merge_tol ) {
        found = true;
        break;
      }
      if ( dist < site_distinct_tol )
        throwNearSpecial( site_dbl, dist );
    }
    if ( !found )
      positions.push_back( img );
  }
  nc_assert_always( positions.size() * nstab == ops.size() );
  res.positions.clear();
  res.positions.reserve( positions.size() );
  for ( auto& p : positions )
    res.positions.push_back( toVector( p, den2 ) );
  res.siteSymmetryOrder = static_cast<unsigned>( nstab );
  return true;
}
