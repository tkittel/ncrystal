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
  auto invalid = [a,b,c,alpha,beta,gamma]() {
    NCRYSTAL_THROW2(BadInput,"Invalid unit cell parameters: a,b,c="<<fmt(a)
                    <<","<<fmt(b)<<","<<fmt(c)<<" Aa, alpha,beta,gamma="
                    <<fmt(alpha)<<","<<fmt(beta)<<","<<fmt(gamma)<<" degrees");
  };
  auto angleOK = []( double x ) { return x > 0.0 && x < 180.0; };
  if ( !( a > 0.0 && b > 0.0 && c > 0.0 && std::isfinite( a + b + c )
          && angleOK( alpha ) && angleOK( beta ) && angleOK( gamma ) ) )
    invalid();
  //(reproducible, and exact for common angles like 90 and 120 degrees):
  const double ca = stable_sincos_deg( alpha ).cos;
  const double cb = stable_sincos_deg( beta ).cos;
  const double cg = stable_sincos_deg( gamma ).cos;
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
  if ( !( det > 0.0 ) )
    invalid();//(e.g. alpha=beta=gamma=150)
  m.recilen = { std::sqrt( c00 / det ), std::sqrt( c11 / det ),
                std::sqrt( c22 / det ) };
  return m;
}

NCU::Uij NCU::convertUstdToUstar( const Uij& u, const CellMetric& m )
{
  const auto& r = m.recilen;
  return { u[0]*r[0]*r[0], u[1]*r[1]*r[1], u[2]*r[2]*r[2],
           u[3]*r[0]*r[1], u[4]*r[0]*r[2], u[5]*r[1]*r[2] };
}

NCU::Uij NCU::convertUstarToUstd( const Uij& u, const CellMetric& m )
{
  const auto& r = m.recilen;
  return { u[0]/(r[0]*r[0]), u[1]/(r[1]*r[1]), u[2]/(r[2]*r[2]),
           u[3]/(r[0]*r[1]), u[4]/(r[0]*r[2]), u[5]/(r[1]*r[2]) };
}

NCU::Uij NCU::applySymOpToUstar( const SymOp& op, const Uij& u )
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

double NCU::calcUeqFromUstar( const Uij& u, const CellMetric& m )
{
  StableSum s;
  for ( unsigned i = 0; i < 3; ++i )
    for ( unsigned j = 0; j < 3; ++j )
      s.add( u[uijIdx( i, j )] * m.g[3*i+j] );
  return s.sum() / 3.0;
}

NC::PairDD NCU::calcUeqWithRoundingError( const Uij& u,
                                          const CellParameters& c,
                                          const Uij& uerr,
                                          const CellParameters& cerr )
{
  auto ueqIn = [&u]( const CellParameters& x ) {
    const auto m = cellMetric( x.a, x.b, x.c, x.alpha, x.beta, x.gamma );
    return calcUeqFromUstar( convertUstdToUstar( u, m ), m );
  };
  const auto m = cellMetric( c.a, c.b, c.c, c.alpha, c.beta, c.gamma );
  const double ueq = calcUeqFromUstar( convertUstdToUstar( u, m ), m );
  StableSum err;
  for ( auto i : ncrange( 6 ) ) {
    if ( !uerr[i] )
      continue;
    Uij unit{};
    unit[i] = 1.0;
    err.add( ncabs( calcUeqFromUstar( convertUstdToUstar( unit, m ), m ) )
             * uerr[i] );
  }
  double CellParameters::* const params[6]
    = { &CellParameters::a, &CellParameters::b, &CellParameters::c,
        &CellParameters::alpha, &CellParameters::beta,
        &CellParameters::gamma };
  for ( auto p : params ) {
    if ( !( cerr.*p ) )
      continue;
    CellParameters x = c;
    x.*p += cerr.*p;
    err.add( ncabs( ueqIn( x ) - ueq ) );
  }
  return { ueq, err.sum() };
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

std::array<double,3> NCU::calcUcartEigenvalues( const Uij& u,
                                                const CellMetric& m )
{
  //Eigenvalues of Ustar*G (similar to Ucart), from the characteristic
  //polynomial (all roots real and positive for positive definite Ustar):
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
  const auto ev = calcUcartEigenvalues( u, m );
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

namespace NCRYSTAL_NAMESPACE {
  namespace UijUtils {
    namespace {
      bool symOpMaps( const SymOp& op, const Vector& a, const Vector& b )
      {
        const Vector img = op.apply( a );
        for ( auto k : ncrange( 3 ) ) {
          const double d = img[k] - b[k];
          if ( !( ncabs( d - std::round( d ) ) < 1e-4 ) )
            return false;
        }
        return true;
      }

      //Structure of the tensors invariant under a site symmetry, from the
      //integer matrix P (6x6, row-major) summing the actions of its
      //operations on the six components: each component is either zero
      //(pivot -1), free (pivot itself), or num/den times an earlier free
      //component (its pivot). Simple if the free components are independent,
      //so these relations describe all invariant tensors:
      std::int64_t gcd64( std::int64_t a, std::int64_t b )
      {
        a = ( a < 0 ? -a : a );
        b = ( b < 0 ? -b : b );
        while ( b ) {
          const std::int64_t t = a % b;
          a = b;
          b = t;
        }
        return a;
      }

      struct UijStructure {
        std::array<int,6> pivot;
        std::array<std::int64_t,6> num, den;
        bool simple;
      };

      UijStructure analyseUijStructure( const std::array<std::int64_t,36>& P )
      {
        UijStructure res;
        unsigned nfree = 0;
        for ( auto i : ncrange( 6 ) ) {
          res.num[i] = res.den[i] = 1;
          const std::int64_t* pi = &P[6*i];
          bool zero = true;
          for ( auto l : ncrange( 6 ) )
            if ( pi[l] )
              zero = false;
          if ( zero ) {
            res.pivot[i] = -1;
            continue;
          }
          res.pivot[i] = i;
          for ( auto j : ncrange( i ) ) {
            if ( res.pivot[j] != j )
              continue;
            const std::int64_t* pj = &P[6*j];
            //Row i = (n/d) * row j, with n/d from the first nonzero entry:
            std::int64_t n = 0, d = 0;
            for ( auto l : ncrange( 6 ) ) {
              if ( pj[l] ) {
                n = pi[l];
                d = pj[l];
                break;
              }
            }
            bool prop = ( n != 0 );
            for ( auto l : ncrange( 6 ) )
              if ( pi[l] * d != pj[l] * n )
                prop = false;
            if ( !prop )
              continue;
            if ( d < 0 ) {
              n = -n;
              d = -d;
            }
            const std::int64_t a = gcd64( n, d );
            res.pivot[i] = j;
            res.num[i] = n / a;
            res.den[i] = d / a;
            break;
          }
          if ( res.pivot[i] == i )
            ++nfree;
        }
        //Rank of P (integer Gaussian elimination, each row divided by the gcd
        //of its entries after each step to keep them small):
        constexpr std::int64_t lim = std::int64_t(1) << 40;
        std::array<std::int64_t,36> m = P;
        unsigned rank = 0;
        for ( unsigned col = 0; col < 6 && rank < 6; ++col ) {
          unsigned piv = rank;
          while ( piv < 6 && m[6*piv+col] == 0 )
            ++piv;
          if ( piv == 6 )
            continue;
          for ( auto l : ncrange( 6u ) )
            std::swap( m[6*rank+l], m[6*piv+l] );
          for ( unsigned r = rank + 1; r < 6; ++r ) {
            const std::int64_t a = m[6*rank+col], b = m[6*r+col];
            if ( !b )
              continue;
            std::int64_t g = 0;
            for ( auto l : ncrange( 6u ) ) {
              m[6*r+l] = a * m[6*r+l] - b * m[6*rank+l];
              nc_assert_always( m[6*r+l] < lim && m[6*r+l] > -lim );
              g = gcd64( g, m[6*r+l] );
            }
            if ( g > 1 )
              for ( auto l : ncrange( 6u ) )
                m[6*r+l] /= g;
          }
          ++rank;
        }
        res.simple = ( rank == nfree );
        return res;
      }

      //Sum of the actions U -> R U R^T of the rotations on the six
      //components, as an integer matrix (row-major):
      std::array<std::int64_t,36>
      sumOfUijActions( Span<const SymOp* const> ops )
      {
        static constexpr unsigned ri[6] = { 0, 1, 2, 0, 0, 1 };
        static constexpr unsigned rj[6] = { 0, 1, 2, 1, 2, 2 };
        std::array<std::int64_t,36> P = {};
        for ( auto op : ops ) {
          //(R U R^T)_ij = sum_kl R_ik U_kl R_jl:
          for ( auto n : ncrange( 6u ) ) {
            for ( auto k : ncrange( 3u ) ) {
              const int rik = op->rot( ri[n], k );
              if ( rik )
                for ( auto l : ncrange( 3u ) )
                  P[6*n+uijIdx( k, l )] += rik * op->rot( rj[n], l );
            }
          }
        }
        return P;
      }

      double uijRelated( const UijStructure& st, unsigned i, double x )
      {
        return ( st.num[i] == st.den[i] ? x
                 : ( x * static_cast<double>( st.num[i] ) )
                 / static_cast<double>( st.den[i] ) );
      }
    }
  }
}

NCU::Uij NCU::mapUstdToPosition( const SGSymmetry& sym, const CellMetric& m,
                                 const Vector& from, const Vector& to,
                                 const Uij& ustd )
{
  for ( auto& op : sym.operations() ) {
    if ( !symOpMaps( op, from, to ) )
      continue;
    //Directly on Ustd if the rotation only mixes axes of equal reciprocal
    //lengths (exact for sign changes and permutations, and consistent with
    //snapUstdToSiteSymmetry):
    bool direct = true;
    for ( auto i : ncrange( 3u ) )
      for ( auto j : ncrange( 3u ) )
        if ( op.rot( i, j ) && m.recilen[i] != m.recilen[j] )
          direct = false;
    if ( direct )
      return applySymOpToUstar( op, ustd );
    return convertUstarToUstd( applySymOpToUstar( op, convertUstdToUstar( ustd,
                                                                         m ) ),
                               m );
  }
  NCRYSTAL_THROW2(BadInput,"No symmetry operation of space group "
                  <<sym.spaceGroup()<<" maps ("<<fmt(from[0],"%.6g")<<", "
                  <<fmt(from[1],"%.6g")<<", "<<fmt(from[2],"%.6g")<<") onto ("
                  <<fmt(to[0],"%.6g")<<", "<<fmt(to[1],"%.6g")<<", "
                  <<fmt(to[2],"%.6g")<<")");
}

bool NCU::detail::uijStructureIsExact( Span<const SymOp* const> ops )
{
  const auto st = analyseUijStructure( sumOfUijActions( ops ) );
  if ( !st.simple )
    return false;
  for ( auto i : ncrange( 6u ) ) {
    //(only factors which are powers of two are applied exactly):
    const std::int64_t n = ( st.num[i] < 0 ? -st.num[i] : st.num[i] );
    const std::int64_t d = st.den[i];
    if ( ( n & ( n - 1 ) ) || ( d & ( d - 1 ) ) )
      return false;
  }
  return true;
}

namespace NCRYSTAL_NAMESPACE {
  namespace UijUtils {
    namespace {
      UijResult uijError( std::string msg )
      {
        UijResult r;
        r.errorMsg = std::move( msg );
        return r;
      }
      Uij uijOrThrow( UijResult&& r )
      {
        if ( !r.uij.has_value() )
          NCRYSTAL_THROW( BadInput, r.errorMsg.value() );
        return r.uij.value();
      }
    }
  }
}

NCU::Uij NCU::snapUstdToSiteSymmetry( const SGSymmetry& sym,
                                      const CellMetric& m, const Vector& site,
                                      const Uij& ustd,
                                      const std::string& descr,
                                      double relTol, const Uij& err )
{
  return uijOrThrow( trySnapUstdToSiteSymmetry( sym, m, site, ustd, descr,
                                                relTol, err ) );
}

NCU::UijResult NCU::trySnapUstdToSiteSymmetry( const SGSymmetry& sym,
                                               const CellMetric& m,
                                               const Vector& site,
                                               const Uij& ustd,
                                               const std::string& descr,
                                               double relTol, const Uij& err )
{
  //The operations leaving the site invariant (at least the identity):
  std::array<const SymOp*,192> stabops;
  unsigned nstab = 0;
  for ( auto& op : sym.operations() ) {
    if ( symOpMaps( op, site, site ) ) {
      nc_assert_always( nstab < stabops.size() );
      stabops[nstab++] = &op;
    }
  }
  nc_assert_always( nstab > 0 );
  const Span<const SymOp* const> stab( stabops.data(),
                                       stabops.data() + nstab );
  //If the rotations only mix axes of equal reciprocal lengths (as for all
  //tabulated settings in cells obeying their constraints), they act on Ustd
  //exactly as on Ustar. The invariant tensors then have a structure of exact
  //zeros and fixed ratios of components (cf. UijStructure), imposed on the
  //average so the result is an exact fixed point (unchanged if snapped again):
  bool direct = true;
  for ( auto op : stab )
    for ( auto i : ncrange( 3u ) )
      for ( auto j : ncrange( 3u ) )
        if ( op->rot( i, j ) && m.recilen[i] != m.recilen[j] )
          direct = false;
  Uij snapped;
  bool done = false;
  if ( direct ) {
    const auto st = analyseUijStructure( sumOfUijActions( stab ) );
    if ( st.simple ) {
      bool invariant = true;
      for ( auto i : ncrange( 6u ) ) {
        const int p = st.pivot[i];
        if ( p < 0 ? ustd[i] != 0.0
             : ( p != static_cast<int>( i )
                 && ustd[i] != uijRelated( st, i, ustd[p] ) ) )
          invariant = false;
      }
      Uij avg = ustd;
      if ( !invariant ) {
        std::array<StableSum,6> sums;
        for ( auto op : stab ) {
          const auto t = applySymOpToUstar( *op, ustd );
          for ( auto i : ncrange( 6u ) )
            sums[i].add( t[i] );
        }
        for ( auto i : ncrange( 6u ) )
          avg[i] = sums[i].sum() / nstab;
      }
      for ( auto i : ncrange( 6u ) ) {
        const int p = st.pivot[i];
        snapped[i] = ( p < 0 ? 0.0
                       : ( p == static_cast<int>( i ) ? avg[i]
                           : uijRelated( st, i, avg[p] ) ) );
      }
      done = true;
    }
  }
  if ( !done ) {
    const auto uf = convertUstdToUstar( ustd, m );
    std::array<StableSum,6> sums;
    for ( auto op : stab ) {
      const auto t = applySymOpToUstar( *op, uf );
      for ( auto i : ncrange( 6u ) )
        sums[i].add( t[i] );
    }
    Uij ufsnap;
    for ( auto i : ncrange( 6u ) )
      ufsnap[i] = sums[i].sum() / nstab;
    snapped = convertUstarToUstd( ufsnap, m );
  }
  Uij diff;
  for ( auto i : ncrange( 6 ) )
    diff[i] = snapped[i] - ustd[i];
  if ( !( maxAbs( diff ) <= relTol * maxAbs( ustd ) + 2.0 * maxAbs( err ) ) )
    return uijError( "The uij values of the site " + descr
                     + " do not obey the site symmetry" );
  if ( !isPositiveDefinite( snapped ) )
    return uijError( "The uij values of the site " + descr
                     + " do not form a positive definite matrix" );
  UijResult res;
  res.uij = snapped;
  return res;
}

NCU::Uij NCU::ustdOfOrbit( const SGSymmetry& sym, const CellMetric& m,
                           const std::vector<Vector>& positions,
                           const std::vector<Uij>& ustds,
                           const std::string& descr, double relTol,
                           const Uij& err )
{
  return uijOrThrow( tryUstdOfOrbit( sym, m, positions, ustds, descr,
                                     relTol, err ) );
}

NCU::UijResult NCU::tryUstdOfOrbit( const SGSymmetry& sym, const CellMetric& m,
                                    const std::vector<Vector>& positions,
                                    const std::vector<Uij>& ustds,
                                    const std::string& descr, double relTol,
                                    const Uij& err )
{
  nc_assert_always( !positions.empty() && positions.size() == ustds.size() );
  auto r = trySnapUstdToSiteSymmetry( sym, m, positions.front(),
                                      ustds.front(), descr, relTol, err );
  if ( !r.uij.has_value() )
    return r;
  const Uij rep = r.uij.value();
  for ( auto i : ncrange( std::size_t(1), positions.size() ) ) {
    const Uij expected = mapUstdToPosition( sym, m, positions.front(),
                                            vectAt( positions, i ), rep );
    Uij diff;
    for ( auto k : ncrange( 6 ) )
      diff[k] = expected[k] - vectAt( ustds, i )[k];
    if ( !( maxAbs( diff ) <= relTol * maxAbs( rep ) + 2.0 * maxAbs( err ) ) )
      return uijError( "The uij values of the symmetry-equivalent positions"
                       " of " + descr + " are not symmetry images of each"
                       " other" );
  }
  return r;
}
