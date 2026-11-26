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

#include "NCrystal/internal/vdos/NCVDOSDebyeXS.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCGaussLegendre.hh"
#include "NCVDOSDebyeXSData.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    namespace D = VDOSDebyeXSData;

    double mapv( double q, double lo, double hi )
    {
      return 2.0 * ( q - lo ) / ( hi - lo ) - 1.0;
    }

    //C2 smoothstep, 0 for t<=0 and 1 for t>=1:
    double smoothStep( double t )
    {
      t = ncmin( 1.0, ncmax( 0.0, t ) );
      return t * t * t * ( 10.0 + t * ( -15.0 + 6.0 * t ) );
    }

    //Legendre polynomials P_0..P_n at t:
    VectD legendreValues( double t, unsigned n )
    {
      VectD p( n + 1 );
      p[0] = 1.0;
      if ( n > 0 )
        p[1] = t;
      for ( unsigned k = 1; k < n; ++k )
        p[k+1] = ( ( 2*k + 1 ) * t * p[k] - k * p[k-1] ) / ( k + 1 );
      return p;
    }

    //sum_i c_i P_i(t), by Clenshaw recurrence:
    double legendreSeries( const VectD& c, double t )
    {
      double b1 = 0.0, b2 = 0.0;
      for ( std::size_t k = c.size(); k-- > 1; ) {
        const double b0 = c[k] + ( 2*k + 1 ) * t * b1 / ( k + 1 ) - ( k + 1 ) * b2 / ( k + 2 );
        b2 = b1;
        b1 = b0;
      }
      return c.front() + t * b1 - 0.5 * b2;
    }

    //Debye-Waller factor per unit alpha: lambda = int_0^x 3/x^3 b coth(b/2) db.
    double debyeLambda( double x )
    {
      auto f = []( double b )
      {
        return b > 1e-8 ? b / std::tanh( 0.5 * b ) : 2.0;
      };
      StableSum s;
      constexpr unsigned n = 16;
      for ( auto i : ncrange( n ) )
        s.add( integrateGaussLegendre8( f, x * i / n, x * ( i + 1 ) / n ) );
      return 3.0 * s.sum() / ( x * x * x );
    }

    //P(2,z) = 1-exp(-z)*(1+z) (regularised lower incomplete gamma):
    double p2( double z )
    {
      if ( z > 0.5 )
        return 1.0 - std::exp( -z ) * ( 1.0 + z );
      //series sum_{k>=2} (-1)^k (k-1) z^k / k!:
      double term = 0.5 * z * z, sum = term;
      for ( unsigned k = 3; k < 20; ++k ) {
        term *= -z / k;
        sum += ( k - 1 ) * term;
      }
      return sum;
    }

    //P(2,zp)-P(2,zm), zm<=zp:
    double p2diff( double zm, double zp )
    {
      if ( zm > 2.0 )
        return std::exp( -zm ) * ( 1.0 + zm ) - std::exp( -zp ) * ( 1.0 + zp );
      return p2( zp ) - p2( zm );
    }
  }
}

std::pair<NC::AtomMass,NC::AtomMass> NC::VDOS::VDOSDebyeXS::massRange()
{
  return { AtomMass{ D::Amin / const_inv_neutron_mass_amu },
           AtomMass{ D::Amax / const_inv_neutron_mass_amu } };
}

double NC::VDOS::VDOSDebyeXS::maxTOverTDebye()
{
  return 1.0 / D::xlow;
}

NC::VDOS::VDOSDebyeXS::VDOSDebyeXS( DebyeTemperature td,
                                    Temperature temp,
                                    AtomMass mass )
{
  td.validate();
  temp.validate();
  mass.validate();
  m_x = td.dbl() / temp.dbl();
  m_A = mass.dbl() * const_inv_neutron_mass_amu;
  if ( !( m_x >= D::xlow ) )
    NCRYSTAL_THROW2( CalcError, "Debye model cross sections are not supported"
                     " for T/TDebye="<<fmt(1.0/m_x)<<" (must be at most "
                     <<fmt(1.0/D::xlow)<<")." );
  if ( !( m_A >= D::Amin && m_A <= D::Amax ) )
    NCRYSTAL_THROW2( CalcError, "Debye model cross sections are not supported"
                     " for atomic mass "<<mass<<"." );
  m_kT = constant_boltzmann * temp.dbl();
  m_ktheta = constant_boltzmann * td.dbl();
  m_lam = debyeLambda( m_x );
  m_sigfree = ncsquare( m_A / ( m_A + 1.0 ) );

  //Fitted multi-phonon part: weights and u-series for fixed x and A (the fit
  //is constant to O(1/x^2) beyond xmax, the cold limit). The pieces for
  //x<xmin blend with those above over [xcblend_lo,xmin]:
  static_assert( sizeof(D::pieces) / sizeof(D::pieces[0]) == 6, "" );
  m_umax = D::umax * ncmax( 1.0, D::xcblend_lo / m_x );
  const double xc = ncmin( m_x, D::xmax );
  const double sx = smoothStep( std::log( xc / D::xblend_lo )
                                / std::log( D::xblend_hi / D::xblend_lo ) );
  const double sc = smoothStep( std::log( xc / D::xcblend_lo )
                                / std::log( D::xmin / D::xcblend_lo ) );
  const double ta = mapv( std::log( m_A + 1.0 ), std::log( D::Amin + 1.0 ),
                          std::log( D::Amax + 1.0 ) );
  for ( auto ip : ncrange( 6 ) ) {
    const auto& p = D::pieces[ip];
    auto& ser = m_series[ip];
    ser.wx = ( p.xlo == D::xlow
               ? 1.0 - sc
               : ( p.xlo == D::xcblend_lo ? sc * ( 1.0 - sx ) : sx ) );
    if ( !( ser.wx > 0.0 ) )
      continue;
    const auto lx = legendreValues( mapv( std::log( xc ), std::log( p.xlo ),
                                          std::log( p.xhi ) ), p.dx );
    const auto la = legendreValues( ta, p.da );
    ser.c.assign( p.du + 1, 0.0 );
    const double * c = &D::coefs[p.offset];
    for ( auto i : ncrange( p.du + 1 ) ) {
      double si = 0.0;
      for ( auto j : ncrange( p.dx + 1 ) )
        for ( auto k : ncrange( p.da + 1 ) )
          si += *c++ * lx[j] * la[k];
      ser.c[i] = si;
    }
  }
  m_sigUmax = sigma1( m_umax * m_ktheta / m_kT ) + sigmaMP( m_umax );
}

double NC::VDOS::VDOSDebyeXS::sigma1( double e ) const
{
  //One-phonon part, e=E/kT: A/(4e) int_{max(-x,-e)}^x G1(b) W1(b) db with
  //G1(b) = 3/(x^3 lam) b/(exp(b)-1), lam*W1 = P(2,lam*alpha+)-P(2,lam*alpha-)
  //and alpha+- = (sqrt(e)+-sqrt(e+b))^2/A, integrated in s=sqrt(e+b) (split
  //at b=0, i.e. at the energy loss/gain boundary):
  const double se = std::sqrt( e );
  const double lamA = m_lam / m_A;
  auto f = [se,e,lamA]( double s )
  {
    const double b = s * s - e;
    const double g = ncabs( b ) > 1e-12 ? b / std::expm1( b ) : 1.0 - 0.5 * b;
    return 2.0 * s * g * p2diff( lamA * ncsquare( se - s ), lamA * ncsquare( se + s ) );
  };
  //Breakpoints at b = 0, +-x*2^-k, down to the smallest feature scale: the
  //exp(-b) fall-off (width ~1) and the alpha band near b=0, where
  //lam*alpha- is O(1) for |b| ~ 2*sqrt(e*A/lam):
  const double blo = -ncmin( m_x, e );
  const double bmin = ncmin( 0.5, 2.0 * std::sqrt( e / lamA ) + 1e-300 );
  const int kmax = ncmin( 20, ncmax( 2, static_cast<int>( std::ceil( std::log2( m_x / bmin ) ) ) ) );
  //Energy gain beyond b=40 is negligible (exp(-40) relative):
  const double bhi = ncmin( m_x, 40.0 );
  VectD bpts = { blo, 0.0, bhi };
  for ( int k = 1; k <= kmax; ++k ) {
    const double b = std::ldexp( m_x, -k );
    if ( b < bhi )
      bpts.push_back( b );
    if ( -b > blo )
      bpts.push_back( -b );
  }
  std::sort( bpts.begin(), bpts.end() );
  StableSum s;
  double sprev = std::sqrt( ncmax( 0.0, e + bpts.front() ) );
  for ( auto it = std::next( bpts.begin() ); it != bpts.end(); ++it ) {
    const double snext = std::sqrt( e + *it );
    s.add( integrateGaussLegendre8( f, sprev, snext ) );
    sprev = snext;
  }
  return m_A / ( 4.0 * e ) * 3.0 / ( m_x * m_x * m_x * m_lam * m_lam ) * s.sum();
}

double NC::VDOS::VDOSDebyeXS::sigmaMP( double u ) const
{
  //Blend of the fitted pieces, giving log(sigma_ge2/sigma_free), for
  //0<u<=m_umax (the sqrt(u) pieces start at u=0):
  const double su = smoothStep( std::log( u / D::ublend_lo )
                                / std::log( D::ublend_hi / D::ublend_lo ) );
  nc_assert( u <= m_umax * ( 1.0 + 1e-12 ) );
  double y = 0.0;
  for ( auto ip : ncrange( 6 ) ) {
    const auto& p = D::pieces[ip];
    const auto& ser = m_series[ip];
    const double w = ser.wx * ( p.sqrtu ? 1.0 - su : su );
    if ( !( w > 0.0 ) )
      continue;
    if ( p.sqrtu ) {
      y += w * ( legendreSeries( ser.c, mapv( std::sqrt( u ), 0.0, std::sqrt( p.uhi ) ) )
                 - 0.5 * std::log( u ) );
    } else {
      y += w * legendreSeries( ser.c, mapv( std::log( u ), std::log( p.ulo ),
                                            std::log( p.uhi ) ) );
    }
  }
  return m_sigfree * std::exp( y );
}

NC::PairDD NC::VDOS::VDOSDebyeXS::sigmaParts( NeutronEnergy ekin ) const
{
  if ( !( ekin.dbl() > 0.0 ) )
    return { 0.0, 0.0 };
  const double s1 = sigma1( ekin.dbl() / m_kT );
  const double u = ekin.dbl() / m_ktheta;
  if ( u <= m_umax )
    return { s1, sigmaMP( u ) };
  //Above umax the total approaches the exact limit as sigma_free + c/E, with
  //c matching at umax:
  return { s1, m_sigfree + ( m_sigUmax - m_sigfree ) * ( m_umax / u ) - s1 };
}

double NC::VDOS::VDOSDebyeXS::sigma( NeutronEnergy ekin ) const
{
  auto p = sigmaParts( ekin );
  return p.first + p.second;
}

NC::IsoXSTable NC::VDOS::VDOSDebyeXS::createXSTabulation() const
{
  //Cubic Hermite pieces of g=sqrt(E)*sigma (the q coefficients of the table),
  //nodes log-spaced in u=E/(k theta) from E<=1e-4*kT (deep in the 1/v
  //regime) to umax (40 per decade from min(0.05,0.1/x) to 30, fewer
  //where sigma is smooth), plus u=1 and nodes crowding towards it from
  //both sides (the one-phonon part starts like (u-1)^1.5 there), with
  //three-point slopes (one-sided at u=1 and the ends). Tails: t0/sqrt(E)+t1
  //below (matching value and slope), and above the same sigma_free+c/E form
  //as used by sigma(..):
  const double ulo = ncmin( 1e-5, 1e-4 / m_x );
  auto logspaced = []( double a, double b, double perdecade )
  {
    return geomspace( a, b, static_cast<unsigned>( std::ceil( perdecade * std::log10( b / a ) ) ) + 1 );
  };
  const double udense = ncmin( 0.05, 0.1 / m_x );//below: E < 0.1 kT, smooth
  VectD us = logspaced( ulo, udense, 12.0 );
  for ( auto v : logspaced( udense, 30.0, 40.0 ) )
    us.push_back( v );
  for ( auto v : logspaced( 30.0, m_umax, 16.0 ) )
    us.push_back( v );
  us.push_back( 1.0 );
  for ( int k = 1; k <= 12; ++k ) {
    us.push_back( 1.0 + std::ldexp( 0.4, -k ) );
    us.push_back( 1.0 - std::ldexp( 0.4, -k ) );
  }
  std::sort( us.begin(), us.end() );
  us.erase( std::unique( us.begin(), us.end() ), us.end() );
  const std::size_t n = us.size();
  const std::size_t ikink = std::lower_bound( us.begin(), us.end(), 1.0 ) - us.begin();
  VectD E( n ), g( n );
  for ( auto i : ncrange( n ) ) {
    E[i] = us[i] * m_ktheta;
    g[i] = std::sqrt( E[i] ) * sigma( NeutronEnergy{ E[i] } );
  }
  //Slope at node i from nodes a<b<c (one of which is i):
  auto slope = [&E,&g]( std::size_t i, std::size_t a, std::size_t b, std::size_t c )
  {
    const double h1 = E[b] - E[a], h2 = E[c] - E[b];
    const double d1 = ( g[b] - g[a] ) / h1, d2 = ( g[c] - g[b] ) / h2;
    if ( i == a )
      return d1 - h1 * ( d2 - d1 ) / ( h1 + h2 );
    if ( i == c )
      return d2 + h2 * ( d2 - d1 ) / ( h1 + h2 );
    return ( h1 * d2 + h2 * d1 ) / ( h1 + h2 );
  };
  auto slopeFrom = [&slope,n]( std::size_t i, bool rightside )
  {
    //Neighbours on one side of the kink (or of the ends) only:
    if ( rightside || i == 0 )
      return slope( i, i, i + 1, i + 2 );
    if ( i + 1 == n )
      return slope( i, i - 2, i - 1, i );
    return slope( i, i - 1, i, i + 1 );
  };
  IsoXSTable res;
  res.E = E;
  res.q.reserve( 4 * ( n - 1 ) );
  for ( auto i : ncrange( n - 1 ) ) {
    const double h = E[i+1] - E[i];
    const double s0 = ( i == ikink ? slopeFrom( i, true ) : slopeFrom( i, false ) );
    const double s1 = ( i + 1 == ikink ? slope( i + 1, i - 1, i, i + 1 ) : slopeFrom( i + 1, false ) );
    const double d0 = h * s0, d1 = h * s1;
    res.q.insert( res.q.end(), { g[i], d0, 3.0 * ( g[i+1] - g[i] ) - 2.0 * d0 - d1,
                                 2.0 * ( g[i] - g[i+1] ) + d0 + d1 } );
  }
  //Below the grid: g = t0 + t1*sqrt(E), matching g and dg/dE at E[0]:
  const double se0 = std::sqrt( E.front() );
  res.tail[1] = 2.0 * se0 * slopeFrom( 0, true );
  res.tail[0] = g.front() - res.tail[1] * se0;
  //Above the grid: the same form as used by sigma(..):
  res.tail[3] = m_sigfree;
  res.tail[2] = ( m_sigUmax - m_sigfree ) * m_umax * m_ktheta;
  return res;
}
