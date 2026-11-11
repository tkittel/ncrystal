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

#include "NCrystal/internal/xstable/NCIsoMacroXSTable.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCIsoXSUtils.hh"
#include "NCPolyBoundHelper.hh"
#include <array>
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"

//FIXME: Get back to reproducibility for the IsoMacroXSTable tables: the
//selected table points depend on floating point details, which might differ
//between platforms and compilers (e.g. due to FMA contractions or
//auto-vectorisation).

namespace NC = NCrystal;
namespace NCX = NCrystal::XSTable;

namespace NCRYSTAL_NAMESPACE {
  namespace XSTable {
    namespace {

      //Relative distance from a discontinuity, at which the cross sections
      //just below and above it are evaluated. Discontinuities closer than
      //twice this are treated as a single one.
      constexpr double edge_eps = 1e-10;

      //The table is built such that it meets a smaller tolerance than
      //requested at all dense sample points and at the midpoints between
      //neighbouring ones. Between neighbouring dense sample points there are no
      //discontinuities, so the error of the table is smooth there and (at the
      //dense spacing) close to a quadratic. A quadratic within a tolerance T at
      //both ends and the midpoint of an interval is within 1.25*T in all of
      //it, so a factor below 0.8 ensures the requested tolerance everywhere
      //(0.7 leaves a margin):
      constexpr double tol_build_factor = 0.7;

      //Relative jump below which a known discontinuity is treated as a kink:
      constexpr double kink_tol = 1e-6;

      //Neighbouring dense sample points closer than this (relative) are not
      //refined further. If the tolerance is not met between such points, the
      //cross section has a discontinuity (or an extremely sharp feature) which
      //was not known in advance, and an exception is thrown:
      constexpr double min_rel_spacing = 1e-9;

      //Max number of refinement iterations (each iteration at least halves
      //the spacing of the points where the tolerance is not met):
      constexpr unsigned max_refine_iterations = 100;

      //Leaf processes without a table (createIsoXSTable), whose cross
      //sections are known to have no features which the table creation could
      //miss: they are smooth, apart from the boundaries of their energy domains
      //(which are included explicitly). Null processes and processes with
      //tables (smooth between their breakpts) need no such validation. Any
      //other process (e.g. a process with absorption resonances) gives a
      //warning, and must be studied before being added here.
      bool isValidatedUntabulated( const std::string& name )
      {
        static const std::set<std::string> validated = { "SABScatterLegacy",
                                                        "AbsOOV",
                                                        "SANSSphereScatter" };
        return validated.count( name ) > 0;
      }

      //Collect the names of the leaf processes (of ProcComposition trees), and
      //of those among them which need a warning (see isValidatedUntabulated):
      void inspectProcess( const ProcImpl::Process& proc,
                           std::set<std::string>& names,
                           std::set<std::string>& unvalidated )
      {
        auto comp = dynamic_cast<const ProcImpl::ProcComposition*>( &proc );
        if ( comp ) {
          for ( auto& c : comp->components() )
            inspectProcess( *c.process, names, unvalidated );
          return;
        }
        names.insert( proc.name() );
        if ( proc.isNull() )
          return;
        if ( proc.createIsoXSTable().has_value() )
          return;
        if ( !isValidatedUntabulated( proc.name() ) )
          unvalidated.insert( proc.name() );
      }

      //Collect the leaf processes (of ProcComposition trees, with their scale
      //factors), separating those with a table (see
      //Process::createIsoXSTable):
      struct UntabulatedLeaf {
        double scale;
        const ProcImpl::Process* process;
        CachePtr cache;
      };
      void collectLeaves( const ProcImpl::Process& proc, double scale,
                          std::vector<std::pair<double,IsoXSTable>>& tabulated,
                          std::vector<UntabulatedLeaf>& untabulated )
      {
        if ( auto comp = dynamic_cast<const ProcImpl::ProcComposition*>( &proc ) ) {
          for ( auto& c : comp->components() )
            collectLeaves( *c.process, scale * c.scale, tabulated, untabulated );
          return;
        }
        if ( proc.isNull() )
          return;
        auto t = proc.createIsoXSTable();
        if ( t.has_value() ) {
          tabulated.emplace_back( scale, std::move( t.value() ) );
          return;
        }
        untabulated.push_back( UntabulatedLeaf{ scale, &proc, nullptr } );
      }

      //Numerical derivative of a function which is smooth between the given
      //(sorted) discontinuities, at a wavelength which is not exactly at one
      //of them. The differences do not cross the nearest discontinuities:
      template<class TFct>
      double numericalSlope( const TFct& f, double wl, const VectD& disc )
      {
        auto it = std::upper_bound( disc.begin(), disc.end(), wl );
        const double lo = ( it == disc.begin() ? 0.0 : *std::prev( it ) );
        const double hi = ( it == disc.end() ? kInfinity : *it );
        double h = 1e-7 * wl;
        if ( wl - 2.0 * h > lo && wl + 2.0 * h < hi ) {
          const double w1 = wl - h;
          const double w2 = wl + h;
          return ( f( w2 ) - f( w1 ) ) / ( w2 - w1 );
        }
        //Second order one-sided differences, on the side with most room:
        const bool above = ( hi - wl > wl - lo );
        h = ncmin( h, 0.25 * ( above ? hi - wl : wl - lo ) );
        if ( !above )
          h = -h;
        return ( ncfma( 4.0, f( wl + h ),
                        ncfma( -3.0, f( wl ), -f( ncfma( 2.0, h, wl ) ) ) )
                 / ( 2.0 * h ) );
      }

      void validateParams( const IsoMacroXSTableParams& params )
      {
        if ( !( params.wlmax > 10.0 * denseWavelengthMin()
                && std::isfinite( params.wlmax ) ) )
          NCRYSTAL_THROW2( BadInput, "IsoMacroXSTable: invalid wlmax "
                           << params.wlmax << " (must be larger than "
                           << 10.0 * denseWavelengthMin() << " Aa)" );
        if ( !( params.tol > 0.0 && params.tol < 1.0 ) )
          NCRYSTAL_THROW2( BadInput, "IsoMacroXSTable: invalid tolerance "
                           << params.tol << " (must be in the range (0,1))" );
      }

      //For materials where all processes have a table,
      //the table is built directly from the merged table of the processes,
      //with the exact errors of the interpolation (see tabulateExact). Even so,
      //a smaller tolerance than requested is used, as a margin against rounding
      //errors:
      constexpr double tol_exact_factor = 0.9;

      struct TablePoints {
        VectD wl, xs, slope;
        std::size_t ncandidates = 0;
        std::size_t ndiscontinuities = 0;
      };

      //The cross section (scaled) of a tabulation and its derivative in a
      //given region: -1 below the table, i for the interval [E[i],E[i+1]] and
      //n-1 above the table (n=E.size()), also at the boundaries of the region:
      class TabulationRegions final {
      public:
        TabulationRegions( const IsoXSTable& t, double scale )
          : m_t( t ), m_scale( scale ), m_n( static_cast<long>( t.E.size() ) ) {}
        long n() const { return m_n; }
        const IsoXSTable& tab() const { return m_t; }
        double scale() const { return m_scale; }
        long region( double e ) const
        {
          if ( e < m_t.E.front() )
            return -1;
          if ( !( e < m_t.E.back() ) )
            return m_n - 1;
          return static_cast<long>( std::upper_bound( m_t.E.begin(), m_t.E.end(), e )
                                    - m_t.E.begin() ) - 1;
        }
        //Value and slope (per wavelength) at a wavelength:
        void eval( long k, double wl, double& y, double& s ) const
        {
          const double e = wl2ekin( wl );
          const auto& t = m_t;
          double xs, dxs;
          if ( k < 0 ) {
            const double se = std::sqrt( e );
            xs = t.tail[0] / se + t.tail[1];
            dxs = -0.5 * t.tail[0] / ( e * se );
          } else if ( k >= m_n - 1 ) {
            xs = t.tail[3] + ( t.tail[2] + t.tail[4] / e ) / e;
            dxs = -( t.tail[2] + 2.0 * t.tail[4] / e ) / ( e * e );
          } else {
            //e*xs = F = P(u) + sqrt(e)*Q(u) + r/e:
            const double w = t.E[k+1] - t.E[k];
            const double u = ( e - t.E[k] ) / w;
            double F = 0.0, dF = 0.0;
            if ( !t.p.empty() ) {
              const double * c = t.p.data() + 4 * k;
              F = evalCubicU( c, u );
              dF = ncfma( u, ncfma( 3.0 * u, c[3], 2.0 * c[2] ), c[1] ) / w;
            }
            if ( !t.q.empty() ) {
              const double * d = t.q.data() + 4 * k;
              const double se = std::sqrt( e );
              const double Q = evalCubicU( d, u );
              const double dQ = ncfma( u, ncfma( 3.0 * u, d[3], 2.0 * d[2] ),
                                       d[1] );
              F = ncfma( se, Q, F );
              dF += ncfma( se, dQ / w, 0.5 * Q / se );
            }
            if ( !t.r.empty() ) {
              F += t.r[k] / e;
              dF -= t.r[k] / ( e * e );
            }
            xs = F / e;
            dxs = ( dF - xs ) / e;
          }
          y = m_scale * xs;
          s = m_scale * dxs * ( -2.0 * e / wl );//de/dwl = -2e/wl
        }
      private:
        const IsoXSTable& m_t;
        double m_scale;
        long m_n;
      };

      //The table built directly from a tabulation (scaled by scale). The
      //candidate table points are the energies of the tabulation in the range,
      //where the cross section can be discontinuous or have kinks. An interval
      //between candidate points is accepted if the interpolation is within the
      //tolerance everywhere, which is decided exactly on each segment between
      //neighbouring candidate points: with t=1/wl=ta+D*z (z in [0,1]) and
      //e=K*t^2, both e^2*xs and e^2*H (H the interpolation) are polynomials
      //in z of degree at most 9, and so are the conditions
      //tol*e^2*xs -/+ e^2*(H-xs) >= 0, which are checked with their Bernstein
      //coefficients. Segments which can not be covered by a single interval
      //are split:
      TablePoints tabulateExact( const IsoXSTable& t, double scale,
                                 double wlmax, double tol )
      {
        const TabulationRegions tab( t, scale );
        const double wlmin = denseWavelengthMin();
        const double K = wl2ekin( 1.0 );
        const double sqrtK = std::sqrt( K );
        constexpr int nb = Poly9::maxDegree + 1;
        using BArr = std::array<double,nb>;

        //Candidate points, and the region of the segments between them:
        VectD lam = { wlmin, wlmax };
        const double ehi = wl2ekin( wlmin ), elo = wl2ekin( wlmax );
        for ( auto e : t.E )
          if ( e > elo && e < ehi )
            lam.push_back( ekin2wl( e ) );
        std::sort( lam.begin(), lam.end() );
        lam.erase( std::unique( lam.begin(), lam.end() ), lam.end() );
        std::vector<long> reg;
        for ( std::size_t m = 0; m + 1 < lam.size(); ++m )
          reg.push_back( tab.region( wl2ekin( 0.5 * ( lam[m] + lam[m+1] ) ) ) );

        //Bernstein coefficients of e^2*xs on segment m (z=0 at the long
        //wavelength end):
        auto computeBG = [&]( std::size_t m )
        {
          const double x0 = lam[m], x1 = lam[m+1];
          const long k = reg[m];
          const double ta = 1.0 / x1;
          const double D = ( x1 - x0 ) / ( x0 * x1 );
          const Poly9 tz( ta, D );
          const Poly9 ez = tz * tz * K;
          const Poly9 e15z = ez * tz * sqrtK;//e^(3/2)
          //e^2*xs = P*e + Q*e^(3/2) + R:
          Poly9 G;
          if ( k < 0 ) {
            G = ez * ez * t.tail[1] + e15z * t.tail[0];
          } else if ( k >= tab.n() - 1 ) {
            G = ez * ez * t.tail[3] + ez * t.tail[2] + Poly9( t.tail[4] );
          } else {
            const double Ek = t.E[k], w = t.E[k+1] - Ek;
            //u=(e-Ek)/w, with the constant term computed directly (not as a
            //difference of the coefficients):
            const Poly9 uz = ( ( ez - Poly9( ez.coeff( 0 ) ) ) * ( 1.0 / w )
                                   + Poly9( ncfma( K * ta, ta, -Ek ) / w ) );
            if ( !t.p.empty() )
              G += composeCubic( t.p.data() + 4 * k, uz ) * ez;
            if ( !t.q.empty() )
              G += composeCubic( t.q.data() + 4 * k, uz ) * e15z;
            if ( !t.r.empty() )
              G += Poly9( t.r[k] );
          }
          BArr r;
          bernsteinCoefficients( G * scale, r.data() );
          return r;
        };
        std::vector<BArr> BG;
        BG.reserve( lam.size() );
        for ( std::size_t m = 0; m + 1 < lam.size(); ++m )
          BG.push_back( computeBG( m ) );

        //Bernstein coefficients of K^2*(1-xa*t)^j*t^(4-j) (j=0..3) on segment
        //m, so e^2*H is their combination with the coefficients of H in powers
        //of wl-xa (cached for the current start point xa):
        struct Basis { BArr b[4]; };
        std::vector<Basis> basisCache;
        std::size_t basisStart = lam.size();
        auto basis = [&]( std::size_t a, std::size_t m ) -> const Basis&
        {
          if ( basisStart != a ) {
            basisCache.clear();
            basisStart = a;
          }
          while ( basisCache.size() <= m - a ) {
            const std::size_t mm = a + basisCache.size();
            const double x0 = lam[mm], x1 = lam[mm+1], xa = lam[a];
            const double ta = 1.0 / x1;
            const double D = ( x1 - x0 ) / ( x0 * x1 );
            const Poly9 tz( ta, D );
            const Poly9 oneminus( ( x1 - xa ) / x1, -xa * D );//1-xa*t
            Poly9 tpow[5];
            tpow[0] = Poly9( 1.0 );
            for ( int i = 1; i <= 4; ++i )
              tpow[i] = tpow[i-1] * tz;
            Poly9 ompow( 1.0 );
            Basis bs;
            for ( int j = 0; j <= 3; ++j ) {
              bernsteinCoefficients( ompow * tpow[4-j] * ( K * K ), bs.b[j].data() );
              ompow = ompow * oneminus;
            }
            basisCache.push_back( bs );
          }
          return basisCache[m - a];
        };

        //Whether the interpolation of the interval [lam[a],xb] is within the
        //tolerance on segment m:
        auto segmentOK = [&]( std::size_t a, std::size_t m, double xb,
                              double ya, double yb, double sa, double sb )
        {
          const double h = xb - lam[a];
          const double dy = ( yb - ya ) / h;
          const double A[4] = { ya, sa,
                                ncfma( 3.0, dy, -ncfma( 2.0, sa, sb ) ) / h,
                                ncfma( -2.0, dy, sa + sb ) / ( h * h ) };
          const Basis& bs = basis( a, m );
          const BArr& bg = BG[m];
          BArr bplus, bminus;
          bool plusOK = true, minusOK = true;
          for ( int i = 0; i < nb; ++i ) {
            const double eh = ncfma( A[3], bs.b[3][i],
                                     ncfma( A[2], bs.b[2][i],
                                            ncfma( A[1], bs.b[1][i],
                                                   A[0] * bs.b[0][i] ) ) );
            bplus[i] = ncfma( 1.0 + tol, bg[i], -eh );
            bminus[i] = ncfma( tol - 1.0, bg[i], eh );
            plusOK = plusOK && bplus[i] >= 0.0;
            minusOK = minusOK && bminus[i] >= 0.0;
          }
          return ( ( plusOK || bernsteinNonNegative( bplus.data() ) )
                   && ( minusOK || bernsteinNonNegative( bminus.data() ) ) );
        };

        //Values and slopes on each side of candidate point m. Differences at
        //the level of rounding errors (e.g. where the tables of the processes
        //are continuous) are not treated as discontinuities or kinks (the
        //latter would add a repeated table point, depending on rounding). The
        //slopes are compared by their effect on the value over a wavelength
        //range of the size of the wavelength:
        auto below = [&]( std::size_t m, double& y, double& s ) { tab.eval( reg[m-1], lam[m], y, s ); };
        auto above = [&]( std::size_t m, double& y, double& s )
        {
          tab.eval( reg[m], lam[m], y, s );
          if ( m > 0 ) {
            double yb, sb;
            below( m, yb, sb );
            if ( ncabs( y - yb ) <= 1e-12 * ncmax( y, yb ) ) {
              y = yb;
              if ( ncabs( s - sb ) * lam[m] <= 1e-12 * ncmax( y, yb ) )
                s = sb;
            }
          }
        };

        //Check the interval from a to b, with the slopes limited such that the
        //interpolation can not become negative:
        auto check = [&]( std::size_t a, std::size_t b, double& s0, double& s1 )
        {
          double ya, sa, yb, sb;
          above( a, ya, sa );
          below( b, yb, sb );
          const double h = lam[b] - lam[a];
          s0 = ncmax( sa, -3.0 * ya / h );
          s1 = ncmin( sb, 3.0 * yb / h );
          for ( std::size_t m = a; m < b; ++m )
            if ( !segmentOK( a, m, lam[b], ya, yb, s0, s1 ) )
              return false;
          return true;
        };

        TablePoints res;
        auto addPoint = [&res]( double x, double y, double s )
        {
          res.wl.push_back( x );
          res.xs.push_back( y );
          res.slope.push_back( s );
        };
        //Wavelength 0 is the limit at high energies (from the tail above the
        //table), with the chord towards wlmin as the slope:
        double y0 = scale * t.tail[3], ymin, smin;
        above( 0, ymin, smin );
        const double chord = ( ymin - y0 ) / wlmin;
        addPoint( 0.0, y0, chord );
        addPoint( wlmin, ymin, chord );

        //Greedy selection of the intervals, with an exponential and then a
        //binary search for the furthest end point:
        std::size_t a = 0;
        while ( a + 1 < lam.size() ) {
          double s0, s1;
          while ( !check( a, a + 1, s0, s1 ) ) {
            const double xm = 0.5 * ( lam[a] + lam[a+1] );
            if ( !( xm > lam[a] && xm < lam[a+1] ) )
              NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: the tolerance can"
                               " not be met at wavelength " << xm << " Aa" );
            lam.insert( lam.begin() + a + 1, xm );
            reg.insert( reg.begin() + a + 1, reg[a] );
            BG[a] = computeBG( a );
            BG.insert( BG.begin() + a + 1, computeBG( a + 1 ) );
            basisStart = lam.size();
          }
          std::size_t good = a + 1, bad = lam.size();
          double good_s0 = s0, good_s1 = s1;
          std::size_t step = 2;
          while ( true ) {
            const std::size_t j = std::min( a + step, lam.size() - 1 );
            if ( j <= good )
              break;
            if ( check( a, j, s0, s1 ) ) {
              good = j;
              good_s0 = s0;
              good_s1 = s1;
              if ( j == lam.size() - 1 )
                break;
              step *= 2;
            } else {
              bad = j;
              break;
            }
          }
          while ( bad < lam.size() && bad - good > 1 ) {
            const std::size_t mid = good + ( bad - good ) / 2;
            if ( check( a, mid, s0, s1 ) ) {
              good = mid;
              good_s0 = s0;
              good_s1 = s1;
            } else {
              bad = mid;
            }
          }
          double ya, sa, yb, sb;
          above( a, ya, sa );
          if ( !( res.wl.back() == lam[a] && res.xs.back() == ya
                  && res.slope.back() == good_s0 ) )
            addPoint( lam[a], ya, good_s0 );
          below( good, yb, sb );
          addPoint( lam[good], yb, good_s1 );
          a = good;
        }
        //Known discontinuities (jumps, not kinks) at the candidate points:
        for ( std::size_t m = 1; m + 1 < lam.size(); ++m ) {
          double yl, sl, yr, sr;
          below( m, yl, sl );
          above( m, yr, sr );
          if ( ncabs( yr - yl ) > kink_tol * ncmax( yl, yr ) )
            ++res.ndiscontinuities;
        }
        res.ncandidates = lam.size();
        return res;
      }

      //The wavelengths of the known discontinuities and kinks of the cross
      //section of a process (see Process::createIsoXSBreakpoints,
      //e.g. the Bragg edges or the boundaries of its energy domain):
      void collectBreakpoints( const ProcImpl::Process& proc, VectD& wl )
      {
        const auto bp = proc.createIsoXSBreakpoints();
        wl.reserve( wl.size() + bp.size() );
        for ( auto e : bp )
          if ( e > 0.0 && std::isfinite( e ) )
            wl.push_back( ekin2wl( e ) );
      }

      //Sorted discontinuities within (0,wlmax), with clusters of nearby ones
      //represented by their largest member (for Bragg edges, the cross section
      //drops to the value above the whole cluster there).
      VectD mergeDiscontinuities( VectD d, double wlmax )
      {
        std::sort( d.begin(), d.end() );
        VectD res;
        res.reserve( d.size() );
        for ( auto v : d ) {
          if ( !( v > 0.0 && v < wlmax*(1.0-4*edge_eps) ) )
            continue;
          if ( !res.empty() && v <= res.back() * ( 1.0 + 2*edge_eps ) )
            res.back() = v;//merge with previous
          else
            res.push_back( v );
        }
        return res;
      }

      //The cubic Hermite interpolation between (x0,y0) and (x1,y1), with the
      //slopes s0 and s1 there:
      inline double hermite( double x0, double x1, double y0, double y1,
                             double s0, double s1, double x )
      {
        const double h = x1 - x0;
        const double u = ( x - x0 ) / h;
        const double v = 1.0 - u;
        const double t0 = ncfma( ncfma( 2.0, u, 1.0 ), y0, u * h * s0 );
        const double t1 = ncfma( ncfma( -2.0, u, 3.0 ), y1, -( v * h * s1 ) );
        return ncfma( v * v, t0, u * u * t1 );
      }

      //Try to cover the sample points from i to j with one interval, given the
      //slope s0 at i. The slope s1 at j is the slope of the cross section
      //there, limited such that the interval and the following one (at least
      //up to the next sample point) can not become negative: h*s0 >= -3*y0 and
      //h*s1 <= 3*y1 for an interval of width h are sufficient for that
      //(Schmidt and Hess 1988). Returns whether all interior points are
      //reproduced within tol:
      bool fitInterval( const VectD& x, const VectD& y, const VectD& slopes,
                        std::size_t i, std::size_t j, double s0, double tol,
                        double& s1 )
      {
        nc_assert( j > i && x[j] > x[i] );
        const double h = x[j] - x[i];
        s1 = ncmin( slopes[j], 3.0 * y[j] / h );
        if ( j + 1 < x.size() && x[j+1] > x[j] )
          s1 = ncmax( s1, -3.0 * y[j] / ( x[j+1] - x[j] ) );
        if ( h * s0 < -3.0 * y[i] * ( 1.0 + 1e-12 ) )
          return false;
        for ( std::size_t k = i + 1; k < j; ++k ) {
          const double yk = hermite( x[i], x[j], y[i], y[j], s0, s1, x[k] );
          if ( !( ncabs( yk - y[k] ) <= tol * y[k] ) )
            return false;
        }
        return true;
      }

    }
  }
}

NCX::ReducedHermitePoints NCX::reduceHermitePoints( const VectD& x,
                                                    const VectD& y,
                                                    const VectD& slopes,
                                                    double tol )
{
  //From the current point (with its slope fixed), find the most distant point
  //which can be reached with an interval satisfying the tolerance. An
  //exponential search followed by a binary search is used, which assumes that
  //feasibility is (approximately) monotonic in the end point. Every selected
  //interval is explicitly verified, so the tolerance is always respected at all
  //the sample points. An interval can span small discontinuities (both values
  //must then be within the tolerance), otherwise a new point is started at
  //the discontinuity, with a fresh slope.
  const std::size_t n = x.size();
  nc_assert_always( n >= 2 && y.size() == n && slopes.size() == n );
  auto startSlope = [&x,&y,&slopes,n]( std::size_t i )
  {
    //The slope of the cross section, limited such that the interval up to
    //the next sample point can not become negative:
    return ( i + 1 < n && x[i+1] > x[i]
             ? ncmax( slopes[i], -3.0 * y[i] / ( x[i+1] - x[i] ) )
             : slopes[i] );
  };
  ReducedHermitePoints res;
  res.idx.push_back( 0 );
  res.slope.push_back( startSlope( 0 ) );
  std::size_t i = 0;
  while ( i + 1 < n ) {
    if ( !( x[i+1] > x[i] ) ) {
      //Discontinuity:
      ++i;
      res.idx.push_back( i );
      res.slope.push_back( startSlope( i ) );
      continue;
    }
    auto reach = [&]( double s0, double& s_good )
    {
      std::size_t good = i + 1;
      (void)fitInterval( x, y, slopes, i, good, s0, tol, s_good );
      std::size_t bad = n;//first known infeasible end point (n: none)
      std::size_t step = 2;
      double s;
      while ( true ) {
        std::size_t j = ncmin( i + step, n - 1 );
        if ( j <= good )
          break;
        if ( x[j] > x[i] && fitInterval( x, y, slopes, i, j, s0, tol, s ) ) {
          good = j;
          s_good = s;
          if ( j == n - 1 )
            break;
          step *= 2;
        } else {
          bad = j;
          break;
        }
      }
      if ( bad < n ) {
        while ( bad - good > 1 ) {
          std::size_t mid = good + ( bad - good ) / 2;
          if ( fitInterval( x, y, slopes, i, mid, s0, tol, s ) ) {
            good = mid;
            s_good = s;
          } else {
            bad = mid;
          }
        }
      }
      return good;
    };
    double s_good;
    std::size_t good = reach( res.slope.back(), s_good );
    if ( good == i + 1 ) {
      //The slope inherited from the previous interval (where it might have
      //been limited to avoid negative values) does not reach beyond the next
      //sample point. Restart with the slope of the cross section instead, by
      //repeating the point (with the same value):
      const double s0 = startSlope( i );
      if ( s0 != res.slope.back() ) {
        res.idx.push_back( i );
        res.slope.push_back( s0 );
        good = reach( s0, s_good );
      }
    }
    res.idx.push_back( good );
    res.slope.push_back( s_good );
    i = good;
  }
  return res;
}

unsigned NCX::autoNDense( double tol, double wlmax )
{
  //10000 points per decade for tol=1e-3. The segments of the table get
  //shorter for smaller tolerances (roughly as sqrt(tol) in smooth regions), so
  //the dense sample must get denser:
  nc_assert_always( tol > 0.0 && wlmax > denseWavelengthMin() );
  const double ndecades = std::log10( wlmax / denseWavelengthMin() );
  const double n = ndecades * 10000.0 * std::sqrt( 1e-3 / tol );
  return static_cast<unsigned>( ncclamp( n, 1000.0, 1e7 ) );
}

double NCX::evalIsoMacroXSTable( const double* wl, const double* xs,
                                 const double* slope, std::size_t n,
                                 double wavelength )
{
  //Evaluate the table. This is both used by NCrystal itself, but it is also a
  //reference implementation for people using our tables elsewhere.
  //
  //Note for NCrystal developers: Do not change this reference implementation in
  //                              a way that changes its behaviour. It might
  //                              have been copied elsewhere outside our
  //                              code-base!
  if ( n == 0 )
    return 0.0;
  if ( !( wavelength > wl[0] ) || n == 1 )
    return xs[0];
  if ( wavelength >= wl[n-1] ) {
    //Linear extrapolation of the value and slope of the last point, clamped
    //at 0:
    if ( slope[n-1] == 0.0 )
      return xs[n-1];
    const double v = ncfma( wavelength - wl[n-1], slope[n-1], xs[n-1] );
    return v > 0.0 ? v : 0.0;
  }
  //Last point with wl[i] <= wavelength (the second of a pair at a
  //discontinuity):
  const std::size_t i
    = static_cast<std::size_t>( std::upper_bound( wl, wl + n,
                                                  wavelength ) - wl ) - 1;
  nc_assert( i + 1 < n && wl[i+1] > wl[i] );
  return hermite( wl[i], wl[i+1], xs[i], xs[i+1], slope[i], slope[i+1],
                  wavelength );
}

NCX::IsoMacroXSTable
NCX::createIsoMacroXSTable( const std::function<double(double)>& xs_of_wl,
                            const std::function<double(double)>& dxs_of_wl,
                            VectD discontinuities_in,
                            const VectD& samplePoints,
                            const IsoMacroXSTableParams& params,
                            IsoMacroXSTableDiagnostics* diagnostics )
{
  validateParams( params );
  if ( params.ndense == 1 )
    NCRYSTAL_THROW( BadInput, "IsoMacroXSTable: ndense must be at least 2" );
  const unsigned ndense = ( params.ndense
                            ? params.ndense
                            : autoNDense( params.tol, params.wlmax ) );
  const double tol = params.tol;
  const double tol_build = tol_build_factor * tol;

  auto evalXS = [&xs_of_wl]( double wl )
  {
    const double v = xs_of_wl( wl );
    if ( !( std::isfinite( v ) && v >= 0.0 ) )
      NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: invalid cross section " << v
                       << " at wavelength " << wl << " Aa" );
    return v;
  };
  auto evalSlope = [&dxs_of_wl]( double wl )
  {
    const double v = dxs_of_wl( wl );
    if ( !std::isfinite( v ) )
      NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: invalid slope " << v
                       << " of the cross section at wavelength " << wl
                       << " Aa" );
    return v;
  };

  const VectD discontinuities = ( params.discontinuities
                                  ? mergeDiscontinuities( std::move( discontinuities_in ),
                                                          params.wlmax )
                                  : VectD() );

  //Dense sample: wavelength 0, log-spaced points, and points just below and
  //above each discontinuity. The x-values are the wavelengths, where a
  //discontinuity contributes two identical x-values (with the cross sections
  //below and above it).
  VectD wl_eval;//wavelengths at which to evaluate
  VectD x;//the corresponding x values of the dense sample
  wl_eval.reserve( 1 + ndense + 2 * discontinuities.size() );
  x.reserve( 1 + ndense + 2 * discontinuities.size() );
  wl_eval.push_back( 0.0 );//NB: the limit is calculated below
  x.push_back( 0.0 );
  {
    VectD grid = geomspace( denseWavelengthMin(), params.wlmax, ndense );
    for ( auto wl : samplePoints )
      if ( wl > denseWavelengthMin() && wl < params.wlmax )
        grid.push_back( wl );
    std::sort( grid.begin(), grid.end() );
    grid.erase( std::unique( grid.begin(), grid.end() ), grid.end() );
    auto itD = discontinuities.begin();
    for ( auto wl : grid ) {
      for ( ; itD != discontinuities.end() && *itD <= wl*(1.0+2*edge_eps); ++itD ) {
        const double e = *itD;
        if ( x.size() > 1 && x.back() >= e*(1.0-2*edge_eps) ) {
          //previous grid point too close to the discontinuity, remove it:
          x.pop_back();
          wl_eval.pop_back();
        }
        wl_eval.push_back( e * ( 1.0 - edge_eps ) );
        x.push_back( e );
        wl_eval.push_back( e * ( 1.0 + edge_eps ) );
        x.push_back( e );
      }
      if ( x.size() > 1 && wl <= x.back()*(1.0+2*edge_eps) )
        continue;//grid point too close to the previous discontinuity
      wl_eval.push_back( wl );
      x.push_back( wl );
    }
    nc_assert_always( itD == discontinuities.end() );
  }

  //The cross sections and their slopes (the slope at wavelength 0 is set
  //below):
  VectD y, slopes;
  y.reserve( wl_eval.size() );
  slopes.reserve( wl_eval.size() );
  y.push_back( 0.0 );
  slopes.push_back( 0.0 );
  for ( std::size_t i = 1; i < wl_eval.size(); ++i ) {
    y.push_back( evalXS( wl_eval[i] ) );
    slopes.push_back( evalSlope( wl_eval[i] ) );
  }
  wl_eval.clear();

  //Known discontinuities where the cross section turns out to be continuous
  //(kinks, where the values just below and above agree to 1e-6) are not
  //counted:
  std::size_t ndisc = 0;
  for ( std::size_t i = 0; i + 1 < x.size(); ++i )
    if ( x[i+1] == x[i]
         && ncabs( y[i+1] - y[i] ) > kink_tol * ncmax( y[i], y[i+1] ) )
      ++ndisc;

  //The limit for wavelength -> 0, from linear extrapolation to 0 from very
  //short wavelengths (which is exact for any cross section of the form a+b*wl
  //there), at three scales. When the limit exists, the differences between
  //the estimates at consecutive scales decrease rapidly (e.g. by a factor of
  //100 for a term c*wl^2), while they do not for cross sections without a
  //limit (e.g. with 1/wl, log(wl) or sqrt(wl) behaviour):
  {
    auto extrapTo0 = [&evalXS]( double h )
    {
      return ncfma( 2.0, evalXS( h ), -evalXS( 2.0 * h ) );
    };
    const double lim_a = extrapTo0( 1e-7 );
    const double lim_b = extrapTo0( 1e-8 );
    const double lim_c = extrapTo0( 1e-9 );
    const double d_ab = ncabs( lim_a - lim_b );
    const double d_bc = ncabs( lim_b - lim_c );
    //Either the estimates agree, or they converge (geometrically). The limit
    //is only used below denseWavelengthMin(), where the tolerance is not
    //required, so a small negative estimate (e.g. -2*c*h^2 for a cross
    //section c*wl^2) simply means a limit of 0:
    const double scale = ncmax( ncabs( lim_c ),
                                evalXS( denseWavelengthMin() ) );
    const bool converges = ( d_bc <= 1e-3 * tol * scale
                             || d_bc <= 0.05 * d_ab );
    if ( !converges )
      NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: the cross section does not"
                       " converge for wavelength -> 0 (estimated limits "
                       << lim_a << ", " << lim_b << " and " << lim_c
                       << " /cm)" );
    y.front() = ( lim_c > 0.0 ? lim_c : 0.0 );
    //Wavelength 0 is a limit, without a slope of its own. The chord towards
    //the next sample point is used:
    slopes.front() = ( y[1] - y[0] ) / ( x[1] - x[0] );
  }

  //Reduce the dense sample, and verify the result at the midpoints between
  //all neighbouring sample points. Midpoints where the tolerance is not met
  //are added to the sample, and the procedure is repeated.
  ReducedHermitePoints red;
  std::size_t nrefine = 0;
  while ( true ) {
    red = reduceHermitePoints( x, y, slopes, tol_build );
    const auto& idx = red.idx;
    VectD add_x, add_y, add_slopes;
    add_x.reserve( x.size() );//at most one point per interval
    add_y.reserve( x.size() );
    add_slopes.reserve( x.size() );
    std::size_t seg = 0;//table segment is from x[idx[seg]] to x[idx[seg+1]]
    for ( std::size_t k = 0; k + 1 < x.size(); ++k ) {
      while ( idx[seg+1] <= k )
        ++seg;
      if ( !( x[k+1] > x[k] ) )
        continue;//known discontinuity
      const std::size_t i = idx[seg];
      const std::size_t j = idx[seg+1];
      const double xm = 0.5 * ( x[k] + x[k+1] );
      const double ym = evalXS( xm );
      const double yc = hermite( x[i], x[j], y[i], y[j],
                                 red.slope[seg], red.slope[seg+1], xm );
      if ( xm < denseWavelengthMin()
           || ncabs( yc - ym ) <= tol_build * ym )
        continue;
      if ( !( x[k+1] - x[k] > min_rel_spacing * x[k+1] ) )
        NCRYSTAL_THROW2( CalcError, "IsoMacroXSTable: the cross section has an"
                         " unexpected discontinuity (or an extremely sharp"
                         " feature) at wavelength " << xm << " Aa" );
      add_x.push_back( xm );
      add_y.push_back( ym );
      add_slopes.push_back( evalSlope( xm ) );
    }
    if ( add_x.empty() )
      break;
    if ( ++nrefine > max_refine_iterations )
      NCRYSTAL_THROW( CalcError, "IsoMacroXSTable: the refinement of the table"
                      " did not converge" );
    //Merge the new points into the (sorted) sample:
    VectD new_x, new_y, new_slopes;
    new_x.reserve( x.size() + add_x.size() );
    new_y.reserve( x.size() + add_x.size() );
    new_slopes.reserve( x.size() + add_x.size() );
    std::size_t ia = 0;
    for ( std::size_t k = 0; k < x.size(); ++k ) {
      for ( ; ia < add_x.size() && add_x[ia] < x[k]; ++ia ) {
        new_x.push_back( add_x[ia] );
        new_y.push_back( add_y[ia] );
        new_slopes.push_back( add_slopes[ia] );
      }
      new_x.push_back( x[k] );
      new_y.push_back( y[k] );
      new_slopes.push_back( slopes[k] );
    }
    nc_assert_always( ia == add_x.size() );
    x.swap( new_x );
    y.swap( new_y );
    slopes.swap( new_slopes );
  }

  VectD res_wl, res_xs;
  res_wl.reserve( red.idx.size() );
  res_xs.reserve( red.idx.size() );
  for ( auto i : red.idx ) {
    res_wl.push_back( x[i] );
    res_xs.push_back( y[i] );
  }
  auto res = IsoMacroXSTable::fromPoints( res_wl, res_xs, red.slope );

  if ( diagnostics ) {
    diagnostics->nDiscontinuities = ndisc;
    diagnostics->nDense = x.size();
    diagnostics->nRefine = nrefine;
    diagnostics->processes.clear();
  }
  return res;
}

NCX::IsoMacroXSTable
NCX::createIsoMacroXSTable( const MatCfg& cfg,
                            const IsoMacroXSTableParams& params,
                            IsoMacroXSTableDiagnostics* diagnostics )
{
  auto info = FactImpl::createInfo( cfg );
  auto scat = FactImpl::createScatter( cfg );
  auto absn = FactImpl::createAbsorption( cfg );
  if ( scat->isOriented() || absn->isOriented() )
    NCRYSTAL_THROW( BadInput, "IsoMacroXSTable: only isotropic materials are"
                    " supported" );

  std::set<std::string> names, unvalidated;
  inspectProcess( *scat, names, unvalidated );
  inspectProcess( *absn, names, unvalidated );
  for ( auto& n : unvalidated )
    NCRYSTAL_WARN( "IsoMacroXSTable: the physics process \"" << n << "\" has"
                   " not been validated for IsoMacroXSTable (features of its"
                   " cross section might be missed)" );
  VectD discontinuities;
  collectBreakpoints( *scat, discontinuities );
  collectBreakpoints( *absn, discontinuities );

  //barn/atom * atoms/Aa^3 = 1/cm:
  const double numberDensity = info->getNumberDensity().dbl();
  CachePtr cp_scat, cp_abs;
  auto xs_of_wl = [&scat,&absn,&cp_scat,&cp_abs,numberDensity]( double wl )
  {
    const NeutronEnergy ekin{ wl2ekin( wl ) };
    return numberDensity
      * ( scat->crossSectionIsotropic( cp_scat, ekin ).dbl()
          + absn->crossSectionIsotropic( cp_abs, ekin ).dbl() );
  };

  //The slope of the cross section: from the merged tables of the processes
  //which have a table, and from numerical differences
  //for the other processes (whose cross sections are smooth between the known
  //discontinuities):
  std::vector<std::pair<double,IsoXSTable>> tabulated;
  std::vector<UntabulatedLeaf> untabulated;
  collectLeaves( *scat, 1.0, tabulated, untabulated );
  collectLeaves( *absn, 1.0, tabulated, untabulated );
  Optional<IsoXSTable> merged;
  {
    std::vector<std::pair<double,const IsoXSTable*>> inputs;
    for ( auto& t : tabulated )
      if ( t.first != 0.0 )
        inputs.emplace_back( t.first, &t.second );
    if ( !inputs.empty() )
      merged = mergeIsoXSTables( inputs );
  }
  VectD sorted_disc = discontinuities;
  std::sort( sorted_disc.begin(), sorted_disc.end() );
  auto dxs_of_wl = [&merged,&untabulated,&sorted_disc,
                    numberDensity]( double wl )
  {
    const double e = wl2ekin( wl );
    //de/dwl = -2e/wl:
    double res = 0.0;
    if ( merged.has_value() )
      res = ( evalIsoXSDerivative( merged.value(), e )
              * ( -2.0 * e / wl ) );
    for ( auto& u : untabulated ) {
      auto f = [&u]( double w )
      {
        const NeutronEnergy ekin{ wl2ekin( w ) };
        return u.process->crossSectionIsotropic( u.cache, ekin ).dbl();
      };
      res += u.scale * numericalSlope( f, wl, sorted_disc );
    }
    return numberDensity * res;
  };

  if ( untabulated.empty() && merged.has_value() ) {
    validateParams( params );
    auto pts = tabulateExact( merged.value(), numberDensity, params.wlmax,
                              tol_exact_factor * params.tol );
    if ( diagnostics ) {
      diagnostics->nDiscontinuities = pts.ndiscontinuities;
      diagnostics->nDense = pts.ncandidates;
      diagnostics->nRefine = 0;
      diagnostics->processes.assign( names.begin(), names.end() );
    }
    return IsoMacroXSTable::fromPoints( pts.wl, pts.xs, pts.slope );
  }

  //Processes with tables are smooth between the energies of their tables,
  //so these are sampled too:
  VectD samplePoints;
  if ( merged.has_value() )
    for ( auto e : merged.value().E )
      samplePoints.push_back( ekin2wl( e ) );
  auto res = createIsoMacroXSTable( xs_of_wl, dxs_of_wl,
                                    std::move( discontinuities ), samplePoints,
                                    params, diagnostics );
  if ( diagnostics )
    diagnostics->processes.assign( names.begin(), names.end() );
  return res;
}

NCX::IsoMacroXSTable::IsoMacroXSTable( std::unique_ptr<double[]> data )
  : m_data( std::move( data ) )
{
}

NCX::IsoMacroXSTable NCX::IsoMacroXSTable::fromPoints( const VectD& wl,
                                                       const VectD& xs,
                                                       const VectD& slope )
{
  //The block in the format documented in examples/ncrystal_example_filter.c:
  const std::size_t n = wl.size();
  nc_assert_always( n >= 2 && xs.size() == n && slope.size() == n
                    && wl.front() == 0.0 && wl[n-2] < wl[n-1] );
  auto data = ncmake_unique_array<double>( 3 * n + 4 );
  data[0] = magic();
  data[1] = formatVersion();
  data[2] = static_cast<double>( n );
  std::copy( wl.begin(), wl.end(), data.get() + 3 );
  std::copy( xs.begin(), xs.end(), data.get() + 3 + n );
  std::copy( slope.begin(), slope.end(), data.get() + 3 + 2 * n );
  data[3 + 3 * n] = magic();
  return IsoMacroXSTable( std::move( data ) );
}
