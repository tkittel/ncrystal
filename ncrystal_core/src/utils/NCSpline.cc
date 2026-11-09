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

#include "NCrystal/internal/utils/NCSpline.hh"
#include "NCrystal/internal/utils/NCFastSearch.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCSpline_FMA.hh"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace NC = NCrystal;

//TODO: For convenience, CubicSpline derivatives should default to
//special values (e.g. inf), in which case derivatives will be approximated,
//either linearly from the edge points or using a numerical derivative
//estimation.

void NC::CubicSpline::set( const VectD& y,
                           double derivative_y_left,
                           double derivative_y_right )
{
  const std::size_t n = y.size();
  nc_assert_always(n>3);
  const std::size_t nm1 = n-1;
  VectD u(nm1,0.0);
  VectD y2(n,0.0);
  y2[0] = -0.5;
  double dy_left = y[1]-y[0];
  u[0] = 3.0*( dy_left-derivative_y_left );
  for (std::size_t i = 1; i<nm1;++i) {
    double p = 2.0 + 0.5 * y2[i-1];
    nc_assert(p!=0.0);
    p = 1.0/p;
    y2[i] = -0.5*p;
    double uu = y[i+1] - 2.0*y[i] + y[i-1];
    u[i] = (3.0*uu -0.5*u[i-1])*p;
  }
  double qn = 0.5;
  double un = 3.0 * (derivative_y_right -(y[nm1]-y[n-2]));
  nc_assert((qn*y2[n-2]+1.0) != 0.0);
  y2[nm1] = (un-qn*u[n-2]) / (qn*y2[n-2]+1.0);
  for (std::size_t k = nm1; k; --k) {
    std::size_t km1 = k-1;
    y2[km1] *= y2[k];
    y2[km1] += u[km1];
  }
  std::vector<PairDD > data;
  data.reserve(y.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    data.push_back(std::make_pair(y[i],y2[i]));
    nc_assert(!ncisnan(data.back().first));
    nc_assert(!ncisnan(data.back().second));
  }
  //all good, set:
  std::swap(m_data,data);
  m_nm2 = n-2;
}

void NC::SplinedLookupTable::set( const VectD& fvals,
                                  double a,double b,
                                  double fprime_a, double fprime_b,
                                  const std::string& name,
                                  const std::string& description  )
{
  nc_assert(b>a);
  nc_assert(fvals.size()>3);
  m_a = a;
  m_b = b;
  double delta = (b-a)/(fvals.size()-1);
  m_spline.set( fvals, fprime_a*delta, fprime_b*delta );
  nc_assert(delta>0.0);
  m_invdelta = 1.0 / delta;
  if (ncgetenv_bool("DEBUG_SPLINES"))
    producefile( 0/*true function not known*/,fprime_a, fprime_b,name,description );

}

void NC::SplinedLookupTable::set( const Fct1D* thefct,
                                  double a,double b,
                                  double fprime_a, double fprime_b,
                                  unsigned npts,
                                  const std::string& name,
                                  const std::string& description )
{
  nc_assert(!ncisnan(fprime_a));
  nc_assert(!ncisnan(fprime_b));
  nc_assert(thefct);
  nc_assert(b>a);
  nc_assert(npts>3);
  m_a = a;
  m_b = b;
  VectD fvals;
  fvals.reserve(npts);
  double delta = (b-a)/(npts-1);
  std::size_t nm1 = npts-1;
  for (unsigned i = 0; i < nm1; ++i) {
    fvals.push_back(thefct->eval(a+i*delta));
    nc_assert(!ncisnan(fvals.back()));
  }
  fvals.push_back(thefct->eval(b));
  nc_assert(!ncisnan(fvals.back()));
  m_spline.set( fvals, fprime_a*delta, fprime_b*delta );
  nc_assert(delta>0.0);
  m_invdelta = 1.0 / delta;


  if (ncgetenv_bool("DEBUG_SPLINES"))
    producefile( thefct,fprime_a, fprime_b,name,description );
}

void NC::SplinedLookupTable::producefile( const Fct1D* thefct,
                                          double fprime_a, double fprime_b,
                                          const std::string& username,
                                          const std::string& userdesc ) const
{

  std::string name = (username.empty()?std::string("unknownspline"):username);
  std::string description = (userdesc.empty()?std::string("no description"):userdesc);
  std::string filename;
  for (unsigned ifile = 0;ifile<=99;++ifile) {
    std::stringstream s;
    s << "ncrystal_splinedlt_"<<name;
    if (ifile>0)
      s<<"_"<<ifile;
    s<<".txt";
    filename = s.str();
    if (!file_exists(filename))
      break;
  }

  if (file_exists(filename)) {
    NCRYSTAL_WARN("Aborted writing of spline data file ncrystal_splinedlt_"
                  <<username<<"_??.txt - too many files in current dir.");
    return;
  }

  std::ofstream ofs (filename.c_str(), std::ofstream::out);
  ofs << std::setprecision(20);
  ofs << "#ncrystal_splined_lookuptable\n";
  ofs << "#name = "<<name<<"\n";
  ofs << "#description = "<<description<<"\n";
  ofs << "#a = "<<getLower()<<"\n";
  ofs << "#b = "<<getUpper()<<"\n";
  ofs << "#fprime_a = "<<fprime_a<<"\n";
  ofs << "#fprime_b = "<<fprime_b<<"\n";
  ofs << "#input_fvals = ";
  for (std::size_t i = 0; i < m_spline.m_data.size(); ++i)
    ofs<<" "<<m_spline.m_data[i].first;
  ofs << "\n#data_colums = x,spline_of_x";
  if (thefct)
    ofs <<",truefct_of_x";
  ofs<<"\n";
  std::size_t numpts = 100*m_spline.m_data.size();
  if (numpts>1000000)
    numpts = std::max<std::size_t>(numpts/10,1000000);
  double delta = (m_b-m_a)/(numpts-1.0);
  for (std::size_t i = 0; i<=numpts;++i) {
    double x = (i==numpts?m_b:m_a+i*delta);
    ofs << x << " "<<eval(x);
    if (thefct)
      ofs<<" "<<thefct->eval(x);
    ofs<<"\n";
  }
  NCRYSTAL_MSG("Wrote "<<filename<<" (since NCRYSTAL_DEBUG_SPLINE is set).");
}

void NC::PiecewiseLinearFct1D::dumpToFile( const std::string& filename ) const {
  std::ofstream ofs (filename.c_str(), std::ofstream::out);
  ofs << std::setprecision(20);
  ofs << "#colnames=x,y\n";
  ofs << "#plotstyle=*-\n";
  ofs << "#overflow=";
  if  ( m_ofVals.overflowYValue.has_value() )
    ofs << "none\n";
  else
    ofs << m_ofVals.overflowYValue.value()<< "\n";
  ofs << "#underflow=";
  if  ( m_ofVals.underflowYValue.has_value() )
    ofs << "none\n";
  else
    ofs << m_ofVals.underflowYValue.value()<< "\n";

  for ( auto i : ncrange(m_x.size()) ) {
    ofs << m_x.at(i) << " "<<m_y.at(i) << "\n";
  }
  ofs << std::flush;
}

double NC::PiecewiseLinearFct1D::evalEdgeCase( VectD::const_iterator it,
                                               double x ) const
{
  if ( it == m_x.end() ) {
    if ( !m_ofVals.overflowYValue.has_value() )
      NCRYSTAL_THROW2(CalcError,"PiecewiseLinearFct1D: Out of bounds: x>xmax"
                      " and no overflow value supplied (x="<<x<<", xmax="<<m_x.back()<<").");
    return m_ofVals.overflowYValue.value();
  };
  if ( x >= m_x.front() )
    return m_y.front();
  if ( !m_ofVals.underflowYValue.has_value() )
    NCRYSTAL_THROW2(CalcError,"PiecewiseLinearFct1D: Out of bounds: x<xmin"
                    " and no underflow value supplied (x="<<x<<", xmin="<<m_x.front()<<").");
  return m_ofVals.underflowYValue.value();
}

NC::PCHIPInterp::PCHIPInterp( Span<const double> x, Span<const double> y,
                              Optional<double> opt_ls, Optional<double> opt_rs,
                              OptimiseForLogSpacedX optlogspacedx )
{
#ifndef NDEBUG
  nc_assert(x.size() == y.size());
  nc_assert(nc_is_grid(x));
  nc_assert(std::isfinite(opt_ls.value_or(1.0)));
  nc_assert(std::isfinite(opt_rs.value_or(1.0)));
  for ( auto e : y )
    nc_assert( std::isfinite(e) );
#endif
  m_data.resize(5u * x.size() - 2u);
  m_n = x.size();

  if ( x.front() > 0.0 && optlogspacedx == OptimiseForLogSpacedX::YES )
    m_guide.emplace( x.data(), x.size() );

  const std::size_t oh = 2 * m_n;
  const std::size_t oi = 3 * m_n - 1;
  const std::size_t om = 4 * m_n - 2;

  for (std::size_t i = 0; i < m_n; ++i) {
    vectAt(m_data, i) = x[i];
    vectAt(m_data, m_n + i) = y[i];
  }

  for (std::size_t i = 0; i + 1 < m_n; ++i) {
    const double h = x[i + 1] - x[i];
    nc_assert( std::isfinite(h) );
    vectAt(m_data, oh + i) = h;
    vectAt(m_data, oi + i) = 1.0 / h;
  }

  auto sec = [this,oi](std::size_t i) {
    return ( ( vectAt(m_data, m_n + i + 1) - vectAt(m_data, m_n + i) )
             * vectAt(m_data, oi + i) );
  };

  auto limit = [](double a, double b) {
    if (b == 0.0 || a * b <= 0.0)
      return 0.0;
    if (ncabs(a) > 3.0 * ncabs(b))
      return 3.0 * b;
    return a;
  };

auto edgeSlope = [this,&sec](std::size_t i0,
                               std::size_t i1) {
    if (m_n == 2)
      return sec(0);

    const double h0 = vectAt(m_data, 2 * m_n + i0);
    const double h1 = vectAt(m_data, 2 * m_n + i1);
    const double d0 = sec(i0);
    const double d1 = sec(i1);

    double a = ((2.0 * h0 + h1) * d0 - h0 * d1)
             / (h0 + h1);

    if (a * d0 <= 0.0)
      a = 0.0;
    else if (d0 * d1 < 0.0 && ncabs(a) > 3.0 * ncabs(d0))
      a = 3.0 * d0;

    return a;
  };

  const double ls = ( opt_ls.has_value() ?
                      opt_ls.value() :
                      edgeSlope(0, 1) );
  const double rs = ( opt_rs.has_value() ?
                      opt_rs.value() :
                      edgeSlope(m_n - 2, m_n - 3) );

  vectAt(m_data, om) = limit(ls, sec(0));

  for (std::size_t i = 1; i + 1 < m_n; ++i) {
    const double a = sec(i - 1);
    const double b = sec(i);
    if (a * b <= 0.0) {
      vectAt(m_data, om + i) = 0.0;
    } else {
      const double w0
        = 2.0 * vectAt(m_data, oh + i) + vectAt(m_data, oh + i - 1);
      const double w1
        = vectAt(m_data, oh + i) + 2.0 * vectAt(m_data, oh + i - 1);
      vectAt(m_data, om + i) = (w0 + w1) / (w0 / a + w1 / b);
    }
  }

  vectAt(m_data, om + m_n - 1) = limit(rs, sec(m_n - 2));
}

double NC::PCHIPInterp::eval(double x) const
{
#ifndef NDEBUG
  nc_assert_always( m_n >= 2 );
  nc_assert(std::isfinite(x));
#endif
  const double* p = m_data.data();
  const double* xp = p;
  const double* yp = p + m_n;
  const double* hp = p + 2 * m_n;
  const double* ip = p + 3 * m_n - 1;
  const double* mp = p + 4 * m_n - 2;
#ifndef NDEBUG
  nc_assert(x >= xp[0]);
  nc_assert(x <= xp[m_n - 1]);
#endif
  if (x == xp[m_n - 1])
    return yp[m_n - 1];
  const std::size_t i = ( m_guide.has_value()
                          ? m_guide.value().upperBoundIdx(xp, m_n, x)
                          : NC::fastUpperBoundIdx(xp, m_n, x) ) - 1;
  const double t = (x - xp[i]) * ip[i];
  const double t2 = t * t;
  const double u = 1.0 - t;
  const double u2 = u * u;
  const double h00 = (1.0 + 2.0 * t) * u2;
  const double h10 = t * u2;
  const double h01 = t2 * (3.0 - 2.0 * t);
  const double h11 = t2 * (t - 1.0);
  return ( h00*yp[i] + h10*hp[i]*mp[i] + h01*yp[i + 1] + h11*hp[i]*mp[i + 1] );
}

NC::PairDD NC::PCHIPInterp::endPointSlopes() const
{
  nc_assert_always(m_n >= 2);
  return { vectAt(m_data, 4 * m_n - 2), vectAt(m_data, 5 * m_n - 3 ) };
}

NC::HymanInterp::HymanInterp( Span<const double> x, Span<const double> y )
  : m_x( x.begin(), x.end() ),
    m_y( y.begin(), y.end() )
{
#ifndef NDEBUG
  nc_assert_always( x.size() == y.size() && x.size() >= 4 );
  nc_assert_always( nc_is_grid( x ) );
  for ( auto e : y )
    nc_assert_always( std::isfinite( e ) );
#endif
  //The slope at a node is the derivative there of the cubic through 4
  //surrounding nodes (the average of the two such cubics with the node in
  //their middle, where available). The Hyman filter then keeps its sign equal
  //to that of the data where the data is monotone, limits it to 3 times the
  //smaller adjacent secant slope (the Fritsch-Carlson bound), and zeroes it at
  //local extrema of the data, which keeps each interval within its node values.
  const std::size_t n = m_x.size();
  //Derivative at x[i] of the cubic through nodes j..j+3 (Lagrange form):
  auto derivCubic = [this]( std::size_t j, std::size_t i )
  {
    double res = 0.0;
    for ( auto a : ncrange( j, j + 4 ) ) {
      double da = 0.0;//derivative of the basis polynomial of node a at x[i]
      for ( auto b : ncrange( j, j + 4 ) ) {
        if ( b == a )
          continue;
        double prod = 1.0 / ( m_x[a] - m_x[b] );
        for ( auto c : ncrange( j, j + 4 ) )
          if ( c != a && c != b )
            prod *= ( m_x[i] - m_x[c] ) / ( m_x[a] - m_x[c] );
        da += prod;
      }
      res += m_y[a] * da;
    }
    return res;
  };
  m_m.resize( n );
  for ( auto i : ncrange( n ) ) {
    //Stencils i-1..i+2 and i-2..i+1, or the 4 nodes nearest the ends:
    const bool s1 = ( i >= 1 && i + 2 < n );
    const bool s2 = ( i >= 2 && i + 1 < n );
    double m;
    if ( s1 && s2 )
      m = 0.5 * ( derivCubic( i - 1, i ) + derivCubic( i - 2, i ) );
    else if ( s1 )
      m = derivCubic( i - 1, i );
    else if ( s2 )
      m = derivCubic( i - 2, i );
    else
      m = derivCubic( i == 0 ? 0 : n - 4, i );
    //Hyman filter (secant slopes on both sides, the one-sided one at ends):
    const double dl = ( i > 0 ? ( m_y[i] - m_y[i-1] ) / ( m_x[i] - m_x[i-1] )
                        : ( m_y[1] - m_y[0] ) / ( m_x[1] - m_x[0] ) );
    const double dr = ( i + 1 < n ? ( m_y[i+1] - m_y[i] ) / ( m_x[i+1] - m_x[i] )
                        : dl );
    m_m[i] = hymanFilter( m, dl, dr );
  }
}

NC::VectD NC::HymanInterp::uniformSlopes( Span<const double> y, double h )
{
#ifndef NDEBUG
  nc_assert_always( y.size() >= 4 );
  nc_assert_always( h > 0.0 && std::isfinite( h ) );
  for ( auto e : y )
    nc_assert_always( std::isfinite( e ) );
#endif
  VectD res( y.size() );
  NCRYSTAL_APPLY_C_NAMESPACE(detail_hymanSlopesUniform)( y.data(), y.size(),
                                                        h, res.data() );
  return res;
}

NC::VectD NC::HymanInterp::uniformCubics( Span<const double> y, double h )
{
#ifndef NDEBUG
  nc_assert_always( y.size() >= 4 );
  nc_assert_always( h > 0.0 && std::isfinite( h ) );
  for ( auto e : y )
    nc_assert_always( std::isfinite( e ) );
#endif
  VectD res( 4 * ( y.size() - 1 ) );
  NCRYSTAL_APPLY_C_NAMESPACE(detail_hymanCubicsUniform)( y.data(), y.size(),
                                                        h, res.data() );
  return res;
}

double NC::HymanInterp::eval( double x ) const
{
  const std::size_t n = m_x.size();
  nc_assert( x >= m_x.front() && x <= m_x.back() );
  if ( x == m_x.back() )
    return m_y.back();
  const std::size_t i = fastUpperBoundIdx( m_x.data(), n, x ) - 1;
  const double h = m_x[i+1] - m_x[i];
  const double t = ( x - m_x[i] ) / h;
  const double t2 = t * t;
  const double u = 1.0 - t;
  const double u2 = u * u;
  return ( ( 1.0 + 2.0 * t ) * u2 * m_y[i] + t * u2 * h * m_m[i]
           + t2 * ( 3.0 - 2.0 * t ) * m_y[i+1] + t2 * ( t - 1.0 ) * h * m_m[i+1] );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    void fitPCInterval( const std::function<double(double)>& f,
                        double a, double fa, double b, double fb,
                        double tol, unsigned depth, PiecewiseCubicFit& res )
    {
      //Fit the cubic in u=(x-a)/(b-a) through f at u=0,1/3,2/3,1 (the Lagrange
      //interpolation of those 4 values, written below in monomial form with
      //the lowest order first), and check it against f at u=1/6,1/2,5/6, the
      //midpoints between those points, close to where its error is largest. If
      //any check fails, bisect geometrically (at sqrt(a*b)), at most 60 levels
      //deep:
      nc_assert( a < b && res.x.back() == a );
      const double h = b - a;
      const double y1 = f( a + h / 3.0 ), y2 = f( a + 2.0 * h / 3.0 );
      const double c[4] = { fa,
                            0.5 * ( -11.0 * fa + 18.0 * y1 - 9.0 * y2 + 2.0 * fb ),
                            0.5 * ( 18.0 * fa - 45.0 * y1 + 36.0 * y2 - 9.0 * fb ),
                            0.5 * ( -9.0 * fa + 27.0 * y1 - 27.0 * y2 + 9.0 * fb ) };
      bool ok = true;
      for ( double u : { 1.0 / 6.0, 0.5, 5.0 / 6.0 } ) {
        const double v = evalCubicU( c, u );
        const double ref = f( a + u * h );
        if ( !( ncabs( v - ref ) <= 0.5 * tol * ncabs( ref ) ) )
          ok = false;
      }
      if ( !ok && depth < 60 ) {
        const double mid = std::sqrt( a * b );
        if ( mid > a && mid < b ) {
          const double fmid = f( mid );
          fitPCInterval( f, a, fa, mid, fmid, tol, depth + 1, res );
          fitPCInterval( f, mid, fmid, b, fb, tol, depth + 1, res );
          return;
        }
      }
      res.x.push_back( b );
      res.c.insert( res.c.end(), c, c + 4 );
    }
  }
}

NC::PiecewiseCubicFit NC::fitPiecewiseCubic( const std::function<double(double)>& f,
                                             double a, double b, double tol,
                                             double maxratio )
{
  nc_assert( a > 0.0 && b > a && std::isfinite( b ) );
  nc_assert( tol > 0.0 && maxratio > 1.0 );
  //Initial intervals span at most a factor maxratio. Each gets the cubic
  //through f at u=0,1/3,2/3,1, and is bisected (geometrically) until it agrees
  //with f within tol/2 at u=1/6,1/2,5/6 (5 evaluations per interval). The
  //check is not a strict bound, hence the factor 2.
  const std::size_t n = 1 + static_cast<std::size_t>
    ( std::ceil( std::log( b / a ) / std::log( maxratio ) ) );
  VectD pts = geomspace( a, b, std::max<std::size_t>( n, 2 ) );
  pts.front() = a;
  pts.back() = b;
  PiecewiseCubicFit res;
  res.x.push_back( a );
  double fa = f( a );
  for ( auto i : ncrange( pts.size() - 1 ) ) {
    const double fb = f( pts[i+1] );
    fitPCInterval( f, pts[i], fa, pts[i+1], fb, tol, 0, res );
    fa = fb;
  }
  return res;
}
