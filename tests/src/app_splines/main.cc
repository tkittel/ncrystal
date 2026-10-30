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
#include <iostream>

namespace NC=NCrystal;

namespace {

  using NC::VectD;
  using NC::floateq;
  using NC::PCHIPInterp;
  using NC::ncmin;
  using NC::ncmax;
  using NC::NullOpt;

  void testPCHIPInterp1()
  {
    VectD x{1.0, 2.0, 4.0, 8.0};
    VectD y{3.0, 5.0, 9.0, 17.0};
    PCHIPInterp p(x,y,2.0, 2.0);
    nc_assert_always(floateq(p.eval(1.0),3.0));
    nc_assert_always(floateq(p.eval(8.0),17.0));
    nc_assert_always(floateq(p.eval(3.0),7.0));
    nc_assert_always(floateq(p.eval(6.0),13.0));
  }

  void testPCHIPInterp2()
  {
    auto checkClose = [](double a, double b, double rtol = 1.0e-12,
                         double atol = 1.0e-14) {
      nc_assert_always(floateq(a, b, rtol, atol));
    };
    auto checkKnotValues = [&](const PCHIPInterp& interp,
                               const VectD& x,
                               const VectD& y) {
      nc_assert_always(x.size() == y.size());
      for (std::size_t i = 0; i < x.size(); ++i)
        checkClose(interp.eval(x.at(i)),y.at(i));
    };

    auto checkIntervalRange = [&](const PCHIPInterp& interp,
                                  const VectD& x,
                                  const VectD& y) {
      nc_assert_always(x.size() == y.size());
      for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        const double x0 = x.at(i);
        const double x1 = x.at(i + 1);
        const double y0 = y.at(i);
        const double y1 = y.at(i + 1);
        const double ymin = ncmin(y0, y1);
        const double ymax = ncmax(y0, y1);
        for (std::size_t j = 1; j < 10; ++j) {
          const double t = static_cast<double>(j) / 10.0;
          const double xx = x0 + t * (x1 - x0);
          const double yy = interp.eval(xx);
          nc_assert_always(yy >= ymin || floateq(yy, ymin, 1.0e-12,1.0e-14));
          nc_assert_always(yy <= ymax || floateq(yy, ymax, 1.0e-12,1.0e-14));
        }
      }
    };

    // An affine function should be reproduced exactly, including its
    // derivative, even on a nonuniform grid.
    {
      const VectD x{0.0, 0.25, 1.0, 2.5, 7.0};
      const VectD y{7.0, 7.5, 9.0, 12.0, 21.0};
      PCHIPInterp interp(x, y);
      checkKnotValues(interp, x, y);
      const auto slopes = interp.endPointSlopes();
      checkClose(slopes.first, 2.0);
      checkClose(slopes.second, 2.0);
      checkClose(interp.eval(0.125), 7.25);
      checkClose(interp.eval(0.625), 8.25);
      checkClose(interp.eval(3.75), 14.5);
    }

    // Constant data should produce zero slopes and a constant result.
    {
      const VectD x{0.0, 0.1, 1.0, 3.0, 10.0};
      const VectD y{4.5, 4.5, 4.5, 4.5, 4.5};
      PCHIPInterp interp(x, y);
      const auto slopes = interp.endPointSlopes();
      checkClose(slopes.first, 0.0);
      checkClose(slopes.second, 0.0);
      checkKnotValues(interp, x, y);
      checkClose(interp.eval(0.37), 4.5);
      checkClose(interp.eval(6.0), 4.5);
    }

    // Check the standard one-sided PCHIP endpoint formula using y = x^2
    // on a unit-spaced grid. The endpoint secants are 1 and 3 on the
    // left, and 3 and 5 on the right.
    //
    // left slope  = ((2 + 1) * 1 - 1 * 3) / 2 = 0.0
    // right slope = ((2 + 1) * 5 - 1 * 3) / 2 = 6.0
    {
      const VectD x{0.0, 1.0, 2.0, 3.0};
      const VectD y{0.0, 1.0, 4.0, 9.0};
      PCHIPInterp interp(x, y);
      const auto slopes = interp.endPointSlopes();
      checkClose(slopes.first, 0.0);
      checkClose(slopes.second, 6.0);
      checkKnotValues(interp, x, y);
    }

    // A monotone data set must not overshoot between neighboring knots.
    {
      const VectD x{0.0, 1.0, 2.0, 4.0, 7.0};
      const VectD y{1.0, 2.0, 4.0, 5.0, 8.0};
      PCHIPInterp interp(x, y);
      checkKnotValues(interp, x, y);
      checkIntervalRange(interp, x, y);
    }

    // A decreasing data set must also remain within the range of each
    // pair of neighboring data values.
    {
      const VectD x{0.0, 0.5, 1.5, 3.0, 6.0};
      const VectD y{8.0, 7.0, 5.0, 4.0, 1.0};
      PCHIPInterp interp(x, y, NullOpt, NullOpt);
      checkKnotValues(interp, x, y);
      checkIntervalRange(interp, x, y);
    }

    // At a local extremum, the interior PCHIP slope must be zero.
    // The interval-range test also checks that the interpolant does
    // not overshoot around the extrema.
    {
      const VectD x{0.0, 1.0, 2.0, 3.0};
      const VectD y{0.0, 1.0, 0.0, 1.0};
      PCHIPInterp interp(x, y, NullOpt, NullOpt);
      checkKnotValues(interp, x, y);
      checkIntervalRange(interp, x, y);
    }

    // With only two points, automatically calculated endpoint slopes
    // must both equal the single secant.
    {
      const VectD x{2.0, 5.0};
      const VectD y{11.0, 20.0};
      PCHIPInterp interp(x, y);
      const auto slopes = interp.endPointSlopes();
      checkClose(slopes.first, 3.0);
      checkClose(slopes.second, 3.0);
      checkKnotValues(interp, x, y);
      checkClose(interp.eval(3.5), 15.5);
    }

    // Supplied slopes are currently treated as preferred slopes and
    // passed through the endpoint limiter. With one interval, the
    // secant is 1.0:
    //
    //   left supplied slope  = 4.0 -> limited to 3.0
    //   right supplied slope = -2.0 -> limited to 0.0
    {
      const VectD x{0.0, 1.0};
      const VectD y{0.0, 1.0};
      const double leftSlope(4.0);
      const double rightSlope(-2.0);
      PCHIPInterp interp(x, y,leftSlope, rightSlope);
      const auto slopes = interp.endPointSlopes();
      checkClose(slopes.first, 3.0);
      checkClose(slopes.second, 0.0);
      checkClose(interp.eval(0.0), 0.0);
      checkClose(interp.eval(1.0), 1.0);
      // With m0 = 3 and m1 = 0, p(0.5) = 0.875.
      checkClose(interp.eval(0.5), 0.875);
    }

    // A supplied slope that already has the correct sign and magnitude
    // should be retained.
    {
      const VectD x{0.0, 1.0};
      const VectD y{0.0, 1.0};
      const double leftSlope(1.5);
      const double rightSlope(0.5);
      PCHIPInterp interp(x, y,leftSlope, rightSlope);
      const auto slopes = interp.endPointSlopes();
      checkClose(slopes.first, 1.5);
      checkClose(slopes.second, 0.5);
    }

    // A supplied slope with the wrong sign is changed to zero by the
    // current limiter.
    {
      const VectD x{0.0, 1.0, 2.0};
      const VectD y{0.0, 1.0, 2.0};
      const double leftSlope(-1.0);
      const double rightSlope(-1.0);
      PCHIPInterp interp(x, y,leftSlope, rightSlope);
      const auto slopes = interp.endPointSlopes();
      checkClose(slopes.first, 0.0);
      checkClose(slopes.second, 0.0);
    }
  }

  void testHymanInterp()
  {
    using NC::HymanInterp;
    //Cubic data with the filter inactive (monotone, slopes within bounds) are
    //reproduced exactly, on an irregular grid:
    {
      auto f = []( double x ) { return 1.0 + x + 0.1 * x * x + 0.01 * x * x * x; };
      NC::VectD x{ 0.0, 0.3, 0.35, 1.0, 1.8, 2.5, 4.0 }, y;
      for ( auto v : x )
        y.push_back( f( v ) );
      HymanInterp interp( x, y );
      for ( double v : NC::linspace( 0.0, 4.0, 401 ) )
        nc_assert_always( NC::floateq( interp.eval( v ), f( v ), 1e-13, 1e-13 ) );
    }
    //O(h^4) convergence for smooth data (exp on a geometric grid):
    {
      double prev = 0.0;
      for ( std::size_t n : { 20, 40, 80 } ) {
        NC::VectD x = NC::geomspace( 0.1, 10.0, n ), y;
        x.front() = 0.1;
        x.back() = 10.0;
        for ( auto v : x )
          y.push_back( std::exp( -v ) );
        HymanInterp interp( x, y );
        double maxerr = 0.0;
        for ( double v : NC::geomspace( 0.1, 9.999, 3001 ) )
          maxerr = ncmax( maxerr, NC::ncabs( interp.eval( v ) / std::exp( -v ) - 1.0 ) );
        if ( prev > 0.0 )
          nc_assert_always( prev / maxerr > 10.0 );//2^4=16 ideally
        prev = maxerr;
      }
    }
    //Monotone data with a sharp step: monotone interpolation, no overshoot;
    //zero slope at extrema of the data; nodeSlopes consistent with eval:
    {
      NC::VectD x{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
      NC::VectD y{ 0.0, 0.0, 0.01, 0.02, 5.0, 5.01, 5.02, 5.02 };
      HymanInterp interp( x, y );
      double last = -1.0;
      for ( double v : NC::linspace( 0.0, 7.0, 7001 ) ) {
        const double e = interp.eval( v );
        nc_assert_always( e >= last - 1e-12 && e >= -1e-12 && e <= 5.02 + 1e-12 );
        last = e;
      }
      NC::VectD xe{ 0.0, 1.0, 2.0, 3.0, 4.0 }, ye{ 1.0, 2.0, 3.0, 2.5, 1.0 };
      HymanInterp ie( xe, ye );
      nc_assert_always( ie.nodeSlopes()[2] == 0.0 );
      for ( double v : NC::linspace( 0.0, 4.0, 401 ) )
        nc_assert_always( ie.eval( v ) <= 3.0 + 1e-12 );
      auto m = ie.nodeSlopes();
      for ( auto i : NC::ncrange( xe.size() - 1 ) ) {
        const double h = xe[i+1] - xe[i];
        const double t = 0.3, t2 = t*t, u = 1.0 - t, u2 = u*u;
        const double ref = ( ( 1.0 + 2.0*t ) * u2 * ye[i] + t * u2 * h * m[i]
                             + t2 * ( 3.0 - 2.0*t ) * ye[i+1]
                             + t2 * ( t - 1.0 ) * h * m[i+1] );
        nc_assert_always( NC::floateq( ie.eval( xe[i] + t * h ), ref, 1e-14, 1e-14 ) );
      }
    }
  }

  void testFitPiecewiseCubic()
  {
    //fitPiecewiseCubic: within tol of smooth functions (checked densely),
    //exact at the nodes, and exact for cubics:
    auto check = []( const std::function<double(double)>& f, double a, double b, double tol )
    {
      auto fit = NC::fitPiecewiseCubic( f, a, b, tol );
      nc_assert_always( fit.x.front() == a && fit.x.back() == b );
      nc_assert_always( fit.c.size() == 4 * ( fit.x.size() - 1 ) );
      double worst = 0.0;
      for ( auto i : NC::ncrange( fit.x.size() - 1 ) ) {
        const double lo = fit.x[i], hi = fit.x[i+1];
        nc_assert_always( lo < hi );
        const double * c = &fit.c[4*i];
        nc_assert_always( c[0] == f( lo ) );
        for ( double u : NC::linspace( 0.0, 1.0, 41 ) ) {
          const double v = c[0] + u * ( c[1] + u * ( c[2] + u * c[3] ) );
          worst = std::max( worst, std::fabs( v / f( lo + u * ( hi - lo ) ) - 1.0 ) );
        }
      }
      nc_assert_always( worst < tol );
      return fit.x.size() - 1;
    };
    check( []( double x ) { return std::exp( -x ) * ( 1.0 + x ); }, 1e-3, 30.0, 1e-10 );
    check( []( double x ) { return std::sqrt( x ) + 1.0; }, 1e-6, 1e3, 1e-8 );
    check( []( double x ) { return 1.0 / ( 1.0 + x * x ); }, 0.1, 10.0, 1e-12 );
    const auto ncubic = check( []( double x ) { return 1.0 + x * ( 2.0 + x * ( 0.5 + x ) ); }, 1.0, 2.0, 1e-14 );
    nc_assert_always( ncubic == 4 );//no bisection needed
  }
}

int main() {
  std::cout<<"Launching PCHIPInterp tests"<<std::endl;
  //Fixme: also test PCHIPInterp::OptimiseForLogSpacedX::YES
  testPCHIPInterp1();
  testPCHIPInterp2();
  testHymanInterp();
  testFitPiecewiseCubic();
  std::cout<<"All OK"<<std::endl;
  return 0;
}
