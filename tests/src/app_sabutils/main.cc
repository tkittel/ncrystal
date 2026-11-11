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

#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/interfaces/NCInfoTypes.hh"
#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/atomdb/NCAtomDB.hh"
#include "NCrystal/core/NCFmt.hh"
#include <iostream>
namespace NC = NCrystal;
namespace NS = NC::SABUtils;

void printrange(NC::Span<double> v,bool eol = true) {
  //TODO: we should think about numpy-like adapters to print span's of values in
  //ostreams. Perhaps even overloading adapters: std::cout << NC::stream(v)
  //<<std::endl; std::cout << NC::stream_shortform(v) <<std::endl;
  std::cout<<"[";
  if (!v.empty()) {
    for (auto e: NC::Span<const double>(v).first(v.size()-1))
      std::cout<<e<<", ";
    std::cout<<v.back();
  }
  std::cout<<"]";
  if (eol)
    std::cout<<std::endl;
}

void test_loglin( long double a, long double fa, long double b, long double fb, long double x, long double calculated_value ) {
  double expected;
  std::ostringstream ss;
  ss << "loglin(a="<<a<<", fa="<<fa<<", b="<<b<<", fb="<<fb<<", x="<<x<<")";
  bool linlin(false);
  if (fa==fb) {
    expected = fa;//trivial
  } else if ( fa!=0 && fb!=0 ) {
    if ( x == a ) {
      expected = fa;
    } else if ( x == b ) {
      expected = fb;
    } else {
      expected = fa*std::pow(fb/fa,(x-a)/(b-a));//efficient rewrite of exp(log(fa)+(logfb-logfa)*r)
      nc_assert_always( NC::floateq( expected,
                                     std::exp(std::log(fa)+(std::log(fb/fa))*(x-a)/(b-a) ),
                                     1.0e-8 ) );
    }
  } else {
    linlin = true;
    if ( x == a ) {
      expected = fa;
    } else if ( x == b ) {
      expected = fb;
    } else {
#if 1
      expected = fa+(fb-fa)*(x-a)/(b-a);//linlin
#else
      NC::StableSum sum;
      sum.add(fa*b);
      sum.add(fb*x);
      sum.add(-fb*a);
      sum.add(-fa*x);
      expected = sum.sum() / ( b-a );
#endif
    }
  }
  if (!NC::floateq(expected, calculated_value,1e-12,1e-12))
    NCRYSTAL_THROW2(CalcError,ss.str()<<" gave "<<NC::fmt(calculated_value)
                    <<" not expected "<<NC::fmt(expected)<<" (absdiff "<<
                    calculated_value-expected<<")"
                    <<(linlin?"linlin":"loglin")<<")");
  //std::cout<<"Testing "<<ss.str()<<" [expected "<<expected<<", got "<<calculated_value<<"] ("<<(linlin?"linlin":"loglin")<<")"<<std::endl;

}

void test_alpha_integrals() {
  //Reference values created with:
  //
  //#!/usr/bin/env python3
  //import numpy as np
  //import mpmath
  //mpmath.mp.dps = 1000
  //mpf=mpmath.mpf
  //def exact_eval_mpf(y):
  //    if y==0.0:
  //        return mpf(0.5)
  //    s1=(mpf(1)-mpf(y))/mpf(2)
  //    s2=(mpf(1)+mpf(y))/mpf(2)
  //    return (s2-s1)/mpmath.mp.log(s2/s1)
  //yvals = np.asarray([1e-200,1e-20,1e-10,1e-4,0.005,0.007,0.09,0.11,0.2,0.5,0.9,0.99,1.0-1e-14])
  //yvals = np.hstack((-np.flip(yvals),np.asarray([0.0]),yvals))
  //for y in yvals:
  //    print('   ref.emplace_back(',y,',',mpmath.nstr(exact_eval_mpf(y),30),');')
  //
  std::vector<NC::PairDD> ref;
  ref.emplace_back( -0.99999999999999 , 0.0303673187635422435892942080719 );
  ref.emplace_back( -0.99 , 0.18702871509983897278285159318 );
  ref.emplace_back( -0.9 , 0.305660944705597718728739460966 );
  ref.emplace_back( -0.5 , 0.455119613313418696807120082868 );
  ref.emplace_back( -0.2 , 0.493260692475286336445129433371 );
  ref.emplace_back( -0.11 , 0.497976784653061267607378554863 );
  ref.emplace_back( -0.09 , 0.49864707156245853926222023281 );
  ref.emplace_back( -0.007 , 0.499991833226619483216329399507 );
  ref.emplace_back( -0.005 , 0.499995833305555191792870951536 );
  ref.emplace_back( -0.0001 , 0.499999998333333328888888705869 );
  ref.emplace_back( -1e-10 , 0.499999999999999999998333333333 );
  ref.emplace_back( -1e-20 , 0.5 );
  ref.emplace_back( -1e-200 , 0.5 );
  ref.emplace_back( 0.0 , 0.5 );
  ref.emplace_back( 1e-200 , 0.5 );
  ref.emplace_back( 1e-20 , 0.5 );
  ref.emplace_back( 1e-10 , 0.499999999999999999998333333333 );
  ref.emplace_back( 0.0001 , 0.499999998333333328888888705869 );
  ref.emplace_back( 0.005 , 0.499995833305555191792870951536 );
  ref.emplace_back( 0.007 , 0.499991833226619483216329399507 );
  ref.emplace_back( 0.09 , 0.49864707156245853926222023281 );
  ref.emplace_back( 0.11 , 0.497976784653061267607378554863 );
  ref.emplace_back( 0.2 , 0.493260692475286336445129433371 );
  ref.emplace_back( 0.5 , 0.455119613313418696807120082868 );
  ref.emplace_back( 0.9 , 0.305660944705597718728739460966 );
  ref.emplace_back( 0.99 , 0.18702871509983897278285159318 );
  ref.emplace_back( 0.99999999999999 , 0.0303673187635422435892942080719 );

  //for (auto& [ y, refval ] : ref ) {
  for (auto& e : ref ) {
    // auto& y = e.first;
    // auto& refval = e.second;
    const double y = e.first;
    const double refval = e.second;
    const double s1 = (1.0-y)/2.0;
    const double s2 = (1.0+y)/2.0;
    const double res = NS::integrateLogLinOrLin(0.0, s1, 1.0, s2);//claims precision better than 1e-15
    const double resfast = NS::integrateLogLinOrLinWCL(0.0, s1, 1.0, s2, std::log(s1), std::log(s2) );//claims precision better than 1e-14
    if (!NC::floateq(res,refval,1e-15,0.0))
      NCRYSTAL_THROW2(LogicError,"problem with integrateLogLinOrLin for y="
                      <<NC::fmt(y)<<" : "<<res<<" differs too much from expected "<<refval);
    if (!NC::floateq(resfast,refval,1e-14,0.0))
      NCRYSTAL_THROW2(LogicError,"problem with integrateLogLinOrLinWCL for y="
                      <<NC::fmt(y)<<" : "<<resfast<<" differs too much from expected "<<refval);
  }
}

namespace {
  //Test that DI_ScatKnlDirect objects whose buildSAB() method provides a
  //kernel at the wrong temperature, consistently fail (previously only the
  //first call failed, while subsequent calls returned the invalid kernel).
  class TestDI final : public NC::DI_ScatKnlDirect {
    NC::Temperature m_tknl;
  public:
    TestDI( NC::Temperature t, NC::Temperature tknl )
      : NC::DI_ScatKnlDirect( 1.0,
                              NC::IndexedAtomData{ NC::AtomDB::getNaturalElement("Al"),
                                                   NC::AtomIndex{0} },
                              t ),
        m_tknl(tknl)
    {
    }
    EGridShPtr energyGrid() const override { return nullptr; }
  protected:
    std::shared_ptr<const NC::SABData> buildSAB() const override
    {
      return std::make_shared<const NC::SABData>( NC::VectD{ 0.1, 0.2 },
                                                  NC::VectD{ -0.1, 0.1 },
                                                  NC::VectD{ 1.0, 1.0, 1.0, 1.0 },
                                                  m_tknl,
                                                  NC::SigmaBound{ 1.0 },
                                                  NC::AtomMass{ 27.0 } );
    }
  };

  void testScatKnlDirect( NC::Temperature t, NC::Temperature tknl )
  {
    TestDI di( t, tknl );
    for ( unsigned i = 0; i < 3; ++i ) {
      nc_assert_always( di.hasBuiltSAB() == ( i > 0 && t == tknl ) );
      try {
        auto sab = di.ensureBuildThenReturnSAB();
        nc_assert_always( sab != nullptr );
        nc_assert_always( t == tknl );
        std::cout << "DI_ScatKnlDirect: Got kernel at T="
                  << sab->temperature() << std::endl;
      } catch ( NC::Error::BadInput& ) {
        nc_assert_always( t != tknl );
        std::cout << "DI_ScatKnlDirect: Got BadInput error" << std::endl;
      }
    }
  }
}

namespace {
  NC::SABData makeSAB( NC::VectD sab )
  {
    return NC::SABData( NC::VectD{ 0.5, 1.0, 2.0, 4.0 },
                        NC::VectD{ -1.0, 0.0, 1.5 },
                        std::move( sab ), NC::Temperature{ 300.0 },
                        NC::SigmaBound{ 1.0 }, NC::AtomMass{ 1.0 } );
  }

  //Expected value: log-linear in alpha (linear if a corner is 0), then
  //linear in beta, within cell (ia,ib):
  double expectedSAB( const NC::SABData& sd, std::size_t ia, std::size_t ib,
                      double alpha, double beta )
  {
    const auto& A = sd.alphaGrid();
    const auto& B = sd.betaGrid();
    auto s = [&sd,&A]( std::size_t i, std::size_t j )
    { return sd.sab().at( j * A.size() + i ); };
    auto inA = [&]( std::size_t j ) {
      const double fa = s( ia, j ), fb = s( ia + 1, j );
      const double r = ( alpha - A.at(ia) ) / ( A.at(ia+1) - A.at(ia) );
      return ( fa > 0.0 && fb > 0.0
               ? std::exp( std::log(fa) + r * std::log( fb / fa ) )
               : fa + r * ( fb - fa ) );
    };
    const double u = ( beta - B.at(ib) ) / ( B.at(ib+1) - B.at(ib) );
    return inA( ib ) + u * ( inA( ib + 1 ) - inA( ib ) );
  }

  void test_evalSABOneShot()
  {
    //4x3 grid, one cell with zero corners (alpha 2..4, beta 0..1.5):
    NC::VectD sab = { 0.3, 0.5, 0.2, 0.1,
                      1.0, 2.0, 0.7, 0.0,
                      0.4, 0.9, 0.0, 0.0 };
    const auto sd = makeSAB( sab );
    const auto& A = sd.alphaGrid();
    const auto& B = sd.betaGrid();
    //Exact at all grid points (incl. the upper edges):
    for ( auto ib : NC::ncrange( B.size() ) )
      for ( auto ia : NC::ncrange( A.size() ) )
        nc_assert_always( NS::evalSABOneShot( sd, A.at(ia), B.at(ib) )
                          == sab.at( ib * A.size() + ia ) );
    //Inside each cell:
    for ( auto ib : NC::ncrange( B.size() - 1 ) )
      for ( auto ia : NC::ncrange( A.size() - 1 ) )
        for ( double ra : { 0.1, 0.5, 0.77 } )
          for ( double rb : { 0.2, 0.5, 0.9 } ) {
            const double a = A.at(ia) + ra * ( A.at(ia+1) - A.at(ia) );
            const double b = B.at(ib) + rb * ( B.at(ib+1) - B.at(ib) );
            const double v = NS::evalSABOneShot( sd, a, b );
            nc_assert_always( NC::floateq( v, expectedSAB( sd, ia, ib, a, b ),
                                           1e-13, 0.0 ) );
          }
    //Outside the grid:
    for ( auto ab : { NC::PairDD{ 0.49, 0.0 }, NC::PairDD{ 4.01, 0.0 },
                      NC::PairDD{ 1.0, -1.01 }, NC::PairDD{ 1.0, 1.51 } } )
      nc_assert_always( NS::evalSABOneShot( sd, ab.first, ab.second ) == 0.0 );
    //Subnormal S: same as the table scaled up by 2^1040 (exact), scaled down
    //again (to within the precision of the subnormal result):
    NC::VectD sab_sub, sab_up;
    for ( auto v : sab ) {
      sab_sub.push_back( std::ldexp( v, -1040 ) );
      sab_up.push_back( std::ldexp( sab_sub.back(), 1040 ) );
    }
    const auto sd_sub = makeSAB( sab_sub );
    const auto sd_up = makeSAB( sab_up );
    for ( double a : { 0.5, 0.7, 1.9, 3.0 } )
      for ( double b : { -1.0, -0.3, 0.9, 1.5 } ) {
        const double v = NS::evalSABOneShot( sd_sub, a, b );
        nc_assert_always( v == 0.0 || v < std::numeric_limits<double>::min() );
        nc_assert_always( NC::floateq( std::ldexp( v, 1040 ),
                                       NS::evalSABOneShot( sd_up, a, b ),
                                       1e-8, 0.0 ) );
      }
    std::cout << "evalSABOneShot tests OK" << std::endl;
  }

  //interpLogLinWCL for extreme ratios fb/fa (where exp(dl*r) alone would
  //overflow or become subnormal): finite, exact at the nodes, within
  //[fa,fb], and accurate (long double reference):
  void test_interpLogLinWCL_extreme()
  {
    const double vals[] = { 4.9406564584124654e-324, 1e-310, 2.5e-308, 1e-300,
                            1.0, 1e250, 1e300, 1.7e308 };
    for ( double fa : vals ) {
      for ( double fb : vals ) {
        const double la = std::log( fa );
        const double lb = std::log( fb );
        for ( double r : { 0.0, 1e-9, 0.1, 0.3, 0.4999, 0.5, 0.7, 1.0 } ) {
          const double a = 2.0;
          const double b = 3.0;
          const double x = ( r == 1.0 ? b : a + r * ( b - a ) );
          const auto res = NS::interpLogLinWCL( a, fa, b, fb, x, la, lb );
          nc_assert_always( std::isfinite( res.first ) );
          nc_assert_always( std::isfinite( res.second ) );
          nc_assert_always( res.first >= NC::ncmin( fa, fb ) );
          nc_assert_always( res.first <= NC::ncmax( fa, fb ) );
          if ( x == a )
            nc_assert_always( res.first == fa );
          if ( x == b )
            nc_assert_always( res.first == fb );
          const long double t = ( (long double)x - a ) / ( (long double)b - a );
          const long double ref = std::exp( (long double)la
                                            + ( (long double)lb - la ) * t );
          const double refd = static_cast<double>( ref );
          if ( refd >= std::numeric_limits<double>::min() )
            nc_assert_always( NC::floateq( res.first, refd, 1e-12, 0.0 ) );
          else//subnormal results: absolute precision of the subnormal range
            nc_assert_always( NC::ncabs( res.first - refd ) <= 1e-12
                              * std::numeric_limits<double>::min() );
        }
      }
    }
    std::cout << "interpLogLinWCL extreme tests OK" << std::endl;
  }
}

int main()
{
  test_evalSABOneShot();
  test_interpLogLinWCL_extreme();
  testScatKnlDirect( NC::Temperature{ 300.0 }, NC::Temperature{ 300.0 } );
  testScatKnlDirect( NC::Temperature{ 300.0 }, NC::Temperature{ 200.0 } );

  {
    //Unit test alpha limits taylor expansion (really for a dedicated KinUtils test)
    const double ekin_div_kT = 6.345703125e-08;//probably doesn't matter (using value where we once saw a failure)
    const double betalim=0.01*ekin_div_kT;//NB: 0.01 should be same factor as used in alphaMinusNeedsTaylor!
    const double b0 = betalim*(1+1e-15);
    const double b1 = betalim*(1-1e-15);
    const NC::EDivKT e{ ekin_div_kT };
    nc_assert_always( NC::detail::alphaMinusNeedsTaylor(ekin_div_kT,b0) != NC::detail::alphaMinusNeedsTaylor(ekin_div_kT,b1) );
    nc_assert_always( NC::detail::alphaMinusNeedsTaylor(ekin_div_kT,-b0) != NC::detail::alphaMinusNeedsTaylor(ekin_div_kT,-b1) );
    nc_assert_always( NC::floateq( NC::getAlphaLimits(e,b0).first, NC::getAlphaLimits(e,b1).first, 1e-11, 1e-299 ) );
    nc_assert_always( NC::floateq( NC::getAlphaLimits(e,-b0).first, NC::getAlphaLimits(e,-b1).first, 1e-11, 1e-299 ) );
  }

  ////////////////////////

  using V = NCrystal::VectD;

  ////////////////////////

  V outb, outsab;
  auto insab = V{  0.111, 0.333, 0.777,
                   1.111, 1.333, 1.777,
                   4.111, 4.333, 4.777 };
  NS::expandBetaAndSABToAllBetas( NC::Span<const double>(V{0,1,4}),
                                  V{ 11.1, 33.3, 77.7 },
                                  insab,
                                  outb,outsab);

  nc_assert_always(outb==V({-4, -1, 0, 1, 4}));
  nc_assert_always(outsab==V({4.111, 4.333, 4.777, 1.111, 1.333, 1.777, 0.111, 0.333, 0.777, 1.111, 1.333, 1.777, 4.111, 4.333, 4.777}));

  nc_assert_always( outsab.size() == insab.size() + 6);
  for ( std::size_t i = 0; i < insab.size(); ++ i) {
    nc_assert_always( insab.at(i) == outsab.at(i+6) );
  }

  printrange(outb);
  printrange(outsab);

  //////////////////////////

  for (auto& a :  { -1e5, -1e-5, 0.0, 0.5, 10.0 ,1e5 } ) {
    for (auto& brel :  { 0.001, 0.9, 1000.0 } ) {
      double b = a + std::abs(a?a:1.0)*brel;
      for (auto& fa :  { 0.0, 1e-5, 1.0, 1e5 } ) {
        for (auto& fb :  { 0.0, 1e-5, 1.0, 1e5 } ) {
          for (auto& xrel : { 0.0, 0.1, 0.5, 0.8, 1.0 } ) {
            double x = (xrel == 1.0 ? b : a + xrel * (b-a) );
            test_loglin( a,fa,b,fb,x, NS::interpolate_loglin_fallbacklinlin( a,fa,b,fb,x ) );
            double logfa(fa>0.0?std::log(fa):-1000.0);
            double logfb(fb>0.0?std::log(fb):-1000.0);
            test_loglin( a,fa,b,fb,x, NS::interpLogLinOrLinWCL( a,fa,b,fb,x, logfa, logfb ).first );
          }
        }
      }
    }
  }

  //////////////////////////
  test_alpha_integrals();

}
