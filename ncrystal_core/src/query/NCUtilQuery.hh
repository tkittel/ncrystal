#ifndef NCrystal_UtilQuery_hh
#define NCrystal_UtilQuery_hh

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

#include "NCrystal/internal/query/NCQuery.hh"
#include "NCBrowseQuery.hh"
#include "NCrystal/threads/NCFactThreads.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCUtilQuery_FMA.hh"
#include <chrono>

namespace NCRYSTAL_NAMESPACE {
  namespace {

    //queryimpl_fmadiagnose_cmulDispatch: see NCUtilQuery_FMA.hh (included
    //above).

    //queryimpl_fmadiagnose_cmulXXX: Version without NCRYSTAL_FMADISPATCH_ATTR
    void queryimpl_fmadiagnose_cmulReference( double* re, double* im,
                                              const double* cre,
                                              const double* cim,
                                              std::size_t n )
    {
      for ( std::size_t i = 0; i < n; ++i ) {
        double a = re[i], b = im[i], c = cre[i], d = cim[i];
        re[i] = ncfma( a, c, -(b*d) );
        im[i] = ncfma( a, d, b*c );
      }
    }

    void queryimpl_fmadiagnose( std::ostream& os )
    {
      const std::size_t n = 32768;
      VectD re0(n), im0(n), cre(n), cim(n);
      for ( auto i : ncrange(n) ) {
        const double x = static_cast<double>(i);
        re0[i] = 1.0/(1.0+0.001*x);
        im0[i] = 0.3*std::sin(0.01*x);
        cre[i] = std::cos(0.02*x);
        cim[i] = std::sin(0.02*x);
      }

      VectD reA(re0), imA(im0), reB(re0), imB(im0);
      queryimpl_fmadiagnose_cmulDispatch( reA.data(), imA.data(),
                                          cre.data(), cim.data(), n );
      queryimpl_fmadiagnose_cmulReference( reB.data(), imB.data(),
                                           cre.data(), cim.data(), n );
      const bool correct = ( reA == reB && imA == imB );

      //Best-of-ntrials timing, in ns/element (kept fast: this is a diagnostic
      //query, not a benchmark):
      auto timeIt = []( void(*f)(double*,double*,const double*,
                                 const double*,std::size_t),
                        const VectD& re0_, const VectD& im0_,
                        const VectD& cre_, const VectD& cim_ )
      {
        const std::size_t nn = re0_.size();
        const int ntrials = 5, nrep = 20;
        double best = 1e300;
        for ( auto trial : ncrange(ntrials) ) {
          (void)trial;
          VectD re(re0_), im(im0_);
          const auto t0 = std::chrono::steady_clock::now();
          for ( auto r : ncrange(nrep) ) {
            (void)r;
            f( re.data(), im.data(), cre_.data(), cim_.data(), nn );
          }
          const auto t1 = std::chrono::steady_clock::now();
          volatile double sink = re[0] + im[0];
          (void)sink;
          const double ns
            = std::chrono::duration<double,std::nano>(t1-t0).count();
          best = ncmin( best, ns / ( static_cast<double>(nrep)
                                    * static_cast<double>(nn) ) );
        }
        return best;
      };
      const double tDispatch = timeIt( queryimpl_fmadiagnose_cmulDispatch,
                                       re0, im0, cre, cim );
      const double tReference = timeIt( queryimpl_fmadiagnose_cmulReference,
                                        re0, im0, cre, cim );
      const double speedup = tReference / tDispatch;

      //Whether this is an optimised (NDEBUG) build: in an unoptimised
      //(e.g. simplebuild debug cache) build, target_clones clones do not
      //get the optimisation needed to show their benefit, so callers
      //should not enforce a speedup requirement unless this is true:
#ifdef NDEBUG
      constexpr bool optimised = true;
#else
      constexpr bool optimised = false;
#endif
      streamJSONDictEntry( os, "enabled", bool(NCRYSTAL_FMADISPATCH_ENABLED),
                           JSONDictPos::FIRST );
      streamJSONDictEntry( os, "correct", correct );
      streamJSONDictEntry( os, "optimised", optimised );
      streamJSONDictEntry( os, "speedup", speedup, JSONDictPos::LAST );
    }

    void queryimpl_version( std::ostream& os, const Query& query )
    {
      if ( query.size() != 1 || query.front() != "version" )
        NCRYSTAL_THROW2(BadInput, "Invalid JSON query (no arguments"
                        " should follow key \"version\")");
      int vmajor = NCRYSTAL_VERSION_MAJOR;
      int vminor = NCRYSTAL_VERSION_MINOR;
      int vpatch = NCRYSTAL_VERSION_PATCH;
      int vint = NCRYSTAL_VERSION;
      if ( ncgetenv_bool("FIX_QUERY_VERSION_FOR_TESTS") ) {
        vmajor = 123;
        vminor = 45;
        vpatch = 6;
        vint = 1000000*vmajor+1000*vminor+vpatch;
      }
      std::array<int,3> vlist = { vmajor, vminor, vpatch };
      streamJSONDictEntry( os, "major", vmajor, JSONDictPos::FIRST );
      streamJSONDictEntry( os, "minor", vminor );
      streamJSONDictEntry( os, "patch", vpatch );
      streamJSONDictEntry( os, "int", vint );
      streamJSONDictEntry( os, "list", vlist );
      std::ostringstream vstr;
      vstr << vmajor << '.' << vminor << '.' << vpatch;
      streamJSONDictEntry( os, "str", vstr.str(), JSONDictPos::LAST );
    }

    void queryimpl_kinutils( std::ostream& os )
    {
      //Dumps getAlphaMinus/getAlphaPlus/getBetaMinus/getBetaPlus on a grid of
      //ekin_div_kT values, with offsets landing densely around the two singular
      //points, for validation against an mpmath reference:
      const VectD evals = { 1e-6, 1e-3, 0.1, 1.0, 1e2, 1e4, 1e8 };
      const VectD offsets =
        { -2.0, -1.0, -0.5, -0.2, -0.1, -0.06, -0.055,
          -0.051, -0.05, -0.04, -0.02, -0.01, -0.001,
          -1e-4, -1e-6, -1e-9, -1e-12, -1e-15, 0.0,
          1e-15, 1e-12, 1e-9, 1e-6, 1e-4, 0.001, 0.01,
          0.02, 0.04, 0.05, 0.051, 0.055, 0.06, 0.1,
          0.2, 0.5, 1.0, 2.0, 5.0 };
      VectD ab_e, ab_beta, ab_aminus, ab_aplus;
      VectD bb_e, bb_alpha, bb_bminus, bb_bplus;
      for ( auto e_val : evals ) {
        const EDivKT e{ e_val };
        for ( auto x : offsets ) {
          //alpha_minus/alpha_plus(e,beta), beta=x*e, singular at x=0.
          //Valid domain needs beta>=-e (kk=e+beta>=0), i.e. x>=-1:
          if ( x >= -1.0 ) {
            const double beta = x*e.dbl();
            ab_e.push_back(e.dbl());
            ab_beta.push_back(beta);
            ab_aminus.push_back( getAlphaMinus(e,beta) );
            ab_aplus.push_back( getAlphaPlus(e,beta) );
          }
          //beta_minus/beta_plus(e,alpha), alpha=(4+x)*e, singular at
          //x=0. Valid domain needs alpha>=0, i.e. x>=-4 (always true
          //here since offsets>=-2):
          {
            const double alpha = (4.0+x)*e.dbl();
            bb_e.push_back(e.dbl());
            bb_alpha.push_back(alpha);
            bb_bminus.push_back( getBetaMinus(e,alpha) );
            bb_bplus.push_back( getBetaPlus(e,alpha) );
          }
        }
      }
      os << "{\"alphabounds\":{\"e\":";
      streamJSON(os,ab_e);
      os << ",\"beta\":";
      streamJSON(os,ab_beta);
      os << ",\"aminus\":";
      streamJSON(os,ab_aminus);
      os << ",\"aplus\":";
      streamJSON(os,ab_aplus);
      os << "},\"betabounds\":{\"e\":";
      streamJSON(os,bb_e);
      os << ",\"alpha\":";
      streamJSON(os,bb_alpha);
      os << ",\"bminus\":";
      streamJSON(os,bb_bminus);
      os << ",\"bplus\":";
      streamJSON(os,bb_bplus);
      os << "}}";
    }

    void queryimpl_util( std::ostream& os, const Query& query )
    {
      auto invalid = [&query](const char * reason){
        std::ostringstream ss;
        ss << "Invalid util query: ";
        streamJSON( ss, query );//trick: use streamJSON for easy format.
        if ( reason )
          ss<< " ("<<reason<<')';
        NCRYSTAL_THROW( BadInput, ss.str() );
      };
      constexpr auto sv_util = StrView::make("util");
      if ( query.size() < 2 || query.front() != sv_util ) {
        invalid(nullptr);
        return;
      }
      //shift off the "util" and key entries:
      const auto& key = query.at(1).trimmed();
      auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
      //auto argstr = [&arg]( std::size_t i ) { return arg(i).to_string(); };
      const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

      constexpr auto sv_list = StrView::make("list");
      constexpr auto sv_wl2ekin = StrView::make("wl2ekin");
      constexpr auto sv_ekin2wl = StrView::make("ekin2wl");
      constexpr auto sv_mathval = StrView::make("mathval");
      constexpr auto sv_sincosdeg = StrView::make("sincosdeg");
      constexpr auto sv_fmadiagnose = StrView::make("fmadiagnose");
      constexpr auto sv_kinutils = StrView::make("kinutils");
      constexpr auto sv_browsedb = StrView::make("browsedb");
      constexpr auto sv_browsefacts = StrView::make("browsefactories");
      constexpr auto sv_factthreads = StrView::make("factorythreads");
      constexpr auto sv_atomdb = StrView::make("atomdb");
      constexpr auto sv_eqrefl = StrView::make("eqrefl");
      constexpr auto sv_jsontest = StrView::make("jsontest");
      if ( key == sv_list ) {
        if ( nargs != 0 )
          invalid("no arguments should come after: [\"util\",\"list\"]");
        os<<"[\"wl2ekin\", \"ekin2wl\", \"mathval\", \"sincosdeg\","
          " \"fmadiagnose\","
          " \"kinutils\", \"browsedb\", \"browsefactories\","
          " \"factorythreads\", \"atomdb\", \"eqrefl\", \"jsontest\"]";
      } else if ( key == sv_eqrefl ) {
        //["util","eqrefl",MODE,SG,"h1,k1,l1","h2,k2,l2",...]: For each
        //requested hkl, return the list of reflections considered
        //symmetry-equivalent by NCrystal for the space group setting SG
        //(parsed in MODE "strict" or "relaxed", cf. SpaceGroupParseMode; e.g.
        //"166:r" for rhombohedral axes), one of each Friedel pair (-h,-k,-l
        //is implied).
        const auto usage = ( "correct usage: [\"util\",\"eqrefl\","
                             "\"strict|relaxed\",SG,"
                             "\"h0,k0,l0\",\"h1,k1,l1\",...]" );

        if ( nargs < 3 || ( arg(0) != "strict" && arg(0) != "relaxed" ) )
          invalid( usage );
        const SpaceGroup sg( arg(1), ( arg(0) == "strict"
                                       ? SpaceGroupParseMode::Strict
                                       : SpaceGroupParseMode::Relaxed ) );
        EqRefl eqrefl( sg );
        os << '[';
        for ( auto i : ncrange( std::size_t(2), nargs ) ) {
          auto parts = arg(i).splitTrimmed(',');
          int hkl[3];
          for ( auto j : ncrange( 3 ) ) {
            auto x = ( parts.size() == 3 ? parts.at(j).toInt32() : NullOpt );
            if ( !x.has_value() )
              invalid( usage );
            hkl[j] = x.value();
          }
          os << ( i>2 ? ",[" : "[" );
          bool first = true;
          for ( auto& e : eqrefl.getEquivalentReflections( hkl[0], hkl[1],
                                                           hkl[2] ) ) {
            os << ( first ? "[" : ",[" )
               << e.h << ',' << e.k << ',' << e.l << ']';
            first = false;
          }
          os << ']';
        }
        os << ']';
      } else if ( key == sv_jsontest ) {
        //["util","jsontest",BITS1,BITS2,...]: doubles given by their IEEE 754
        //bit patterns (as decimal integers), written by streamJSON one at a
        //time and as an array by streamJSONHugeDblVect (which only uses its
        //fallback unless huge arrays are enabled): {"scalar":[..],"array":[..]}
        static_assert( sizeof(double) == sizeof(std::uint64_t), "" );
        VectD v;
        v.reserve( nargs );
        for ( auto i : ncrange( nargs ) ) {
          auto bits = arg(i).toUInt64();
          if ( !bits.has_value() )
            invalid( "correct usage: [\"util\",\"jsontest\",BITS1,...] with"
                     " IEEE 754 bit patterns of doubles as decimal integers" );
          double x;
          std::memcpy( &x, &bits.value(), sizeof(x) );
          v.push_back( x );
        }
        os << "{\"scalar\":[";
        for ( auto i : ncrange( v.size() ) ) {
          if ( i )
            os << ',';
          streamJSON( os, v[i] );
        }
        os << "],\"array\":";
        streamJSONHugeDblVect( os, std::move( v ) );
        os << '}';
      } else if ( key == sv_atomdb ) {
        if ( nargs != 0 )
          invalid("no arguments should come after: [\"util\",\"atomdb\"]");
        BrowseQuery::atomDB( os );
      } else if ( key == sv_factthreads ) {
        if ( nargs != 0 )
          invalid("no arguments should come after:"
                  " [\"util\",\"factorythreads\"]");
        //NB: These also apply any NCRYSTAL_FACTORY_THREADS setting:
        streamJSONDictEntry( os, "threads_available",
                             FactoryThreadPool::threadsAvailable(),
                             JSONDictPos::FIRST );
        streamJSONDictEntry( os, "user_configured",
                             FactoryThreadPool::userConfigured() );
        const auto nt = FactoryThreadPool::currentThreadCount().get();
        streamJSONDictEntry( os, "nthreads",
                             static_cast<unsigned>( nt ),
                             JSONDictPos::LAST );
      } else if ( key == sv_browsedb ) {
        std::vector<StrView> args;
        for ( auto i : ncrange( nargs ) )
          args.push_back( arg(i) );
        BrowseQuery::browseDB( os, args );
      } else if ( key == sv_browsefacts ) {
        if ( nargs != 0 )
          invalid("no arguments should come after:"
                  " [\"util\",\"browsefactories\"]");
        BrowseQuery::browseFactories( os );
      } else if ( key == sv_fmadiagnose ) {
        if ( nargs != 0 )
          invalid("no arguments should come after:"
                  " [\"util\",\"fmadiagnose\"]");
        queryimpl_fmadiagnose( os );
      } else if ( key == sv_kinutils ) {
        if ( nargs != 0 )
          invalid("no arguments should come after:"
                  " [\"util\",\"kinutils\"]");
        queryimpl_kinutils( os );
      } else if ( isOneOf(key,sv_wl2ekin,sv_ekin2wl) ) {
        double val = ( nargs == 1
                       ? arg(0).toDbl().value_or(-1.0)
                       : -1.0 );
        if ( !(val>=0.0) )
          invalid( key==sv_wl2ekin
                   ? "correct usage: [\"util\",\"wl2ekin\",VALUE_ANGSTROM]"
                   : "correct usage: [\"util\",\"ekin2wl\",VALUE_EV]" );
        streamJSON( os, ( key==sv_wl2ekin
                          ? wl2ekin( val )
                          : ekin2wl( val ) ) );
      } else if ( key == sv_sincosdeg ) {
        //["util","sincosdeg",x1,x2,...]: [[sin,cos],...] of the angles in
        //degrees, from stable_sincos_deg (cross platform reproducible). Angles
        //can be prefixed with '@' (needed for negative values):
        if ( nargs == 0 )
          invalid("correct usage: [\"util\",\"sincosdeg\",x1,x2,...]");
        os << '[';
        for ( auto i : ncrange( nargs ) ) {
          auto x = ( arg(i).startswith('@') ? arg(i).substr(1)
                                            : arg(i) ).toDbl();
          if ( !x.has_value() || !std::isfinite( x.value() ) )
            invalid("angles must be finite numbers");
          const auto sc = stable_sincos_deg( x.value() );
          os << ( i ? ",[" : "[" );
          streamJSON( os, sc.sin );
          os << ',';
          streamJSON( os, sc.cos );
          os << ']';
        }
        os << ']';
      } else if ( key == sv_mathval ) {
        if ( nargs != 1 || arg(0) != "kpowx" )
          invalid("Only supported mathval arg is for now [\"kpowx\"]");
        VectD testsamplek = { 1e-300, 1e-200, 1e-30, 1e-10, 1e-3, 0.5,
          0.95, 0.99,
          1.0-1.1e-4, 1.0-1e-4, 1.0-0.9e-4, 1.0,
          1.0+1.1e-4, 1.0+1e-4, 1.0+0.9e-4,
          1.01, 2.0, 1e2, 1e5, 1e30, 1e200, 1e300 };
        VectD testsampleR = { 1e-300,1e-20,1e-5,1e-2,0.1,0.5,0.9,
                              0.9899,0.9901,1.0-1e-5,1.0-1e-10,1.0 };
        VectD testsample;
        testsample.reserve( testsamplek.size() * testsampleR.size() );
        for ( auto& k : testsamplek ) {
          for ( auto& R : testsampleR )
            testsample.push_back(randKPowX( k, std::log(k), R ));
        }

        os << "{\"sample_kpowx\":{";
        os << "\"k\":";
        streamJSON(os,testsamplek);
        os << ",\"R\":";
        streamJSON(os,testsampleR);
        os << ",\"samples\":";
        streamJSON(os,testsample);
        os<<"}}";
      } else {
        invalid(nullptr);
      }

    }
  }
}

#endif
