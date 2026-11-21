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
#include "NCrystal/internal/utils/NCNumFmt.hh"
#include "NCUtilQuery_FMA.hh"
#include "NCQueryImpl.hh"
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
      constexpr auto sv_sincosdeg = StrView::make("sincosdeg");
      constexpr auto sv_fmadiagnose = StrView::make("fmadiagnose");
      constexpr auto sv_browsedb = StrView::make("browsedb");
      constexpr auto sv_factthreads = StrView::make("factorythreads");
      constexpr auto sv_atomdb = StrView::make("atomdb");
      constexpr auto sv_numfmt = StrView::make("numfmt");
      if ( key == sv_list ) {
        if ( nargs != 0 )
          invalid("no arguments should come after: [\"util\",\"list\"]");
        os<<"[\"sincosdeg\", \"fmadiagnose\", \"browsedb\","
          " \"factorythreads\", \"atomdb\", \"numfmt\"]";
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
      } else if ( key == sv_fmadiagnose ) {
        if ( nargs != 0 )
          invalid("no arguments should come after:"
                  " [\"util\",\"fmadiagnose\"]");
        queryimpl_fmadiagnose( os );
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
      } else if ( key == sv_numfmt ) {
        //["util","numfmt",RTOL,GRID,STRIP,MAXDEN,HEXVALUES]: texts of the
        //values from fmtRelTolVector (GRID and STRIP 0 or 1, values encoded
        //as for QueryImpl::decodeHexDoubles), as a single JSON string with
        //the texts separated by spaces (cheap to decode for long vectors):
        auto isBool = []( const StrView& e ) { return e == "0" || e == "1"; };
        if ( nargs != 5 || !isBool( arg(1) ) || !isBool( arg(2) ) )
          invalid("correct usage: [\"util\",\"numfmt\",RTOL,GRID,STRIP,"
                  "MAXDEN,HEXVALUES]");
        auto r = arg(0).toDbl();
        if ( !r.has_value() || !( r.value() >= 0.0 && r.value() <= 1e-3 ) )
          invalid("RTOL must be in [0,1e-3]");
        auto maxden = arg(3).toInt();
        if ( !maxden.has_value() || maxden.value() < 0
             || maxden.value() > 1000 )
          invalid("MAXDEN must be an integer in [0,1000]");
        const VectD x = QueryImpl::decodeHexDoubles( arg(4) );
        for ( auto v : x )
          if ( !std::isfinite( v ) )
            invalid("values must be finite");
        const auto texts
          = fmtRelTolVector( x, r.value(), arg(1) == "1", arg(2) == "1",
                             static_cast<unsigned>( maxden.value() ) );
        std::string res;
        for ( auto& t : texts ) {
          if ( !res.empty() )
            res += ' ';
          res += t;
        }
        //(no characters needing JSON escapes, so no need for streamJSON,
        //which is slow for long strings):
        os << '"' << res << '"';
      } else {
        invalid(nullptr);
      }

    }
  }
}

#endif
