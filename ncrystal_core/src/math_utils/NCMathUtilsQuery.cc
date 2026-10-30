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

#include "NCrystal/internal/math_utils/NCMathUtilsQuery.hh"
#include "NCrystal/internal/math_utils/NCFourierInt.hh"
#include "NCrystal/internal/math_utils/NCConvolve.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <chrono>

namespace NC = NCrystal;
namespace NCMU = NCrystal::MathUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace MathUtils {
    namespace {

      [[noreturn]] void invalid( const Query& query,
                                 const std::string& reason )
      {
        std::ostringstream ss;
        ss << "Invalid mathutils query: ";
        streamJSON( ss, query );
        ss << " (" << reason << ")";
        NCRYSTAL_THROW( BadInput, ss.str() );
      }

      //KEY=VALUE options from query entries [ifirst,end):
      std::map<std::string,StrView> parseOptions( const Query& query,
                                                  std::size_t ifirst )
      {
        std::map<std::string,StrView> res;
        for ( std::size_t i = ifirst; i < query.size(); ++i ) {
          auto parts = query.at(i).trimmed().split<2>('=');
          if ( parts.size() != 2 )
            invalid( query, "options must have the form KEY=VALUE" );
          auto key = parts.at(0).trimmed().to_string();
          if ( res.count(key) )
            invalid( query, "repeated option "+key );
          res[key] = parts.at(1).trimmed();
        }
        return res;
      }

      void queryFourierInt( std::ostream& os, const Query& query )
      {
        if ( query.size() < 5 )
          invalid( query, "usage: [\"mathutils\",\"fourierint\",FR,FI,W,"
                   "OPTIONS...]" );
        VectD fr = takeJSONInputDblVect( query.at(2) );
        VectD fi;
        if ( query.at(3) != "none" )
          fi = takeJSONInputDblVect( query.at(3) );
        VectD w = takeJSONInputDblVect( query.at(4) );
        auto opts = parseOptions( query, 5 );

        Optional<double> h;
        FourierScheme scheme = FourierScheme::PiecewiseLinear;
        FourierEval evalmode = FourierEval::Direct;
        unsigned oversample = 8;
        FourierTail tail = FourierTail::None;
        for ( auto& kv : opts ) {
          const auto& key = kv.first;
          const auto& val = kv.second;
          if ( key == "h" ) {
            h = val.toDbl();
            if ( !h.has_value() || !(h.value()>0.0) || !std::isfinite(h.value()) )
              invalid( query, "invalid value for h" );
          } else if ( key == "scheme" ) {
            if ( val == "trap" ) scheme = FourierScheme::Trapezoid;
            else if ( val == "trapec" ) scheme = FourierScheme::TrapezoidEC;
            else if ( val == "lin" ) scheme = FourierScheme::PiecewiseLinear;
            else if ( val == "cubic" ) scheme = FourierScheme::CubicSpline;
            else invalid( query, "invalid scheme" );
          } else if ( key == "eval" ) {
            if ( val == "direct" ) evalmode = FourierEval::Direct;
            else if ( val == "fftlin" ) evalmode = FourierEval::FFTLinear;
            else if ( val == "fftcubic" ) evalmode = FourierEval::FFTCubic;
            else invalid( query, "invalid eval" );
          } else if ( key == "oversample" ) {
            auto v = val.toInt();
            if ( !v.has_value() || v.value() < 1 || v.value() > 1024 )
              invalid( query, "invalid value for oversample" );
            oversample = static_cast<unsigned>( v.value() );
          } else if ( key == "tail" ) {
            if ( val == "none" ) tail = FourierTail::None;
            else if ( val == "geometric" ) tail = FourierTail::Geometric;
            else invalid( query, "invalid tail" );
          } else {
            invalid( query, "unknown option "+key );
          }
        }
        if ( !h.has_value() )
          invalid( query, "missing option h" );
        if ( fr.size() < 2 )
          invalid( query, "FR must have at least two entries" );
        if ( !fi.empty() && fi.size() != fr.size() )
          invalid( query, "FI and FR must have the same size" );

        auto t0 = std::chrono::steady_clock::now();
        VectD F = fourierIntegral( h.value(), fr, fi, w, scheme, evalmode,
                                   oversample, tail );
        std::chrono::duration<double> dt = std::chrono::steady_clock::now() - t0;
        os << "{\"F\":";
        streamJSONHugeDblVect( os, std::move(F) );
        os << ",\"cpu_seconds\":";
        streamJSON( os, dt.count() );
        os << '}';
      }

      void queryConvolve( std::ostream& os, const Query& query )
      {
        if ( query.size() < 7 )
          invalid( query, "usage: [\"mathutils\",\"convolve\",XF,F,XG,G,Z,"
                   "OPTIONS...]" );
        VectD xf = takeJSONInputDblVect( query.at(2) );
        VectD f = takeJSONInputDblVect( query.at(3) );
        VectD xg = takeJSONInputDblVect( query.at(4) );
        VectD g = takeJSONInputDblVect( query.at(5) );
        VectD z = takeJSONInputDblVect( query.at(6) );
        auto opts = parseOptions( query, 7 );
        ConvolveScheme scheme = ConvolveScheme::ExactLinear;
        double step = 0.0;
        for ( auto& kv : opts ) {
          const auto& key = kv.first;
          const auto& val = kv.second;
          if ( key == "scheme" ) {
            if ( val == "direct" ) scheme = ConvolveScheme::DirectTrapezoid;
            else if ( val == "exactlin" ) scheme = ConvolveScheme::ExactLinear;
            else if ( val == "fft" ) scheme = ConvolveScheme::FFT;
            else invalid( query, "invalid scheme" );
          } else if ( key == "step" ) {
            auto v = val.toDbl();
            if ( !v.has_value() || !(v.value()>=0.0) || !std::isfinite(v.value()) )
              invalid( query, "invalid value for step" );
            step = v.value();
          } else {
            invalid( query, "unknown option "+key );
          }
        }
        auto checkGrid = [&query]( const VectD& x, const VectD& y,
                                   const char * name )
        {
          if ( x.size() < 2 || x.size() != y.size() )
            invalid( query, std::string(name) + " grid must have at least two"
                     " points and the same size as its values" );
          for ( std::size_t i = 1; i < x.size(); ++i )
            if ( !( x[i] > x[i-1] ) )
              invalid( query, std::string(name) + " grid must be strictly"
                       " increasing" );
        };
        checkGrid( xf, f, "XF" );
        checkGrid( xg, g, "XG" );

        auto t0 = std::chrono::steady_clock::now();
        VectD h = convolve( xf, f, xg, g, z, scheme, step );
        std::chrono::duration<double> dt = std::chrono::steady_clock::now() - t0;
        os << "{\"h\":";
        streamJSONHugeDblVect( os, std::move(h) );
        os << ",\"cpu_seconds\":";
        streamJSON( os, dt.count() );
        os << '}';
      }
    }
  }
}

void NCMU::JSONQuery( std::ostream& os, const Query& query )
{
  if ( query.size() < 2 || query.front() != "mathutils" )
    invalid( query, "usage: [\"mathutils\",SUBCOMMAND,...]" );
  const auto& sub = query.at(1);
  if ( sub == "fourierint" )
    queryFourierInt( os, query );
  else if ( sub == "convolve" )
    queryConvolve( os, query );
  else
    invalid( query, "unknown subcommand (available: \"fourierint\","
             " \"convolve\")" );
}
