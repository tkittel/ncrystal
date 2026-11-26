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

#include "NCrystal/internal/cif/NCCIFTemperature.hh"
#include "NCrystal/internal/cif/NCCIFUtils.hh"

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    namespace {
      constexpr double tmin = 0.1;
      constexpr double tmax = 3000.0;
      constexpr double tfallback = 293.0;

      bool consistent( double a, double b )
      {
        return ncabs( a - b ) < 0.05 * ncmax( a, b );
      }

      struct Field {
        const char * name;
        std::string value;//(printed)
        Optional<double> number;
      };
    }
  }
}

NCC::CIFUisoTemperature
NCC::resolveUisoTemperature( StrView ambient, StrView measurement,
                             Optional<double> chosen, CIFPermissive perm,
                             const CIFWarnFct& warnfct )
{
#ifndef NDEBUG
  nc_assert_always( bool(warnfct) );
#endif
  const bool permissive = ( perm == CIFPermissive::YES );
  std::vector<Field> present;
  auto add = [&present]( const char * name, StrView v ) {
    if ( !v.has_value() )
      return;
    Field f{ name, v.to_string(), NullOpt };
    const auto x = parseCIFNumber( v );
    if ( x.has_value() )
      f.number = x.value().value;
    present.push_back( std::move( f ) );
  };
  add( "_diffrn_ambient_temperature", ambient );
  add( "_cell_measurement_temperature", measurement );
  const std::string head = ( "When converting CIF to NCMAT, the Uiso (msd)"
                             " values in the CIF data were" );
  auto quoted = []( const Field& f ) {
    return std::string( f.name ) + " (\"" + f.value + "\")";
  };

  if ( chosen.has_value() ) {
    const double T = chosen.value();
    if ( !( T >= tmin && T <= tmax ) )
      NCRYSTAL_THROW2( BadInput, "Invalid temperature chosen for the Uiso"
                       " (msd) values: " << fmt( T ) << "K (must be at least "
                       << fmt( tmin ) << "K and at most " << fmt( tmax )
                       << "K)." );
    std::string s = ( head + " explicitly chosen to be associated with a"
                      " temperature of " + dbl2shortstr( T ).to_string()
                      + "K" );
    if ( present.size() == 1 ) {
      s += ( ", overriding the " + std::string( present.front().name )
             + " CIF field (\"" + present.front().value + "\")" );
    } else if ( present.size() == 2 ) {
      const auto& a = present[0];
      const auto& b = present[1];
      const bool ok = ( a.number.has_value() && b.number.has_value()
                        && consistent( a.number.value(), b.number.value() ) );
      s += ( std::string( ", overriding the " ) + ( ok ? "" : "inconsistent " )
             + quoted( a ) + " and " + quoted( b ) + " CIF fields" );
    }
    return { T, s + ".", false };
  }

  VectS notes;
  bool isWarning = false;
  std::vector<const Field*> valid;
  for ( auto& f : present ) {
    if ( f.number.has_value() && f.number.value() >= tmin
         && f.number.value() <= tmax ) {
      valid.push_back( &f );
      continue;
    }
    if ( !permissive )
      NCRYSTAL_THROW2( BadInput, "The " << f.name << " CIF field (\""
                       << f.value << "\") is not a valid temperature. Specify"
                       " the temperature of the Uiso (msd) values explicitly,"
                       " or try permissive mode to ignore it." );
    warnfct( "The " + std::string( f.name ) + " CIF field (\"" + f.value
             + "\") is not a valid temperature and is ignored (accepted due"
             " to permissive mode)." );
    notes.push_back( "The " + std::string( f.name ) + " CIF field (\""
                     + f.value + "\") was not a valid temperature and was"
                     " ignored, which was only accepted because the"
                     " conversion was done in permissive mode." );
    isWarning = true;
  }
  auto joined = [&notes]( std::string s ) {
    for ( auto& n : notes )
      s += " " + n;
    return s;
  };
  if ( valid.empty() ) {
    const char * v = ( notes.empty() ? "" : "valid " );
    const std::string tf = dbl2shortstr( tfallback ).to_string();
    warnfct( std::string( "The CIF data specifies no " ) + v + "temperature,"
             " so its Uiso (msd) values are assumed to be associated with a"
             " fallback temperature of " + tf + "K. Specify the"
             " temperature explicitly, as this is likely wrong." );
    return { tfallback,
             joined( head + " assumed to be associated with a fallback"
                     " temperature of " + tf + "K, since the CIF"
                     " data specified no " + v + "temperature. This is likely"
                     " wrong." ),
             true };
  }
  const Field& f0 = *valid.front();
  const double T = f0.number.value();
  const std::string Ts = dbl2shortstr( T ).to_string();
  const std::string s = ( head + " assumed to be associated with a"
                          " temperature of " + Ts + "K, based on the "
                          + f0.name + " CIF field (\"" + f0.value + "\")." );
  for ( auto i : ncrange( std::size_t(1), valid.size() ) ) {
    const Field& f = *valid[i];
    if ( consistent( T, f.number.value() ) ) {
      notes.insert( notes.begin(), "The " + std::string( f.name )
                    + " CIF field (\"" + f.value + "\") was consistent with"
                    " this (less than 5% different)." );
      continue;
    }
    if ( !permissive )
      NCRYSTAL_THROW2( BadInput, "The " << quoted( f0 ) << " and "
                       << quoted( f ) << " CIF fields are inconsistent, so it"
                       " is unclear with which temperature the Uiso (msd)"
                       " values are associated. Specify the temperature"
                       " explicitly, or try permissive mode to use the "
                       << f0.name << " value." );
    warnfct( "The " + quoted( f0 ) + " and " + quoted( f ) + " CIF fields"
             " are inconsistent, so the Uiso (msd) values are assumed to be"
             " associated with a temperature of " + Ts + "K, based on"
             " the former (accepted due to permissive mode)." );
    notes.insert( notes.begin(), "The " + std::string( f.name ) + " CIF"
                  " field differed (\"" + f.value + "\"), an inconsistency"
                  " which was only accepted because the conversion was done"
                  " in permissive mode." );
    isWarning = true;
  }
  return { T, joined( s ), isWarning };
}
