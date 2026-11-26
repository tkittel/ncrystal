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

#include "NCrystal/internal/cif/NCCIFDescribe.hh"

namespace NC = NCrystal;
namespace NCC = NCrystal::CIF;

namespace NCRYSTAL_NAMESPACE {
  namespace CIF {
    namespace {

      //Trimmed text of an optional value (empty if absent or "?" or "."):
      std::string text( const Optional<Value>& v )
      {
        if ( !v.has_value() || isNull( v.value() ) )
          return {};
        return v.value().text.trimmed().to_string();
      }
      std::string text( const Optional<TaggedValue>& v )
      {
        return v.has_value() ? text( Optional<Value>( v.value().value ) )
          : std::string();
      }

      //Non-empty trimmed texts of values:
      VectS texts( const std::vector<Value>& v )
      {
        VectS res;
        for ( auto& e : v ) {
          auto t = text( Optional<Value>( e ) );
          if ( !t.empty() )
            res.push_back( std::move( t ) );
        }
        return res;
      }

      //Digits only, as number (0 if not):
      std::uint64_t idOf( const std::string& s )
      {
        if ( s.empty() || !StrView( s ).contains_only( "0123456789" ) )
          return 0;
        auto i = StrView( s ).toInt();
        return ( i.has_value() && i.value() > 0
                 ? static_cast<std::uint64_t>( i.value() ) : 0 );
      }
    }
  }
}

NCC::CIFDescription NCC::getDescription( const StructureItems& s,
                                         const Optional<ProvenanceItems>& prov,
                                         std::uint64_t requestedCodID,
                                         std::uint64_t requestedMPID,
                                         const CIFWarnFct& warnfct )
{
#ifndef NDEBUG
  nc_assert_always( bool(warnfct) );
#endif
  CIFDescription res;
  if ( s.formulaSum.has_value() && !isNull( s.formulaSum.value() ) )
    res.formula = s.formulaSum.value().text.to_string();
  const ProvenanceItems* p = ( prov.has_value() ? &prov.value() : nullptr );
  std::uint64_t codEmbedded = 0, mpEmbedded = 0;
  if ( p ) {
    if ( res.formula.empty() )
      res.formula = text( p->codOriginalFormulaSum );
    codEmbedded = idOf( text( p->codDatabaseCode ) );
    const auto mp = text( p->mpDatabaseCode );
    if ( StrView( mp ).startswith( "mp-" ) )
      mpEmbedded = idOf( mp.substr( 3 ) );
  }
  if ( codEmbedded && requestedCodID && codEmbedded != requestedCodID )
    warnfct( "Embedded CODID " + std::to_string( codEmbedded ) + " is"
             " different than requested " + std::to_string( requestedCodID )
             + "!" );
  if ( mpEmbedded && requestedMPID && mpEmbedded != requestedMPID )
    warnfct( "Embedded MPID " + std::to_string( mpEmbedded ) + " is"
             " different than requested " + std::to_string( requestedMPID )
             + "!" );
  res.codID = requestedCodID ? requestedCodID : codEmbedded;
  res.mpID = requestedMPID ? requestedMPID : mpEmbedded;

  VectS ll;
  const auto title = p ? text( p->title ) : std::string();
  const VectS authors = ( p && p->authors.has_value()
                          ? texts( p->authors.value().values ) : VectS() );
  const auto journal = p ? text( p->journal ) : std::string();
  const auto year = p ? text( p->year ) : std::string();
  const auto doi = p ? text( p->doi ) : std::string();
  std::string a;
  if ( authors.size() == 1 )
    a = authors.front();
  else if ( authors.size() == 2 )
    a = authors.front() + " and " + authors.back();
  else if ( !authors.empty() )
    a = authors.front() + ", et al.";
  if ( !title.empty() )
    ll.push_back( "\"" + title + "\"" );
  if ( !journal.empty() && !year.empty() && a.empty() )
    a = "<unknown authors>";
  const std::string jy = ( !journal.empty() && !year.empty()
                           ? journal + ", " + year
                           : ( !journal.empty() ? journal : year ) );
  if ( !a.empty() )
    ll.push_back( !jy.empty() ? a + " [" + jy + "]"
                  : ( authors.size() == 1 ? "Author: " : "Authors: " ) + a );
  if ( !doi.empty() )
    ll.push_back( "DOI: https://dx.doi.org/" + doi );
  const auto method = p ? text( p->auditCreationMethod ) : std::string();
  const auto date = p ? text( p->auditCreationDate ) : std::string();
  std::string auditAuthors;
  for ( auto& e : ( p ? texts( p->auditAuthors ) : VectS() ) )
    auditAuthors += ( auditAuthors.empty() ? "" : ", " ) + e;
  if ( !method.empty() )
    ll.push_back( "CIF creation method: " + method );
  if ( !date.empty() )
    ll.push_back( "CIF creation date: " + date );
  if ( !auditAuthors.empty() )
    ll.push_back( "CIF created by: " + auditAuthors );
  if ( res.codID ) {
    ll.push_back( "Crystallography Open Database entry "
                  + std::to_string( res.codID ) );
    ll.push_back( "https://www.crystallography.net/cod/"
                  + std::to_string( res.codID ) + ".html" );
  }
  if ( res.mpID ) {
    ll.push_back( "The Materials Project" );
    ll.push_back( "https://www.materialsproject.org/materials/mp-"
                  + std::to_string( res.mpID ) );
  }
  //Multi-line values split into lines (each right-trimmed):
  for ( auto& l : ll ) {
    std::size_t b = 0;
    while ( true ) {
      const auto e = l.find_first_of( "\r\n", b );
      res.lines.push_back(
        StrView( l ).substr( b, e == std::string::npos ? std::string::npos
                             : e - b ).rtrimmed().to_string() );
      if ( e == std::string::npos )
        break;
      b = e + ( ( l[e] == '\r' && e + 1 < l.size() && l[e+1] == '\n' )
                ? 2 : 1 );
    }
  }
  return res;
}
