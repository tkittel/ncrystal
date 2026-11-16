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

#include "NCrystal/internal/sabscatter/NCVDOSDebyeGroupScatter.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

NC::VDOSDebyeGroupScatter::VDOSDebyeGroupScatter
( shared_obj<const DynInfoUtils::VDOSDebyeGroupData> data, SigmaBound sb )
  : m_data( std::move( data ) ), m_sb( sb )
{
  if ( !( sb.dbl() > 0.0 && std::isfinite( sb.dbl() ) ) )
    NCRYSTAL_THROW2( BadInput, "VDOSDebyeGroupScatter: invalid SigmaBound"
                     " (must be positive and finite): " << sb );
}

NC::CrossSect
NC::VDOSDebyeGroupScatter::crossSectionIsotropic( CachePtr&,
                                                  NeutronEnergy ekin ) const
{
  return CrossSect{ m_data->xs().evalScaled( ekin.dbl(), 1.0, m_sb.dbl() ) };
}

NC::ScatterOutcomeIsotropic
NC::VDOSDebyeGroupScatter::sampleScatterIsotropic( CachePtr&, RNG& rng,
                                                   NeutronEnergy ekin ) const
{
  return m_data->sampleScatter( rng, ekin );
}

NC::Optional<NC::IsoXSTable>
NC::VDOSDebyeGroupScatter::createIsoXSTable() const
{
  return scaleIsoXSTable( m_data->xs().table(), 1.0, m_sb.dbl() );
}

NC::Optional<std::string>
NC::VDOSDebyeGroupScatter::specificJSONDescription() const
{
  const auto& d = *m_data;
  const auto sigfree = m_sb.free( d.mass() );
  std::ostringstream ss;
  {
    std::ostringstream tmp;
    tmp << "n=" << d.nEntries() << ";nclusters=" << d.clusters().size()
        << ";TDebye=" << DebyeTemperature{ d.debyeTempRange().first };
    if ( d.debyeTempRange().second != d.debyeTempRange().first )
      tmp << ".." << DebyeTemperature{ d.debyeTempRange().second };
    tmp << ";T=" << d.temperature() << ";M=" << d.mass()
        << ";sigma_free=" << sigfree;
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "n", d.nEntries() );
  streamJSONDictEntry( ss, "nclusters", d.clusters().size() );
  streamJSONDictEntry( ss, "TDebye_min", d.debyeTempRange().first );
  streamJSONDictEntry( ss, "TDebye_max", d.debyeTempRange().second );
  streamJSONDictEntry( ss, "T", d.temperature().dbl() );
  streamJSONDictEntry( ss, "M", d.mass().dbl() );
  streamJSONDictEntry( ss, "sigma_bound", m_sb.dbl() );
  streamJSONDictEntry( ss, "sigma_free", sigfree.dbl() );
  std::vector<std::string> clusters;
  for ( auto& s : d.clusters() ) {
    const auto& p = s.sabext->processor();
    const auto& sab = *p.sabDataPtr();
    std::ostringstream tmp;
    tmp << "n=" << s.n << ";scale=" << fmt( s.scale )
        << ";TDebye=" << s.debyeTemperature
        << ";nalpha=" << sab.alphaGrid().size()
        << ";nbeta=" << sab.betaGrid().size()
        << ";Emax=" << p.getEMaxInfo().ekin;
    clusters.push_back( tmp.str() );
  }
  streamJSONDictEntry( ss, "clusters", clusters, JSONDictPos::LAST );
  return ss.str();
}
