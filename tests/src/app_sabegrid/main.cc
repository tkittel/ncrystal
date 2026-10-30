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


#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/sab/NCSABCfg.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

int main()
{
  //Phase-space integrals of the adaptive energy grid (knllux=3) against a
  //dense geometric grid (2500 points) over the same range, for a difficult
  //(cryogenic) case:
  const char * cfgstr = "Al_sg225.ncmat;temp=20K";
  auto info = NC::FactImpl::createInfo( NC::MatCfg( cfgstr ) );
  nc_assert_always( info->getDynamicInfoList().size() == 1 );
  auto di = dynamic_cast<const NC::DI_VDOS*>( info->getDynamicInfoList().front().get() );
  nc_assert_always( di != nullptr );
  auto sab = NC::extractSABDataFromDynInfo( di );
  const auto cfg = NC::SABCfg::createConfig( 3 );
  NCS::SABProcessor p( cfg, sab, nullptr, NCS::SABProcessor::SampleSupport::NO );
  const auto& x = p.getEDivKTGrid();
  const double kT = p.kT();
  auto dense = std::make_shared<NC::VectD>( NC::geomspace( x.front() * kT, x.back() * kT, 2500 ) );
  dense->front() = x.front() * kT;
  dense->back() = x.back() * kT;
  NCS::SABProcessor pd( cfg, sab, dense, NCS::SABProcessor::SampleSupport::NO );
  double maxdev = 0.0;
  for ( auto t : NC::geomspace( x.front(), x.back(), 20001 ) ) {
    const NC::NeutronEnergy e{ t * kT };
    const double ref = pd.phaseSpaceIntegral( e );
    maxdev = std::max( maxdev, std::fabs( p.phaseSpaceIntegral( e ) / ref - 1.0 ) );
  }
  const bool ok = maxdev < cfg.egrid_tolerance;
  std::cout << cfgstr << ": adaptive egrid " << ( ok ? "within " : "NOT within " )
            << cfg.egrid_tolerance << " of a dense egrid" << std::endl;
  nc_assert_always( ok );
  return 0;
}
