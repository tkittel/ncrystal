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

#include "NCrystal/internal/vdos/NCVDOSLux.hh"

namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace VDOS {
    namespace {
      //Bits: 0-6 lvl, 7 legacy, 8 no grouping (raw 100+lvl):
      constexpr unsigned noGroupingBit = 256u;

      int vdosluxEncodeRaw( unsigned lvl, bool isLegacy, bool noGrouping )
      {
        int v = static_cast<int>(lvl);
        if ( isLegacy )
          v += 1000;
        if ( noGrouping )
          v += 100;
        return v;
      }

      unsigned vdosluxDecodeRaw( int rawval )
      {
        //0 + val : nextgen
        //100 + val : nextgen, no grouping of VDOSDebye entries
        //1000 + val : legacy
        //2000 + val : nextgen
        if ( rawval >= 0 && rawval <= 6 )
          rawval += 2000;
        unsigned res;
        if ( rawval >= 100 && rawval <= 106 ) {
          res = static_cast<unsigned>(rawval-100) | noGroupingBit;
        } else if ( rawval >= 2000 && rawval <= 2006 ) {
          res = static_cast<unsigned>(rawval-2000);
          nc_assert_always( res <= 127 );
        } else if ( rawval >= 1000 && rawval <= 1005 ) {
          res = static_cast<unsigned>(rawval-1000);
          nc_assert_always( res <= 127u );
          res |= 128u;
        } else {
          NCRYSTAL_THROW2(BadInput,"Invalid vdoslux value: "<<rawval);
        }
        return res;
      }
    }
  }
}

NC::VDOS::VDOSLux::VDOSLux( int rawval )
  : m_v(vdosluxDecodeRaw(rawval))
{
  nc_assert( lvl() <= (isLegacy()?5u:6u) );
  nc_assert( raw() == rawval || m_v == vdosluxDecodeRaw(raw()) );
}

//NB: Keep this carefully synchronised with the vdoslux implementation in
//NCCfgVars.hh.

NC::VDOS::VDOSLux::VDOSLux()
  : m_v( 3u )//nextgen 3
{
  nc_assert( *this == VDOSLux(3) );
}

int NC::VDOS::VDOSLux::raw() const
{
  return vdosluxEncodeRaw( lvl(),
                           isLegacy(),
                           static_cast<bool>( m_v & noGroupingBit ) );
}

NC::VDOS::VDOSLux NC::VDOS::VDOSLux::reduceLvl( unsigned v ) const
{
  return VDOSLux( vdosluxEncodeRaw( ( lvl() > v ? lvl() - v : 0u ),
                                    isLegacy(),
                                    static_cast<bool>( m_v & noGroupingBit ) ) );
}

NC::VDOS::VDOSLux NC::VDOS::VDOSLux::expansionLux() const
{
  return VDOSLux( vdosluxEncodeRaw( lvl(), isLegacy(), false ) );
}

