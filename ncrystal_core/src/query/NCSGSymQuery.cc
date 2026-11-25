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

#include "NCSGSymQuery.hh"
#include "NCQueryImpl.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/sgsym/NCSGIdentify.hh"
#include "NCrystal/internal/cif/NCCIFUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

void NC::SGSymQuery::JSONQuery( std::ostream& os, const Query& rawquery )
{
  const Query query = QueryImpl::expandRawData( rawquery );
  auto invalid = [&query](const char * reason){
    std::ostringstream ss;
    ss << "Invalid sgsym query: ";
    streamJSON( ss, query );//trick: use streamJSON for easy format.
    if ( reason )
      ss<< " ("<<reason<<')';
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  constexpr auto sv_sgsym = StrView::make("sgsym");
  if ( query.size() < 2 || query.front() != sv_sgsym ) {
    invalid(nullptr);
    return;
  }
  //shift off the "sgsym" and key entries:
  const auto& key = query.at(1).trimmed();
  auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
  //Arguments which may start with "-" are given with a "@" prefix:
  auto argNoAt = [&arg]( std::size_t i ) {
    return arg(i).startswith('@') ? arg(i).substr(1) : arg(i);
  };
  const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_info = StrView::make("info");
  constexpr auto sv_sites = StrView::make("sites");
  constexpr auto sv_positions = StrView::make("positions");
  constexpr auto sv_fromops = StrView::make("fromops");
  constexpr auto sv_fromhall = StrView::make("fromhall");
  constexpr auto sv_fromhm = StrView::make("fromhm");
  constexpr auto sv_cell = StrView::make("cell");
  //Settings as a JSON list like ["14:b1","14:b2"]:
  auto streamSettings = [&os]( const std::vector<SpaceGroup>& v ) {
    os << '[';
    bool first = true;
    for ( auto& sg : v ) {
      if ( !first )
        os << ',';
      streamJSON( os, sg.toString() );
      first = false;
    }
    os << ']';
  };
  if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"sgsym\",\"list\"]");
    os<<"[\"info\", \"sites\", \"positions\", \"fromops\","
      " \"fromhall\", \"fromhm\", \"cell\"]";
  } else if ( key == sv_fromops ) {
    //["sgsym","fromops",OP,...]: the settings with exactly these symmetry
    //operations (like "x,-y,z+1/2", cf. SymOp, with a "@" prefix if starting
    //with "-"), as ["14:b1",...]:
    std::vector<SymOp> ops;
    ops.reserve( nargs );
    for ( auto i : ncrange( nargs ) )
      ops.emplace_back( argNoAt(i) );
    streamSettings( settingsFromOperations( ops ) );
  } else if ( key == sv_cell ) {
    //["sgsym","cell",TOL,DOTHROW,SG,a,b,c,alpha,beta,gamma]: the cell
    //parameters snapped to the constraints of the setting SG (parsed in strict
    //mode) with TOL "strict" or "lenient" (cf. SGCellConstraints::snap), as
    //{"cell":{"a":..,...,"gamma":..},"error_msg":null}. If they do not fit, a
    //BadInput is thrown if DOTHROW is 1, and if it is 0 the result is
    //{"cell":null,"error_msg":".."}. The values may be printed CIF numbers
    //(e.g. "4.123(2)"), whose rounding errors are then used with lenient
    //tolerances:
    const auto usage = "correct usage: [\"sgsym\",\"cell\","
      "\"strict|lenient\",0|1,SG,a,b,c,alpha,beta,gamma]";
    if ( nargs != 9 || ( arg(0) != "strict" && arg(0) != "lenient" )
         || ( arg(1) != "0" && arg(1) != "1" ) )
      invalid( usage );
    const bool lenient = ( arg(0) == "lenient" );
    const bool dothrow = ( arg(1) == "1" );
    double v[6], e[6];
    for ( auto i : ncrange( 6 ) ) {
      auto x = parseCIFNumber( arg( 3 + i ) );
      if ( !x.has_value() )
        invalid( usage );
      v[i] = x.value().value;
      e[i] = lenient ? x.value().roundingError : 0.0;
    }
    using Tol = SGCellConstraints::Tolerances;
    const auto rr = SGSymmetry::get( SpaceGroup( arg(2),
                                                 SpaceGroupParseMode::Strict ) )
      .cellConstraints()
      .trySnap( { v[0], v[1], v[2], v[3], v[4], v[5] },
                lenient ? Tol::lenient() : Tol::strict(),
                { e[0], e[1], e[2], e[3], e[4], e[5] } );
    if ( !rr.cell.has_value() ) {
      if ( dothrow )
        NCRYSTAL_THROW(BadInput,rr.errorMsg.value());
      os << "{\"cell\":null";
      streamJSONDictEntry( os, "error_msg", rr.errorMsg.value(),
                           JSONDictPos::LAST );
      return;
    }
    const auto& r = rr.cell.value();
    os << "{\"cell\":";
    streamJSONDictEntry( os, "a", r.a, JSONDictPos::FIRST );
    streamJSONDictEntry( os, "b", r.b );
    streamJSONDictEntry( os, "c", r.c );
    streamJSONDictEntry( os, "alpha", r.alpha );
    streamJSONDictEntry( os, "beta", r.beta );
    streamJSONDictEntry( os, "gamma", r.gamma, JSONDictPos::LAST );
    os << ",\"error_msg\":null}";
  } else if ( key == sv_fromhall || key == sv_fromhm ) {
    //["sgsym","fromhall",SYMBOL] and ["sgsym","fromhm",SYMBOL]: the settings
    //matching a Hall or Hermann-Mauguin symbol (cf. NCSGIdentify.hh, with a
    //"@" prefix if starting with "-"), as ["14:b1",...]:
    if ( nargs != 1 )
      invalid( key == sv_fromhall
               ? "correct usage: [\"sgsym\",\"fromhall\",SYMBOL]"
               : "correct usage: [\"sgsym\",\"fromhm\",SYMBOL]" );
    streamSettings( key == sv_fromhall ? settingsFromHallSymbol( argNoAt(0) )
                    : settingsFromHMSymbol( argNoAt(0) ) );
  } else if ( key == sv_info ) {
    //["sgsym","info",MODE,SG]: Symmetry operations of a space group
    //setting (SG as e.g. "225", "166:r", "227:2", "14:c1" or "62:cab",
    //parsed in MODE "strict" or "relaxed", cf. SpaceGroupParseMode), both
    //as CIF-like strings ("ops", incl. centring) and as the distinct
    //rotation parts (row-major integer matrices acting on fractional
    //coordinates):
    const auto usage
      = "correct usage: [\"sgsym\",\"info\",\"strict|relaxed\",SG]";
    if ( nargs != 2 || ( arg(0) != "strict" && arg(0) != "relaxed" ) )
      invalid( usage );
    const SpaceGroup sg( arg(1), ( arg(0) == "strict"
                                   ? SpaceGroupParseMode::Strict
                                   : SpaceGroupParseMode::Relaxed ) );
    const SGSymmetry& sym = SGSymmetry::get( sg );
    streamJSONDictEntry( os, "spacegroup", sg.toString(),
                         JSONDictPos::FIRST );
    streamJSONDictEntry( os, "settingidx",
                         static_cast<unsigned>( sg.settingIdx().get() ) );
    streamJSONDictEntry( os, "order", sym.order() );
    streamJSONDictEntry( os, "latticesymbol",
                         std::string( 1, sym.latticeSymbol() ) );
    streamJSONDictEntry( os, "centrosymmetric", sym.isCentrosymmetric() );
    os << ",\"ops\":[";
    bool first = true;
    for ( auto& op : sym.operations() ) {
      if ( !first )
        os << ',';
      streamJSON( os, op.toString() );
      first = false;
    }
    //Representatives have one operation per distinct rotation part:
    os << "],\"rotations\":[";
    first = true;
    for ( auto& op : sym.representatives() ) {
      os << ( first ? "[" : ",[" );
      for ( auto i : ncrange( 9u ) )
        os << ( i ? "," : "" ) << op.rot( i / 3, i % 3 );
      os << ']';
      first = false;
    }
    os << "]}";
  } else if ( key == sv_sites || key == sv_positions ) {
    //["sgsym","sites",MODE,SG,TOL,CELL,"lbl,x,y,z[,occupancy]",...]: group
    //unit cell positions (fractional) into crystallographic sites, one per
    //orbit of the space group, in the cell CELL ("a,b,c,alpha,beta,gamma"),
    //with the site tolerances TOL ("strict" or "lenient", cf.
    //SGSiteTolerances) and occupancies (default 1, the first given for a
    //label applies to all of its positions): {"spacegroup_setting":"227:2",
    //"is_default_setting":b,"sites":[[lbl,[x,y,z],multiplicity],...]}.
    //SG is parsed in MODE "strict" or "relaxed" (cf. SpaceGroupParseMode),
    //and in relaxed mode the setting is moreover deduced from the
    //positions among all settings of the number with the same axes (so a
    //bare number finds e.g. an origin choice 2 structure), while in strict
    //mode the positions must have the symmetry of exactly the given
    //setting (an error otherwise).
    //["sgsym","positions",MODE,SG,TOL,CELL,"lbl,x,y,z[,occupancy]",...]:
    //the inverse,
    //expanding sites in the setting SG into all unit cell positions:
    //[[lbl,[x,y,z]],...].
    const bool tosites = ( key == sv_sites );
    const auto usage = ( tosites
                         ? "correct usage: [\"sgsym\",\"sites\","
                         "\"strict|relaxed\",SG,\"strict|lenient\","
                         "\"a,b,c,alpha,beta,gamma\","
                         "\"lbl,x,y,z[,occupancy]\",...]"
                         : "correct usage: [\"sgsym\",\"positions\","
                         "\"strict|relaxed\",SG,\"strict|lenient\","
                         "\"a,b,c,alpha,beta,gamma\","
                         "\"lbl,x,y,z[,occupancy]\",...]" );
    if ( nargs < 5 || ( arg(0) != "strict" && arg(0) != "relaxed" )
         || ( arg(2) != "strict" && arg(2) != "lenient" ) )
      invalid( usage );
    const auto tol = ( arg(2) == "strict" ? SGSiteTolerances::strict()
                       : SGSiteTolerances::lenient() );
    CellParameters cell;
    {
      auto parts = arg(3).splitTrimmed(',');
      double v[6];
      if ( parts.size() != 6 )
        invalid( usage );
      for ( auto j : ncrange( 6 ) ) {
        auto x = parts.at(j).toDbl();
        if ( !x.has_value() || !( x.value() > 0.0 ) )
          invalid( usage );
        v[j] = x.value();
      }
      cell = { v[0], v[1], v[2], v[3], v[4], v[5] };
    }
    const bool strict = ( arg(0) == "strict" );
    const SpaceGroup given( arg(1), ( strict ? SpaceGroupParseMode::Strict
                                      : SpaceGroupParseMode::Relaxed ) );
    auto wrap = []( double x ) {
      double r = x - std::floor( x );
      return ( r >= 1.0 || r == 0.0 ) ? 0.0 : r;
    };
    //Positions grouped by label, in order of first appearance:
    std::vector<std::string> labels;
    std::vector<std::vector<Vector>> positions;
    std::vector<double> occupancies;//(per label)
    for ( auto i : ncrange( std::size_t(4), nargs ) ) {
      auto parts = arg(i).splitTrimmed(',');
      if ( parts.size() != 4 && parts.size() != 5 )
        invalid( usage );
      double occ = 1.0;
      if ( parts.size() == 5 ) {
        auto x = parts.at(4).toDbl();
        if ( !x.has_value() || !( x.value() > 0.0 && x.value() <= 1.0 ) )
          invalid( usage );
        occ = x.value();
      }
      double v[3];
      for ( auto j : ncrange( 3 ) ) {
        auto x = parts.at(j+1).toDbl();
        if ( !x.has_value() || !std::isfinite( x.value() ) )
          invalid( usage );
        v[j] = wrap( x.value() );
      }
      const std::string lbl = parts.at(0).to_string();
      std::size_t ilbl = 0;
      while ( ilbl < labels.size() && labels[ilbl] != lbl )
        ++ilbl;
      if ( ilbl == labels.size() ) {
        labels.push_back( lbl );
        positions.emplace_back();
        occupancies.push_back( occ );
      }
      positions[ilbl].emplace_back( v[0], v[1], v[2] );
    }
    auto streamPos = [&os]( const Vector& p ) {
      os << '[';
      streamJSON( os, p.x() );
      os << ',';
      streamJSON( os, p.y() );
      os << ',';
      streamJSON( os, p.z() );
      os << ']';
    };
    if ( tosites ) {
      Optional<SpaceGroup> setting;
      SGStructureCheck res;
      if ( strict ) {
        setting = given;
        res = checkStructureSymmetry( given, positions );
      } else {
        const bool rhombo = ( given.settingInfo().hexFamilyAxes
                              == SGHexFamilyAxes::Rhombohedral );
        res = checkStructureSymmetry( given.number(), rhombo, positions );
        if ( res.settingIdx.has_value() )
          setting.emplace( res.settingIdx.value() );
      }
      if ( !setting.has_value() || !( res.maxDeviation <= 1e-2 ) )
        NCRYSTAL_THROW2(BadInput,"Atom positions are not consistent with"
                        " space group "<<arg(1)<<" (largest deviation of"
                        " a symmetry image from an atom of the same type: "
                        <<fmt(res.maxDeviation,"%.2g")<<" in fractional"
                        " coordinates)");
      const SGSymmetry& sym = SGSymmetry::get( setting.value() );
      streamJSONDictEntry( os, "spacegroup_setting",
                           setting.value().toString(), JSONDictPos::FIRST );
      streamJSONDictEntry( os, "is_default_setting",
                           setting.value().isDefaultSetting() );
      os << ",\"sites\":[";
      bool first = true;
      for ( auto ilbl : ncrange( labels.size() ) ) {
        //(sorted, so the representative is the first of its orbit):
        auto& plist = positions[ilbl];
        std::sort( plist.begin(), plist.end(),
                   []( const Vector& a, const Vector& b )
                   { return std::make_tuple( a.x(), a.y(), a.z() )
                            < std::make_tuple( b.x(), b.y(), b.z() ); } );
        for ( auto& orbit : groupPositionsIntoOrbits( sym, cell, plist,
                                                      occupancies[ilbl],
                                                      tol ) ) {
          os << ( first ? "[" : ",[" );
          first = false;
          streamJSON( os, labels[ilbl] );
          os << ',';
          streamPos( vectAt( plist, orbit.front() ) );
          os << ',' << orbit.size() << ']';
        }
      }
      os << "]}";
    } else {
      const SGSymmetry& sym = SGSymmetry::get( given );
      os << '[';
      bool first = true;
      for ( auto ilbl : ncrange( labels.size() ) ) {
        for ( auto& site : positions[ilbl] ) {
          const auto orbit = expandSiteToOrbit( sym, cell, site,
                                                occupancies[ilbl], tol );
          for ( auto& p : orbit.positions ) {
            os << ( first ? "[" : ",[" );
            first = false;
            streamJSON( os, labels[ilbl] );
            os << ',';
            streamPos( p );
            os << ']';
          }
        }
      }
      os << ']';
    }
  } else {
    invalid("unknown key");
  }
}
