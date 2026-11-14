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
#include "NCrystal/internal/phys_utils/NCUijUtils.hh"
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
  constexpr auto sv_fromnumber = StrView::make("fromnumber");
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
      " \"fromhall\", \"fromhm\", \"fromnumber\", \"cell\"]";
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
    //["sgsym","cell","strict|lenient",SG,a,b,c,alpha,beta,gamma]: the cell
    //parameters snapped to the constraints of the setting SG (parsed in strict
    //mode), cf. SGCellConstraints::snap, as {"a":..,...,"gamma":..}. The
    //values may be printed CIF numbers (e.g. "4.123(2)"), whose rounding
    //errors are then used with lenient tolerances:
    const auto usage = "correct usage: [\"sgsym\",\"cell\","
      "\"strict|lenient\",SG,a,b,c,alpha,beta,gamma]";
    if ( nargs != 8 || ( arg(0) != "strict" && arg(0) != "lenient" ) )
      invalid( usage );
    const bool lenient = ( arg(0) == "lenient" );
    double v[6], e[6];
    for ( auto i : ncrange( 6 ) ) {
      auto x = parseCIFNumber( arg( 2 + i ) );
      if ( !x.has_value() )
        invalid( usage );
      v[i] = x.value().value;
      e[i] = lenient ? x.value().roundingError : 0.0;
    }
    using Tol = SGCellConstraints::Tolerances;
    const auto r = SGSymmetry::get( SpaceGroup( arg(1),
                                                SpaceGroupParseMode::Strict ) )
      .cellConstraints()
      .snap( { v[0], v[1], v[2], v[3], v[4], v[5] },
             lenient ? Tol::lenient() : Tol::strict(),
             { e[0], e[1], e[2], e[3], e[4], e[5] } );
    streamJSONDictEntry( os, "a", r.a, JSONDictPos::FIRST );
    streamJSONDictEntry( os, "b", r.b );
    streamJSONDictEntry( os, "c", r.c );
    streamJSONDictEntry( os, "alpha", r.alpha );
    streamJSONDictEntry( os, "beta", r.beta );
    streamJSONDictEntry( os, "gamma", r.gamma, JSONDictPos::LAST );
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
  } else if ( key == sv_fromnumber ) {
    //["sgsym","fromnumber",N]: all settings of the space group number N
    //(1-230), in table order, as ["227:1","227:2"]:
    auto n = ( nargs == 1 ? arg(0).toInt() : NullOpt );
    if ( !n.has_value() || n.value() < 1 || n.value() > 230 )
      invalid( "correct usage: [\"sgsym\",\"fromnumber\",N],"
               " with N in 1..230" );
    streamSettings( SpaceGroup::settingsOfNumber(
                      static_cast<unsigned>( n.value() ) ) );
  } else if ( key == sv_info ) {
    //["sgsym","info",MODE,SG]: H-M symbol ("hm", ITVB computer-entry
    //format) and symmetry operations of a space group setting (SG as e.g. "225", "166:r", "227:2", "14:c1" or "62:cab",
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
    streamJSONDictEntry( os, "hm", sg.hmSymbol() );
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
    //["sgsym","positions",MODE,SG,TOL,CELL,ENTRY,...]: expand sites into all
    //their positions in the setting SG (parsed in MODE "strict" or
    //"relaxed", cf. SpaceGroupParseMode), in the cell CELL
    //("a,b,c,alpha,beta,gamma"), with the site tolerances TOL ("strict" or
    //"lenient", cf. SGSiteTolerances). Each ENTRY is a site,
    //"lbl,x,y,z[,occupancy][,u11,u22,u33,u12,u13,u23]" (occupancy default 1,
    //the first given for a label applies to all of its entries; uij as Ustd,
    //cf. NCUijUtils.hh). Numbers may be printed CIF numbers (e.g.
    //"0.0123(4)"), whose rounding errors are used with lenient tolerances
    //(angles printed as 90 or 120 without uncertainty are exact).
    //Result: {"spacegroup_setting":"227:2","sites":[{"label":..,
    //"multiplicity":n,"split":null|{"positions":k,"distance":d},
    //"uij":null|[6 values, snapped to the site symmetry],"ueq":null|[ueq,
    //rounding error],"positions":[[x,y,z],...],"uijs":null|[[6 values],...]},
    //...]}, where "split" describes a split site (cf. expandSiteToOrbit),
    //"ueq" is that of the given uij (cf. calcUeqWithRoundingError) and "uijs"
    //are the tensors at the positions.
    //["sgsym","sites",MODE,SG,TOL,CELL,ENTRY,...]: the inverse, grouping unit
    //cell positions (entries as above, with uij at each position) into sites,
    //one per orbit: {"spacegroup_setting":"227:2","is_default_setting":b,
    //"sites":[{"label":..,"position":[x,y,z],"multiplicity":n,"split":..,
    //"uij":..},...]}, with the representative position (the first of the
    //orbit, sorted) and its uij (snapped, and the uij at the other positions
    //must be its symmetry images). In relaxed mode the setting is moreover
    //deduced from the positions among all settings of the number with the
    //same axes (so a bare number finds e.g. an origin choice 2 structure),
    //while in strict mode the positions must have the symmetry of exactly the
    //given setting (an error otherwise).
    const bool tosites = ( key == sv_sites );
    const auto usage = ( tosites
                         ? "correct usage: [\"sgsym\",\"sites\","
                         "\"strict|relaxed\",SG,\"strict|lenient\","
                         "\"a,b,c,alpha,beta,gamma\","
                         "\"lbl,x,y,z[,occupancy][,u11,u22,u33,u12,u13,u23]\","
                         "...]"
                         : "correct usage: [\"sgsym\",\"positions\","
                         "\"strict|relaxed\",SG,\"strict|lenient\","
                         "\"a,b,c,alpha,beta,gamma\","
                         "\"lbl,x,y,z[,occupancy][,u11,u22,u33,u12,u13,u23]\","
                         "...]" );
    if ( nargs < 5 || ( arg(0) != "strict" && arg(0) != "relaxed" )
         || ( arg(2) != "strict" && arg(2) != "lenient" ) )
      invalid( usage );
    const bool strict = ( arg(0) == "strict" );
    const bool lenient = ( arg(2) == "lenient" );
    const SpaceGroup given( arg(1), ( strict ? SpaceGroupParseMode::Strict
                                      : SpaceGroupParseMode::Relaxed ) );
    const auto tol = ( lenient ? SGSiteTolerances::lenient()
                       : SGSiteTolerances::strict() );
    auto number = [&invalid,&usage]( StrView s ) {
      auto x = parseCIFNumber( s );
      if ( !x.has_value() )
        invalid( usage );
      return x.value();
    };
    CellParameters cell, cellerr{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    {
      auto parts = arg(3).splitTrimmed(',');
      if ( parts.size() != 6 )
        invalid( usage );
      double v[6], e[6];
      for ( auto j : ncrange( 6 ) ) {
        const auto x = number( parts.at(j) );
        v[j] = x.value;
        //(angles printed as 90 or 120 without uncertainty are exact):
        const bool exact = ( j >= 3 && ( v[j] == 90.0 || v[j] == 120.0 )
                             && !parts.at(j).contains( '(' ) );
        e[j] = ( lenient && !exact ) ? x.roundingError : 0.0;
        if ( !( v[j] > 0.0 ) )
          invalid( usage );
      }
      cell = { v[0], v[1], v[2], v[3], v[4], v[5] };
      cellerr = { e[0], e[1], e[2], e[3], e[4], e[5] };
    }
    const auto metric = UijUtils::cellMetric( cell.a, cell.b, cell.c,
                                              cell.alpha, cell.beta,
                                              cell.gamma );
    auto wrap = []( double x ) {
      double r = x - std::floor( x );
      return ( r >= 1.0 || r == 0.0 ) ? 0.0 : r;
    };
    //Entries grouped by label, in order of first appearance:
    struct Entry {
      Vector pos;
      Optional<UijUtils::Uij> uij;
      UijUtils::Uij uijerr = UijUtils::Uij{};
    };
    std::vector<std::string> labels;
    std::vector<std::vector<Entry>> entries;
    std::vector<double> occupancies;//(per label)
    for ( auto i : ncrange( std::size_t(4), nargs ) ) {
      auto parts = arg(i).splitTrimmed(',');
      const std::size_t np = parts.size();
      if ( np != 4 && np != 5 && np != 10 && np != 11 )
        invalid( usage );
      Entry e;
      double v[3];
      for ( auto j : ncrange( 3 ) ) {
        v[j] = number( parts.at(j+1) ).value;
        if ( !std::isfinite( v[j] ) )
          invalid( usage );
        v[j] = wrap( v[j] );
      }
      e.pos = Vector( v[0], v[1], v[2] );
      double occ = 1.0;
      if ( np == 5 || np == 11 ) {
        occ = number( parts.at(4) ).value;
        if ( !( occ > 0.0 && occ <= 1.0 ) )
          invalid( usage );
      }
      if ( np >= 10 ) {
        UijUtils::Uij u;
        for ( auto j : ncrange( 6 ) ) {
          const auto x = number( parts.at( np - 6 + j ) );
          u[j] = x.value;
          if ( lenient )
            e.uijerr[j] = x.roundingError;
        }
        e.uij = u;
      }
      const std::string lbl = parts.at(0).to_string();
      std::size_t ilbl = 0;
      while ( ilbl < labels.size() && labels[ilbl] != lbl )
        ++ilbl;
      if ( ilbl == labels.size() ) {
        labels.push_back( lbl );
        entries.emplace_back();
        occupancies.push_back( occ );
      }
      if ( !entries[ilbl].empty()
           && entries[ilbl].front().uij.has_value() != e.uij.has_value() )
        NCRYSTAL_THROW2(BadInput,"Anisotropic displacements (uij) must be"
                        " given for all or none of the entries of "<<lbl);
      entries[ilbl].push_back( std::move( e ) );
    }
    auto streamVals = [&os]( const double* v, std::size_t n ) {
      os << '[';
      for ( auto j : ncrange( n ) ) {
        if ( j )
          os << ',';
        streamJSON( os, v[j] );
      }
      os << ']';
    };
    auto streamPos = [&streamVals]( const Vector& p ) {
      const double v[3] = { p.x(), p.y(), p.z() };
      streamVals( v, 3 );
    };
    auto streamUij = [&os,&streamVals]( const Optional<UijUtils::Uij>& u ) {
      if ( u.has_value() )
        streamVals( u.value().data(), 6 );
      else
        os << "null";
    };
    auto streamSplit = [&os]( const SGSiteOrbit& orbit ) {
      if ( !orbit.nearCopies ) {
        os << "null";
        return;
      }
      streamJSONDictEntry( os, "positions", orbit.nearCopies + 1,
                           JSONDictPos::FIRST );
      streamJSONDictEntry( os, "distance", orbit.nearestCopyDistance,
                           JSONDictPos::LAST );
    };
    auto maxErr = []( const std::vector<Entry>& v ) {
      UijUtils::Uij m = UijUtils::Uij{};
      for ( auto& e : v )
        for ( auto j : ncrange( 6 ) )
          m[j] = ncmax( m[j], e.uijerr[j] );
      return m;
    };
    auto descr = [&labels]( std::size_t ilbl, const Vector& p ) {
      std::ostringstream ss;
      ss << labels[ilbl] << " (" << fmt( p.x(), "%.6g" ) << ", "
         << fmt( p.y(), "%.6g" ) << ", " << fmt( p.z(), "%.6g" ) << ")";
      return ss.str();
    };
    if ( tosites ) {
      std::vector<std::vector<Vector>> positions;
      for ( auto& v : entries ) {
        positions.emplace_back();
        for ( auto& e : v )
          positions.back().push_back( e.pos );
      }
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
        auto& elist = entries[ilbl];
        std::stable_sort( elist.begin(), elist.end(),
                          []( const Entry& a, const Entry& b )
                          { return std::make_tuple( a.pos.x(), a.pos.y(),
                                                    a.pos.z() )
                              < std::make_tuple( b.pos.x(), b.pos.y(),
                                                 b.pos.z() ); } );
        std::vector<Vector> plist;
        for ( auto& e : elist )
          plist.push_back( e.pos );
        for ( auto& orbit : groupPositionsIntoOrbits( sym, cell, plist,
                                                      occupancies[ilbl],
                                                      tol ) ) {
          const Vector& rep = vectAt( plist, orbit.front() );
          const auto expanded = expandSiteToOrbit( sym, cell, rep,
                                                   occupancies[ilbl], tol );
          Optional<UijUtils::Uij> uij;
          if ( vectAt( elist, orbit.front() ).uij.has_value() ) {
            std::vector<Vector> opos;
            std::vector<UijUtils::Uij> ou;
            std::vector<Entry> oe;
            for ( auto i : orbit ) {
              opos.push_back( vectAt( plist, i ) );
              ou.push_back( vectAt( elist, i ).uij.value() );
              oe.push_back( vectAt( elist, i ) );
            }
            uij = UijUtils::ustdOfOrbit( sym, metric, opos, ou,
                                         descr( ilbl, rep ), 0.05,
                                         maxErr( oe ) );
          }
          if ( !first )
            os << ',';
          first = false;
          streamJSONDictEntry( os, "label", labels[ilbl], JSONDictPos::FIRST );
          os << ",\"position\":";
          streamPos( rep );
          streamJSONDictEntry( os, "multiplicity",
                               static_cast<unsigned>( orbit.size() ) );
          os << ",\"split\":";
          streamSplit( expanded );
          os << ",\"uij\":";
          streamUij( uij );
          os << '}';
        }
      }
      os << "]}";
    } else {
      const SGSymmetry& sym = SGSymmetry::get( given );
      streamJSONDictEntry( os, "spacegroup_setting", given.toString(),
                           JSONDictPos::FIRST );
      os << ",\"sites\":[";
      bool first = true;
      for ( auto ilbl : ncrange( labels.size() ) ) {
        for ( auto& e : entries[ilbl] ) {
          const auto orbit = expandSiteToOrbit( sym, cell, e.pos,
                                                occupancies[ilbl], tol );
          Optional<UijUtils::Uij> uij;
          if ( e.uij.has_value() )
            uij = UijUtils::snapUstdToSiteSymmetry( sym, metric,
                                                    orbit.positions.front(),
                                                    e.uij.value(),
                                                    descr( ilbl, e.pos ),
                                                    0.05, e.uijerr );
          if ( !first )
            os << ',';
          first = false;
          streamJSONDictEntry( os, "label", labels[ilbl], JSONDictPos::FIRST );
          streamJSONDictEntry( os, "multiplicity",
                               static_cast<unsigned>( orbit.positions.size() ) );
          os << ",\"split\":";
          streamSplit( orbit );
          os << ",\"uij\":";
          streamUij( uij );
          os << ",\"ueq\":";
          if ( e.uij.has_value() ) {
            const auto ueq
              = UijUtils::calcUeqWithRoundingError( e.uij.value(), cell,
                                                    e.uijerr, cellerr );
            const double v[2] = { ueq.first, ueq.second };
            streamVals( v, 2 );
          } else {
            os << "null";
          }
          os << ",\"positions\":[";
          for ( auto& p : orbit.positions ) {
            if ( &p != &orbit.positions.front() )
              os << ',';
            streamPos( p );
          }
          os << "],\"uijs\":";
          if ( uij.has_value() ) {
            os << '[';
            for ( auto& p : orbit.positions ) {
              if ( &p != &orbit.positions.front() )
                os << ',';
              const auto u = UijUtils::mapUstdToPosition( sym, metric,
                                                           orbit.positions.front(),
                                                           p, uij.value() );
              streamVals( u.data(), 6 );
            }
            os << ']';
          } else {
            os << "null";
          }
          os << '}';
        }
      }
      os << "]}";
    }
  } else {
    invalid("unknown key");
  }
}
