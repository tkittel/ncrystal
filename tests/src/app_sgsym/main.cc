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

// Tests of the sgsym component (space groups and their symmetries).

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include <iostream>
#include <sstream>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    //fwd declare internal fct:
    std::vector<SymOp> rawSymOpsFromHallSymbol( StrView );
  }
}

////////////////////////////////////////////////////////////////////////////////
// Space group table
////////////////////////////////////////////////////////////////////////////////

//Tests the SpaceGroup and SpaceGroupSettingIdx classes. The complete table of
//space group settings is printed, so the reference log also guards against
//accidental changes to it.

namespace test_table {

namespace {
  void expectBad( const char * s, NC::SpaceGroupParseMode mode
                  = NC::SpaceGroupParseMode::Relaxed )
  {
    try {
      NC::SpaceGroup sg( s, mode );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  \"" << s << "\" -> BadInput: " << e.what() << std::endl;
      return;
    }
    nc_assert_always( false );
  }
  void expectBadHall( std::uint16_t hn )
  {
    try {
      NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  setting index " << hn << " -> BadInput: " << e.what()
                << std::endl;
      return;
    }
    nc_assert_always( false );
  }
}

void run()
{
  static_assert( sizeof( NC::SpaceGroup ) == 2, "" );
  static_assert( std::is_trivially_copyable<NC::SpaceGroup>::value, "" );
  //Constructors are explicit (no implicit conversions):
  static_assert( !std::is_convertible<const char*,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<std::string,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<NC::StrView,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<NC::SpaceGroupSettingIdx,
                                      NC::SpaceGroup>::value, "" );
  static_assert( std::is_constructible<NC::SpaceGroup,std::string&&,
                                       NC::SpaceGroupParseMode>::value, "" );
  static_assert( !std::is_constructible<NC::SpaceGroup,std::string&&>::value,
                 "" );//(the parse mode is never implicit)

  std::cout << "All settings (setting index, setting, Hermann-Mauguin symbol,"
               " Hall symbol, flags):" << std::endl;
  std::set<unsigned> numbers, numbers_multi;
  unsigned last_number = 0;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
    REQUIRE( sg.settingIdx().get() == hn );
    REQUIRE( sg.number() >= 1 && sg.number() <= 230 );
    REQUIRE( sg.number() >= last_number );//sorted
    REQUIRE( sg.isDefaultSetting() == ( sg.number() != last_number ) );
    last_number = sg.number();
    numbers.insert( sg.number() );
    if ( sg.hasMultipleSettings() )
      numbers_multi.insert( sg.number() );
    REQUIRE( *sg.hallSymbol() );
    REQUIRE( *sg.hmSymbol() );
    //Empty choice codes exactly for numbers with a single setting:
    REQUIRE( sg.choice() != nullptr );
    REQUIRE( ( *sg.choice() != 0 ) == sg.hasMultipleSettings() );
    //String round trip and stream output:
    const std::string str = sg.toString();
    constexpr auto relaxed = NC::SpaceGroupParseMode::Relaxed;
    constexpr auto strict = NC::SpaceGroupParseMode::Strict;
    REQUIRE( NC::SpaceGroup( str, relaxed ) == sg );
    REQUIRE( NC::SpaceGroup( str, strict ) == sg );
    REQUIRE( NC::SpaceGroup( std::string( str ), relaxed ) == sg );//temporary
    REQUIRE( NC::SpaceGroup( NC::StrView( str ), relaxed ) == sg );
    REQUIRE( NC::SpaceGroup( str.c_str(), relaxed ) == sg );
    //setting index -> SpaceGroup -> string -> SpaceGroup -> setting index:
    REQUIRE( NC::SpaceGroup( str, relaxed ).settingIdx().get() == hn );
    std::ostringstream ss;
    ss << sg;
    REQUIRE( ss.str() == str );
    //Bare number gives default setting:
    const std::string nstr = std::to_string( sg.number() );
    REQUIRE( ( NC::SpaceGroup( nstr, relaxed ) == sg ) == sg.isDefaultSetting() );
    REQUIRE( ( NC::SpaceGroup( nstr, relaxed ) < sg ) == !sg.isDefaultSetting() );
    //Strictly, only for numbers with a single setting:
    if ( sg.hasMultipleSettings() ) {
      try {
        NC::SpaceGroup sg2( nstr, strict );
        REQUIRE( false );
      } catch ( NC::Error::BadInput& ) {
      }
    } else {
      REQUIRE( NC::SpaceGroup( nstr, strict ) == sg );
    }
    const auto settings = NC::SpaceGroup::settingsOfNumber( sg.number() );
    REQUIRE( settings.front().isDefaultSetting() );
    REQUIRE( ( settings.size() > 1 ) == sg.hasMultipleSettings() );
    REQUIRE( std::find( settings.begin(), settings.end(), sg ) != settings.end() );
    std::cout << "  " << hn << " " << sg << " \"" << sg.hmSymbol() << "\" \""
              << sg.hallSymbol() << "\""
              << ( sg.isDefaultSetting() ? " [default]" : "" )
              << ( sg.hasMultipleSettings() ? "" : " [single]" ) << std::endl;
  }
  REQUIRE( numbers.size() == 230 );
  std::cout << "Numbers with multiple settings: " << numbers_multi.size()
            << std::endl;
  REQUIRE( numbers_multi.size() == 90 );

  std::cout << "Examples:" << std::endl;
  for ( auto s : { "227", "227:1", "227:2", "166", "166:h", "166:r", "14",
                   "14:b1", "14:c1", "62", "62:abc", "62:cab", "68:1cab", "70:2",
                   "15:-b2", "1", "225" } ) {
    const NC::SpaceGroup sg( s, NC::SpaceGroupParseMode::Relaxed );
    std::cout << "  \"" << s << "\" -> " << sg << " (setting index "
              << sg.settingIdx() << ", Hall symbol \"" << sg.hallSymbol()
              << "\")" << std::endl;
  }
  REQUIRE( NC::SpaceGroup( "227", NC::SpaceGroupParseMode::Relaxed ).settingIdx().get() == 525 );
  REQUIRE( NC::SpaceGroup( "227:2", NC::SpaceGroupParseMode::Strict ).settingIdx().get() == 526 );
  std::cout << "Settings of 227: " << NC::SpaceGroup::listSettings( 227 )
            << "; of 62: " << NC::SpaceGroup::listSettings( 62 )
            << "; of 225: " << NC::SpaceGroup::listSettings( 225 ) << std::endl;

  std::cout << "Invalid input (relaxed parsing):" << std::endl;
  for ( auto s : { "", "0", "231", "999", "1000", "-1", "+1", " 227", "227 ",
                   "abc", "227:", "227:3", "227:1:2", "166:H", "62:bac",
                   "225:1", "14:b", "3:b1", ":2", "22a" } )
    expectBad( s );
  std::cout << "Invalid input (strict parsing):" << std::endl;
  for ( auto s : { "227", "166", "14", "62", "227:3", "167:H", "225:1" } )
    expectBad( s, NC::SpaceGroupParseMode::Strict );
  expectBadHall( 0 );
  expectBadHall( 531 );
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Symmetry operations
////////////////////////////////////////////////////////////////////////////////

//Tests the SymOp class (symmetry operations).

namespace test_symop {

namespace {

  using Mat = std::array<int,9>;

  Mat matMult( const Mat& a, const Mat& b )
  {
    Mat q;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        q[3*i+j] = a[3*i]*b[j] + a[3*i+1]*b[3+j] + a[3*i+2]*b[6+j];
    return q;
  }

  //Independent construction of the 64 allowed rotations: the 48 signed
  //permutation matrices, and the 24 rotations of 6/mmm in hexagonal axes
  //(generated here by 6-fold rotations about z, 2-fold rotations about the
  //a-axis, and inversion):
  std::set<Mat> allowedRotations()
  {
    std::set<Mat> res;
    const int perms[6][3] = { {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1},
                              {2,1,0} };
    for ( auto& pm : perms ) {
      for ( int signs = 0; signs < 8; ++signs ) {
        Mat m{};
        for ( int i = 0; i < 3; ++i )
          m[3*i+pm[i]] = ( signs & ( 1 << i ) ) ? -1 : 1;
        res.insert( m );
      }
    }
    std::set<Mat> hex{ Mat{ { 1,0,0, 0,1,0, 0,0,1 } } };
    const Mat six{ { 1,-1,0, 1,0,0, 0,0,1 } };
    const Mat twofold_a{ { 1,-1,0, 0,-1,0, 0,0,-1 } };
    const Mat inv{ { -1,0,0, 0,-1,0, 0,0,-1 } };
    for ( bool changed = true; changed; ) {
      changed = false;
      for ( auto m : std::vector<Mat>( hex.begin(), hex.end() ) )
        for ( auto& g : { six, twofold_a, inv } )
          if ( hex.insert( matMult( m, g ) ).second )
            changed = true;
    }
    nc_assert_always( res.size() == 48 && hex.size() == 24 );
    res.insert( hex.begin(), hex.end() );
    return res;
  }

  bool entriesOK( const NC::SymOp& op )
  {
    for ( unsigned i = 0; i < 3; ++i ) {
      for ( unsigned j = 0; j < 3; ++j )
        if ( !( op.rot(i,j) >= -1 && op.rot(i,j) <= 1 ) )
          return false;
      if ( !( op.trans(i) >= 0 && op.trans(i) < NC::SymOp::tdenom ) )
        return false;
    }
    return true;
  }

  //Group generated by the given operations (only for testing):
  std::vector<NC::SymOp> closure( std::vector<NC::SymOp> gens )
  {
    std::set<NC::SymOp> g{ NC::SymOp() };
    for ( bool changed = true; changed; ) {
      changed = false;
      std::vector<NC::SymOp> cur( g.begin(), g.end() );
      for ( auto& a : cur )
        for ( auto& b : gens )
          if ( g.insert( a * b ).second )
            changed = true;
    }
    return { g.begin(), g.end() };
  }

  void expectBad( const char * s )
  {
    try {
      NC::SymOp op( s );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  \"" << s << "\" -> BadInput: " << e.what() << std::endl;
      return;
    }
    nc_assert_always( false );
  }
}

void run()
{
  static_assert( sizeof( NC::SymOp ) == 12, "" );
  static_assert( std::is_trivially_copyable<NC::SymOp>::value, "" );
  static_assert( !std::is_convertible<const char*,NC::SymOp>::value, "" );
  static_assert( !std::is_convertible<std::string,NC::SymOp>::value, "" );
  static_assert( !std::is_convertible<NC::StrView,NC::SymOp>::value, "" );
  static_assert( std::is_constructible<NC::SymOp,std::string&&>::value, "" );

  REQUIRE( NC::SymOp().isIdentity() );
  REQUIRE( NC::SymOp().toString() == "x,y,z" );

  //All 3^9 matrices with entries in {-1,0,1}:
  const std::array<std::array<int,3>,5> translations = { {
      { { 0, 0, 0 } }, { { 12, 0, 0 } }, { { 1, 5, 23 } }, { { 8, 16, 3 } },
      { { -3, 33, 48 } } } };
  const std::set<Mat> allowed = allowedRotations();
  REQUIRE( allowed.size() == 64 );
  unsigned nvalid = 0;
  std::map<int,unsigned> ntypes;
  for ( int code = 0; code < 19683; ++code ) {
    Mat m;
    for ( int i = 0, c = code; i < 9; ++i, c /= 3 )
      m[i] = c % 3 - 1;
    const bool valid = allowed.count( m ) > 0;
    REQUIRE( valid == NC::SymOp::isAllowedRotation( m ) );
    bool accepted = true;
    try {
      NC::SymOp tmp( m, { { 0, 0, 0 } } );
    } catch ( NC::Error::BadInput& ) {
      accepted = false;
    }
    REQUIRE( valid == accepted );
    if ( !valid )
      continue;
    ++nvalid;
    for ( auto& t : translations ) {
      const NC::SymOp op( m, t );
      REQUIRE( entriesOK( op ) );
      for ( unsigned i = 0; i < 3; ++i ) {
        REQUIRE( op.trans(i) == ( ( t[i] % 24 ) + 24 ) % 24 );
        for ( unsigned j = 0; j < 3; ++j )
          REQUIRE( op.rot(i,j) == m[3*i+j] );
      }
      //Rotation type, determinant, order, powers:
      const int rt = op.rotationType();
      REQUIRE( ( rt > 0 ) == ( op.determinant() == 1 ) );
      if ( &t == &translations.front() )
        ++ntypes[rt];
      const unsigned order = op.rotationOrder();
      NC::SymOp p = op;
      for ( unsigned k = 1; k < order; ++k ) {
        REQUIRE( entriesOK( p ) );
        REQUIRE( !( p.isIdentity() ) || op.hasTranslation() );
        p = p * op;
      }
      REQUIRE( entriesOK( p ) );
      //R^order=I (translations need not vanish, e.g. for screw axes):
      REQUIRE( NC::SymOp( { { p.rot(0,0), p.rot(0,1), p.rot(0,2),
                              p.rot(1,0), p.rot(1,1), p.rot(1,2),
                              p.rot(2,0), p.rot(2,1), p.rot(2,2) } },
                          { { 0, 0, 0 } } ).isIdentity() );
      //Inverse:
      const NC::SymOp inv = op.inverse();
      REQUIRE( entriesOK( inv ) );
      REQUIRE( ( op * inv ).isIdentity() );
      REQUIRE( ( inv * op ).isIdentity() );
      REQUIRE( inv.inverse() == op );
      //String round trips and stream output:
      const std::string s = op.toString();
      REQUIRE( NC::SymOp( s ) == op );
      REQUIRE( NC::SymOp( NC::StrView( s ) ) == op );
      REQUIRE( NC::SymOp( s.c_str() ) == op );
      std::ostringstream ss;
      ss << op;
      REQUIRE( ss.str() == s );
    }
  }
  std::cout << "Accepted rotation matrices: " << nvalid << std::endl;
  REQUIRE( nvalid == 64 );

  //All products and inverses within each of the two families are valid, but
  //mixing hexagonal-only and orthogonal-only rotations is not:
  unsigned nmixed_invalid = 0;
  for ( auto& a : allowed ) {
    for ( auto& b : allowed ) {
      const NC::SymOp opa( a, { { 0, 0, 0 } } );
      const NC::SymOp opb( b, { { 0, 0, 0 } } );
      bool ok = true;
      try {
        NC::SymOp tmp = opa * opb;
        REQUIRE( entriesOK( tmp ) );
      } catch ( NC::Error::LogicError& ) {
        ok = false;
      }
      REQUIRE( ok == ( allowed.count( matMult( a, b ) ) > 0 ) );
      if ( !ok )
        ++nmixed_invalid;
    }
  }
  std::cout << "Products of two accepted rotations which are invalid: "
            << nmixed_invalid << " (of " << 64*64 << ")" << std::endl;
  try {
    NC::SymOp( "x-y,x,z" ) * NC::SymOp( "z,x,y" );
    REQUIRE( false );
  } catch ( NC::Error::LogicError& e ) {
    std::cout << "  Example: " << e.what() << std::endl;
  }
  std::cout << "Numbers per rotation type:";
  for ( auto& e : ntypes )
    std::cout << " " << e.first << ":" << e.second;
  std::cout << std::endl;

  //Associativity and group properties within cubic and hexagonal groups
  //(incl. translations):
  for ( auto gens : { std::vector<NC::SymOp>{ NC::SymOp( "z,x,y" ),
                                               NC::SymOp( "-y+1/4,x+3/4,z+1/4" ),
                                               NC::SymOp( "-x,-y,-z" ) },
                      std::vector<NC::SymOp>{ NC::SymOp( "x-y,x,z+1/6" ),
                                               NC::SymOp( "y,x,-z" ),
                                               NC::SymOp( "-x,-y,-z" ) } } ) {
    const auto g = closure( gens );
    std::cout << "Closure of";
    for ( auto& e : gens )
      std::cout << " " << e;
    std::cout << ": " << g.size() << " operations (modulo 1)" << std::endl;
    for ( auto& a : g ) {
      REQUIRE( entriesOK( a ) );
      for ( auto& b : g ) {
        REQUIRE( entriesOK( a * b ) );
        for ( std::size_t k = 0; k < g.size(); k += 7 ) {
          const auto& c = g[k];
          REQUIRE( ( a * b ) * c == a * ( b * c ) );
        }
      }
    }
  }

  //Apply:
  {
    const NC::SymOp op( "-y+1/4,x+3/4,z+1/4" );
    const NC::Vector v = op.apply( NC::Vector( 0.1, 0.2, 0.3 ) );
    REQUIRE( NC::floateq( v.x(), -0.2 + 0.25, 1e-15, 1e-15 ) );
    REQUIRE( NC::floateq( v.y(), 0.1 + 0.75, 1e-15, 1e-15 ) );
    REQUIRE( NC::floateq( v.z(), 0.3 + 0.25, 1e-15, 1e-15 ) );
    //Inverse gives back original position modulo 1 (since translations are
    //reduced modulo 1):
    const NC::Vector w = op.inverse().apply( v );
    auto mod1 = []( double x ) { return x - std::floor( x ); };
    REQUIRE( NC::floateq( mod1( w.x() ), 0.1, 1e-14, 1e-14 ) );
    REQUIRE( NC::floateq( mod1( w.y() ), 0.2, 1e-14, 1e-14 ) );
    REQUIRE( NC::floateq( mod1( w.z() ), 0.3, 1e-14, 1e-14 ) );
  }

  std::cout << "Accepted input:" << std::endl;
  for ( auto s : { "x,y,z", " X , Y , Z ", "1/2+x,1/2-y,z", "-y+x,x,z+1/6",
                   "+x,-y,+z", "x+1/2+1/4,y,z", "x-1/2,y,z", "x+3/2,y,z",
                   "x+2/4,y,z", "x+1/8,y,z", "x,y,z+1", "x,y,z-1/24",
                   "-x+y,-x,z+2/3", "y,x,-z+5/12" } )
    std::cout << "  \"" << s << "\" -> " << NC::SymOp( s ) << std::endl;

  std::cout << "Invalid input:" << std::endl;
  for ( auto s : { "", "x,y", "x,y,z,x", ",y,z", "2x,y,z", "x,x,z",
                   "x+y+z,y,z", "x+1/5,y,z", "x+0.5,y,z", "x,y,z+.5",
                   "a,b,c", "x y,y,z", "x+,y,z", "x+1/0,y,z", "x-x,y,z",
                   "x+1//2,y,z", "x,y,-z-", "--x,y,z", "x,y,1/2",
                   "x+1234567890,y,z" } )
    expectBad( s );
  try {
    NC::SymOp op( { { 2,0,0, 0,1,0, 0,0,1 } }, { { 0, 0, 0 } } );
    REQUIRE( false );
  } catch ( NC::Error::BadInput& e ) {
    std::cout << "  (2,0,0,0,1,0,0,0,1) -> BadInput: " << e.what()
              << std::endl;
  }
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Symmetry of space groups
////////////////////////////////////////////////////////////////////////////////

//Tests SGSymmetry (symmetry operations derived from the Hall symbols of the
//530 space group settings): group properties, canonical ordering, total
//counts, consistency of the Laue orbits with the EqRefl class (which derives
//its tables from SGSymmetry, so this checks its encoding and sign handling
//against a direct orbit calculation), and errors from the Hall symbol parser.

namespace test_symmetry {

namespace {
  using Mat = std::array<int,9>;
  Mat rotOf( const NC::SymOp& op )
  {
    Mat m;
    for ( unsigned i = 0; i < 9; ++i )
      m[i] = op.rot( i / 3, i % 3 );
    return m;
  }
  constexpr Mat identity = { { 1,0,0, 0,1,0, 0,0,1 } };
  constexpr Mat inversion = { { -1,0,0, 0,-1,0, 0,0,-1 } };

  //Does EqRefl support the setting (only default orientations of Laue
  //classes, i.e. unique axis b for monoclinic, and hexagonal or rhombohedral
  //axes for trigonal)?
  bool eqReflApplies( const NC::SpaceGroup& sg )
  {
    if ( sg.number() >= 3 && sg.number() <= 15 ) {
      const std::string c = sg.choice();
      return !c.empty() && ( c[0] == 'b' || c.substr( 0, 2 ) == "-b" );
    }
    return true;
  }

  //Compare Laue orbits of hkl (from rotation parts plus inversion) with
  //EqRefl families. Returns number of hkl points checked:
  unsigned checkEqRefl( const NC::SpaceGroup& sg,
                        const NC::SGSymmetry& sym )
  {
    const NC::EqRefl eqrefl( static_cast<int>( sg.number() ),
                             std::string( sg.choice() ) == "r" );
    unsigned n = 0;
    for ( int h = -4; h <= 4; ++h ) {
      for ( int k = -4; k <= 4; ++k ) {
        for ( int l = -4; l <= 4; ++l ) {
          if ( !h && !k && !l )
            continue;
          //Reflections transform with the transposed rotation matrix:
          std::set<std::array<int,3>> orbit;
          for ( auto& op : sym.representatives() ) {
            std::array<int,3> v;
            for ( unsigned j = 0; j < 3; ++j )
              v[j] = op.rot(0,j) * h + op.rot(1,j) * k + op.rot(2,j) * l;
            orbit.insert( v );
            orbit.insert( { { -v[0], -v[1], -v[2] } } );
          }
          std::set<std::array<int,3>> fam;
          for ( auto& e : eqrefl.getEquivalentReflections( h, k, l ) ) {
            fam.insert( { { e.h, e.k, e.l } } );
            fam.insert( { { -e.h, -e.k, -e.l } } );
          }
          REQUIRE( orbit == fam );
          ++n;
        }
      }
    }
    return n;
  }

  void expectBadHall( const char * s )
  {
    try {
      NC::detail::rawSymOpsFromHallSymbol( s );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  \"" << s << "\" -> BadInput: " << e.what() << std::endl;
      return;
    }
    nc_assert_always( false );
  }
}

void run()
{
  std::cout << "setting index, setting, order, number of representatives,"
            << " lattice symbol, centrosymmetric:" << std::endl;
  unsigned total_ops = 0, total_reps = 0, n_eqrefl = 0, n_eqrefl_hkl = 0;
  std::set<Mat> all_rots;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
    const NC::SGSymmetry& sym = NC::SGSymmetry::get( sg );
    REQUIRE( &sym == &NC::SGSymmetry::get( sg ) );//cached
    REQUIRE( sym.spaceGroup() == sg );
    const auto& ops = sym.operations();
    const auto& reps = sym.representatives();
    const auto& cent = sym.centringVectors();
    REQUIRE( sym.order() == ops.size() );
    REQUIRE( ops.size() == reps.size() * cent.size() );
    total_ops += sym.order();
    total_reps += static_cast<unsigned>( reps.size() );
    //Canonical order: identity first, then blocks per centring vector:
    REQUIRE( ops.front().isIdentity() && reps.front().isIdentity()
             && cent.front().isIdentity() );
    for ( std::size_t ic = 0; ic < cent.size(); ++ic )
      for ( std::size_t ir = 0; ir < reps.size(); ++ir )
        REQUIRE( ops.at( ic * reps.size() + ir ) == cent.at( ic ) * reps.at( ir ) );
    REQUIRE( std::is_sorted( reps.begin() + 1, reps.end() ) );
    REQUIRE( std::is_sorted( cent.begin(), cent.end() ) );
    //Centring vectors are pure translations, representatives have distinct
    //rotations, and each has the smallest translation among operations with
    //that rotation:
    for ( auto& c : cent )
      REQUIRE( rotOf( c ) == identity );
    std::set<Mat> rots;
    for ( auto& r : reps ) {
      REQUIRE( rots.insert( rotOf( r ) ).second );
      for ( auto& op : ops )
        if ( rotOf( op ) == rotOf( r ) )
          REQUIRE( !( op < r ) );
    }
    all_rots.insert( rots.begin(), rots.end() );
    //Group properties (closure, inverses):
    const std::set<NC::SymOp> opset( ops.begin(), ops.end() );
    REQUIRE( opset.size() == ops.size() );
    for ( auto& a : ops ) {
      REQUIRE( opset.count( a.inverse() ) );
      for ( auto& b : ops )
        REQUIRE( opset.count( a * b ) );
    }
    //Lattice symbol and centrosymmetry:
    const char lat = sym.latticeSymbol();
    const std::size_t ncent_expected = ( lat == 'p' ? 1 : lat == 'r' ? 3
                                         : lat == 'f' ? 4 : 2 );
    REQUIRE( std::string( "pabcirf" ).find( lat ) != std::string::npos );
    REQUIRE( cent.size() == ncent_expected );
    REQUIRE( sym.isCentrosymmetric() == ( rots.count( inversion ) > 0 ) );
    //Laue classes consistent with EqRefl:
    if ( eqReflApplies( sg ) ) {
      ++n_eqrefl;
      n_eqrefl_hkl += checkEqRefl( sg, sym );
    }
    std::cout << "  " << hn << " " << sg << " " << sym.order() << " "
              << reps.size() << " " << lat << " "
              << ( sym.isCentrosymmetric() ? 1 : 0 ) << std::endl;
  }
  std::cout << "Total number of operations: " << total_ops
            << ", representatives: " << total_reps
            << ", distinct rotations: " << all_rots.size() << std::endl;
  REQUIRE( total_ops == 7388 && total_reps == 4462 && all_rots.size() == 64 );
  std::cout << "Laue classes consistent with EqRefl for " << n_eqrefl
            << " settings (" << n_eqrefl_hkl << " hkl points in total)"
            << std::endl;

  std::cout << "Examples:" << std::endl;
  for ( auto s : { "14:b1", "15:b1", "166:h", "166:r", "194", "227:2" } ) {
    const NC::SpaceGroup sg( s, NC::SpaceGroupParseMode::Strict );
    const auto& sym = NC::SGSymmetry::get( sg );
    std::cout << "  " << sg << " (Hall symbol \"" << sg.hallSymbol()
              << "\"):" << std::endl << "    representatives:";
    for ( auto& op : sym.representatives() )
      std::cout << " " << op;
    std::cout << std::endl << "    centring vectors:";
    for ( auto& op : sym.centringVectors() )
      std::cout << " " << op;
    std::cout << std::endl;
  }

  std::cout << "Invalid Hall symbols:" << std::endl;
  for ( auto s : { "", "X 1", "-", "P 5", "P 2q", "P 22", "P 2 2 2 2",
                   "P 3* 2\"", "P 6 3*", "P 4 (0 0", "P 1 (0 0 1 1)",
                   "P 2 (a b c)", "P 1 (0 0 1) x", "P 32z*", "P 3 3" } )
    expectBadHall( s );
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Crystal systems, settings and cell constraints
////////////////////////////////////////////////////////////////////////////////

//Tests crystal systems, setting information and cell constraints of the 530
//space group settings. The cell constraints (derived from the symmetry
//operations) are checked against the expected patterns, and against the
//existing checkAndCompleteLattice[Angles] functions for all settings those
//support.

namespace test_cell {

namespace {
  using LC = NC::SGCellConstraints::LengthConstraint;
  using AC = NC::SGCellConstraints::AngleConstraint;
  using CS = NC::SGCrystalSystem;

  //Expected constraints, by hand from the crystal system and setting:
  NC::SGCellConstraints expectedConstraints( const NC::SpaceGroup& sg )
  {
    const auto si = sg.settingInfo();
    switch ( sg.crystalSystem() ) {
    case CS::Triclinic:
      return { LC::Free, LC::Free, AC::Free, AC::Free, AC::Free };
    case CS::Monoclinic:
      switch ( si.uniqueAxis ) {
      case NC::SGUniqueAxis::b: case NC::SGUniqueAxis::minus_b:
        return { LC::Free, LC::Free, AC::Is90, AC::Free, AC::Is90 };
      case NC::SGUniqueAxis::c: case NC::SGUniqueAxis::minus_c:
        return { LC::Free, LC::Free, AC::Is90, AC::Is90, AC::Free };
      case NC::SGUniqueAxis::a: case NC::SGUniqueAxis::minus_a:
        return { LC::Free, LC::Free, AC::Free, AC::Is90, AC::Is90 };
      default:
        REQUIRE( false );
      }
      break;
    case CS::Orthorhombic:
      return { LC::Free, LC::Free, AC::Is90, AC::Is90, AC::Is90 };
    case CS::Tetragonal:
      return { LC::EqualToA, LC::Free, AC::Is90, AC::Is90, AC::Is90 };
    case CS::Trigonal:
    case CS::Hexagonal:
      if ( si.hexFamilyAxes == NC::SGHexFamilyAxes::Rhombohedral )
        return { LC::EqualToA, LC::EqualToA, AC::Free, AC::EqualToAlpha,
                 AC::EqualToAlpha };
      return { LC::EqualToA, LC::Free, AC::Is90, AC::Is90, AC::Is120 };
    case CS::Cubic:
      return { LC::EqualToA, LC::EqualToA, AC::Is90, AC::Is90, AC::Is90 };
    }
    REQUIRE( false );
    return {};
  }

  bool operator==( const NC::SGCellConstraints& x,
                   const NC::SGCellConstraints& y )
  {
    return ( x.b == y.b && x.c == y.c && x.alpha == y.alpha
             && x.beta == y.beta && x.gamma == y.gamma );
  }

}

void run()
{
  std::cout << "setting index, setting, crystal system, setting info, cell"
            << " constraints:" << std::endl;
  std::set<std::string> patterns;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
    const auto si = sg.settingInfo();
    const auto cs = sg.crystalSystem();
    const auto& sym = NC::SGSymmetry::get( sg );
    const auto& cc = sym.cellConstraints();

    //Setting info fields apply exactly where expected:
    REQUIRE( ( si.uniqueAxis != NC::SGUniqueAxis::NotApplicable )
             == ( cs == CS::Monoclinic ) );
    REQUIRE( si.cellChoice == NC::SGCellChoice::NotApplicable
             || cs == CS::Monoclinic );
    REQUIRE( ( si.axisPermutation != NC::SGAxisPermutation::NotApplicable )
             == ( cs == CS::Orthorhombic ) );
    REQUIRE( ( si.hexFamilyAxes != NC::SGHexFamilyAxes::NotApplicable )
             == ( cs == CS::Trigonal || cs == CS::Hexagonal ) );
    REQUIRE( si.originChoice == NC::SGOriginChoice::NotApplicable
             || cs == CS::Orthorhombic || cs == CS::Tetragonal
             || cs == CS::Cubic );
    //Hexagonal axes iff some rotation part is not a signed permutation:
    bool hexrot = false;
    for ( auto& op : sym.representatives() )
      for ( unsigned i = 0; i < 3; ++i ) {
        int nnonzero = 0;
        for ( unsigned j = 0; j < 3; ++j )
          if ( op.rot( i, j ) )
            ++nnonzero;
        if ( nnonzero != 1 )
          hexrot = true;
      }
    REQUIRE( hexrot == ( si.hexFamilyAxes == NC::SGHexFamilyAxes::Hexagonal ) );

    //Cell constraints as expected:
    REQUIRE( cc == expectedConstraints( sg ) );
    std::ostringstream ss;
    ss << cc;
    patterns.insert( ss.str() );

    std::cout << "  " << hn << " " << sg << " " << cs << " " << si << " ["
              << cc << "]" << std::endl;
  }
  std::cout << "Distinct cell constraint patterns: " << patterns.size()
            << std::endl;
  for ( auto& p : patterns )
    std::cout << "  " << p << std::endl;
  REQUIRE( patterns.size() == 9 );
  std::cout << "Examples of errors:" << std::endl;
  auto showErr = []( const char * sgstr, NC::CellParameters cp, bool complete )
  {
    const auto& cc = NC::SGSymmetry::get( NC::SpaceGroup( sgstr, NC::SpaceGroupParseMode::Strict ) )
      .cellConstraints();
    try {
      if ( complete )
        cc.complete( cp );
      else
        cc.check( cp );
      REQUIRE( false );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  " << sgstr << ( complete ? " complete: " : " check: " )
                << e.what() << std::endl;
    }
  };
  showErr( "225", { 4.0, 0.0, 0.0, 90.0, 90.0, 90.0 }, false );
  showErr( "14:b1", { 4.0, 5.0, 6.0, 90.0, 0.0, 90.0 }, true );
  showErr( "14:c1", { 4.0, 5.0, 6.0, 90.0, 100.0, 90.0 }, false );
  showErr( "194", { 4.0, 4.1, 6.0, 90.0, 90.0, 120.0 }, false );
  showErr( "194", { 4.0, 4.0, 6.0, 90.0, 90.0, 119.9 }, false );
  showErr( "166:r", { 4.0, 0.0, 0.0, 120.0, 0.0, 0.0 }, true );
  showErr( "166:r", { 4.0, 4.0, 4.0, 60.0, 60.0, 61.0 }, false );
  showErr( "1", { 4.0, 5.0, 6.0, 90.0, 180.0, 90.0 }, false );
  showErr( "1", { -4.0, 5.0, 6.0, 90.0, 90.0, 90.0 }, false );
  NC::CellParameters cp{ 4.0, 0.0, 0.0, 70.0, 0.0, 0.0 };
  NC::SGSymmetry::get( NC::SpaceGroup( "166:r", NC::SpaceGroupParseMode::Strict ) ).cellConstraints()
    .complete( cp );
  REQUIRE( cp.b == 4.0 && cp.c == 4.0 && cp.beta == 70.0 && cp.gamma == 70.0 );
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Expansion of sites into orbits
////////////////////////////////////////////////////////////////////////////////

//Tests expandSiteToOrbit: general and special positions in all 530 space
//group settings, symmetrisation of rounded special positions, rejection of
//sites near (but not on) special positions, and wrapping of coordinates.

namespace test_orbit {

namespace {
  double pdist( const NC::Vector& a, const NC::Vector& b )
  {
    double m = 0.0;
    const double d[3] = { a.x() - b.x(), a.y() - b.y(), a.z() - b.z() };
    for ( double x : d )
      m = NC::ncmax( m, NC::ncabs( x - std::round( x ) ) );
    return m;
  }

  bool inUnitCell( const NC::Vector& v )
  {
    return ( v.x() >= 0.0 && v.x() < 1.0 && v.y() >= 0.0 && v.y() < 1.0
             && v.z() >= 0.0 && v.z() < 1.0 );
  }

  NC::SGSiteOrbit expand( const char * sg, const NC::Vector& site )
  {
    return NC::expandSiteToOrbit( NC::SGSymmetry::get( NC::SpaceGroup( sg, NC::SpaceGroupParseMode::Relaxed ) ),
                                  site );
  }

  void show( const char * sg, const NC::Vector& site, unsigned expected_mult )
  {
    const auto orbit = expand( sg, site );
    std::cout << "  " << sg << " site (" << NC::fmt( site.x(), "%.6g" )
              << ", " << NC::fmt( site.y(), "%.6g" ) << ", "
              << NC::fmt( site.z(), "%.6g" ) << "): "
              << orbit.positions.size() << " positions, site symmetry order "
              << orbit.siteSymmetryOrder << std::endl;
    REQUIRE( orbit.positions.size() == expected_mult );
  }

  void expectBad( const char * sg, const NC::Vector& site )
  {
    try {
      expand( sg, site );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  " << sg << " -> BadInput: " << e.what() << std::endl;
      return;
    }
    REQUIRE( false );
  }
}

void run()
{
  //General position in all settings:
  const NC::Vector generic( 0.1234567, 0.2345678, 0.3456789 );
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
    const auto& sym = NC::SGSymmetry::get( sg );
    const auto orbit = NC::expandSiteToOrbit( sym, generic );
    REQUIRE( orbit.positions.size() == sym.order() );
    REQUIRE( orbit.siteSymmetryOrder == 1 );
    REQUIRE( pdist( orbit.positions.front(), generic ) == 0.0 );
    for ( std::size_t i = 0; i < orbit.positions.size(); ++i ) {
      REQUIRE( inUnitCell( orbit.positions[i] ) );
      for ( std::size_t j = 0; j < i; ++j )
        REQUIRE( pdist( orbit.positions[i], orbit.positions[j] ) > 1e-2 );
    }
  }
  std::cout << "General positions OK in all 530 settings" << std::endl;

  std::cout << "Special positions:" << std::endl;
  show( "225", NC::Vector( 0.0, 0.0, 0.0 ), 4 );
  show( "225", NC::Vector( 0.25, 0.25, 0.25 ), 8 );
  show( "227:1", NC::Vector( 0.0, 0.0, 0.0 ), 8 );
  show( "227:2", NC::Vector( 0.125, 0.125, 0.125 ), 8 );
  show( "194", NC::Vector( 1.0/3, 2.0/3, 0.25 ), 2 );
  show( "166:h", NC::Vector( 0.0, 0.0, 0.0 ), 3 );
  show( "166:r", NC::Vector( 0.0, 0.0, 0.0 ), 1 );
  show( "14:b1", NC::Vector( 0.0, 0.0, 0.0 ), 2 );
  show( "62", NC::Vector( 0.1, 0.25, 0.2 ), 4 );
  show( "221", NC::Vector( 0.5, 0.5, 0.5 ), 1 );

  std::cout << "Rounded special positions:" << std::endl;
  //4 digits are symmetrised to the exact special position:
  for ( auto s : { NC::Vector( 0.3333, 0.6667, 0.25 ),
                   NC::Vector( 0.3333, 0.6666, 0.25 ),
                   NC::Vector( 0.33333, 0.66667, 0.25 ) } ) {
    const auto orbit = expand( "194", s );
    REQUIRE( orbit.positions.size() == 2 );
    REQUIRE( orbit.siteSymmetryOrder == 12 );
    REQUIRE( pdist( orbit.positions.front(),
                    NC::Vector( 1.0/3, 2.0/3, 0.25 ) ) < 1e-14 );
  }
  std::cout << "  194 site (0.3333, 0.6667, 0.25) etc.: symmetrised to"
            << " (1/3, 2/3, 1/4)" << std::endl;
  //3 digits are rejected:
  expectBad( "194", NC::Vector( 0.333, 0.667, 0.25 ) );
  expectBad( "227:2", NC::Vector( 0.125, 0.125, 0.12 ) );

  std::cout << "Sites near a mirror plane (Pm with y->-y):" << std::endl;
  {
    //Merged and symmetrised (distance 4e-4):
    const auto o1 = expand( "6:b", NC::Vector( 0.1, 0.0002, 0.2 ) );
    REQUIRE( o1.positions.size() == 1 && o1.siteSymmetryOrder == 2 );
    REQUIRE( o1.positions.front().y() == 0.0 );
    //Distinct (distance 1.2e-2):
    const auto o2 = expand( "6:b", NC::Vector( 0.1, 0.006, 0.2 ) );
    REQUIRE( o2.positions.size() == 2 && o2.siteSymmetryOrder == 1 );
    std::cout << "  y=0.0002: merged (y=0), y=0.006: two positions" << std::endl;
    //Rejected (distance 8e-3):
    expectBad( "6:b", NC::Vector( 0.1, 0.004, 0.2 ) );
    //Across the periodic boundary (y=0.9998 is close to y=0.0002):
    const auto o3 = expand( "6:b", NC::Vector( 0.1, 0.9999, 0.2 ) );
    REQUIRE( o3.positions.size() == 1 && o3.positions.front().y() == 0.0 );
  }

  std::cout << "Wrapping and invalid input:" << std::endl;
  {
    const auto o = expand( "1", NC::Vector( 1.25, -0.5, 2.0 ) );
    REQUIRE( o.positions.size() == 1 );
    REQUIRE( o.positions.front().x() == 0.25 && o.positions.front().y() == 0.5
             && o.positions.front().z() == 0.0 );
    const auto o2 = expand( "1", NC::Vector( -0.0, -1e-30, 0.0 ) );
    REQUIRE( inUnitCell( o2.positions.front() ) );
    std::cout << "  (1.25, -0.5, 2.0) -> (0.25, 0.5, 0)" << std::endl;
  }
  expectBad( "1", NC::Vector( std::numeric_limits<double>::quiet_NaN(),
                              0.0, 0.0 ) );
  expectBad( "1", NC::Vector( 0.0, NC::kInfinity, 0.0 ) );
  std::cout << "All OK" << std::endl;
}

}

////////////////////////////////////////////////////////////////////////////////
// Structure symmetry check
////////////////////////////////////////////////////////////////////////////////

//Tests checkStructureSymmetry: structures expanded from a generic and a
//special site in each of the 530 settings must be matched to that setting
//with zero deviation, while perturbed or wrongly declared structures give the
//deviation. Also haveSameRotationParts, which tells settings apart only by
//their symmetry-equivalent reflections.

namespace test_structcheck {

namespace {
  //Positions by species: one generic orbit (species 0) and the orbit of the
  //origin (species 1):
  std::vector<std::vector<NC::Vector>> structure( const NC::SGSymmetry& sym )
  {
    std::vector<std::vector<NC::Vector>> res;
    res.push_back( NC::expandSiteToOrbit( sym, NC::Vector( 0.1234, 0.2345,
                                                           0.3456 ) ).positions );
    res.push_back( NC::expandSiteToOrbit( sym, NC::Vector( 0.0, 0.0,
                                                           0.0 ) ).positions );
    return res;
  }
  bool sameRots( const char * a, const char * b )
  {
    return NC::haveSameRotationParts( NC::SGSymmetry::get( NC::SpaceGroup( a, NC::SpaceGroupParseMode::Relaxed ) ),
                                      NC::SGSymmetry::get( NC::SpaceGroup( b, NC::SpaceGroupParseMode::Relaxed ) ) );
  }
}

void run()
{
  unsigned nok = 0;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupSettingIdx{ hn } };
    const bool rhombo = ( std::string( sg.choice() ) == "r" );
    const auto& sym = NC::SGSymmetry::get( sg );
    const auto res = NC::checkStructureSymmetry( sg.number(), rhombo,
                                                 structure( sym ) );
    REQUIRE( res.maxDeviation < 1e-12 );//(rounding only)
    if ( res.settingIdx != hn ) {
      //A few settings have identical operations (e.g. 68:1 and 68:1ba-c),
      //in which case the first one is returned:
      const NC::SpaceGroup other{ NC::SpaceGroupSettingIdx{ res.settingIdx } };
      REQUIRE( res.settingIdx < hn );
      auto ops1 = sym.operations();
      auto ops2 = NC::SGSymmetry::get( other ).operations();
      std::sort( ops1.begin(), ops1.end() );
      std::sort( ops2.begin(), ops2.end() );
      REQUIRE( ops1 == ops2 );
      std::cout << "  " << sg << " has the same operations as " << other
                << std::endl;
      continue;
    }
    ++nok;
  }
  std::cout << "Structures in all other " << nok << " settings matched to"
            << " their setting (deviation < 1e-12)" << std::endl;
  {
    auto pos = structure( NC::SGSymmetry::get( NC::SpaceGroup( "19", NC::SpaceGroupParseMode::Relaxed ) ) );
    pos.at(0).at(0) += NC::Vector( 0.01, 0.0, 0.0 );
    const auto res = NC::checkStructureSymmetry( 19, false, pos );
    const NC::SpaceGroup setting{ NC::SpaceGroupSettingIdx{ res.settingIdx } };
    std::cout << "SG-19 structure with one atom moved by 0.01: setting "
              << setting << ", deviation " << NC::fmt( res.maxDeviation, "%.3g" )
              << std::endl;
    REQUIRE( NC::floateq( res.maxDeviation, 0.01 ) );
  }
  for ( unsigned sgtrue : { 149u, 150u } ) {
    const unsigned sgclaim = ( sgtrue == 149 ? 150 : 149 );
    auto pos = structure( NC::SGSymmetry::get( NC::SpaceGroup( std::to_string( sgtrue ), NC::SpaceGroupParseMode::Relaxed ) ) );
    const auto res = NC::checkStructureSymmetry( sgclaim, false, pos );
    std::cout << "SG-" << sgtrue << " structure checked as SG-" << sgclaim
              << ": deviation " << NC::fmt( res.maxDeviation, "%.3g" ) << std::endl;
    REQUIRE( res.maxDeviation > 0.05 );
  }
  //Rhombohedral axes only match R settings:
  {
    auto pos = structure( NC::SGSymmetry::get( NC::SpaceGroup( "166:r", NC::SpaceGroupParseMode::Strict ) ) );
    REQUIRE( NC::checkStructureSymmetry( 166, true, pos ).settingIdx
             == NC::SpaceGroup( "166:r", NC::SpaceGroupParseMode::Strict ).settingIdx().get() );
    REQUIRE( NC::checkStructureSymmetry( 3, true, pos ).settingIdx == 0 );
  }
  std::cout << "Same symmetry-equivalent reflections: 62 vs 62:cab "
            << sameRots( "62", "62:cab" ) << ", 227 vs 227:2 "
            << sameRots( "227", "227:2" ) << ", 14 vs 14:b2 "
            << sameRots( "14", "14:b2" ) << ", 14 vs 14:c1 "
            << sameRots( "14", "14:c1" ) << ", 166 vs 166:r "
            << sameRots( "166", "166:r" ) << std::endl;
  REQUIRE( sameRots( "62", "62:cab" ) && sameRots( "227", "227:2" )
           && sameRots( "14", "14:b2" ) && !sameRots( "14", "14:c1" )
           && !sameRots( "166", "166:r" ) );
}

}

int main()
{
  std::cout << "==== Space group table ====" << std::endl;
  test_table::run();
  std::cout << "==== Symmetry operations ====" << std::endl;
  test_symop::run();
  std::cout << "==== Symmetry of space groups ====" << std::endl;
  test_symmetry::run();
  std::cout << "==== Crystal systems, settings and cell constraints ====" << std::endl;
  test_cell::run();
  std::cout << "==== Expansion of sites into orbits ====" << std::endl;
  test_orbit::run();
  std::cout << "==== Structure symmetry check ====" << std::endl;
  test_structcheck::run();
  return 0;
}
