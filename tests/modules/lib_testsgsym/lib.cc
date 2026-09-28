
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

//Gives Python tests access to the SpaceGroup class and its hardwired table of
//settings (fixme: replace by a JSON query).

#include "NCTestUtils/NCTestModUtils.hh"
#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    //Not declared in any header (only for testing):
    std::vector<SymOp> rawSymOpsFromHallSymbol( StrView );
  }
}

NCTEST_CTYPE_DICTIONARY
{
  return
    "int nctest_sgsym_number( int );"
    "const char * nctest_sgsym_choice( int );"
    "const char * nctest_sgsym_hallsymbol( int );"
    "const char * nctest_sgsym_hmsymbol( int );"
    "const char * nctest_sgsym_tostring( int );"
    "int nctest_sgsym_parse( const char * );"
    "int nctest_sgsym_parse_strict( const char * );"
    "const char * nctest_sgsym_parse_error( const char * );"
    "const char * nctest_sgsym_symop( const char * );"
    "const char * nctest_sgsym_symmetry( int );"
    "const char * nctest_sgsym_rawhall( const char * );"
    "const char * nctest_sgsym_settinginfo( int );"
    ;
}

namespace {
  NC::SpaceGroup sgFromSettingIdx( int hn )
  {
    nc_assert_always( hn >= 1 && hn <= 530 );
    return NC::SpaceGroup(
      NC::SpaceGroupSettingIdx{ static_cast<std::uint16_t>( hn ) } );
  }
}

NCTEST_CTYPES int nctest_sgsym_number( int hn )
{
  int res = -1;
  try {
    res = static_cast<int>( sgFromSettingIdx( hn ).number() );
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_choice( int hn )
{
  const char * res = nullptr;
  try {
    res = sgFromSettingIdx( hn ).choice();//static storage
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_hallsymbol( int hn )
{
  const char * res = nullptr;
  try {
    res = sgFromSettingIdx( hn ).hallSymbol();//static storage
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_hmsymbol( int hn )
{
  const char * res = nullptr;
  try {
    res = sgFromSettingIdx( hn ).hmSymbol();//static storage
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_tostring( int hn )
{
  static std::string s_buf;//NB: Not thread-safe, but fine for this test
  try {
    s_buf = sgFromSettingIdx( hn ).toString();
  } NCCATCH;
  return s_buf.c_str();
}

NCTEST_CTYPES int nctest_sgsym_parse_strict( const char * str )
{
  //Returns setting index, or 0 in case of BadInput:
  int res = -1;
  try {
    nc_assert_always( str );
    try {
      res = NC::SpaceGroup( str, NC::SpaceGroupParseMode::Strict ).settingIdx().get();
    } catch ( NC::Error::BadInput& ) {
      res = 0;
    }
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_parse_error( const char * str )
{
  //The BadInput message of strict parsing (empty if it succeeds):
  static std::string s_buf;//NB: Not thread-safe, but fine for this test
  s_buf.clear();
  try {
    nc_assert_always( str );
    try {
      NC::SpaceGroup sg( str, NC::SpaceGroupParseMode::Strict );
    } catch ( NC::Error::BadInput& e ) {
      s_buf = e.what();
    }
  } NCCATCH;
  return s_buf.c_str();
}

NCTEST_CTYPES int nctest_sgsym_parse( const char * str )
{
  //Returns setting index, or 0 in case of BadInput:
  int res = -1;
  try {
    nc_assert_always( str );
    try {
      res = NC::SpaceGroup( str, NC::SpaceGroupParseMode::Relaxed ).settingIdx().get();
    } catch ( NC::Error::BadInput& ) {
      res = 0;
    }
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_symop( const char * str )
{
  //Parses symmetry operation and returns its string representation (in
  //NCrystal's format) followed by the 9 rotation matrix entries and 3
  //translations (units of 1/24), separated by spaces. Returns "" in case of
  //BadInput.
  static std::string s_buf;//NB: Not thread-safe, but fine for this test
  try {
    nc_assert_always( str );
    try {
      const NC::SymOp op( str );
      std::ostringstream ss;
      ss << op;
      for ( unsigned i = 0; i < 3; ++i )
        for ( unsigned j = 0; j < 3; ++j )
          ss << ' ' << op.rot( i, j );
      for ( unsigned i = 0; i < 3; ++i )
        ss << ' ' << op.trans( i );
      s_buf = ss.str();
    } catch ( NC::Error::BadInput& ) {
      s_buf.clear();
    }
  } NCCATCH;
  return s_buf.c_str();
}

namespace {
  std::string joinOps( const std::vector<NC::SymOp>& ops )
  {
    std::ostringstream ss;
    for ( auto& op : ops )
      ss << ( &op == &ops.front() ? "" : ";" ) << op;
    return ss.str();
  }
}

NCTEST_CTYPES const char * nctest_sgsym_symmetry( int hn )
{
  //Returns "<lattice> <centrosymmetric> <order> | <operations> |
  //<representatives> | <centring vectors>", with operations separated by
  //semicolons.
  static std::string s_buf;//NB: Not thread-safe, but fine for this test
  try {
    const auto& sym = NC::SGSymmetry::get( sgFromSettingIdx( hn ) );
    std::ostringstream ss;
    ss << sym.latticeSymbol() << ' ' << ( sym.isCentrosymmetric() ? 1 : 0 )
       << ' ' << sym.order() << '|' << joinOps( sym.operations() ) << '|'
       << joinOps( sym.representatives() ) << '|'
       << joinOps( sym.centringVectors() );
    s_buf = ss.str();
  } NCCATCH;
  return s_buf.c_str();
}

NCTEST_CTYPES const char * nctest_sgsym_rawhall( const char * hall )
{
  //Operations from a Hall symbol (separated by semicolons), or "ERROR: <msg>"
  //in case of BadInput:
  static std::string s_buf;//NB: Not thread-safe, but fine for this test
  try {
    nc_assert_always( hall );
    try {
      s_buf = joinOps( NC::detail::rawSymOpsFromHallSymbol( hall ) );
    } catch ( NC::Error::BadInput& e ) {
      s_buf = std::string( "ERROR: " ) + e.what();
    }
  } NCCATCH;
  return s_buf.c_str();
}

NCTEST_CTYPES const char * nctest_sgsym_settinginfo( int hn )
{
  //Returns "<crystal system> <hex family axes> <unique axis> <cell choice>
  //<origin choice> <axis permutation>":
  static std::string s_buf;//NB: Not thread-safe, but fine for this test
  try {
    const auto sg = sgFromSettingIdx( hn );
    const auto si = sg.settingInfo();
    std::ostringstream ss;
    ss << sg.crystalSystem() << ' ' << si.hexFamilyAxes << ' '
       << si.uniqueAxis << ' ' << si.cellChoice << ' ' << si.originChoice
       << ' ' << si.axisPermutation;
    s_buf = ss.str();
  } NCCATCH;
  return s_buf.c_str();
}
