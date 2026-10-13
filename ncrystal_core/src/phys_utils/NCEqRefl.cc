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

#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/sgsym/NCSGSymmetry.hh"

namespace NC = NCrystal;

NC::EqRefl::EqRefl( int sg, bool rhombohedral_axes )
{
  if (sg<1||sg>230)
    NCRYSTAL_THROW(BadInput,"Space group number is not in the range 1 to 230");
  if ( rhombohedral_axes && !isRhombohedralSpaceGroup( sg ) )
    NCRYSTAL_THROW2(BadInput,"Rhombohedral axes requested for space group "
                    <<sg<<" which is not rhombohedral.");
  //The symmetry-equivalent reflections are those of the Laue class, i.e. the
  //rotation parts of the operations combined with inversion. As (h,k,l) and
  //(-h,-k,-l) are trivially equivalent here, we keep the distinct rotation
  //parts modulo sign (the coset representatives have one operation per
  //distinct rotation part). Among the settings of a number, only the
  //monoclinic unique axis and the R/hexagonal axes change the Laue class in
  //the hkl basis, so the default setting (or the R setting) is used:
  std::string spec = std::to_string( sg );
  if ( rhombohedral_axes )
    spec += ":r";
  init( SGSymmetry::get( SpaceGroup( spec, SpaceGroupParseMode::Relaxed ) ) );
}

NC::EqRefl::EqRefl( const SpaceGroup& sg )
{
  init( SGSymmetry::get( sg ) );
}

void NC::EqRefl::init( const SGSymmetry& sym )
{
  //Encoding of a column (coefficients of h,k,l) as an index into the value
  //table of calcRaw, i.e. {h,k,l,-h,-k,-l,h+k,-h-k,h-k,k-h}:
  auto encodeColumn = []( int ch, int ck, int cl ) -> unsigned
  {
    static const int tbl[10][3] = { {1,0,0}, {0,1,0}, {0,0,1}, {-1,0,0},
                                    {0,-1,0}, {0,0,-1}, {1,1,0}, {-1,-1,0},
                                    {1,-1,0}, {-1,1,0} };
    for ( unsigned i = 0; i < 10; ++i )
      if ( tbl[i][0] == ch && tbl[i][1] == ck && tbl[i][2] == cl )
        return i;
    NCRYSTAL_THROW(LogicError,"Unexpected rotation part in space group"
                   " symmetry operation");
  };
  //Sign flips of the value table (index i -> index of the negated value):
  static const std::uint8_t negate[10] = { 3, 4, 5, 0, 1, 2, 7, 6, 9, 8 };
  for ( auto& op : sym.representatives() ) {
    RowIdx r;
    for ( unsigned j = 0; j < 3; ++j )
      r[j] = static_cast<std::uint8_t>( encodeColumn( op.rot(0,j),
                                                      op.rot(1,j),
                                                      op.rot(2,j) ) );
    //Modulo sign: of R and -R keep the one with the smallest encoding:
    RowIdx rneg = { { negate[r[0]], negate[r[1]], negate[r[2]] } };
    if ( rneg < r )
      r = rneg;
    bool seen = false;
    for ( unsigned i = 0; i < m_nrot; ++i )
      if ( m_rows[i] == r )
        seen = true;
    if ( seen )
      continue;
    nc_assert_always( m_nrot < 24 );
    m_rows[m_nrot++] = r;
  }
#ifndef NDEBUG
  const RowIdx identity = { { 0, 1, 2 } };
  nc_assert( m_nrot >= 1 && m_rows[0] == identity );
#endif
}
