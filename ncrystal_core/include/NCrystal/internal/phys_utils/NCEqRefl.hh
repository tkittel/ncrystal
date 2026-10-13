#ifndef NCrystal_EqRefl_hh
#define NCrystal_EqRefl_hh

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

#include "NCrystal/core/NCTypes.hh"

//Class EqRefl provides symmetry-equivalent reflections for a given space group
//number, by providing a list of all (h,k,l) indices symmetry-equivalent to a
//given (h,k,l) index. Or more precisely, half of the list, since the symmetry
//(h,k,l)->(-h,-k,-l) is trivial and would be a waste of memory to expand
//explicitly.
//
//The equivalences are those of the Laue class of the space group, derived at
//construction from the symmetry operations (the rotation parts of the
//operations, modulo sign).
//
//Note that while two equivalent (h,k,l) indices will always have identical
//values of d-spacing and structure factors (F^2), the opposite is not in
//general true. So (h,k,l) groupings based on (d,F^2) values might contain more
//than one symmetry-equivalency group in each (d,F^2) group.

namespace NCRYSTAL_NAMESPACE {

  class SpaceGroup;
  class SGSymmetry;

  class EqRefl final {
  public:

    //Initialise with space group number in 1..230. For the rhombohedral
    //(R) space groups, rhombohedral_axes indicates that the unit cell is
    //described in rhombohedral rather than hexagonal axes (see also
    //usesRhombohedralAxes in NCLatticeUtils.hh). Apart from that, the
    //default setting of the space group number is assumed (which only
    //matters for monoclinic groups, where unique axis b is assumed):
    EqRefl( int spacegroup, bool rhombohedral_axes = false );

    //Initialise from a specific space group setting (sgsym):
    explicit EqRefl( const SpaceGroup& );

    class EquivReflList {
      //List with half the HKL points in a given family of symmetric equivalent
      //HKL points (only half, since (h,k,l) and (-h,-k,-l) are trivially in the
      //same group). This uses non-dynamic storage to keep malloc overhead low,
      //and thus takes up roughly 100B of storage (it is thus mostly suitable
      //for short-lived usage on the stack).
    public:
      const HKL* begin() const noexcept;
      const HKL& front() const noexcept { return *begin(); }
      const HKL* end() const noexcept;
      std::size_t size() const noexcept;
      bool contains( HKL ) const noexcept;//is (h,k,l) or (-h,-k,-l) in the list?

      //Empty (for later assignment):
      EquivReflList();

      //Copy only:
      EquivReflList( const EquivReflList& );
      EquivReflList& operator=( const EquivReflList& );

    private:
      HKL m_data[24];
      HKL * m_end;
      friend class EqRefl;
      void add(int, int, int);
    };

    EquivReflList getEquivalentReflections(int h, int k, int l) const;
    EquivReflList getEquivalentReflections( const HKL& hkl ) const
    {
      return getEquivalentReflections( hkl.h, hkl.k, hkl.l );
    }

    //Same as getEquivalentReflections(..).front() but slightly faster (for when
    //only that value is needed):
    HKL getEquivalentReflectionsRepresentativeValue( const HKL& hkl ) const;
    HKL getEquivalentReflectionsRepresentativeValue( int h, int k, int l ) const;

    //Number of distinct rotation parts modulo sign (i.e. half the order of
    //the Laue class, at most 24):
    unsigned nRotations() const noexcept { return m_nrot; }

  private:
    //Distinct rotation parts of the operations modulo sign, the first being
    //the identity. Reflections transform with the transposed matrix, h'_j =
    //sum_i h_i R_ij, and every column of the allowed rotation parts is one of
    //just ten combinations of (h,k,l) (see calcRaw), so each is stored as
    //three indices into a table of those combinations:
    using RowIdx = std::array<std::uint8_t,3>;
    unsigned m_nrot = 0;
    RowIdx m_rows[24];
    EquivReflList calcRaw(int h, int k, int l) const;
    void init( const SGSymmetry& );
  };
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  inline const HKL* EqRefl::EquivReflList::begin() const noexcept { return &m_data[0]; }
  inline const HKL* EqRefl::EquivReflList::end() const noexcept { return m_end; }
  inline std::size_t EqRefl::EquivReflList::size() const noexcept { return std::distance(begin(),end()); }
  inline EqRefl::EquivReflList::EquivReflList() { m_end = &m_data[0]; }

  inline EqRefl::EquivReflList::EquivReflList( const EquivReflList& o )
  {
    m_end = &m_data[0];
    for ( auto& e : o )
      *m_end++ = e;
  }

  inline EqRefl::EquivReflList& EqRefl::EquivReflList::operator=( const EquivReflList& o )
  {
    m_end = &m_data[0];
    for ( auto& e : o )
      *m_end++ = e;
    return *this;
  }

  inline void EqRefl::EquivReflList::add(int h, int k, int l)
  {
    //Only ever consider the form of (h,k,l) and (-h,-k,-l) which has first
    //non-zero coordinate positive:
    HKL a(h,k,l);
    auto am = a.flipped();
    *m_end++ =  am < a ? am : a;
  }

  inline EqRefl::EquivReflList EqRefl::calcRaw(int h, int k, int l) const
  {
    //The possible values of a transformed index (keep in sync with the
    //encoding in the constructor):
    const int v[10] = { h, k, l, -h, -k, -l, h+k, -h-k, h-k, k-h };
    EquivReflList e;
    for ( unsigned i = 0; i < m_nrot; ++i ) {
      const RowIdx& r = m_rows[i];
      e.add( v[r[0]], v[r[1]], v[r[2]] );
    }
    return e;
  }

  inline EqRefl::EquivReflList EqRefl::getEquivalentReflections(int h, int k, int l) const
  {
    auto raw = calcRaw(h,k,l);
    //Sort and discard duplicates:
    std::sort((HKL*)&raw.m_data[0],raw.m_end);
    raw.m_end = std::unique((HKL*)&raw.m_data[0],raw.m_end);
    return raw;
  }

  inline HKL EqRefl::getEquivalentReflectionsRepresentativeValue( int h, int k, int l ) const
  {
    auto raw = calcRaw(h,k,l);
    auto itLow = raw.begin();
    auto itE = raw.end();
    for ( auto it = std::next(itLow); it!=itE; ++it ) {
      if ( *it < *itLow )
        itLow = it;
    }
    return *itLow;
  }

  inline HKL EqRefl::getEquivalentReflectionsRepresentativeValue( const HKL& e ) const
  {
    return getEquivalentReflectionsRepresentativeValue( e.h, e.k, e.l );
  }

  inline bool EqRefl::EquivReflList::contains( HKL a ) const noexcept
  {
    //At worst 24 entries, often just a few handfuls => use simple linear search.
    HKL am(-a.h,-a.k,-a.l);
    if ( am < a )
      a = am;
    for ( auto& e : *this ) {
      if ( e == a )
        return true;
    }
    return false;
  }

}

#endif
