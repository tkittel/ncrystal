#ifndef NCrystal_SGSymmetry_hh
#define NCrystal_SGSymmetry_hh

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

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"

namespace NCRYSTAL_NAMESPACE {

  //Unit cell parameters (lengths in Aa, angles in degrees):
  struct CellParameters { double a, b, c, alpha, beta, gamma; };

  //Constraints on the unit cell parameters imposed by the symmetry of a space
  //group setting (i.e. the cell metric must be invariant under all rotation
  //parts of the operations), which are derived from the operations. The
  //length a and the angle alpha are never constrained to other parameters.
  //
  //TODO: Replace checkAndCompleteLattice[Angles]() in NCLatticeUtils.hh with
  //this.
  struct SGCellConstraints {
    enum class LengthConstraint { Free, EqualToA };
    enum class AngleConstraint { Free, Is90, Is120, EqualToAlpha };
    LengthConstraint b, c;
    AngleConstraint alpha, beta, gamma;

    //Check fully specified parameters, with exact equality required for
    //constrained values. Also requires positive lengths, angles in (0,180),
    //and alpha<120 if all angles must be equal. Throws BadInput if invalid
    //(including for parameters which are 0, cf. complete(..) below):
    void check( const CellParameters& ) const;

    //Fill in any parameters given as 0 whose values are implied by the
    //constraints (e.g. b and c for cubic space groups, or gamma=120 for
    //hexagonal axes), and then check(..) the result:
    void complete( CellParameters& ) const;
  };

  //Output like "a, b=a, c, alpha=90, beta, gamma=90":
  std::ostream& operator<<( std::ostream&, const SGCellConstraints& );

  //Immutable symmetry information of one of the 530 space group settings,
  //derived from its Hall symbol. Instances are created on first use and live
  //until the end of the program, so references to them remain valid. Obtain
  //them via SGSymmetry::get(spacegroup).
  //
  //The operations are listed in a fixed canonical order: first one operation
  //per distinct rotation part (the one with the smallest translation; for
  //centred lattices, the other operations with the same rotation part differ
  //from it only by a centring vector), in sorted order but with the identity
  //first, followed by those same operations shifted by each additional centring
  //vector (in sorted order).

  class SGSymmetry final : private NoCopyMove {
  public:

    static const SGSymmetry& get( SpaceGroup );

    SpaceGroup spaceGroup() const noexcept { return m_sg; }

    //All operations (incl. the centring copies), the representatives (one
    //operation per distinct rotation part, see above), and the centring
    //vectors (as pure translations, starting with the zero translation):
    const std::vector<SymOp>& operations() const noexcept { return m_ops; }
    const std::vector<SymOp>& representatives() const noexcept;
    const std::vector<SymOp>& centringVectors() const noexcept;

    unsigned order() const noexcept;//number of operations
    char latticeSymbol() const noexcept;//'p','a','b','c','i','r' or 'f' (as
                                        //in the ITVB Hall entry)
    bool isCentrosymmetric() const noexcept;

    //Constraints on the unit cell parameters:
    const SGCellConstraints& cellConstraints() const noexcept;

  private:
    struct internal_t {};
    SGSymmetry( internal_t, SpaceGroup );
    SpaceGroup m_sg;
    char m_lattice;
    bool m_centrosymmetric;
    SGCellConstraints m_cellconstraints;
    std::vector<SymOp> m_ops, m_reps, m_centring;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline const std::vector<SymOp>& SGSymmetry::representatives() const noexcept
  {
    return m_reps;
  }
  inline const std::vector<SymOp>& SGSymmetry::centringVectors() const noexcept
  {
    return m_centring;
  }
  inline unsigned SGSymmetry::order() const noexcept
  {
    return static_cast<unsigned>( m_ops.size() );
  }
  inline char SGSymmetry::latticeSymbol() const noexcept
  {
    return m_lattice;
  }
  inline bool SGSymmetry::isCentrosymmetric() const noexcept
  {
    return m_centrosymmetric;
  }
  inline const SGCellConstraints& SGSymmetry::cellConstraints() const noexcept
  {
    return m_cellconstraints;
  }
}

#endif
