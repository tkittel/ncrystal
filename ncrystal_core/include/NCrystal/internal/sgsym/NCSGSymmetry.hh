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

    //The parameters with each constrained one set exactly to its required
    //value (b and c to a, beta and gamma to alpha, or 90 or 120 degrees), if
    //its deviation is at most relTol times its value plus the rounding errors
    //of the parameters involved (e.g. of the data printed in a file, all zero
    //by default), and then check(..)'ed. Throws BadInput as check(..)
    //otherwise:
    struct Tolerances {
      double relTol;
      static constexpr Tolerances strict() noexcept { return { 1e-12 }; }
      static constexpr Tolerances lenient() noexcept { return { 1e-6 }; }
    };
    CellParameters snap( const CellParameters&, const Tolerances&,
                         const CellParameters& roundingErrors
                         = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ) const;
    //As snap(..), but returns the message of the BadInput instead of throwing
    //it (the cell is then NullOpt):
    struct SnapResult {
      Optional<CellParameters> cell;
      Optional<std::string> errorMsg;
    };
    SnapResult trySnap( const CellParameters&, const Tolerances&,
                        const CellParameters& roundingErrors
                        = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ) const;
    //Fill in any parameters given as 0 whose values are implied by the
    //constraints (e.g. b and c for cubic space groups, or gamma=120 for
    //hexagonal axes), and then snap(..) the result with strict tolerances:
    void complete( CellParameters& ) const;
  private:
    //The message check(..) would throw, if any:
    Optional<std::string> checkError( const CellParameters& ) const;
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

  //Expansion of a site (fractional coordinates) into its orbit, i.e. all its
  //symmetry-equivalent positions, in a cell (for distances in Aa). Images of
  //the site closer than averageBelow mean a special position: the site is
  //replaced by their average (e.g. 0.3333 becomes 1/3). Images closer than
  //separateFrom are distinct copies of the atom, allowed only if the occupancy
  //is at most 1/k for k such copies (split sites), else BadInput is thrown
  //(e.g. for 0.333 instead of 1/3). The strict tolerances are for NCMAT data,
  //the lenient ones for data of limited precision (e.g. CIF files):
  struct SGSiteTolerances {
    double averageBelow, separateFrom;
    static constexpr SGSiteTolerances strict() noexcept { return { 1e-3, 0.5 }; }
    static constexpr SGSiteTolerances lenient() noexcept { return { 1e-2, 0.5 }; }
  };
  struct SGSiteOrbit {
    std::vector<Vector> positions;//wrapped into [0,1), symmetrised site first
    unsigned siteSymmetryOrder;//= group order / positions.size()
    unsigned nearCopies = 0;//positions closer than separateFrom to the first
    double nearestCopyDistance = 0.0;//(in Aa, if nearCopies>0)
  };
  SGSiteOrbit expandSiteToOrbit( const SGSymmetry&, const CellParameters&,
                                 const Vector& site, double occupancy,
                                 const SGSiteTolerances& );

  //The inverse grouping of positions into orbits (each matched to the nearest
  //symmetry image within separateFrom/2). Returns the indices of the
  //positions of each orbit, in increasing order, the orbits ordered by their
  //first index. Throws BadInput if the positions do not form complete orbits:
  std::vector<std::vector<unsigned>>
  groupPositionsIntoOrbits( const SGSymmetry&, const CellParameters&,
                            const std::vector<Vector>& positions,
                            double occupancy, const SGSiteTolerances& );

  //Check of a crystal structure against the settings of a space group number:
  //positions (fractional coordinates, grouped by species) are mapped by the
  //coset representatives and centring translations of each setting of the
  //number (only R settings for rhombohedral_axes, and never otherwise), and the
  //deviation of a setting is the largest per-coordinate periodic distance from
  //an image to the nearest position of the same species (0 for an exactly
  //symmetric structure, in which case the structure is invariant under all
  //operations). Returns the setting with the smallest deviation (the default
  //setting on ties), or no setting if the number has none compatible
  //with the axes:
  struct SGStructureCheck {
    Optional<SpaceGroupSettingIdx> settingIdx;
    double maxDeviation = kInfinity;
  };
  SGStructureCheck checkStructureSymmetry( unsigned sgnumber,
                                           bool rhombohedral_axes,
                                           const std::vector<std::vector<Vector>>& positions );

  //Same deviation, but for one given setting only:
  SGStructureCheck checkStructureSymmetry( const SpaceGroup& setting,
                                           const std::vector<std::vector<Vector>>& positions );

  //Whether two settings have the same rotation parts modulo sign, i.e. the
  //same symmetry-equivalent reflections (e.g. true for the orthorhombic
  //axis permutations and the two origin choices, false for monoclinic
  //unique axis b versus c):
  bool haveSameRotationParts( const SGSymmetry&, const SGSymmetry& );

  //The settings of the space group number (default first) which are compatible
  //with the provided cell (which may have implied parameters left as 0):
  std::vector<SpaceGroup> settingsCompatibleWithCell( unsigned number,
                                                      const CellParameters& );

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
