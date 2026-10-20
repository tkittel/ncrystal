# Crystal structure: unit cells, settings, sites and structure checks

## Summary

A crystalline material in NCrystal is defined by its conventional unit cell
(lengths and angles), its space group, and the atoms in the cell, given
either as the full list of positions or, since NCMAT v8, as one site per
crystallographic orbit. From this the library derives the cell volume,
density, number of atoms per cell, the composition fractions, the
symmetry-equivalent reflections needed for the hkl list (chapter 12), and,
with the atomic displacements of chapter 35, everything the elastic models
need. NCrystal has always validated this input carefully (cell parameters
against the crystal system, positions against each other, fractions against
the cell content); 5.0.0 adds the full space-group symmetry machinery
(chapter 11), so that the setting of the space group is known exactly,
positions are verified to have the declared symmetry, sites can be expanded
and positions grouped into sites, and several subtleties (lattice matrix for
triclinic cells, Laue class $\bar31m$, rhombohedral axes, non-default
settings) are handled correctly.

Code: `ncrystal_core/src/infobld/NCInfoBldStructure.cc`,
`ncrystal_core/src/ncmat/NCLoadNCMAT.cc`,
`ncrystal_core/include/NCrystal/internal/utils/NCLatticeUtils.hh`,
`ncrystal_core/include/NCrystal/interfaces/NCInfoTypes.hh`.

## Physics and conventions

- The unit cell is the conventional cell of the space-group setting
  (centred where the lattice symbol says so: an fcc metal is the cubic cell
  with four atoms), with lengths $a, b, c$ in angstrom and angles
  $\alpha, \beta, \gamma$ in degrees; fractional coordinates refer to this
  cell. The lattice is oriented with $\vec a$ along $x$ and $\vec b$ in the
  $xy$ plane, the reciprocal lattice is the inverse transpose times $2\pi$,
  and
  $V = abc\sqrt{1-\cos^2\alpha-\cos^2\beta-\cos^2\gamma+2\cos\alpha\cos\beta\cos\gamma}$.
- The symmetry is the space group of the structure in one of the 530
  settings of ITVB Table A1.4.2.7 (chapter 11). The atoms are the complete
  content of the cell: every position of every orbit, with its element (or
  isotope, or `@ATOMDB` mixture, chapter 43). Occupancies other than 1 and
  anisotropic displacements are not part of the model yet.
- Per element the model needs the number of atoms per cell $n_i$ (hence the
  composition fraction $n_i/\sum n_j$), the positions for the structure
  factors, and one isotropic displacement (chapter 35).
- Density: the mass density follows from the cell content and volume, the
  number density is $n_a/V$; both can be overridden by the `density` cfg
  parameter (chapter 42), and `strain` scales $a, b, c$ by $1 + $`strain`.

## Inputs

NCMAT (`docs/ncmat_doc.md`):

- `@CELL`: lengths and angles (`cubic a` as a shorthand). Parameters implied
  by the symmetry may be left as 0 and are completed (e.g. $b = c = a$ and
  $90^\circ$ angles for cubic groups, $\gamma = 120^\circ$ for hexagonal
  axes); constrained parameters that are given must match exactly.
- `@SPACEGROUP`: the number 1-230 (optional before v8 but recommended;
  without it hkl families are grouped by $(d, |F|^2)$ instead of symmetry).
  Before v8 a bare number is interpreted leniently: the setting is deduced
  from the positions among all settings of the number (so origin-choice-2
  structures keep loading), and the file is rejected if none fits. In v8 the
  setting choice code is appended for the 90 numbers with several settings
  (`227:2`, `14:c1`, `62:cab`, `166:r`, exactly the ITVB codes with `abc`
  for the default orthorhombic axis order), while the 140 single-setting
  numbers are written bare.
- `@ATOMPOSITIONS`: element and fractional coordinates per atom (fractions
  like `1/3` allowed); or, in v8, `@ATOMSITES` with one position per orbit
  (requires `@SPACEGROUP`; never both sections). Sites on special positions
  should be given as exact fractions.
- `@DYNINFO` fractions must agree with the cell content; isotopes and
  mixtures via labels and `@ATOMDB`.
- Other crystal sources: `.laz`/`.lau` files (chapter 52), `Info` objects
  built by plugins or by the Python `NCMATComposer` (chapter 60), including
  its CIF and ASE importers.

Cfg parameters: `dcutoff`, `dcutoffup` (hkl range, chapter 12), `strain`,
`temp`, `density`, `atomdb`.

## Implementation

When the `Info` object is built, the structure is validated and completed in
this order:

1. The setting is resolved: a setting declared by the data source is
   validated (a setting of the number, matching hexagonal versus
   rhombohedral axes); otherwise the number plus the atom positions
   determine it among all settings of the number compatible with the axes
   (chapter 11).
2. The cell is completed and checked with the cell constraints derived from
   the setting's symmetry operations (exact equality for constrained
   parameters). Monoclinic cells must have $\alpha = \gamma = 90^\circ$
   unless another unique axis is declared; hexagonal cells with
   $\gamma \ne 120^\circ$ are an error rather than silently corrected.
3. Sites are expanded to positions, or positions grouped into sites when
   they form complete orbits; a site given twice, two sites that are the
   same position, or two sites that are symmetry-equivalent are errors
   naming the sites.
4. Positions are checked against the symmetry of the setting: every
   symmetry image must coincide with a position of the same species. Small
   deviations (below about 1e-2 in fractional coordinates, the level at
   which structures clearly do not have the symmetry) are tolerated, with a
   warning when they exceed the rounding noise of typical data; the
   symmetric structure is then assumed. Error messages list the settings of
   the number and mark those the positions do fit (e.g. `17:abc` declared
   while the positions are in `17:cab`).
5. Positions are mapped into $[0,1)$, sorted, and checked for duplicates
   (periodically, also across cell boundaries), with messages such as "The
   same atom position is used twice: Al (0, 0, 0) and Al (0, 0, 1)".
6. Debye temperatures and displacements must be present for all atoms or
   none; the volume, the number of atoms per cell and the composition
   fractions follow.

The resulting `Info` object carries, per material, the space group number
and setting, the cell parameters, volume and atom count, and, per element,
the positions, the sites (one per orbit, with multiplicity; absent when
there is no space group or the positions do not form complete orbits), the
Debye temperature and displacement, and a link to the element's dynamics.
Python exposes this as `Info.structure_info` and `Info.getAtomInfo()` (with
`AtomInfo.sites`, `.positions`, `.count`, ...); JSON queries on `Info`
objects give the structure and the sites, and utility queries
(`["util","atomsites",...]`, `["util","atompositions",...]`) convert
between positions and sites for any space group, deducing the setting from
a bare number.

## Validation

- `tests/scripts/latticecells.py` and `long_latticecells.py` (many
  synthetic cells and all stdlib crystals against metric-tensor
  calculations and spglib), `cellvalidation.py`.
- `tests/scripts/ncmatsettings.py`, `ncmatsites.py`, `atomsites.py` (the v8
  settings policy, @ATOMSITES expansion, grouping, error messages),
  `tests/src/app_symchecks` (the position-based symmetry check against the
  structure-factor self-check of chapter 12), `app_cellconstraints` (the
  derived constraints against the pre-5.0.0 rules).
- `tests/scripts/composer.py`, `composer2.py`, `cif2ncmat.py`,
  `cif2ncmatgen.py` and `long_cif2ncmatgen.py` (synthetic CIFs in all 230
  space groups and many styles) exercise the structure pipeline end to end.
- CPC 2020 paper: positions and hkl lists of all data-library files
  cross-checked against nxslib/SgInfo.

## Approximations and limitations

- Full occupancy only: partial occupancy is emulated by `NCMATComposer`
  through an `@ATOMDB` mixture with a sterile pseudo-isotope (correct
  scattering including the disorder term, wrong mass density and
  composition bookkeeping); proper support is designed
  (`docs/devel_site_occupancy_roadmap.md`) and postponed.
- Isotropic displacements only, one value per element.
- Only the 530 standard settings; structures in a non-standard origin or
  basis must be transformed first (`NCMATComposer.refine_crystal_structure`
  via spglib does this).
- No chemical disorder beyond `@ATOMDB` mixtures, no superstructures,
  modulated or magnetic structures.

## New or changed in 5.0.0

NCrystal 4.4.6 already validated cells per crystal system, completed
implied parameters, checked duplicate positions, derived the composition
from the cell, and generated symmetry-correct hkl families for all 230
space groups in their default settings. New in 5.0.0: space-group settings
as first-class information (declared in NCMAT v8, deduced for older files),
rhombohedral axes and non-default settings (monoclinic unique axis c, axis
permutations, origin choice 2) loadable and correct, load-time verification
that the positions have the declared symmetry with informative messages,
`@ATOMSITES` and the sites on `Info`, and fixes to the lattice matrix of
triclinic and rhombohedral-axes cells (issue #383, changing the densities of
two stdlib files), to duplicate-position detection and to the monoclinic and
hexagonal cell checks.

## References

- U. Shmueli, S. R. Hall, R. W. Grosse-Kunstleve, "Space-group symbols for
  numeric and symbolic computations", International Tables for
  Crystallography Vol. B (2010), Appendix A1.4.2, Table A1.4.2.7.
- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, sections 2.2 and 3.1
  (NCMAT v1 structure input and its validation).
- `docs/ncmat_doc.md` (the @CELL, @SPACEGROUP, @ATOMPOSITIONS and
  @ATOMSITES sections), `docs/devel_site_occupancy_roadmap.md`.
