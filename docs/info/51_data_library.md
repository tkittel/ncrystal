# The standard data library

## Summary

NCrystal ships 134 NCMAT files (133 materials plus `void.ncmat`), embedded
in the library by default: 110 crystals, 16 amorphous solids, 5 noble gases
at STP and 2 liquids (light and heavy water). Almost all crystals and
amorphous solids carry phonon densities of states, so their inelastic
physics is modelled at any temperature (chapter 31); the water files are the
only ones with direct kernels. Every file documents its provenance in its
header (structure source, VDOS source, validation remarks), which the
`ncrystal browse` tool exposes for searching. More materials come from
plugins (chapter 53) and the `ncrystal-extra` repository, and new files are
made with the tools of chapter 60.

Location: `data/` in the repository, available as the `stdlib::` data
source, installed on disk or embedded in the library.

## Content

- Crystals: elements (metals, semiconductors, graphite, diamond, ...) and
  compounds (oxides, carbides, nitrides, fluorides, halides, hydrides and
  deuterides, scintillators, minerals, ...), with VDOS curves for all but
  three (`Ca_sg229`, `CaSiO3` and `B4C` use `vdosdebye`). Naming is
  `<Formula>_sg<N>[_<Name>].ncmat`.
- Amorphous solids (16, added in 2.7.0 with the state-of-matter concept):
  acrylic glass, epoxy, Kapton, four nylons, PEEK, PVC, polycarbonate,
  polyester, polyethylene, polylactide, polypropylene, polystyrene, rubber,
  most generated with `ncrystal_hfg2ncmat` from hydrogen functional groups
  (chapter 60) with `@STATEOFMATTER solid` and `@DENSITY`.
- Gases: Ar, He, Kr, Ne, Xe as free gases at STP (densities hard-coded for
  273.15 K and 1 atm; use `gasmix::` for other conditions, chapter 42).
- Liquids: `LiquidWaterH2O_T293.6K.ncmat` (H kernel from the CAB model via
  ENDF/B-VIII.0, O as free gas) and `LiquidHeavyWaterD2O_T293.6K.ncmat`;
  now extended above their kernel range with the SCT model at an
  analyser-estimated effective temperature of about 1100 K (chapter 33).
  More temperatures are in the `WaterData` plugin.
- `void.ncmat`: an empty material for porosity in multiphase materials.

## Provenance of the physics data

- Crystal structures: mostly converted from the Crystallography Open
  Database (COD) or the Materials Project with `ncrystal_onlinedb2ncmat`
  / `cif2ncmat`, with the COD or MP identifiers and literature references in
  the header, and often a note on which other database entries are
  compatible.
- VDOS curves: from first-principles phonon calculations (mostly published
  Phonopy results), from digitised published figures, or from existing
  evaluations (e.g. ENDF/B-VIII.0 for Be, graphite and acrylic glass
  hydrogen), in many cases cross-checked against other sources. Each file
  header states the origin of its curve and any validation remarks.
- Cross sections, masses: the atom database (chapter 43).
- Densities of amorphous solids and gases: literature values quoted in the
  headers.

## Validation

- Every file is loaded and dumped in the test suite, with reference cross
  sections, so any change in physics or in a file shows up as a diff.
- The data-library wiki page carries per-material plots and validation
  remarks (cross sections versus EXFOR where available; e.g. the CPC 2020
  paper compared Mo, Ni, Zn, UO$_2$, Sn with EXFOR, and refined simulated
  powder patterns of 17 materials).
- Structures are checked against their space groups at load time (chapter
  10); the 5.0.0 symmetry work found and corrected one rounded coordinate
  in beta quartz and two wrong densities from the lattice-matrix bug
  (boric acid, wollastonite).

## Approximations and limitations

- One VDOS per element, at a single reference temperature in the header
  but used at any temperature; no anharmonic or temperature-dependent
  lattice parameters (use `strain` or separate files).
- Digitised VDOS curves inherit the resolution of the published figures.
- The library aims at materials relevant to neutron instruments and
  shielding (structural metals, moderators, reflectors, shielding,
  detector and sample-environment materials, common minerals), not at
  completeness.
- Gases are free gases: no molecular dynamics.

## New or changed in 5.0.0

- Boric acid and wollastonite densities corrected (issue #383); beta
  quartz coordinates symmetrised to exact fractions.
- The water files gain SCT extension with estimated effective temperatures
  (physics change above the kernel range), and all VDOS files get the 5.0.0
  expansion and kernel processing (sub-percent changes where the legacy
  processing had artefacts).
- `ncrystal browse` and `NCrystal.browse` for finding materials by name,
  header text and physics properties (`docs/data_browsing.md`).

## References

- The NCrystal wiki "Data-library" page (per-material table with plots).
- A. Togo, I. Tanaka, "First principles phonon calculations in materials
  science", Scripta Materialia 108 (2015) 1 (the phonon database behind
  many VDOS curves).
- ENDF/B-VIII.0 thermal scattering sublibrary (LEAPR inputs used for some
  curves and for the water kernels, the latter via the CAB model of
  Marquez Damian et al., Ann. Nucl. Energy 65 (2014) 280).
