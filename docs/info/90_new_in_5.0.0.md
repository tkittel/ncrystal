# What is new in NCrystal 5.0.0

## Framing

NCrystal 4.4.6 (June 2026) was already a mature library: kinematically
exact Bragg diffraction in powders and mosaic and layered single crystals
with models published and validated to high precision, self-consistent
incoherent elastic scattering, inelastic scattering from on-the-fly
multi-phonon expansions of phonon spectra at any temperature with exact
rejection sampling, free-gas extension, dedicated UCN production, a 134-file
validated data library, multiphase materials, plugins, and integrations in
OpenMC, Geant4 and McStas. Release 5.0.0 is the first ABI-breaking release
since 4.0.0 and uses the opportunity for a generation change in several
places. The two headline features are the full space-group symmetry
support in crystallography (chapters 10-12) and the new generation of
inelastic modelling (the VDOS expansion and kernel processing of chapters
31-32), with the short-collision-time extension (chapter 33) and the NCMAT
v8 format (chapter 50) as the other major additions. The list below is
organised by what a user of the physics gains, with pointers to the
chapters that carry the details. The curated changelog input is in
`docs/changes_for_changelog.txt` and `docs/changes_for_changelog2.txt`.

## Crystallography and Bragg diffraction (chapters 10-12, 20-22)

- Full space-group symmetry in the library: all 530 standard settings,
  exact symmetry operations, orbits and sites, structure validation
  against the declared symmetry with informative messages, cell
  constraints derived from the operations, equivalent reflections from the
  operations (fixing Laue class $\bar31m$), rhombohedral axes and
  non-default settings (monoclinic unique axis c, axis permutations, origin
  choice 2).
- NCMAT v8: explicit settings in `@SPACEGROUP`, `@ATOMSITES` (one site per
  orbit), `NCMATComposer` and `cif2ncmat` writing them.
- Exact, bit-identical d-spacings, structure-factor symmetry self-check,
  all species included in $|F|^2$, corrected lattice matrices for
  triclinic/rhombohedral cells (issue #383), better family merging without
  space group.
- Sites and the space-group setting on the `Info` object, with
  site/position utility queries.

## Inelastic scattering (chapters 30-35)

- A new generation of the VDOS expansion (the default): exact detailed
  balance at low temperatures, robust handling of narrow spectral structure
  and of numerical noise, resource-based limits with a fixed 5 eV target,
  errors instead of silent range reductions.
- A new generation of the kernel processing (the default): exact cell-wise
  phase-space integrals, automatic energy grids with monotone cubic
  interpolation, cell-based rejection sampling, caching per kernel.
- Short-collision-time extension above the kernel range with effective
  temperatures from the VDOS, from the NCMAT v8 `effective_temperature`
  keyword, or estimated from kernels; the free gas remains as the
  $T_{eff} = T$ special case.
- Effective temperature and mean-squared displacement estimated from
  tabulated kernels through the moment sum rules of the Gaussian
  approximation; hence elastic physics for direct-kernel solids (and e.g.
  the stdlib water kernels extended with SCT at $T_{eff} \approx 1100$ K).
- NCMAT v8 dynamics keywords (`msd`, `debye_temp`, `effective_temperature`
  on `scatknl`; `msd ... at_temperature` on `vdosdebye`), kernels in
  crystals, state-of-matter rules; `@DEBYETEMPERATURE` removed.
- Pending for the release: the switch of the default luxury values to the
  new generation (decided); `ucnmode` only with legacy processing and
  scheduled for removal afterwards; revisiting the free-gas sampler against
  the OpenMC/Gelbard scheme; SANS and $S(q,\omega)$ in NCMAT.

## Numerics and performance (chapter 70)

- Cross-platform bit-reproducibility of the VDOS-to-kernel pipeline,
  structure factors and mosaic setup (portable elementary functions,
  explicit fused multiply-adds, noise-robust discrete decisions).
- Runtime FMA dispatch for portably built binaries, branchless searches (up
  to 24% faster Bragg sampling, 15% kernel sampling), faster hkl generation
  (10-25%).
- 32-bit and WebAssembly support, per-thread C API error state, C API
  robustness fixes.

## Formats, data, tools (chapters 50-61)

- NCMAT v8 (chapter 50), `.laz`/`.lau` settings header, composer writing
  v8, synthetic-CIF tests in all space groups, full-precision cells.
- Data library: two densities corrected, beta quartz symmetrised, water
  files SCT-extended.
- `ncrystal browse` / `NCrystal.browse` for data discovery incl. physics
  property selection; `ncrystal query` and the extended JSON query system.
- Filter tables (C++, C API, Python) for attenuation-only GPU components.

## ABI and API

- ABI-breaking: the space-group setting and the sites on the `Info`
  object, the public setting type, process-interface additions
  (multi-energy cross-section evaluation, combined cross-section and
  sampling calls, a pure-elastic flag), and the pending discussion items
  (site occupancy, anisotropic displacements, effective-temperature data on
  direct kernels, origin-shifted cells).
- Obsolete and removed: global Debye temperature accessors, `packfact`,
  `bkgd`; scheduled: `ucnmode`, the legacy luxury numbering.

## Compatibility

All NCMAT v1-v7 files load with identical physics except where fixes
apply (listed above: #383 densities, $\bar31m$ structure factors, damped
species, beta quartz), and except for the inelastic physics, where the
next-gen expansion and kernel processing are the 5.0.0 defaults (changes at
the sub-percent level where the legacy processing had artefacts, documented
by the `sabxs` reference values) and all kernels are extended above their
energy range with the SCT model instead of the free gas at the material
temperature (VDOS-expanded kernels always; direct kernels whenever an
effective temperature is given or can be estimated, which is the case for
the stdlib water files). Files that do not have the symmetry of
their declared space group are now rejected.
