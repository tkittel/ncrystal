# NCrystal physics capabilities: source notes

This directory collects, in one place and in a form that does not depend on
any particular output medium, what NCrystal can model and how. It is written
for NCrystal release 5.0.0 and describes both what is there in absolute terms
and what is new in that release. It is meant as the basis from which other
documentation is produced later (website reference pages, user-facing
summaries, papers, and slides), so it favours completeness and precision over
brevity, and it states where each statement comes from (code, tests, papers).

Status: first complete draft of all chapters, written 2026-10-20 on the `tk_volatile` branch. The
chapters are written one at a time; the table below says which exist. All
commits touching this directory carry the subject prefix
`FIXMECLAUDEINFODOC:` and touch nothing outside `docs/info/`.

## Conventions used in the chapters

- Every chapter has the same shape: a short summary, the physics (with the
  governing formulas), the inputs (NCMAT sections and cfg-string parameters),
  the implementation (files, classes, algorithms, numerical choices), how it
  is validated, the approximations and known limitations, a dated "New or
  changed in 5.0.0" section, and references.
- Subsections marked "(derivation)" give the full argument and can be dropped
  when a shorter text is produced.
- The "Approximations and limitations" sections are kept frank, one per
  chapter; in user-facing renderings (website) they are meant to appear at
  the bottom of the chapter in a highlighted box titled "Possible future
  improvements:" rather than as a list of shortcomings.
- Formulas are written in LaTeX notation between `$` signs, in plain ASCII.
  Units are eV, angstrom (Aa), kelvin and barn unless stated.
- Code locations are given as paths relative to the repository root, with
  class and function names, so that statements can be checked against the
  source. The source is the primary reference: where the papers and the code
  differ, the chapters describe the code.
- Version statements ("since 2.0.0", "new in 5.0.0") come from the CHANGELOG,
  the release announcements on the wiki, and the draft change notes in
  `docs/changes_for_changelog.txt` and `docs/changes_for_changelog2.txt`.

## Chapters

| File | Subject | Status |
|---|---|---|
| `00_overview.md` | Scope of NCrystal, the process model, material classes | written |
| `01_conventions.md` | Kinematic variables, normalisations, units | written |
| `10_crystal_structure.md` | Unit cells, settings, sites, structure checks | written |
| `11_space_group_symmetry.md` | The 530 settings, symmetry operations, orbits | written |
| `12_hkl_and_structure_factors.md` | hkl lists, structure factors, d-spacings | written |
| `20_powder_bragg.md` | Bragg diffraction in powders | written |
| `21_single_crystal_mosaic.md` | Single crystals with Gaussian mosaicity | written |
| `22_layered_crystals.md` | Layered crystals (pyrolytic graphite) | written |
| `23_incoherent_elastic.md` | Incoherent elastic scattering | written |
| `24_amorphous_solids.md` | Elastic scattering in amorphous solids | written |
| `30_dynamics_overview.md` | Atomic dynamics and inelastic scattering: overview | written |
| `31_vdos_to_sab.md` | From a phonon density of states to $S(\alpha,\beta)$ | written |
| `32_sab_cross_sections_and_sampling.md` | Cross sections and sampling from kernels | written |
| `33_kernel_extension_sct.md` | Beyond the kernel: free gas, SCT, kernel analysis | written |
| `34_free_gas.md` | The free-gas model | written |
| `35_msd_debye_teff.md` | Displacements, Debye temperatures, effective temperatures | written |
| `40_absorption.md` | Absorption | written |
| `41_sans.md` | Small-angle scattering | written |
| `42_multiphase_and_density.md` | Multiphase materials, density, temperature | written |
| `43_atom_data.md` | Atomic data | written |
| `50_ncmat_format_physics.md` | The NCMAT format from the physics side | written |
| `51_data_library.md` | The standard data library | written |
| `52_other_inputs.md` | Other inputs: .laz/.lau, quick materials, in-memory data | written |
| `53_plugins.md` | Physics provided by plugins | written |
| `60_material_creation_tools.md` | NCMATComposer, cif2ncmat, hfg2ncmat, PhononDOSAnalyser | written |
| `61_analysis_and_export.md` | ncmat2endf, filter tables, MiniMC | written |
| `70_numerics.md` | Numerics and cross-platform reproducibility | written |
| `90_new_in_5.0.0.md` | What is new in 5.0.0 | written |
| `91_validation_and_tests.md` | Validation evidence and the test suite | written |
| `92_references.md` | Bibliography | written |
| `93_glossary.md` | Glossary | written |

The `survey/` subdirectory holds the raw, unedited inventories the chapters
are built from; see its own README.
