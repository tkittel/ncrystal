# The NCMAT format from the physics side

## Summary

NCMAT is NCrystal's own text format for material data: a short, human-
readable, ASCII file (UTF-8 only in comments) with `@`-sections describing a
crystal structure, the atomic dynamics per element, composition tweaks,
density, temperature, state of matter, secondary phases and plugin-specific
custom data. It was designed so that a crystalline material needs only its
unit cell and one dynamic parameter (v1, 2019), and grew by one version per
physics capability: kernels and VDOS (v2), isotopes and mixtures (v3),
VDOS-derived displacements (v4), amorphous solids and states of matter (v5),
multiphase (v6), temperature declarations (v7) and, for 5.0.0, explicit
displacement/effective-temperature data, kernels in crystals, space-group
settings and sites (v8). Every version remains readable by newer releases,
and files from the standard library can be read by all releases back to the
version that introduced their format. This chapter maps the sections to the
physics they drive; the authoritative specification is `docs/ncmat_doc.md`.

Code: `ncrystal_core/src/ncmat/`.

## Sections and the physics they feed

| Section | Since | Drives |
|---|---|---|
| `@CELL` lengths/angles (`cubic`, `!!` from v4) | v1 | Lattice, volume, density, d-spacings (chapters 10, 12) |
| `@SPACEGROUP` number, `:setting` from v8 | v1 | Symmetry-equivalent reflections, orbits, structure validation (chapter 11) |
| `@ATOMPOSITIONS` | v1 | Structure factors, composition fractions (chapters 10, 12) |
| `@ATOMSITES` (one site per orbit) | v8 | Same, expanded by the setting's symmetry |
| `@DEBYETEMPERATURE` global or per element (removed in v8; not needed for `vdos` elements from v4; global disallowed from v4) | v1 | Debye-Waller factors and the Debye fallback dynamics (chapters 31, 35) |
| `@DYNINFO` `element`, `fraction`, `type` sterile/freegas/scatknl/vdos/vdosdebye, plus type-specific fields; v8 `msd`, `debye_temp`, `effective_temperature` | v2 | Inelastic model per element (chapters 30-35), displacements (35) |
| `@DENSITY` for non-crystals | v2 | Macroscopic conversion (chapter 42) |
| `@ATOMDB` constants, isotopes, mixtures, `nodefaults` | v3 | Scattering lengths, cross sections, masses (chapter 43) |
| `@CUSTOM_<NAME>` free-form lines | v3 | Plugins (chapter 53); the core's `@CUSTOM_HARDSPHERESANS` (chapter 41); `@CUSTOM_UNOFFICIALHACKS` for developer switches |
| `@STATEOFMATTER` solid/liquid/gas (mandatory for non-crystals in v8) | v5 | Elastic physics of amorphous solids (chapter 24) |
| `@OTHERPHASES` fraction + cfg-string lines | v6 | Multiphase materials (chapter 42) |
| `@TEMPERATURE` default or locked value | v7 | Material temperature (chapter 42) |
| `NCRYSTALMATCFG[...]` in a comment | any | Embedded default cfg parameters (e.g. `lcaxis` for graphite) |

Other format rules with physics consequences: element names are symbols,
`El<A>` isotopes, `D`/`T` aliases, or mixture labels; `fraction`s may be
rationals (`2/3`); the first VDOS point must be at least 1e-5 eV; kernel
grids have at most 65534 points; `temperature` of a `scatknl` section locks
the material temperature; `inelas` cfg overrides can replace the dynamics
declared in the file.

## Rules linking sections

- Crystals (v2+): either @DYNINFO for every element or none (then
  `vdosdebye` from the Debye temperatures); from v5 crystals without
  @DEBYETEMPERATURE need all sections `vdos`/`vdosdebye`, from v8 also
  `scatknl` with usable displacement data, and never `freegas`/`sterile`.
- Fractions in @DYNINFO must agree with the cell content (snapped to the
  exact values); non-crystals need @DENSITY; the `element` labels must be
  resolvable through the atom database or @ATOMDB.
- `scatknl` temperatures and @TEMPERATURE must agree; `effective_temperature`
  must be at least the kernel temperature; `msd ... at_temperature` needs
  the temperature unless the section has one.
- The state of matter is implied solid by a cell or VDOS-type dynamics and
  must then not be declared otherwise.
- Version policy: a file may only use what its declared version allows (so
  that older readers fail cleanly), and NCrystal validates the full file
  before building anything.

## What NCMAT cannot (yet) express

- $S(q,\omega)$ kernels (the `sqw` field names are reserved in the
  specification but not implemented), partial site occupancy, anisotropic
  displacements, per-site VDOS, texture or extinction parameters (plugins
  use `@CUSTOM_` sections), $I(Q)$ or $S(Q)$ curves for SANS or coherent
  elastic scattering of amorphous materials, resonance data. The 5.0.0
  review list includes SANS/$I(Q)$ and `sqw` in NCMAT.

## Tooling around the format

`NCMATComposer` writes files programmatically (chapter 60; it writes v8
with sites and settings for crystals, otherwise the lowest version that
holds the content), `ncrystal_cif2ncmat`, `ncrystal_hfg2ncmat`,
`ncrystal_vdos2ncmat`, `ncrystal_endf2ncmat` produce them, the parsed
content is available as JSON through the C and Python APIs,
`ncrystal_ncmat2cpp` embeds files in C++ (the standard library is embedded
this way), `nctool --dump` shows the resulting `Info`, and
`ncrystal browse` lists and searches files (chapter 51). The `NCrystal.ncmat`
and `NCrystal.datasrc` Python modules load text from files, strings or
in-memory registrations (chapter 52).

## Validation

- `tests/scripts/ncmatparser.py` (hundreds of accept/reject cases per
  version incl. the v8 keywords), `ncmatmisc.py`, `ncmatsettings.py`,
  `ncmatsites.py`, `ncmatv8msd.py`, `ncmatphases.py`, `composer*.py`, and
  the `.log`-checked dumps of every stdlib file.
- The `ncrystal-verify` package runs the pure-Python subset of these tests
  on installed NCrystal builds.

## New or changed in 5.0.0

NCMAT v8 (chapter 35 and 10 for details): `msd`/`debye_temp`/
`effective_temperature` keywords, `scatknl` in crystals with automatic or
explicit displacements, mandatory @STATEOFMATTER for non-crystals,
@DEBYETEMPERATURE removed, `freegas`/`sterile` banned in crystals,
@SPACEGROUP setting codes with the strict policy for multi-setting numbers,
@ATOMSITES, and rhombohedral axes for R groups (also in pre-v8 files, where
the setting is deduced). The composer emits v8 for crystals and for `msd`
data. The specification text is in `docs/ncmat_doc.md` ("The NCMAT v8
format").

## References

- `docs/ncmat_doc.md` (the specification, per version).
- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, section 3.1 (NCMAT v1).
- NCrystal 2.0.0 (v2), 2.1.0 (v3), 2.6.0 (v4), 2.7.0 (v5), 3.0.0 (v6),
  3.2.0 (v7) release notes in `CHANGELOG`.
