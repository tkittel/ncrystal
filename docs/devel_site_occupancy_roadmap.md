# Roadmap: partial site occupancy

Status: **in progress for NCrystal 5.0.0** (branch `tk_occu`). This note
records the agreed design, updated after the 2026-11 design discussion.

## Motivation and current state

Neither the NCMAT format nor the C++ core supports partial site occupancy. The
Python `NCMATComposer` (and hence `ncrystal_cif2ncmat`) emulates an occupancy
*o* with an `@ATOMDB` mixture: *o* of the real element plus (1-*o*) of a
hijacked isotope (`Og299`, `Og298`, ...) acting as a "sterile" atom with the
mass of the element and zero cross sections (see
`allow_siteoccu_ncmatv5_hack` in `_ncmatimpl.py`).

- Scattering is correct: the coherent scattering length becomes o*b, and the
  mixture variance term 4pi*o(1-o)*b^2 is exactly the disorder ("Laue
  diffuse") scattering of randomly placed vacancies.
- Bookkeeping is wrong: vacancies count as atoms with full mass (density too
  high), the composition contains a fake element (breaks Geant4/OpenMC and
  other downstream users), and the output is a hack other tools can not
  interpret.

## Agreed design

1. **Occupancy is a property of the atom role (`AtomInfo`)**, not of the
   individual crystallographic sites, nor of `AtomData`. The `AtomData`
   (possibly a mixture of elements/isotopes) keeps meaning "what sits on the
   site when it is occupied". Sites of the same element with different
   occupancies become different roles (via `@ATOMDB` aliases, as the composer
   already does), which is natural since such sites most likely also differ
   in dynamics. A vacancy pseudo-component in `AtomData`/`@ATOMDB` was
   rejected (it would e.g. make the averaged mass o*m, which is wrong for the
   dynamics).
   - C++: new `double AtomInfo::occupancy() const` (1 if not specified). The
     placeholder `AtomInfo::Site::occupancy` (added in the 5.0.0 development
     cycle) is removed again, so there is a single source of truth.
   - NCMAT v8: new section `@ATOMOCCUPANCIES`, listing each atom label at
     most once (e.g. `Al 0.8`, fractions like `2/3` allowed, valid range
     0<o<=1, default 1). Older format versions reject the section, and
     existing files are unaffected.
2. **Mixed sites (substitution) only via `@ATOMDB` mixtures**, which already
   include the disorder term. A partially occupied mixed site combines a
   mixture role with an `@ATOMOCCUPANCIES` entry. Several roles on the same
   positions remain a load error: allowing it would need cross terms between
   roles in the incoherent disorder term (4pi*[sum(o_i*b_i^2)-(sum(o_i*b_i))^2]
   per site). Separate dynamics for the components of a mixture, and
   molecular kernels covering several roles, would instead relax the
   one-to-one `AtomInfo`/`DynamicInfo` relation; that is out of scope for
   5.0.0.
3. **Counting distinguishes positions from (average) atoms:**
   `StructureInfo::n_positions` and `AtomInfo::numberOfPositions()` count
   positions, while `StructureInfo::n_atoms_average` and
   `AtomInfo::averageNumberOfAtoms()` are the occupancy-weighted counts
   (positions*o). They replace `StructureInfo::n_atoms` and
   `AtomInfo::numberPerUnitCell()`, so code assuming that a count of atoms
   relates to e.g. the number density must choose explicitly. Composition,
   density, number density and the consistency check of `@DYNINFO`
   fractions use the average counts.
4. **Anisotropic displacements are not part of NCMAT v8**; they will come with
   a future format version (v9), rather than reserving syntax now. In C++,
   `AtomInfo::Site::uij` stays as a reserved per-site placeholder (the tensor
   refers to the site's representative position), and `InfoBuilder` rejects
   it being set until its meaning is defined. The intention is that it will
   only provide the shape of the displacements, the magnitude coming from the
   role's msd.

## Physics

For a role with occupancy o, coherent scattering length b and Debye-Waller
factor exp(-W):

| Quantity | Treatment |
|---|---|
| Structure factor | Absorb o into the per-role coherent scattering lengths used by `NCFillHKL.cc` (the `csl` arrays), so no per-position weights or special fast paths are needed. |
| Atoms per unit cell | positions*o, used for composition, density and number density. |
| Incoherent elastic | o*sigma_inc*exp(-2W), plus the disorder term 4pi*o(1-o)*b^2*exp(-2W), with b the (mixture-averaged) coherent scattering length of the `AtomData`. |
| Inelastic, absorption | Unchanged per real atom (they follow the occupancy-weighted composition). |

The disorder term assumes randomly placed vacancies (no short-range order),
the standard approximation and what the current workaround gives.

## Interfaces

- Public C++ headers change deliberately (5.0.0 is an ABI break anyway, and
  `NCRYSTAL_ALLOW_ABI_BREAKAGE` is not used): `AtomInfo` gains the
  occupancy, `AtomInfo::Site` loses its occupancy placeholder.
- C API: the new `ncrystal_info_cquery` gives keyed access to e.g.
  `n_positions`, `n_atoms_average` and `cell_numberdensity` (the latter is
  what McStas components need), and `ncrystal_info_getstructure` is
  obsolete. Occupancies are exposed through the JSON query path
  (`ncrystal_info_jsonquery`, `Info._query`), and from there in the Python
  `Info.AtomInfo` objects (`numberOfPositions`, `averageNumberOfAtoms`,
  `occupancy`; `count` is an obsolete alias).

## Plan

Each step its own commit(s), with tests:

1. This design update.
2. C++ data model: `AtomInfo::occupancy()`, removal of
   `Site::occupancy`, `InfoBuilder` validation (range, reserved `uij`), and
   occupancy-weighted composition and density.
3. Physics: o*b in `NCFillHKL.cc`, and o*sigma_inc plus the disorder term in
   the incoherent elastic model.
4. NCMAT v8 parsing of `@ATOMOCCUPANCIES`, with error tests.
5. JSON query and Python exposure.
6. Validation tests (see below).
7. Composer and `cif2ncmat`: write native occupancies (v8 only when needed),
   read them back in `from_ncmat`/`from_info`, and remove the `Og`
   workaround.
8. Documentation: `docs/ncmat_doc.md` (v8 section) and CHANGELOG input.

## Validation

- The same structure written with the new section and with the old
  workaround must give identical coherent and incoherent elastic cross
  sections, while density and composition differ (the new ones correct).
- Per-plane F^2 checked against an independent calculation including
  occupancy (reusing `tests/pypath/NCTestUtils/latticecells.py`).
- Limits: o=1 bit-identical to no section; continuity for small o;
  combination with `@ATOMDB` mixtures.
- `cif2ncmat` round trips of partially occupied CIF data (occupancy support
  to be added to `tests/pypath/NCTestUtils/cifgen.py`): composition and
  density must match the CIF formula, without fake elements.
- Errors: o<=0 or o>1, the section in pre-v8 files, unknown or repeated
  labels, `@DYNINFO` fractions inconsistent with the weighted composition.

## Related items

- Anisotropic (per-site) displacements: see design point 4.
- The single crystal d-spacing spread (`delta_d`) is only a stub which throws
  if used; implementing it needs a new broadening model.
