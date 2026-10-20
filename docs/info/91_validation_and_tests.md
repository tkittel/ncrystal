# Validation evidence and the test suite

## Published validation

- CPC 2020 (framework paper): hkl lists and positions cross-checked
  against nxslib/SgInfo for the whole library; McStas-simulated powder
  patterns refined with GSAS-II (4 materials, R < 1.9%) and idealised
  patterns with FullProf (17 materials, R < 2.7%); measured $|F|^2$ of
  monoclinic SiLu$_2$O$_5$; EXFOR total cross sections of Mo, Ni, Zn, UO$_2$,
  Sn; Geant4 absorption identical to Geant4's own.
- JCP 2019 (sampling paper): kernel-sampled cross sections within 0.05%
  (0.1 eV) and 0.2% (1 eV) of the analytic free gas; sampled
  energy-transfer distributions statistically identical to the analytic
  ones; a free-gas water sphere within a fraction of a percent of MCNP6's
  classical treatment (and exposing MCNP6 artefacts).
- CPC 2021 (elastic paper): mosaic rocking curves against 100-digit
  references; sampled azimuthal distributions; 1e10-step consistency walks;
  isotropic illumination reproducing the powder cross section to 1e-3;
  Ge-511 reflectivity against Sears and Wuttke; layered-crystal model
  against brute-force rotation averages; incoherent elastic against EXFOR
  (V, Ni).
- 2.0.0 announcement: BeO kernels from VDOS at eight temperatures against
  ENDF/B-VIII.0 kernels; Debye-model fallback compared.
- 3.1.0 announcement: UCN production cross sections in superfluid helium.
- Plugin and user publications: beryllium transmission with extinction
  (Xu 2025), nanodiamond SANS (Rizzi 2023), NJOY+NCrystal libraries (Ramic
  2022), MCNP extension (Grammer 2026), graphite criticality benchmarks with
  OpenMC (Jalil 2026), HighNESS moderator studies.

## In-repository validation (5.0.0 development)

- Reference-value tests: `sabxs` (cross sections of many materials over
  all `vdoslux`/`knllux` levels against stored references with
  level-dependent tolerances; `long_sabxs_*` for low temperatures, light and
  heavy atoms, high Debye temperatures, `knllux` 6), `vdosgn.py` and
  `kinutils.py` against mpmath, `app_gnconvnoisefloor`, `app_gnerange`,
  `app_sctext` (SCT against quadrature and the free-gas limit),
  `sabanalyse.py` (kernel analyser against VDOS truth),
  `sabcelleval.py`/`sabcellsample.py`/`sabsample.py` (processor against
  `SABRefEval`), `app_symchecks`, `app_eqrefltables`, `app_cellconstraints`,
  `latticecells.py` (cells against metric tensors and spglib), `sgsym.py`
  (against spglib and gemmi), `cif2ncmatgen.py` (synthetic CIFs in all 230
  groups).
- Reproducibility: the same references on Linux, macOS, Windows, across
  compilers, with forced FMA contraction, in 32-bit and wasm32 developer
  builds; `app_fmaprobe`, `app_fmamath`, `app_ncerf`, `app_stablemathfct`.
- Regression: every stdlib file dumped and cross-section-sampled with
  `.log` references; MiniMC multiple-scattering scenarios (`mmc_*.py`) for
  powders, single crystals, SANS, phases; `ncrystal-verify` runs the
  pure-Python subset on installed builds.

## The test suite in numbers

82 C++ test applications under `tests/src/` and 133 Python test scripts
under `tests/scripts/` (30 of them `long_` tests enabled with `--long`),
with shared utilities in `tests/pypath/NCTestUtils/` (reference samplers,
mpmath contexts, histogram and statistics helpers, the `sabxs` driver, the
MiniMC reference implementation, CIF generators). All C++ builds with
`-Wall -Wextra -pedantic -Werror`; CI covers many compilers and platforms
and Python 3.9-3.14.

## Gaps

- No dedicated publication yet on the inelastic model (VDOS expansion,
  kernel processing, SCT, kernel analysis).
- No measured-data benchmark of the incoherent approximation for
  amorphous solids or of the SANS framework in the repository.
- Displacements estimated from kernels are provisional pending comparisons
  with evaluated reference data.
