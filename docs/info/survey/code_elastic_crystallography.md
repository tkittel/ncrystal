<!-- Unedited AI-generated survey inventory (2026-10-20), kept as raw source notes for docs/info. Not reviewed. To be removed from history before merging to develop. -->

# Inventory of elastic and crystallographic physics in NCrystal (branch tk_volatile, towards 5.0.0)

There is no `v4.4.6` tag. As the baseline I used commit **c4bfa91e9** ("Update version numbers to 4.4.7"). There are 446 commits from there to HEAD.

Commits per area since the baseline (src / headers):

| Area | src commits | header commits | Assessment |
|---|---|---|---|
| powderbragg | 1 | 0 | essentially stable |
| scbragg | 0 | 0 | essentially stable |
| lcbragg | 0 | 0 | essentially stable |
| bkgdextcurve | 0 | 0 | essentially stable |
| elincscatter | 2 | 0 | small change |
| cfgutils | 0 | 3 | small change |
| filter | 1 | 1 | entirely new |
| lazlau | 7 | - | heavily changed |
| infobld | 15 | - | heavily changed |
| sgsym | 12 | 13 | entirely new |
| extd_utils | 22 | 11 | heavily changed |

- The sgsym component is entirely new (about 2.3k lines).
- In infobld, `NCInfoBldStructure`, `NCInfoBldComposition` and `NCInfoBldOverrides` were split out of `NCInfoBuilder.cc`.
- In extd_utils, `NCFillHKL.cc` was largely rewritten (+1204/-... lines). `NCFillHKL_FMA.hh` and `NCFillHKL_WINFMA.cc` are new.

The curated change notes are in `/workspace/docs/changes_for_changelog.txt` (section "Crystallography and Bragg diffraction") and `/workspace/docs/changes_for_changelog2.txt` (sgsym, KSqCalc, sites, NCMAT v8). Developer notes are in `/workspace/docs/devel_sgsym_plans.md` and `/workspace/docs/devel_site_occupancy_roadmap.md`. The NCMAT v8 spec is in `/workspace/docs/ncmat_doc.md`, lines 881-1080.

---

## 1. How the standard scatter factory puts elastic components together
- **File:** `/workspace/ncrystal_core/src/stdscatfactory/NCStdScatFact.cc`, class `StdScatFact` (priority 100).
- **Crystalline materials:**
  - With `incoh_elas` on, `ElIncScatter(info)` is added if `ElIncScatter::hasSufficientInfo`.
  - With `coh_elas` on and no orientation given, `PowderBragg(info)` is used. The comment says "Layered polycrystals get same treatment as unlayered polycrystals".
  - With `coh_elas` on for a single crystal:
    - If `lcaxis` is set, `LCBragg(info, sco, mos, lcaxis, lcmode, 0, pp, mosprec, 0)` is used.
    - Otherwise `SCBragg(info, sco, mos, 0.0, pp, mosprec, 0)` is used.
    - **sccutoff trick:** the `PlaneProviderWCutOff` wrapper holds back planes with d < `sccutoff` (default 0.4 Aa) and returns them as (d, 2*fsq) pairs. These go to an extra `PowderBragg(StructureInfo, VectDFM)`, i.e. they are treated as having infinite mosaicity.
- **Amorphous solids** (non-crystalline, `StateOfMatter::Solid`): `ElIncScatter` with `use_sigma_incoherent = incoh_elas` and `use_sigma_coherent = coh_elas`, i.e. coherent elastic in the incoherent approximation.
  - The `@CUSTOM_UNOFFICIALHACKS` keyword `no_cohelas_via_incohapprox_for_amorphous_solids` turns the coherent part off.
- **`inelas=external`** adds `BkgdExtCurve` (see section 7).
- **No extinction model, strain broadening, texture or grain size** appears anywhere in the core factory.

---

## 2. Powder Bragg diffraction (`PowderBragg`)
- **Files:**
  - `/workspace/ncrystal_core/include/NCrystal/internal/powderbragg/NCPowderBragg.hh`
  - `/workspace/ncrystal_core/src/powderbragg/NCPowderBragg.cc`
- **Physics:** ideal powder / non-textured polycrystal.
  - The per-atom cross section is sigma(lambda) = lambda^2/(2*V_0*N_atoms) * Sum_{2d>=lambda} d*|F|^2*mult.
  - The code stores `xsectfact = 0.5/(V0*natoms)*wl2ekin(1)` so the cross section is a multiplication by 1/E.
  - Scattering is elastic. The scattering angle comes from sintheta_B = lambda/2d, giving mu = 1 - 2*sin^2theta_B.
  - The header states: "Does not account for texture, grain size, dspacing-deviations and other similar effects."
- **Inputs:**
  - `Info::hklList()` (dspacing, fsquared, multiplicity) and `StructureInfo` (volume, n_atoms).
  - Planes are set by `dcutoff`/`dcutoffup` at Info level and switched on/off by `coh_elas`/`bragg`.
  - Specialised constructors take a list of (d, F^2*mult) pairs directly.
- **Algorithms:**
  - Planes are sorted by d (descending) and merged when |Deltad| < 1e-11.
  - Two arrays are built: `m_2dE` (Bragg-edge energies) and `m_fdm_commul`, a cumulative sum computed with `StableSum`.
  - A plane is picked by binary search on the cumulative array. The search now uses `fastUpperBoundIdx`/`fastLowerBoundIdx` from the new `NCFastSearch.hh` (the only change on this branch).
  - `CachePowderBragg` caches the energy and the last valid plane index.
  - `createMerged` merges plane lists (comment: "TODO: can we do better" about the subtraction).
  - Under `NCRYSTAL_ALLOW_ABI_BREAKAGE` there are also `isPureElasticScatter`, `evalXSAndSampleScatter[Isotropic]` and `evalManyXSIsotropic`.
  - `domain()` has a lower bound at the Bragg threshold.

## 3. Single-crystal Bragg diffraction with Gaussian mosaicity (`SCBragg` + `GaussMos` + `GaussOnSphere`)
- **Files:**
  - `/workspace/ncrystal_core/include/NCrystal/internal/scbragg/NCSCBragg.hh`
  - `/workspace/ncrystal_core/src/scbragg/NCSCBragg.cc`
  - `/workspace/ncrystal_core/include/NCrystal/internal/phys_utils/NCGaussMos.hh` and `src/phys_utils/NCGaussMos.cc`
  - `NCGaussOnSphere.hh` and `.cc`
  - `/workspace/ncrystal_core/include/NCrystal/internal/extd_utils/NCOrientUtils.hh` (`getCrystal2LabRot`)
  - `NCPlaneProvider.hh`
- **Physics:** "Bragg diffraction in single crystals with given orientation and mosaicity, assuming a Gaussian distribution of crystallite orientations."
  - The distribution is a truncated 2D Gaussian on the unit sphere, specified by its FWHM or sigma.
  - The cross section for one plane is Q*(integral of the Gaussian density along a circle).
    - The circle is the cone of normals satisfying the Bragg condition.
    - Q' = lambda^3/(sintheta*costheta) and Q = Q'*xsfact/2, where xsfact = |F|^2/(V_0*N_atoms).
  - Reflection families are grouped by (d, |F|^2).
  - The `delta_d` parameter is accepted, but `GaussMos::setDSpacingSpread` throws: "not actually implemented and debugged fully yet". The factory always passes 0.
- **Configuration:** `mos` (FWHM, in (0, pi/2]), `dir1`, `dir2` (syntax `@crys:` or `@crys_hkl:` ... `@lab:`), `dirtol` (default 1e-4; 180 deg gives a monochromator setup), `mosprec` (default 1e-3, range [1e-7, 0.1]), `sccutoff`.
- **Algorithms:**
  - `GaussOnSphere` uses spline lookup tables to avoid calls to exp, acos and erf.
  - Circle integrals use Romberg integration, or approximation formulas where appropriate.
  - The truncation angle ntrunc is estimated from the precision (`estimateNTruncFromPrec`).
  - Points on the circle are sampled exactly by rejection sampling.
  - Normal and anti-normal are pre-checked together.
  - The energy is rounded to 15 decimals for the cache key (for reproducibility).
  - (d, fsq) keys are stored as integers with 1e-10 precision.
  - The environment variable `NCRYSTAL_DEBUG_GAUSSONSPHERE` prints statistics.
- **On this branch:** no source changes.

## 4. Layered crystals such as pyrolytic graphite (`LCBragg`)
- **Files:**
  - `/workspace/ncrystal_core/include/NCrystal/internal/lcbragg/NCLCBragg.hh`
  - `/workspace/ncrystal_core/src/lcbragg/NCLCBragg.cc`
  - `/workspace/ncrystal_core/include/NCrystal/internal/extd_utils/NCLCUtils.hh` and `src/extd_utils/NCLCUtils.cc`
  - `NCLCRefModels.hh` and `.cc`
- **Physics:** single-crystal mosaic model (via `GaussMos`) where crystallites are additionally rotated uniformly at random around `lcaxis`.
- **Configuration:** `lcaxis` (direct-lattice vector, e.g. "0,0,1"); requires `dir1`, `dir2` and `mos`.
- **`lcmode` options:**
  - `0` (default): the fast, accurate `LCHelper`.
  - `N > 0`: `LCBraggRef`, which averages `SCBragg` over N rotations. It is slow and serves as a reference.
  - `N < 0`: `LCBraggRndmRot`, which re-samples N random rotations at every `crossSection` call. Marked "NOT MT-safe".
- **Algorithms used by `LCHelper`:**
  - `LCPlaneSet` groups normals with the same angle to `lcaxis` and the same d. The comment says "alpha ... is theta_n in the upcoming NCrystal Bragg paper".
  - `LCROIFinder` finds the phi-ranges with non-zero cross section.
  - `LCStdFrame` works in a frame where lcaxis = (0,0,1) and the neutron is at (s3, 0, c3).
  - The cross section is integrated over phi with an adaptive Simpson integrator.
  - Sampling uses MC overlays.
  - The cache key discretises (lambda, c3).
  - `mode=0` asserts `delta_d == 0`.
- **On this branch:** source unchanged (only small edits in `NCLCUtils.cc`).

## 5. Incoherent elastic scattering (`ElIncScatter` / `ElIncXS`)
- **Files:**
  - `/workspace/ncrystal_core/include/NCrystal/internal/elincscatter/NCElIncScatter.hh`
  - `/workspace/ncrystal_core/src/elincscatter/NCElIncScatter.cc`
  - `/workspace/ncrystal_core/include/NCrystal/internal/phys_utils/NCElIncXS.hh` and `.cc`
  - `NCDebyeMSD.hh` and `.cc`
- **Physics:** isotropic harmonic Debye-Waller model, following CPC 2019 (doi 10.1016/j.cpc.2019.07.015), section 2.3, eq. 23 and eq. 40.
  - Differential: dsigma/dOmega = sigma_inc/4pi * exp(-2W), with 2W = Q^2<u^2>.
  - Total: sigma = sigma_inc*(1 - e^{-t})/t, with t = 4k^2<u^2>.
  - mu is sampled from a density proportional to exp(a*mu), a = 2k^2<u^2>.
  - `eval_1mexpmtdivt` is a numerically safe evaluation of (1 - e^{-t})/t.
  - Polyatomic materials use a weighted sum over elements of (msd, sigma*scale).
  - `ElIncScatterCfg` has `scale_factor`, `use_sigma_incoherent` and `use_sigma_coherent`. The last gives the incoherent approximation for the coherent part, used for amorphous solids.
- **Where the MSD comes from:**
  1. `AtomInfo::msd()`, falling back to `debyeIsotropicMSD(T_D, T, M)` (Glauber 1955, eqs. 11-12).
  2. Otherwise the DynInfo list, via `extractMSDFromDynInfo` in `NCDynInfoUtils.cc`:
     - VDOS: `VDOSEval::getMSD`.
     - VDOSDebye: Debye-model MSD.
     - scatknl: the explicit v8 `msd`/`debye_temp` value, or else an estimate from `SABAnalyser::estimateTeffMSD` (only for NCMAT v8 and later).
- **On this branch:**
  - Switched from vdos-specific casts to `extractMSDFromDynInfo`.
  - New `NCRYSTAL_WARN` when components lack an MSD (partial contribution).
  - freegas and sterile components are skipped.
  - The header comment still says "for now scatknl ... will be ignored", which is out of date.

## 6. Structure factors and hkl-list generation (`calculateHKLPlanes` / `FillHKLCfg`)
- **Files:**
  - `/workspace/ncrystal_core/include/NCrystal/internal/extd_utils/NCFillHKL.hh`
  - `/workspace/ncrystal_core/src/extd_utils/NCFillHKL.cc`
  - `NCFillHKL_FMA.hh`, `NCFillHKL_WINFMA.cc`
  - Called from `/workspace/ncrystal_core/src/ncmat/NCLoadNCMAT.cc` (around line 855) through a lazy `HKLListGenFct`.
- **Physics:**
  - |F|^2 = |Sum_species b_coh * exp(-W) * Sum_pos exp(2pii h*r)|^2.
  - W = 1/2*k^2*msd (comment cites Sears, Acta Cryst. A53 (1997) 35) with k^2 = (2pi/d)^2. This is an isotropic Debye-Waller factor per species.
  - Every atom must have an msd.
  - Missing so far: occupancy, anisotropic displacements (Uij) and anomalous scattering.
- **`FillHKLCfg` fields:**
  - `dcutoff` (default 0.5; the NCMAT loader auto-selects 0.1 Aa, or 0.2 Aa if the cell has more than 40 atoms, or 0.5*dcutoffup if that conflicts).
  - `dcutoffup`.
  - `fsquarecut` = 1e-5 barn.
  - `merge_tolerance` = 1e-6, for the no-space-group path. The comment explains why 1e-6 was chosen.
  - `max_buffered_points`.
  - `use_unit_debye_waller_factor`. There is a FIXME: "time to remove forceunitdebyewallerfactor again".
- **Two code paths:**
  - **With a space group:** `calculateHKLPlanesWithSymEqRefl`.
    - Brute-force hkl loop over half of the space.
    - Each hkl is mapped to its representative with `EqRefl`, and `SymHKLSeenTracker` is a bit table of representatives already seen.
    - The bit table is sized to the actual hkl range and limited by `NCRYSTAL_FILLHKL_MEMLIM` (default 20 MB).
    - A loose k^2 preselection uses a cheap matrix kernel.
    - Result type: `HKLInfoType::SymEqvGroup`.
    - **New:** a symmetry self-check compares F^2 for equivalent planes with |h|,|k|,|l| <= 4. A difference above 1% (relative to max(F^2, 1e-4*F^2max)) is an error.
  - **Without a space group:** `calculateHKLPlanesNoSym`.
    - Families are built by (d, F^2) using buffered sort-and-sweep merging (chunks of at most 16 MB).
    - Family values are averages computed as first + mean(member - first), which is order-independent (Sterbenz).
    - Result type: `ExplicitHKLs`.
- **Numerical tricks:**
  - A cheap upper bound on F^2 (all terms in phase) gives early rejection.
  - `stable_exp` is used instead of `std::exp` (libm-free, reproducible).
  - Explicit FMA dispatch (`NCRYSTAL_FMADISPATCH_*`) for the k^2 and F^2 kernels.
  - `estimateHKLRange` uses the exact bound |h| <= a/dcutoff with a 0.1% margin.
  - `estimateNHKLPoints` = (2pi/3)*V*(kmax^3 - kmin^3).
  - A guard throws above 1e6 families.
  - Environment variables `FILLHKL_SELECTHKL=h,k,l` and `FILLHKL_IGNOREFSQCUT`.
- **On this branch (from the changelog notes):**
  - Heavily damped species are no longer skipped (this gave up to about 5% bias on weak planes).
  - The family bookkeeping was rewritten.
  - The symmetry self-check was added.
  - KSqCalc is now used (next entry).
  - The Laue class -31m fix (space groups 149, 151, 153, 157, 159, 162, 163) via the new `EqRefl`.

### KSqCalc (exact d-spacings)
- **Files:** `/workspace/ncrystal_core/include/NCrystal/internal/utils/NCLatticeUtils.hh` (`class KSqCalc`) and `src/utils/NCLatticeUtils.cc`.
- **What it does:**
  - Computes |k|^2 = (2pi/d)^2 from the reciprocal metric.
  - Integer parts are formed exactly in 64-bit integers, and terms with identical factors share one integer sum.
  - Exact sines and cosines are used for 90 deg and 120 deg.
  - Only the FMA chain rounds.
  - As a result, symmetry-equal planes (e.g. 333 and 511 in a cubic cell, 700 and 530 in a hexagonal cell) get bit-identical d on all platforms.
- **Users:** FillHKL, the Lazy loader and `Info::dspacingFromHKL`.
- It must stay identical to `detail_fillhkl_ksq` in `NCFillHKL_FMA.hh`.
- **Related:** hkl lists are now sorted strictly by (d, multiplicity, hkl), without the old 1e-6 Aa tolerance.
- **Related bug fix (issue #383):** the wrong real-space lattice matrix when cos alpha != cos beta*cos gamma (triclinic and rhombohedral cells) was fixed. BoricAcid density changes 2.288 -> 1.562 g/cm^3 and Wollastonite 3.096 -> 2.909 g/cm^3.

### EqRefl (symmetry-equivalent reflections)
- **Files:** `/workspace/ncrystal_core/include/NCrystal/internal/phys_utils/NCEqRefl.hh` and `src/phys_utils/NCEqRefl.cc`.
- **What it does:** half-lists of equivalent (h,k,l) for the Laue class.
  - Now derived at construction from the sgsym rotation parts modulo sign; this is new.
  - The old hand-written `EqReflTables` were moved to `tests/src/app_eqrefltables` as a reference.
  - Can be constructed from an SG number plus a `rhombohedral_axes` flag, or from a `SpaceGroup` setting.
  - Distinguishes -31m from -3m1.
  - Supports rhombohedral axes for R groups.

### Plane providers and orientation
- `/workspace/ncrystal_core/include/NCrystal/internal/extd_utils/NCPlaneProvider.hh`:
  - `PlaneProvider` loops over demi-normals.
  - `createStdPlaneProvider` prefers explicit normals, then expanded hkl, then space group plus multiplicity.
  - `ExpandHKLHelper(const Info&)` "digs out the setting" (changed on this branch).
- `HKLInfoType` values: SymEqvGroup, ExplicitHKLs, ExplicitNormals, Minimal (`/workspace/ncrystal_core/include/NCrystal/interfaces/NCInfoTypes.hh`).
- `NCOrientUtils.hh`: `getCrystal2LabRot(SCOrientation, reci_lattice)`.

### NCInfoUtils (setting helpers, new on this branch)
- **Files:** `/workspace/ncrystal_core/include/NCrystal/internal/extd_utils/NCInfoUtils.hh` and `.cc`.
- `getSpaceGroupFromStructureInfo` throws if the setting is unknown.
- `getSpaceGroupCandidateForHKLFromStructureInfo` returns the stored setting, or else the first setting compatible with the cell, provided all such settings have the same equivalent reflections.
- `completeAndCheckStructureInfoCell`.

## 7. Background from an external cross-section curve (`BkgdExtCurve`)
- **Files:** `/workspace/ncrystal_core/include/NCrystal/internal/bkgdextcurve/NCBkgdExtCurve.hh` and `src/bkgdextcurve/NCBkgdExtCurve.cc`.
- **Physics:** the cross section is `Info::xsectScatNonBragg(E)`, taken from the Info "XSectProvider" (`providesNonBraggXSects`). Scattering is isotropic and elastic.
- **Usage:** only with `inelas=external`; the cfg doc says "only currently supported by .nxs files".
- Unchanged on this branch.

## 8. Space group symmetry (sgsym) - entirely new
- **Files:**
  - `/workspace/ncrystal_core/include/NCrystal/internal/sgsym/NCSpaceGroup.hh`, `NCSGSymmetry.hh`, `NCSymOp.hh`
  - `/workspace/ncrystal_core/src/sgsym/NCSGSymTable.cc` and `.hh`, `NCSGSymmetry.cc`, `NCSpaceGroup.cc`, `NCSymOp.cc`
  - Public type `SpaceGroupSettingIdx` (1..530) in `/workspace/ncrystal_core/include/NCrystal/core/NCTypes.hh`
- **Settings table (`NCSGSymTable`):**
  - Hard-wired copy of the 530 rows of ITVB 2010 Table A1.4.2.7 (Shmueli, Hall & Grosse-Kunstleve): number, choice code, Hermann-Mauguin symbol, Hall symbol.
  - Adds the explicit "abc" choice code for the 46 orthorhombic numbers.
  - Cross-checked against SgInfo, spglib, gemmi and ASE.
- **`SpaceGroup` class:**
  - Strings like "227:2", "166:r", "14:c1", "62:cab"; `SpaceGroupParseMode` Strict or Relaxed.
  - Methods: `hmSymbol`, `hallSymbol`, `crystalSystem` (`SGCrystalSystem`), `settingInfo` (`SGHexFamilyAxes`, `SGUniqueAxis`, `SGCellChoice`, `SGOriginChoice`, `SGAxisPermutation`), `settingsOfNumber`, `listSettings`.
  - FIXME: migrate `isRhombohedralSpaceGroup` and `usesRhombohedralAxes` here from NCLatticeUtils.
- **`SymOp` class:**
  - Exact (R, t). R is one of the 64 allowed rotations (48 signed permutations and 24 of 6/mmm in hexagonal axes).
  - t is stored in units of 1/24, modulo 1.
  - CIF-style strings like "x,-y,z+1/2".
  - Composition, inverse, `rotationType`, `rotationOrder`.
- **`SGSymmetry::get(sg)`:**
  - Generates the full group by parsing the Hall symbol, then caches it lazily for the lifetime of the program.
  - Canonical order of operations.
  - Accessors: `representatives`, `centringVectors`, `latticeSymbol`, `isCentrosymmetric`.
  - `cellConstraints()` returns `SGCellConstraints`, derived from the invariance of the metric tensor under the rotation parts (via a row-space test). It has `check()` and `complete()`, which fill in parameters left as 0.
- **`expandSiteToOrbit`:**
  - Images within 5e-4 are treated as the same position, and the site is symmetrised (e.g. 0.3333 becomes 1/3).
  - Images that differ but lie within 1e-2 cause a BadInput error ("close to special position").
  - These thresholds come from an analysis of COD CIF files.
  - Returns the positions and `siteSymmetryOrder`.
- **`groupPositionsIntoOrbits`** is the inverse operation.
- **`checkStructureSymmetry(number, rhombohedral_axes, positions)`** returns the best setting and the largest deviation. A per-setting overload exists.
- **Other functions:** `haveSameRotationParts`, `settingsCompatibleWithCell`.
- **Limitations:** only the 530 standard settings are supported; non-standard origins are rejected.

## 9. Building Info objects (infobld)
- **Files:**
  - `/workspace/ncrystal_core/include/NCrystal/internal/infobld/NCInfoBuilder.hh`
  - `/workspace/ncrystal_core/src/infobld/NCInfoBuilder.cc`
  - `NCInfoBldStructure.{hh,cc}` (new)
  - `NCInfoBldComposition.{hh,cc}` (new; composition, temperature and density validation; fractions snapped to sum exactly to 1)
  - `NCInfoBldOverrides.{hh,cc}` (new; density and cfg overrides reuse a bounded cache of Info objects)
  - `NCInfoBldInternalState.hh`
- **Structure checks (`NCInfoBldStructure.cc`):**
  - The setting comes from the declaration, or is deduced from the atom positions among all settings of the number. The cell is checked and completed with `SGCellConstraints`.
  - Tolerances: `tol_warn = 5e-4` (warning, symmetric structure assumed) and `tol_error = 1e-2` (BadInput).
  - Error messages list the settings of the number and mark those that fit the positions. There is a FIXME to add Hermann-Mauguin symbols.
  - Sites are expanded to positions, or positions grouped into sites.
  - Positions are mapped into [0,1) and sorted.
  - Duplicate positions are detected with `pos_tol = 1e-4`, periodic across the x=1 boundary. The comment says this must match `cifutils.py`.
  - Debye temperatures and msd must be given for all atoms or none.
  - The volume is computed from the cell.
- **HKL lists:** lazy generator (`HKLPlanes::source`, `dspacingRange`) and `validateAndCompleteHKLList`.
  - Info API: `hklListPartialCalc`, `getBraggThreshold` (partial loading), `hklDMinVal`/`hklDMaxVal`, `dspacingFromHKL`.
- **New public data (ABI change):**
  - `StructureInfo::spacegroup_setting`.
  - `AtomInfo::Site` with position, multiplicity, and `occupancy`/`uij` fields that are placeholders, "not yet supported" (`/workspace/ncrystal_core/include/NCrystal/interfaces/NCInfoTypes.hh`, line 161).
  - `AtomInfo::sites()`.
  - Python `AtomInfo.sites`; Info JSON queries "sites" and "structure".
- **Where Debye temperatures and MSDs come from in NCMAT loading** (`/workspace/ncrystal_core/src/ncmat/NCLoadNCMAT.cc`, lines ~560-800):
  - v8 `msd X at_temperature T` is converted to a Debye temperature and then back to the material temperature.
  - The per-element Debye temperature is converted to MSD (CPC paper section 2.5).
  - For crystalline v8 materials with scatknl dynamics, the msd is resolved at load time from the explicit value or concurrently by `SABAnalyser`, with a load error if no confident value results. FIXME: "Revisit for NCrystal 5.0.0 ... cleaner at time of usage".
  - Without dynamics but with Debye temperatures, `DI_VDOSDebye` entries are synthesised.
  - Config parameter `strain` (range [-0.5, 0.5]) scales a, b, c by (1 + strain).

## 10. .laz / .lau files (McStas PowderN / Single_crystal)
- **Files:** `/workspace/ncrystal_core/src/lazlau/NCLazy.hh`, `NCLazy.cc`, `NCLauLazyFact.cc` (the comment calls it a "small test factory ... probably not fully functional").
- **Physics:**
  - |F|^2 values are taken from the file. The hkl list is validated and normalised against the space group (`validateAndNormaliseHKLFsqList`), which checks that F^2 is consistent within each family, that listed d matches the calculated d, and that listed multiplicity is not above the expected value.
  - Dynamics come from `DI_VDOSDebye` with `debye_temperature` (default 300 K) and temperature (default 293.15 K).
  - Configuration: `dcutoff`, `dcutoffup`, `temp`.
- **On this branch:**
  - Optional `spacegroup_setting` header field (strict parse).
  - Files with ncmat2hklversion < 3 are assumed to be in the default setting, and files affected by issue #383 are rejected.
  - Cells are checked with `completeAndCheckStructureInfoCell`.
  - KSqCalc is used.
  - Rhombohedral axes are rejected (only hexagonal axes are allowed).

## 11. Filter tables (`Filter::createTable`) - entirely new; not a physics model
- **Files:**
  - `/workspace/ncrystal_core/include/NCrystal/internal/filter/NCFilterTable.hh`
  - `/workspace/ncrystal_core/src/filter/NCFilterTable.cc`
  - `/workspace/ncrystal_core/src/query/NCFilterQuery.{hh,cc}` (JSON query `["filtertable", CFGSTR, ...]`)
  - C API `ncrystal_filtertable` (`/workspace/ncrystal_core/include/NCrystal/cinterface/ncrystal.h`, line 1310)
  - Python `/workspace/ncrystal_python/src/NCrystal/filter.py` (`NCrystalFilter`)
  - `/workspace/examples/ncrystal_example_filter.c`
  - Design study: `/workspace/devel/filterstudy/` (README, `wiki_filter_tables.md`; GitHub issue #338)
- **What it is:** piecewise-linear tables of total cross section sigma(lambda) (scattering + absorption, barn/atom) for isotropic materials, for wavelengths 0 to `wlmax`. It is aimed at GPU-friendly, attenuation-only McStas filter components.
  - Bragg edges and energy-domain boundaries are represented as duplicate lambda points.
  - `evalTable` is the reference evaluator, with linear extrapolation clamped at 0.
- **Parameters (`TableParams`):** `wlmax` = 500 Aa, `tol` = 1e-3 (relative to max(sigma, 1e-12 b)), `ndense` (0 = automatic via `autoNDense`), `discontinuities`, `ReductionAlgo` (Greedy or DouglasPeucker).
- **Algorithm:**
  - Dense log-spaced sampling from 1e-9 Aa, plus points at +-1e-10 relative around each edge.
  - Reduction via `reducePoints`.
  - Midpoint refinement, at most 100 iterations, built with tolerance x 0.98.
  - Final check at random points.
  - Throws on an unknown discontinuity (relative spacing < 1e-9), on non-finite values, or if sigma does not converge as lambda -> 0.
- **Supported processes:** NullScatter, NullAbsorption, FreeGas, ElIncScatter, SABScatter, PowderBragg, AbsOOV. Any other process triggers a warning. Notably `SABScatterNG` is not in the list.
- **Open item:** FIXME "revisit the name 'filter' ... CrossSectionTable might be more appropriate" and check for duplicated utilities.

## 12. Extinction
- **There is no extinction code in the tracked repository.** A grep for extinction, Sabine, Becker, Coppens and Darwin finds only platform-name ("Darwin") matches.
- Single-crystal "domain" in the code means `EnergyDomain`. Mosaicity exists only as the Gaussian `mos` model (sections 3-4).
- **External only:** the untracked directory `/workspace/nocheck/pluginsurvey/ncplugin-CrysExtn/` holds a plugin (S. Xu et al., J. Appl. Cryst. 58 (2025), doi 10.1107/S1600576725007939).
  - Models: Sabine uncorrelated and correlated block models (rectangular or triangular tilt distribution), Becker & Coppens, Cooper & Rouse, and Kulda's RED model, for primary and secondary extinction in polycrystals.
  - It is configured through `@CUSTOM_CRYSEXTN`.
- `nocheck/pluginsurvey/ncplugin-CrysText` (texture) is also untracked and marked "not ready".

## 13. Configuration-string variables for elastic physics
Defined in `/workspace/ncrystal_core/include/NCrystal/internal/cfgutils/NCCfgVars.hh`; pseudo-variables in `src/cfgutils/NCCfgVars.cc` (docs) and `src/cfgutils/NCCfgManip.cc` (parsing).

| Variable | Default | Effect |
|---|---|---|
| `dcutoff` | 0 (auto; -1 is also treated as 0) | Planes below this d are ignored. Range [1e-3, 1e5] Aa. |
| `dcutoffup` | inf | Planes above this d are ignored. |
| `sccutoff` | 0.4 Aa | Single-crystal planes below this d are treated as powder. |
| `mos` | none | Gaussian mosaic FWHM, in (0, pi/2]. Requires `dir1` and `dir2`. |
| `dir1`, `dir2` | none | Crystal orientation (`@crys:` / `@crys_hkl:` ... `@lab:`). |
| `dirtol` | 1e-4 | Tolerance on the dir1/dir2 opening angle, in (0, pi]. |
| `mosprec` | 1e-3 | Precision of the mosaic model, in [1e-7, 0.1]. |
| `lcaxis` | none | Axis of a layered crystal. |
| `lcmode` | 0 | Choice of layered-crystal model (section 4). |
| `strain` | 0 | Isotropic lattice scaling by (1 + strain), range [-0.5, 0.5]. |
| `incoh_elas` | true | Incoherent elastic component. |
| `coh_elas` | true | Coherent elastic component (Bragg for crystals, incoherent approximation for amorphous solids). |
| `sans` | true | SANS components. |
| `temp` | - | Material temperature; drives the Debye-Waller factors. |
| `inelas=external` | - | Gives `BkgdExtCurve`. |

Pseudo-variables:
- `bragg` is an alias for `coh_elas`.
- `elas` sets `coh_elas`, `incoh_elas` and `sans` together.
- `comp=` takes a list from elas, incoh_elas, coh_elas, bragg, inelas, sans, and only disables what is not listed.
- `bkgd=0|none` is obsolete and equivalent to `inelas=0;incoh_elas=0;sans=0`.

On this branch, only `knllux` and the `vdoslux` extensions are new in cfgutils (both inelastic). No new elastic variables were added.

## 14. Tests
**Compiled test apps** (`/workspace/tests/src/`):
- `app_sgsym` (10 commits on this branch)
- `app_symchecks`: the positions-based `checkStructureSymmetry` versus the FillHKL F^2 self-check
- `app_cellconstraints`: legacy lattice checks versus `SGCellConstraints`
- `app_eqrefltables`: new `EqRefl` versus the legacy tables, for all settings
- `app_eqhkl`
- `app_fillhkl`: Bragg threshold and partial hkl lists
- `app_lc`: LCBragg
- `app_backscat`: SC backscattering
- `app_debyemsd`
- `app_filtertable`, `app_capifilter`, `app_filterexample`
- `app_cfg`, `app_cfgparse`, `app_extracfg`
- `app_sctext` (SCT, inelastic)

**Python test scripts** (`/workspace/tests/scripts/`):
- sgsym (needs spglib and gemmi)
- ncmatsettings, ncmatsites, atomsites (v8 settings and @ATOMSITES)
- ncmatv8msd, multimsd
- cellvalidation, latticecells and long_latticecells (spglib)
- fillhkl, hklobj, nosg (no space group)
- strain
- ncmat2hkl, ncmat2hkl2 (laz/lau round trip)
- cif2ncmat, cif2ncmatgen, long_cif2ncmatgen (random CIFs in all 230 space groups)
- composer, composer2, verifyatompos (ase)
- dkelas (scatknl materials, Info and process level, including elastic)
- pgbe (issue #171, LCBragg)
- mmc_scge, mmc_scpg (single-crystal Ge and PG via MiniMC)
- azimuth (single-crystal Bragg)
- filtertable, filtertable_api

Most of these have `.log` reference files.
