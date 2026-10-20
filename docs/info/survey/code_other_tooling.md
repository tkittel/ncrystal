<!-- Unedited AI-generated survey inventory (2026-10-20), kept as raw source notes for docs/info. Not reviewed. To be removed from history before merging to develop. -->

# NCrystal survey: physics and tooling outside Bragg/elastic and inelastic kernels (branch tk_volatile, towards 5.0.0)

There is no `v4.4.6` tag. I used commit `aee2e1552` ("Update version numbers to 4.4.6") as the baseline: there are 447 commits since then. "NEW" below means the file or directory did not exist at that commit. "Deltan" is the number of commits touching the path since then. The two draft changelogs `docs/changes_for_changelog.txt` and `docs/changes_for_changelog2.txt` are both new and are the authoritative list of what changed on this branch.

## 1. Absorption
- **Files:** `ncrystal_core/src/absfact/NCStdAbsFact.cc` (factory "stdabs", Priority 100, `MultiPhaseCapability::Both`) and `ncrystal_core/src/absoov/NCAbsOOV.cc` (class `AbsOOV`). Delta0 for both.
- **Model:** only the pure 1/v law, `sigma(E) = sigma_2200 * sqrt(E_2200/E)`, implemented as `m_c/sqrt(E)` (infinite at E=0).
- **Where sigma_2200 comes from:** `Info::getXSectAbsorption()` (`src/interfaces/NCInfo.cc:190`) is the composition-weighted sum of `AtomData::captureXS()`. So it comes from atomdb, `@ATOMDB` or the `atomdb=` cfg.
- **Multiphase:** handled automatically because the Info carries the combined composition. `createMerged` adds the scaled coefficients.
- **Limitations:**
  - No tabulated or energy-dependent absorption other than 1/v.
  - No resonances. The only mention is a comment in the filter table code about "a future process with absorption resonances".
  - Absorption factories can be swapped with the cfg var `absnfactory`, and plugins can register their own absorption factories.
- In MiniMC, absorption is applied by weight attenuation (engine option `absorption=0|1`).

## 2. SANS
- **Hard-sphere model:** `ncrystal_core/src/sanshardsphere/NCSANSFact.cc` (factory `hardspheresans`, registered as plugin "stdexpscat" via `register_experimentalscatfact`) and `NCSANSSphScat.cc` (`SANSSphereScatter`). The header calls it "a small test factory ... mainly intended to test and demonstrate the feasibility of SANS models".
  - **Physics:** dilute approximation, S(Q)=1, form factor `P(Q)=[3(sin x - x cos x)/x^3]^2`. The xs is in closed form, `C_sans*(3pi/32)R^3*(Rk)^-6[32(Rk)^4-8(Rk)^2+4Rk sin4Rk+cos4Rk-1]`, with a Taylor expansion at low Rk. Sampling is in x=QR.
  - **Activation:** an NCMAT `@CUSTOM_HARDSPHERESANS` section holding one number, the radius in Aa (0 < R < 1e6). It only works in a multiphase context (`MPOnly`, Priority 200) and requires `sans=true`.
- **Shared utilities:** `src/extd_utils/NCSANSUtils.cc` and `include/.../NCSANSUtils.hh` (Delta1).
  - Generic formula: `sigma_sans(k)=(2pi*phi*drho^2/nd)*(Vp/k^2)Integral_0^{2k} Q P(Q) S(Q) dQ`, with strong type `SANSScaleFactor`.
  - `extractCustomDataForSANSPlugin` / `hasCustomDataForSANSPlugin` scan the phase list. Phase "A" (with the custom section) pairs with the following phase "B". The SLD contrast is `Info::getSLD()`, and phi is the smaller volume fraction of the pair. The scan recurses into nested multiphases and respects `;sans=0` on sub-phases.
- **I(Q) helper:** `src/phys_utils/NCIofQHelper.cc` (`IofQHelper`) integrates and samples Q*I(Q) on [0,2k] from user (Q,I) points, flat below the first point (giving 1/v at low E). It is unused in core; the out-of-repo SANSND plugin uses it.
- **Cfg:**
  - `sans` (bool, default true, described as "primarily added to support future developments").
  - `elas=0` also turns `sans` off.
  - `comp=` accepts `sans`.
- **Python:** `NCMATComposer.add_hard_sphere_sans_model(sphere_radius)`; `nctool --sans`.
- **Not done:** a TODO in `src/stdscatfactory/NCStdMPScatFact.cc:95` says native (non-@CUSTOM) SANS in NCMAT is not implemented. `ncmat_doc` v6 also hints at future "phase-contrast physics (SANS)".

## 3. Atom DB
- **Files:** `src/atomdb/NCAtomDB.cc` (Delta0) and `NCAtomDBExtender.cc` (handles `@ATOMDB` and `atomdb=` lines; dedups identical AtomData instances via a global cache).
- **Content:** 348 entries (natural elements and isotopes), keyed by (Z,A) with A=0 meaning the natural element. Values are mass, b_coh, sigma_inc and sigma_abs. The auto-generated table cites:
  - Scattering lengths: the ATI Vienna table (Scattering_lengths_table_20010419, parsed from the HTML), following Landolt-Bornstein I/16A and ADNDT 49 (1991) 65 (i.e. Sears 1992).
  - Masses and natural abundances: AME2016 (Huang et al., CPC 41 030002). Natural-element masses are computed from isotopes, so they deviate slightly from CIAAW values.
- **Manual fixes:**
  1. Four legacy elements kept from pre-v2.1.
  2. Xe sigma_inc set to 1.16133 b (was 0).
  3. Ca b_coh 4.7->5.795 fm and Ca40 -> 5.9096 fm, per NIMA 2019 doi:10.1016/j.nima.2019.02.072.
- **Abundances:** natural abundances are present in the generator but commented out. Isotope mixtures are only via `@ATOMDB`.
- **Labels:** D and T are aliases for H2 and H3. Generic labels X, X1..X99 may denote mixtures only.
- **`@ATOMDB` syntax** (NCMAT v3+):
  - `Sym mass u b_coh fm sigma_inc b sigma_abs b`
  - `X is f1 A f2 B ...`
  - `nodefaults`
  - Lines are evaluated in order.
- **`atomdb=` cfg:** same syntax, with `:` acting as a space and `@` as a line separator. It appends to the file's `@ATOMDB`.
- **Quick factories** (`solid::`, `gasmix::`) accept atomdb lines with `_` separators, e.g. `B_is_0.95_B10_0.05_B11`.
- **Browsing:**
  - NEW JSON query `["util","atomdb"]`.
  - NEW `ncrystal browse --atomdb` (e.g. `-w "absxs > 1000 and natural"`) and `NCrystal.browse.AtomDBBrowser` / `query_atomdb()`.
  - Python `NCrystal.atomdata`: `atomDB(Z,A)`, `iterateAtomDB`, element name/Z helpers.

## 4. Multiphase, density, temperature, state of matter, composition
- **Multiphase:**
  - Cfg syntax `phases<FRAC1*CFG1&..&FRACN*CFGN>[;COMMONCFG]` (fractions by volume; parsed in `src/factories/NCMatCfg.cc` ~L205). Examples in `NCMatCfg.hh`: Al+voids (`void.ncmat`), Al2O3 powder in water.
  - `phasechoice` is a top-level pseudo-param that selects a child phase at the Info level.
  - `@OTHERPHASES` (NCMAT v6) gives `frac cfgstr` lines; the primary phase becomes phase 0.
  - The std multiphase scatter factory is `stdmpscat`.
- **Density:**
  - Top-level `density=` (`src/cfgutils/NCCfgManip.cc` ~L193) takes units `gcm3`, `kgm3`, `perAa3`, or `x` (scale factor; repeated uses multiply).
  - `packfact` is obsolete (throws; `MatCfg::get_packfact()` returns 1.0).
  - Overrides are implemented in `src/infobld/NCInfoBldOverrides.cc`, which shares data and caches override variants. NEW split of InfoBuilder (Delta15).
  - Density and number density are completed from the unit cell or `@DENSITY` in `NCInfoBldComposition.cc`.
- **Temperature:**
  - `temp=` (default -1 means 293.15 K unless the data fixes it). Units K, C or F (`src/cfgutils/NCCfgTypes.cc:70`).
  - Allowed range is 0.001 K to 1e6 K.
  - NCMAT `@TEMPERATURE` (v7): `default T` sets a default, a bare `T` locks it.
  - `scatknl` temperature fields must agree.
  - `gasmix::` takes T in K, C or F.
- **Other Info cfg vars:**
  - `strain` (isotropic lattice scaling, range [-0.5,0.5], added in 4.4.0).
  - `atomdb`.
  - `infofactory`, `scatfactory`, `absnfactory`.
  - `ucnmode=refine|only|remove[:thr]` (UCN split, default 300 neV).
  - Component switches: `incoh_elas`, `coh_elas`, `sans`, `inelas`; aliases `bragg`, `elas`, `comp=` (list from elas, bragg, coh_elas, incoh_elas, inelas, sans).
  - `bkgd` is obsolete.
  - Embedded cfg via `NCRYSTALMATCFG[...]`.
- **State of matter:**
  - `Info::StateOfMatter` is one of Unknown, Solid, Gas or Liquid.
  - NCMAT v5 introduced `@STATEOFMATTER`.
  - v8 makes it mandatory for non-crystals.
  - Non-crystalline solids get incoherent-approximation elastic scattering from MSD. NEW: this also applies to `scatknl` when solid (v8 or explicit msd).
- **Composition:**
  - `NCInfoBldComposition.hh/.cc` (NEW file split) runs these validation steps in order: validateAtomIndexes, validateAndCompleteComposition, validateTemperatures, calculateAverageAtomMass, validateAndCompleteDensities, validateAndCompleteStateOfMatter.
  - Fractions are snapped to sum to 1, with tolerance 1e-9.
  - `src/misc/NCCompositionUtils.cc` handles flattening and breakdowns.

## 5. NCMAT format (`docs/ncmat_doc.md`, 1113 lines) by version
All files start with `NCMAT vN`, are ASCII or UTF-8 (non-ASCII only in comments), and use LF or CRLF line endings.

- **v1** (all releases): `@CELL` (lengths, angles), optional `@SPACEGROUP` (1-230), `@ATOMPOSITIONS` (fractional), `@DEBYETEMPERATURE` (global or per element).
  - Physics is the Debye model only (CPC 246 (2020) 106851).
  - Comments are allowed only before the first section.
  - Rhombohedral groups: from 5.0.0 rhombohedral axes are accepted, and the setting is deduced from the positions.
- **v2** (NC 2.0.0):
  - `@DYNINFO` per element with `element`, `fraction` (rationals like 2/3 allowed) and `type`:
    - `sterile`.
    - `freegas`.
    - `scatknl`: `alphagrid` / `betagrid` / `sab` or `sab_scaled` (S'=S*e^{-beta/2}, optionally beta>=0 only), mandatory `temperature`, grid sizes 5..65534, `egrid` (1, 3 or >=10 values), and the `<v>r<n>` repeat syntax.
    - `vdos`: `vdos_egrid` (2 points means endpoints), `vdos_density` (>=5 values), E^2 extrapolation below the grid, first point >=1e-5 eV, internally re-gridded.
    - `vdosdebye`.
  - `@DENSITY` in `atoms_per_aa3`, `kg_per_m3` or `g_per_cm3`.
  - The cell sections become optional, which enables non-crystals.
  - "D" element, fractions in positions, comments anywhere.
  - `inelas` overrides.
  - Planned S(Q,omega) (`qgrid`/`omegagrid`/`sqw`) is not implemented.
- **v3** (2.1.0): isotope labels (H2, Li6, ...), D/T aliases, `@ATOMDB` (see Sec. 3), and `@CUSTOM_<A-Z+>` sections (raw lists of strings, a warning is emitted, used by plugins).
- **v4** (2.6.0): `@DEBYETEMPERATURE` not needed for elements with `vdos`; global Debye temperature disallowed; `@CELL` gains `cubic a` and `!!` (repeat previous length).
- **v5** (2.7.0): `debye_temp` inside `vdosdebye` (allows amorphous solids); `@STATEOFMATTER` (solid, liquid or gas). Crystals without `@DEBYETEMPERATURE` need all-vdos/vdosdebye.
- **v6** (3.0.0): `@OTHERPHASES` (volume fraction plus cfg-string per line, whitespace normalised).
- **v7** (3.2.0): `@TEMPERATURE` (`default T` or locked T; range 0 < T <= 1e6 K).
- **v8** (5.0.0, NEW on this branch):
  - `msd <Aa^2> [at_temperature T]` as an alternative to `debye_temp` (converted via the isotropic Debye model). Both are allowed in `scatknl` only for solids.
  - `effective_temperature` in `scatknl` (ENDF-style Teff, used for short-collision-time extension). It is auto-estimated when absent, for any NCMAT version.
  - Crystals may use `scatknl`; the msd comes from the keywords or confident automatic kernel analysis (SABAnalyser), and the automatic extraction is gated to v8+.
  - `@STATEOFMATTER` mandatory for non-crystals.
  - `@DEBYETEMPERATURE` removed; freegas/sterile banned in crystals.
  - `@SPACEGROUP` setting codes from ITVB Table A1.4.2.7 (530 settings, e.g. `227:2`, `14:c1`, `62:cab`, `62:abc`, `166:r`, lower-case only). A bare number is allowed only for the 140 single-setting groups; positions are verified against the symmetry.
  - New `@ATOMSITES` section (one position per orbit, expanded by sgsym; mutually exclusive with `@ATOMPOSITIONS`).
- **Implementation:** `src/ncmat/` (`NCParseNCMAT.cc`, `NCLoadNCMAT.cc`, `NCNCMATData.cc`, `NCFactory_NCMAT.cc`; Delta11). There is a FIXME at `NCLoadNCMAT.cc:616` ("Revisit for NCrystal 5.0.0").

## 6. Other formats, quick factories, virtual API, plugins
- **`.laz`/`.lau`:** `src/lazlau/NCLauLazyFact.cc` and `NCLazy.cc` (Delta7). Info factory "stdlaz", described as a "small test factory... probably not fully functional".
  - Reads the SPCGRP/TITLE headers and the NEW optional `# spacegroup_setting` header (written by ncmat2hkl, `ncmat2hklversion` 3).
  - Old ncmat2hkl (<=4.4) files with cos(alpha) != cos(beta)cos(gamma) are flagged (issue #383).
  - Supports `atomdb`.
  - Writing is done by `ncrystal_ncmat2hkl` / `mcstasutils.cfgstr_2_hkl`.
- **Quick factories** (TextData, `Priority::OnlyOnExplicitRequest`; `src/quickfact/`). They emit NCMAT v7 (deliberately not v8).
  - `freegas::FORMULA/DENSITY[/atomdb]`: free gas, `@STATEOFMATTER gas`.
  - `solid::FORMULA/DENSITY[/TDebye<T>K[_El]][/atomdb]`: amorphous solid with `vdosdebye`, default Debye temperature 300 K.
  - Density units: `perAa3`, `kgm3`, `gcm3`.
  - `gasmix::` (`NCGasMixFact.cc`, `src/gasmix/NCGasMixUtils.cc`):
    - Molar (default) or `massfractions` mixtures, e.g. `0.7xCO2+0.3xAr`; `air` alias.
    - Pressure units Pa, mbar, hPa, kPa, bar, atm (default 1 atm). Temperature in K, C or F. Or a density.
    - `relhumidity` uses the saturated vapour pressure formula from doi:10.1175/JAMC-D-17-0334.1.
    - Ideal gas law; emits `@TEMPERATURE` (locked), freegas `@DYNINFO`.
    - Examples: `He/10bar`, `BF3/2atm/25C/B_is_0.95_B10_0.05_B11`.
- **Virtual API:** `src/virtualapi/` and `include/NCrystal/virtualapi/NCVirtAPI_Type1_v1.hh` (added 4.1.0; used by OpenMC).
  - ABI-stable abstract interface: createScatter(cfgstr), clone, crossSectionUncached, sampling.
  - Example in `examples/virtualapi_project`; test in `tests/src/app_vapit1v1`.
- **Plugin system:** `src/plugins/NCPluginMgmt.cc` (Delta0).
  - Standard built-in plugins: stddatasrc, stdscat, stdmpscat, stdexpscat (hard-sphere SANS), stdlaz, stdabs, stdncmat, stdquick.
  - Dynamic plugins come via `ncrystal-pluginmanager` (pip-installed, `ncrystal_pypluginmgr/`), `NCRYSTAL_PLUGIN_LIST` and `NCRYSTAL_PLUGIN_DATADIRS`, and can be disabled with `NCRYSTAL_DISABLE_DYNLOAD`.
  - Plugins register Info, Scatter, Absorption or TextData factories (`FactImpl`). Typically they wrap `globalCreateScatter` and add a process when a `@CUSTOM_<NAME>` section is present.
  - Data files are named `plugins::<Name>/file.ncmat`.
  - Templates: `examples/plugin` (DummyPlugin: step-function incoherent elastic via `@CUSTOM_DUMMYPLUGIN`) and `examples/plugin_dataonly`.
  - The CHANGELOG names `ncrystal-plugin-WaterData`.
  - A local untracked survey in `/workspace/nocheck/pluginsurvey/` lists:
    - ncrystal-plugin-ClathrateHydrateData and GraphiticData (HighNESS data).
    - CrysExtn (extinction: Sabine, Becker-Coppens, Cooper-Rouse, Kulda RED).
    - CrysText.
    - MagScat (paramagnetic O2-hydrate/O2-C60, Zimmer). The changelog notes its files now fail the new symmetry check.
    - SANSND (nanodiamond SANS using IofQHelper, Teshigawara fit).
    - nxslib (legacy .nxs).
    - BzScope (coherent phonon S(Q,omega), CINEMA/IHEP).
    - LiquidHeliumData.
    - Liquids (`@CUSTOM_LIQUIDS`).
    - ncrystal-extra.
- **Data source namespaces** (`src/factories/NCDataSources.cc`): `stdlib`, `stdpath` (NCRYSTAL_DATA_PATH), `plugins`, `relpath`, `abspath`, `virtual` (in-memory, `registerInMemoryFileData`), `freegas`, `solid`, `gasmix`. Usage is `name::file`.

## 7. Python tooling (`ncrystal_python/src/NCrystal/`)
- **`ncmat.py` NCMATComposer** (impl `_ncmatimpl.py` Delta15, `ncmat.py` Delta3).
  - Sources: `from_cif` (incl. `codid::`/`mpid::` online, `mp_apikey`), `from_ase`, `from_info`, `from_cfgstr`, `from_ncmat`, `from_hfg`.
  - Structure: `set_cellsg[_cubic]`, `set_atompos`, `refine_crystal_structure` / `verify_crystal_structure` (spglib, symprec), `as_spglib_cell`.
  - Composition: `set_composition`, `remap_atom`, `update_atomdb`, `set_fraction`, `set_density`.
  - Dynamics:
    - `set_dyninfo_vdos`, `_vdosdebye`, `_debyetemp`.
    - `_msd(temperature=)` and `_uiso`: NEW, they write the v8 `msd` keyword.
    - `_scatknl`, `_freegas`, `_sterile`, `_from_object`, `transfer_dyninfo_objects`, `allow_fallback_dyninfo`.
  - Other: `set_state_of_matter`, `add_secondary_phase`, `add_hard_sphere_sans_model`, custom sections, `lock_temperature` / `set_default_temperature`.
  - Output and loading: `create_ncmat`, `write`, `register_as`, `load*`, `plot_xsect`, `inspect`.
  - NEW: writes v8 with `@ATOMSITES` and the setting for crystals; full-precision cells.
- **`cifutils.py`:** `CIFSource`, `CIFLoader` (gemmi plus spglib, merging of equivalent atoms, Uiso handling), `produce_validation_plot(s)`.
  - Doc bug: `_cli_cif2ncmat.py:59` refers to a non-existent "CIFAnalyser" class.
- **`vdos.py`:** `createVDOSDebye`, `debyeIsotropicMSD`, `debyeTempFromIsotropicMSD`, `analyseVDOS` (msd, Debye temperature, Teff, gamma0, integral), `extractGn`, `extractKnl`.
  - `PhononDOSAnalyser`: raw, AnyVDOS or QuantumEspresso matdyn.dos input; merge, drop, cutoff, regularise, plot, `apply_to(composer)`.
  - `misc.AnyVDOS`.
- **`mcstasutils.py`:** `cfgstr_2_unioncfg`, `cfgstr_2_union_instrument_code`, `cfgstr_2_hkl` (laz/lau), `cfgstr_detect_components`.
- **`filter.py`** (NEW): `NCrystalFilter(cfgstr)`, the macroscopic total xs (scattering plus absorption) vs wavelength as a piecewise-linear table to 0.1%, for transmission exp(-xs*L).
  - C++ side: `src/filter/NCFilterTable.cc` (NEW, FIXMEs: "revisit the name filter... CrossSectionTable", reproducibility).
  - C API: `ncrystal_filtertable`; example `examples/ncrystal_example_filter.c`; JSON query `["filtertable",...]`.
  - Intended for GPU McStas attenuation-only filter components.
  - Not mentioned in either changelog draft.
- **`browse.py`** (NEW, 1293 lines): `DataBrowser`, `DataEntry`, `PhysicsProps`, `query_data`, `AtomDBBrowser`.
- **Other modules:** `plot.py` (plot_xsect(s), plot_vdos, plot_knl, plot_vdos_Gn), `hist.py` (Hist1D), `datasrc.py`, `cfgstr.py` (generateCfgStrDoc, decodeCfg), `misc.py` (evaluate_query, MaterialSource, AnyTextData), `hfg2ncmat.py`, `ncmat2endf.py`, `ncmat2cpp.py`, `atomdata.py`, `constants.py`, `plugins.py`.
- **CLI tools** (`ncrystal <cmd>` or `ncrystal_<cmd>`; scripts in `pyproject.toml`):
  - `nctool`: main investigation tool (dump, plots, `--browse`/`--extract`/`--plugins`, `--doc`, `--bench`, `--mc`, `--phases`, component flags).
  - `browse` (NEW): search and filter data by name, glob or regex, header text (`-s`, `-E`) and physics `-w` expressions; outputs `--columns`, `--sort`, `--csv`/`--json`/`--html`, `--info`, `--atomdb`, `--plugins`. Documented in NEW `docs/data_browsing.md`; FIXME at `_cli_browse.py:565`.
  - `cif2ncmat`: CIF to NCMAT (gemmi plus spglib); `--uisotemp`, `--debyetemp` (default 300), `--dynamics`, `--via-ase`.
  - `vdos2ncmat`: described as "somewhat outdated" VDOS preparation.
  - `hfg2ncmat`: hydrogen-rich amorphous materials from hydrogen functional groups (CHali, CHaro, CH2, CH3, NH, NH2, NH3, OH, SH; Romanelli et al. JPCM 2021). Default Debye temperature 400 K for non-H atoms.
  - `endf2ncmat`: ENDF TSL to `@DYNINFO` via PyNE, expert only.
  - `ncmat2endf`: ENDF-6 TSL export via IAEA endf-parserpy; `elastic_mode` greater/mixed/scaled, othertemps, isotopic expansion.
  - `ncmat2hkl`: laz/lau for McStas PowderN or Single_crystal.
  - `ncmat2cpp`: C++ code with inline NCMAT data.
  - `mcstasunion`: McStas-Union instrument code.
  - `minimc`: see Sec. 8.
  - `query`: low-level JSON queries.
  - `verifyatompos`: still uses ASE spacegroup, now partly redundant with the sgsym load-time check.
  - `config`.
- **JSON query machinery** (`src/query/`, Delta38, many NEW files). Top-level keys: `version`, `list`, `mmc`, `util`, `sab`, `vdos`, `filtertable`.
  - `util` subkeys: wl2ekin, ekin2wl, mathval, fmadiagnose, kinutils, browsedb, browsefactories, factorythreads, atomdb, eqrefl, atomsites, atompositions, sgsym (most NEW).
  - `sab` subkeys: refeval, samplepb, surveyor, sglcell, integschemes, proc, refsample, analyse.
  - `vdos`: expand.
  - NEW per-object queries via `ncrystal_info_jsonquery`: keys list, structure, sites.

## 8. MiniMC
- **Files:** `ncrystal_core/src/minimc/NCMMC_*.cc` (Delta5) and Python `minimc.py`, `minimc_objects.py`, `_mmc*.py`. Public since 4.3.0.
- **Physics:** a single convex volume. The neutron is split at each step into a transmitted (tallied) part and a forced collision; weight-based absorption; roulette `P,W,N`; `nscatlimit`.
- **Geometry:** `sphere;r=`, `box;dx;dy;dz`, `slab;dz`, `cyl;r;dy` (about the y axis, dy=0 means infinite). Lengths in m.
- **Sources:** `constant` (pencil), `circular` (disk), `isotropic` (point or sphere). Common params `n` (default 1e6), `ekin`/`wl`.
- **Engine options:** ignoremiss, nthreads, seed (64-bit), tallyref (truth or src), absorption, nscatlimit, roulette, tally, tallybins, tallybreakdown.
- **Tallies:** de, e, l, mu, nscat, nscat_uw, q, theta, w (plus highres, lowres, nobreakdown flags), scored as neutrons exit the volume.
- **Scenario strings:**
  - Form: `ENERGY [pencil] on THICKNESS [sphere|slab] [N times]`.
  - Energy units: Aa, eV, meV, keV, MeV, GeV, kT, BT (Bragg threshold).
  - Thickness: a length or `mfp`.
  - Defaults: 0.8BT, 1kT or 1.8Aa; 1e6 neutrons (1e5 if oriented).
- **Entry points:** `minimc.run(cfgstr, geomcfg, srccfg | scenario, enginecfg, callback)`; `nctool --mc`.

## 9. Standard data library (`data/`, Delta1)
- **Size:** 134 `.ncmat` files (133 materials plus `void.ncmat`).
- **Versions:** 38 v4, 94 v5, 2 v7 (B4C, BoricAcid), 0 v8.
- **Crystals (110):** elements and compounds, mostly with `vdos` @DYNINFO (123 files have at least one vdos).
  - Exceptions using `vdosdebye` only: Ca_sg229, CaSiO3, B4C.
  - Issue #383 fixes changed the densities of BoricAcid (1.562 g/cm^3) and Wollastonite.
- **Gases (5, freegas, STP):** Ar, He, Kr, Ne, Xe.
- **Liquids (2, the only `scatknl` files):** LiquidWaterH2O_T293.6K (O freegas) and LiquidHeavyWaterD2O_T293.6K. Now SCT-extended at Teff~1100 K.
- **Amorphous solids (16, `@STATEOFMATTER solid`, `@DENSITY`):** AcrylicGlass, Epoxy_Araldite506, Kapton, Nylon11/12/610/66or6, PEEK, PVC, Polycarbonate, Polyester, Polyethylene, Polylactide, Polypropylene, Polystyrene, Rubber.
- **`void.ncmat`:** the only file with `@ATOMDB` (zeroed H, sterile, density 1e-200).
- **Not used in stdlib:** `@OTHERPHASES`, `@TEMPERATURE`, `@CUSTOM_`.

## 10. Integrations (from README and CHANGELOG)
- **OpenMC:** `openmc.Material.from_ncrystal(cfgstr)`, OpenMC >= 13.3 built with `-DOPENMC_USE_NCRYSTAL`, via the virtual API. CI in `.github/workflows/openmc.yml`; `devel/openmc/basicsim.py`.
- **Geant4:** separate `ncrystal-geant4` package; `.github/workflows/geant4.yml`.
- **McStas 3:** `NCrystal_sample.comp` and Union `NCrystal_process`; McStasScript `add_ncrystal_union_material`; `mcstasutils --union`; laz/lau export; the new filter tables target GPU filter components.
- **ENDF ecosystem:** endf2ncmat (PyNE) and ncmat2endf (endf-parserpy) cover NJOY/MCNP-style users. MCNP and PHITS appear only in the SANSND plugin README.
- **External libraries:** ASE (CIF and structures), spglib, gemmi, COD and Materials Project online CIF.
- **Not mentioned:** McXtrace, Vitess, Mantid.

## Other new or devel material relevant to these areas
- `docs/devel_site_occupancy_roadmap.md`: `AtomInfo::sites()` has placeholders for occupancy and anisotropic displacement, not yet supported.
- `docs/devel_sgsym_plans.md`; `src/sgsym/` (NEW).
- NEW `NCSCTUtils` (short-collision-time) and `NCSABAnalyser` (Teff/msd from kernels): these are inelastic-area items.
