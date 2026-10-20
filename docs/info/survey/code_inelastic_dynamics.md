<!-- Unedited AI-generated survey inventory (2026-10-20), kept as raw source notes for docs/info. Not reviewed. To be removed from history before merging to develop. -->

# NCrystal inelastic / atomic-dynamics inventory (branch tk_volatile, towards 5.0.0)

Paths: `src` = /workspace/ncrystal_core/src, `inc` = /workspace/ncrystal_core/include/NCrystal/internal.

**Change baseline.** There is no `v4.4.6` tag. I used commit `aee2e1552` ("Update version numbers to 4.4.6") as the baseline. Commits since then, per directory: vdos 50, sab 59, sabscatter 2, dyninfoutils 5, phys_utils 10, freegas 0, gasmix 2, stdscatfactory 2, cfgutils 3. The diff is about +11.7k / -1.3k lines.

Draft changelog text exists in /workspace/docs/changes_for_changelog.txt (lines 1-130 cover the inelastic work) and /workspace/docs/changes_for_changelog2.txt.

**Defaults are still legacy.** Both next-gen paths are opt-in for now:
- `vdoslux` defaults to 3, which means legacy VDOS expansion. Next-gen is selected with 2000..2006, a hidden developer-only range.
- `knllux` defaults to -1, which means the legacy SAB integrator/sampler. Next-gen is 0..6, or 100..106 / 200..206.

---

## 1. vdos (VDOS -> S(alpha,beta))

### 1.1 VDOSEval: density evaluation and integrals
- **Files:** inc/vdos/NCVDOSEval.hh, src/vdos/NCVDOSEval.cc. Changed on branch (+137).
- **Physics:**
  - f(E) is linearly interpolated on a regularised grid, assumed to scale as E^2 below the first point (solid-state assumption), and 0 above emax. It is auto-normalised over [0, emax].
  - `calcGamma0()` integrates coth(E/2kT)/E (Sjolander).
  - `getMSD(gamma0)` gives msd = gamma0*hbar^2/(2M) (Sjolander p326, Q*gamma0 = 2W).
  - `calcEffectiveTemperature()` gives kTeff = Integral rho(E)(E/2)coth(E/2kT) (Sjolander I.6). The comment openly worries about a factor 2 or 3 but notes Teff -> T at high T.
  - `evalG1Symmetric` / `evalG1Asymmetric` give the Sjolander G1 (eq. II.29). Above 200kT the asymmetric form must be evaluated directly rather than through a detailed-balance factor (`numericallySafeG1SymmetricELimitInUnitsOfKT`). `evalG1AsymmetricAtEPair` returns both +-E.
- **Free functions:** `regulariseVDOSGrid`, `checkIsRegularVDOSGrid`. The regular grid is equidistant and, if extended downward, hits 0 exactly.
- **Numerics:**
  - `safe_xcothx`, `safe_x3cothx`, `safe_xdivsinhx` use Taylor expansions near 0.
  - The parabolic [0, emin] part uses `integrateRomberg33`.
  - Bin integrals are accumulated with `StableSum`.
  - Teff is range-checked.
- **NCMAT input:** `vdos_egrid`, `vdos_density` (first egrid point must be >= 1e-5 eV).

### 1.2 VDOSGn: Sjolander multi-phonon Gn functions
- **Files:** inc/vdos/NCVDOSGn.hh; src/vdos/NCVDOSGn.cc (+734), NCVDOSGn_FMA.hh, NCVDOSGn_WINFMA.cc (FMA dispatch, new).
- **Physics:** asymmetric Gn functions (Sjolander 1958, eq. II.27). Gn is obtained by self-convolving G1.
  - G{2n} = Gn * Gn and G{2n+1} = Gn * G{n+1}, so only about log2(n) convolutions are needed.
  - Convolutions use FFT on a fixed binning. Orders up to 10000 or more have been tested.
- **Interface:** `enum class Cfg { Default, Legacy, MaxLux }`, `growMaxOrder(Order)` (not MT-safe), `eval`, `evalMany`, `eRange(n, relthr)`, `binWidth`, `getRawSpectrum`. `Order` is a 1-based checked integer.
- **Cfg knobs (`CfgDecoded`):** thinning (minThinOrder=4, thinNBins=1000, aggressive thinning from order 20 at 300 bins), truncation (threshold 1e-13), `g1MaxBinWidthKT=0.25`, `g1InterpTol=3e-3` (`g1InterpThickenFactor`), `g1InterpMaxNBins=1e5`, `minTemperature=0.1K`, `thinSmoothTol=0.1`, `thinSmoothFloor=1e-6`.
  - `Legacy` reproduces NCrystal 2.x-4.x: thin/trunc from order 5, threshold 1e-14, `convolveLegacy`, no G1 thickening.
  - `MaxLux` uses `FastConvolve::convolveDirect` (StableSum, no FFT floor). It is not wired to any vdoslux yet ("fixme: wire MaxLux to vdoslux=2006").
- **New on branch:**
  - Grids are anchored at E=0 (`startIdx`), so thinning and truncation are independent of rounding.
  - FFT noise is snapped to 0 using `estimateFFTConvolutionNoiseFloor`.
  - Truncation edge is found with `estimateSpectrumCrossing` (windowed quadratic fit of ln spec) plus `applyCrossingTaper`.
  - Smoothness-gated thinning.
  - G1 bin width capped at 0.25kT (fixes detailed-balance errors at low T; T < 0.1K is rejected).
  - "FIXME: Revisit all of the above for a dependency on vdoslux".

### 1.3 VDOSExpand: `expandVDOSToGnFcts` (new file)
- **Files:** inc/vdos/NCVDOSExpand.hh, src/vdos/NCVDOSExpand.cc.
- **Returns `GnExpansion`:** Gn, `alpha2x` (= 2W/alpha = 2*m_n/hbar^2 * kT * msd), `suggestedEmax`, `sabRange` (Rectangle), and per-order `abRanges` (the comment says these may be diagnostics only).
- **Algorithm:** add orders until the (alpha, beta) rectangle of order n (from `Gn.eRange` and `rangeXNexpMX`) no longer intersects the kinematic region for targetEmax (`findABExtentWithinKB`).
  - **Target Emax:**
    - legacy: lux-dependent {0.5, 1, 3, 5, 8, 12} eV;
    - next-gen: always 5 eV;
    - an explicit `egrid` Emax overrides both.
  - **relcontriblvl:**
    - legacy: 10^-(3+2*lux);
    - next-gen: 1e-6 .. 1e-13 for lux 0..6.
  - **Resource limits (next-gen):** total Gn points 1e6/2e6/3e6/5e6/1e7/3e7/1e8 for lux 0..6, single convolution 5e6 points, at most 60000 orders. Legacy uses order_limit 1000 (x10 for lux 5 or explicit Emax, /10 for lux 0).
  - If a limit is hit, Emax is reduced in steps of x0.99. Next-gen raises a CalcError if this would go below 0.1 eV ("FIXME: 0.1eV just a first guess"). Next-gen warns on any reduction.
  - FFT caches are trimmed to 10 MB if they exceed 100 MB (`FFTCacheGuard`).

### 1.4 VDOSToScatKnl: `createScatteringKernel`
- **Files:** inc/vdos/NCVDOSToScatKnl.hh, src/vdos/NCVDOSToScatKnl.cc (+FMA/WINFMA files, new). Heavily reworked.
- **Physics:**
  - Incoherent approximation plus Gaussian approximation (harmonic): S(alpha,beta) = Sum_n exp(-x)*x^n/n! * Gn(beta), with x = alpha2x*alpha (Sjolander II.28, without the sigma*kT/4E factor).
  - Above n = 16 the factorial uses Stirling's series (10 terms, `stirling_threshold=16`); below, it is built recursively. Half of exp(-x) is split out to extend the usable range.
  - Next-gen `alpha_factor_floor` = 1e-100.
- **`ScaleGnContributionFct`:** scales the Gn contribution per order. Example use: remove sigma_coh from G1 when coherent one-phonon scattering is modelled elsewhere.
- **Grids:**
  - Legacy: `setupLegacyAlphaGrid` / `setupLegacyBetaGrid` (region schemes, `reducePtsInDistribution` with gap filling, extra points near beta=0).
  - Next-gen: `determineAlphaBetaGridFromGn`.
- **Output:** `ScatKnlData` (KnlType SAB, `betaGridOptimised=true`, suggestedEmax).

### 1.5 VDOSKnlGrid: next-gen (alpha, beta) grid (new)
- **Files:** inc/vdos/NCVDOSKnlGrid.hh, src/vdos/NCVDOSKnlGrid.cc.
- **`gridDimFromLux`:**
  - legacy: 50*2^lux x 100*2^lux;
  - next-gen (nalpha x nbeta): 60x130, 110x240, 170x400, 280x650, 470x1080, 770x1770, 2000x4000 for lux 0..6.
  - Env var `HACK_NBINS` overrides.
- **`getCombinedGnFct`:** weighted Gn sum with weights ~1/n^2, normalised, used to place beta points.
- **`setupE0ABGrid`:** points along alpha=beta for the E->0 limit, weighted by sqrt(beta)*f(x,n).
- **Point placement:**
  - Beta points via `reducePtsByEquidistributionRobust` (noise tolerant, wider curvature stencil; in utils/NCMath.hh alongside `reducePtsByEquidistribution` and `reducePtsInDistribution`, `PtReduceCfg`).
  - Alpha points via powspace p=2 (linear in sqrt(x)), starting from xlow=0.005.
  - Both grids are guaranteed to contain 0 and have the exact requested size (`spaceOutGrid`, `topOffGrid`, `mergeGridsWithTol`). The beta=0 node is snapped exactly to 0.
- **Comments:** "TODO: reduce nc_assert_always ...", "fixme: parallelise", "Tune powspace powers".

### 1.6 VDOSUtils (new)
- **Files:** inc/vdos/NCVDOSUtils.hh, src/vdos/NCVDOSUtils.cc (+FMA).
- **Functions:** `rangeXNexpMX`, `estimateSpectrumCrossing`, `applyCrossingTaper`, `estimateGnErange`, `estimateFFTConvolutionNoiseFloor` (about 8*peak*eps*sqrt(n), calibrated in app_gnconvnoisefloor), `findABExtentWithinKB` (intersection with (alpha-beta)^2 < 4(E/kT)alpha), `spaceOutGrid`, `topOffGrid`, `equidistantGridPoint` (single-rounding std::fma), `makeCommonGrid`, `mergeGridsWithTol`, `coverEquidistantGrids`, `PWLFct`, `pwlNarrowToPos`, `evalPWLSum`, `trimTailByIntegral`, `trimEquidistantGridUpperEdge`.

### 1.7 VDOSLux (new)
- **Files:** inc/vdos/NCVDOSLux.hh, src/vdos/NCVDOSLux.cc.
- **Encoding:** raw values 0..5 and 1000..1005 mean legacy; 2000..2006 mean next-gen. `lvl()`, `isLegacy()` (bit 128), `reduceLvl`. The default is legacy 3.
- **Debye-model reduction:** for vdosdebye, legacy reduces the level by 3; next-gen applies no reduction (`vdoslux_reduced_for_debye` in dyninfoutils: "too confusing and erratic").

### 1.8 VDOSCache
- **Files:** inc/vdos/NCVDOSCache.hh, src/vdos/NCVDOSCache.cc. Not changed on branch.
- **Contents:** `VDOSDataHashPtr` (hash plus content compare), `EnergyGridPtr` (shared egrid plus UID), `getCachedVDOSDataHashPtr` (deduplicates identical VDOSData so caches work), and `DI_VDOSShPtr`.

### 1.9 Debye-model VDOS
- **Where:** `createVDOSDebye` in src/dyninfoutils/NCDynInfoUtils.cc.
- **Model:** rho ~ E^2 on 20 points over [0.5*E_D, E_D]; below that, VDOSEval's E^2 extrapolation is used. Must stay in sync with the Python version.
- **Caching:** results are cached by rounded (vdoslux, M, T, T_D, sigma_b) (`VDOSDebye2SABFactory`, 1e-7 rounding).

---

## 2. sab / sabscatter (kernels, integration, sampling)

### 2.1 Data, conversion and helpers
- **Files:** inc/sab/NCScatKnlData.hh, NCSABUtils.hh, NCSABEval.hh, NCSABIdx.hh (new); src/sab/NCScatKnlData.cc, NCSABUtils.cc.
- **`transformKernelToStdFormat`:** converts SAB, SCALED_SAB, SCALED_SYM_SAB (betagrid starting at 0, expanded with S'(alpha,-beta)=S'(alpha,beta)) and SQW to the standard asymmetric, unscaled S (`validateScatKnlData`). It calls `trimZeroEdgesFromKernel` (with a WARNING).
- **Interpolation:** log-linear in alpha (linear fallback where S=0), linear in beta. `integrateAlphaInterval(_fast)` uses the analytic (a2-a1)(s2-s1)/log(s2/s1).
- **NCSABIdx:** packed uint16 (ialpha, ibeta) indices in a uint32, so grids are limited to 65534 points.

### 2.2 Legacy path (knllux -1..-6; also all vdosdebye and ucnmode)
- **Files:** NCSABIntegrator, NCSABSampler, NCSABSamplerModels, NCSABXSProvider, NCSABScatterHelper, NCSABFactory (`createScatterHelper(WithCache)`), sabscatter/NCSABScatter.
- **Energy grid:** geomspace with 300 points by default. Emax comes from `suggestedEmax` or `determineEMax`; Emin from `determineEMin`. The NCMAT `egrid` field is honoured.
- **Cross section:** linear interpolation on the grid, 1/sqrt(E) ("1/v") extrapolation below the grid, and k/E + extender_xs(E) above Emax (continuity-anchored).
- **Sampling:**
  - `SABSamplerAtE_Alg1` implements Algorithm 1 of Cai & Kittelmann, JCP 380 (2019). Since v3.1.0 the first beta bin gets special treatment ("beta-endpoint fix", shifted by -1/3 of the bin width).
  - `SABSamplerAtE_Gridded` does cell-based rejection sampling.
  - `SABSamplerAtE_NoScatter` handles the null case.
  - Rejection uses an `EGridMargin` overlay.
- **`LegacySABAlgOpts`:** DEFAULT=-1, OVERSAMPLE10/50, NOBETAFIX, and combinations, -6..-1.
- **UCN:** NCSABUCN (`UCNScatter`, `ExcludeUCNScatter`, cfg `ucnmode`) is legacy-only. stdscatfactory throws for ucnmode with knllux != -1.

### 2.3 Next-gen path (knllux 0..6): SABProcessor plus cells (all new)
- **SABCfg** (inc/sab/NCSABCfg.hh, src/sab/NCSABCfg.cc):
  - `createConfig(sablux)` selects integration schemes (Trapez/Simpson/Romberg N, Flex5..Flex65 = adaptive Romberg with target precision 1e-3..1e-12).
  - Per-level settings: `egrid_npts` 100/140/200/300/450/600/2400; `sct_table_npts` 90 by default, 1000 at lux 6; `egrid_emin_accuracy`; `fullCellSamplingARThreshold=0.15`; `bcSamplingLargeSRatioThreshold` (marked "not used yet").
  - 100..106 switches `SIntegralInterpolation` to Linear (default Smooth).
  - 200..206 forces `ExtenderModel::FreeGas` (default SCT).
- **SABSurveyor** (NCSABSurveyor.hh/.cc): orders cells by the E/kT at which the kinematic region touches them and fully covers them. `SABCellSurvey` splits a cell into alpha regions so each is smooth for Romberg.
- **SABCellInteg** (`StdLogLinCellIntegrator`, NCSABCellInteg*):
  - `integrateFullCell`: trapezoid in beta times analytic log-lin integral in alpha, accumulated with `StableSumKahan`.
  - `integrateWithinKB`: cells crossed by the kinematic boundary, using the chosen scheme.
- **SABProcessor** (NCSABProcessor.hh/.cc, 1500 lines):
  - **Energy grid:**
    - Emax is chosen automatically from full-cell integrals of touched cells (trailing near-constant integral discarded), or taken from the egrid input or `suggestedEmax`.
    - Emin is root-found where S-integral/sqrt(E) departs from its low-E constant.
    - Grid layout: 44 geometric points over [1e-4*Emin, Emin] (justified by a 0.9 acceptance-rate argument), then geometric up to Emax; at least 80 points.
  - **`sIntegralAtE`:**
    - fully covered cells are cached incrementally;
    - an "irrelevant cell" shortcut uses a smooth quintic-smootherstep taper (the comment notes code duplicated with VDOSUtils).
  - **Phase-space integral interpolation:** PCHIP by default (Smooth), with left slope checked against the sqrt(E) law and a CalcError otherwise ("FIXME: No check of slope at Emax"). Below the grid it scales as sqrt(E); above, it is clamped. XS = sigma_b/(4*E/kT) x integral.
  - **Sampling:**
    - Look up the first grid energy at or above E (`fastLowerBoundIdx`).
    - If the full-cell acceptance rate there is at least 0.15, use full-cell sampling: pick a cell from the cumulative full-cell integrals (smallest first, for stability), sample with `FullCellSampler`, reject if outside the kinematic region at E.
    - Otherwise use bounded-cell sampling (`BoundedCellSampler`, with packed `BCSData` in a `MixedDataVector`), then the same kinematic rejection.
  - Has JSON diagnostics (`toJSON`, `toJSONProcessInfo`). Many "fixme" notes ("major cleanup needed", "FIXME: Why do we need the binary search?").
- **Cell samplers** (NCSABCellSample.hh/.cc):
  - `LogLinDistSampler` samples alpha under log-linear interpolation (`randKPowX`), with a linear fallback.
  - `FullCellSampler` picks a beta edge by its alpha integral, then samples a linear-in-beta triangle.
  - `BoundedCellSampler` uses overlay values with restricted alpha ranges.
  - `RefCellSampler` is a slow test reference.
  - NCSABKBCellSmpl (uniform sampling in cell intersect kinematic region) is marked "Fixme: Is this actually used ultimately?" and is only referenced by NCSABQuery.
- **SABRefSampler** (NCSABRefSampler): reference multi-sampling, optionally extended with a `RefSampleExtension`.
- **SABExtended** (NCSABExtended.hh/.cc): combines a SABProcessor (E <= Emax) with a SABExtender.
  - **XS above Emax:** k/E + ext_xs(E), with k = (table_xs - ext_xs)*Emax at Emax.
  - **Sampling above Emax:** a two-region mixture. The table is sampled at Emax with weight S_table(Emax); the extender with weight S_ext(E) - S_ext(Emax), rejecting points inside the Emax kinematic curve.
  - For E in (Emax, 1.1*Emax), it samples at 1.1*Emax and rejects, to avoid an O(1/epsilon) efficiency loss.
  - `create()` picks SCT when Teff is available and cfg allows, otherwise free gas. The comment requires SigmaBound = 1 ("fixme").
  - "Fixme: duplicated ... is it technically correct that we end up isotropic??"
- **Factory** (inc/sab/NCSABFactory.hh): `createSABExtendedNoCache` / `WithCache`. Cache key is (sab uid, egrid uid, knllux, teff bits).
- **SABScatterNG** (inc/sabscatter/NCSABScatterNG.hh, .cc, new): Process wrapper scaled by SigmaBound ("fixme: rename NG<->Legacy?").

### 2.4 Extenders (beyond the kernel's Emax)
- **File:** inc/sab/NCSABExtender.hh.
- **Classes:**
  - `SABFGExtender`: single-target free gas at T.
  - `SABSCTExtender` (new): takes (T, Teff, M, sigma_b, sct_table_npts) and wraps `SCTXSProvider` / `SCTSampler`.
  - `SABNullExtender`: for e.g. coherent one-phonon kernels.
- **Comments:** mention possible future multi-target extenders.

### 2.5 SABAnalyser (new; Teff and msd from a kernel)
- **Files:** inc/extd_utils/NCSABAnalyser.hh, src/extd_utils/NCSABAnalyser.cc.
- **Method:** per-alpha-row beta moments using Gaussian-approximation sum rules.
  - msd_row = -log1p(-m0)/alpha / (alpha2x_factor*kT), where m0 is the zeroth moment, i.e. 1 - exp(-2W).
  - Teff_row = (M/m_n)*T*var/(2*alpha).
  - Rows must pass a recoil-centre check: mean ~ -alpha*m_n/M within `center_tol=0.05`.
  - Other gates (hardcoded, justified by a sensitivity scan over 263 kernels): alpha >= 0.1, +-4 sigma wing coverage, 8 points per +-2 sigma, floor-marker masking, m0 in [0.05, 0.95], and an m0 trapezoid-convergence self-test.
- **Result:** the median over rows, plus an IQR/median spread.
- **Acceptance (`estimateTeffMSD`, cached):**
  - Teff: at least 10 rows and spread <= 0.1 (worst error 5.3% over 501 records).
  - msd: at least 20 rows, spread <= 0.03, and msd <= 10 Aa^2.
- **Diagnostics:** `RowStatus`, `Diagnostics`; query `["sab","analyse",cfg,label,"diag"]`.

### 2.6 Queries and filter tables
- **src/query/NCSABQuery.cc subcommands:** analyse, refeval, refsample, samplepb, sglcell, surveyor, fc, proc, list, integschemes, legacy, legacy_oversample. It notes "SABCellEval is not trustworthy".
- **Filter tables** (inc/filter/NCFilterTable.hh, new): total-cross-section-versus-wavelength piecewise-linear tables with duplicated points at discontinuities (Bragg edges), for McStas GPU filters. Not inelastic-specific.
  - "FIXME ... revisit the name 'filter' ... CrossSectionTable".
  - HEAD commit: "Add the study behind the design of the filter tables".
  - Discontinuity handling exists only here (no inelastic-specific code).

---

## 3. phys_utils, freegas, gasmix, dyninfoutils

### 3.1 NCSCTUtils (new): short collision time model
- **Files:** inc/phys_utils/NCSCTUtils.hh (long physics header with references [1]-[8]), src/phys_utils/NCSCTUtils.cc.
- **Physics:** impulse approximation (Lamb), using the GA/ENDF SCT law:
  - S_SCT(alpha_E, -|beta|) = exp(-(alpha_E - |beta|)^2*T/(4*alpha_E*Teff) - |beta|/2) / sqrt(4pi*alpha_E*Teff/T);
  - the upscatter side follows from detailed balance at T (not at Teff);
  - Teff = T reduces it exactly to free gas.
  - The header catalogues published misprints (Sutton et al.).
- **`SCTSampler`:** runs `FreeGasSampler` at Teff, scales the result by r = Teff/T, and accepts upscatter candidates with probability exp(-beta*(r-1)) (rejection sampling). Downscatter is identical event-for-event to the free-gas case.
- **`SCTXSProvider`:** sigma = sigma_FG(Teff) - sigma_b*A/(4c) x Integral_{beta>0} I(c,beta)(1 - e^{-(r-1)beta}).
  - The free-gas beta kernel I has a closed erf form (`evalFGBetaKernel`).
  - The correction is exactly (r-1)/r when saturated (c >= csat).
  - Otherwise a natural cubic spline lookup table in ln c is used (fixed-order Romberg-17 node values; PCHIP was rejected for smoothness). Accuracy scales as npts^-4: 90 points gives about 1e-5.
  - Below the table range it falls back to adaptive quadrature.
  - Teff < T is clamped if the dip is tiny.
- **Stated limitation:** total XS is changed by only O(kTeff/E); the main effect is on secondary energy/angle distributions.

### 3.2 NCFreeGasUtils and freegas
- **Files:** inc/phys_utils/NCFreeGasUtils.hh, src/phys_utils/NCFreeGasUtils.cc; inc/freegas/NCFreeGas.hh, src/freegas/NCFreeGas.cc (unchanged).
- **XS:** sigma_free*[(1 + 1/(2a^2))*erf(a) + exp(-a^2)/(sqrt(pi)*a)], with a^2 = A*E/kT and sigma_free = sigma_b*(A/(1+A))^2 (sampling paper eq. 20).
- **Sampling (`FreeGasSampler`):** quantum free gas.
  - beta is sampled with bounds via an erfc lookup table (lower/upper bounds).
  - alpha given beta samples f(x) = exp(-1/x - c*x)/sqrt(x), choosing between two methods.
  - Near the kinematic endpoint, alpha is sampled uniformly (`muIsotropicAtBeta`).
- **`FreeGas` Process:** used for `@DYNINFO freegas` and `inelas=freegas`.

### 3.3 NCKinUtils (changed)
- **Functions:** `getAlphaLimits(WithDiff)`, `getBetaLimits`, `convertAlphaBetaToDeltaEMu`, `muIsotropicAtBeta` (beta <= -E/kT*(1 - 1e-14)).
- **New:** alpha-minus Taylor series (9 terms, Horner with std::fma) for |beta| < 0.01*E/kT.
- **NCrystal alpha convention:** alpha = (E' + E - 2mu*sqrt(EE'))/kT, without division by A.

### 3.4 NCDebyeMSD (unchanged)
- **Functions:** `debyeIsotropicMSD(T_D, T, M)` (Glauber 1955 eq. 11-12: 3hbar^2/(M*k*T_D) x [1/4 + x^2 Integral_0^{1/x} t/(e^t - 1) dt], x = T/T_D), `calcDebyeMSDShape`, `calcDebyeMSDScale`, and `debyeTempFromIsotropicMSD` (inverse, slow).
- **Usage:** NCMAT v8 msd <-> Debye temperature conversion and temperature transfer.

### 3.5 GasMix
- **Files:** inc/gasmix/NCGasMixUtils.hh, src/gasmix/NCGasMixUtils.cc (2 trivial commits).
- **Contents:** `GasMixRequest` (molar or mass fractions, relative humidity, T, P or density, atomDB lines), `analyseGasMixRequest`, `saturatedVapourPressureOfWater` (JAMC-D-17-0334.1), `flattenComponentList`, and request string round-trip.
- **Dynamics:** gasmix-generated materials rely on free-gas dynamics. The in-memory writers remain at NCMAT v7.

### 3.6 dyninfoutils (changed)
- **Files:** inc/dyninfoutils/NCDynInfoUtils.hh, NCDynInfoWithTeffMSD.hh (new); src/dyninfoutils/NCDynInfoUtils.cc.
- **`extractSABDataFromDynInfo(di, vdoslux, useCache, vdos2sabExcludeFlag)`:**
  - VDOSDebye goes through the rounded-key cache;
  - direct kernels use `ensureBuildThenReturnSAB`;
  - VDOS goes through `VDOS2SABFactory` (keyed by VDOSData uid, requested Emax from egrid, exclude flag, vdoslux).
  - `vdos2sabExcludeFlag` = MODE + 4*LOW + 40000*HIGH removes coh/incoh parts of phonon orders. It is set via the unofficial hack `vdos2sab_ignorecontrib`.
- **`extractTeffFromDynInfo` (new):**
  - DI_VDOS: `VDOSEval::calcEffectiveTemperature`;
  - otherwise the explicit NCMAT v8 `effective_temperature`;
  - otherwise `SABAnalyser` for direct kernels;
  - otherwise NullOpt. VDOSDebye returns NullOpt, but vdosdebye is legacy-only anyway.
- **`extractMSDFromDynInfo` (new):**
  - VDOS: `VDOSEval::getMSD`;
  - VDOSDebye: `debyeIsotropicMSD`;
  - otherwise the explicit v8 msd;
  - for pre-v8 files (`allow_msd_estimation=false`): NullOpt;
  - otherwise `SABAnalyser`.
  - Consumed by ElIncScatter (src/elincscatter/NCElIncScatter.cc:116).
- **`DI_ScatKnlDirectWithTeffMSD::TeffMSDData`:** {explicit_teff, explicit_msd, allow_msd_estimation}.

---

## 4. stdscatfactory: choosing the inelastic model (src/stdscatfactory/NCStdScatFact.cc, +139)

**inelas=auto resolves in this order:**
1. `external` if the Info provides non-Bragg cross sections (.nxs files).
2. `dyninfo` if dynamic info exists (all .ncmat files).
3. If a temperature exists:
   - `vdosdebye` if the state of matter is solid or unknown and a Debye temperature exists;
   - otherwise `0` if HKL info exists, else `freegas`.
4. Otherwise `0`.

Values `none`, `0`, `false` and `sterile` disable inelastic scattering.

**inelas=dyninfo, per DI entry** (entries with zero scattering XS are skipped):
- **DI_ScatKnl (VDOS, VDOSDebye or direct):**
  - With knllux in -6..-1 (legacy): `createScatterHelperWithCache` then `SABScatter`, plus the optional UCN split (Refine, Remove, Only).
  - Otherwise (next-gen): Teff = `extractTeffFromDynInfo`, then `createSABExtendedWithCache(knllux, sab, teff, egrid)`, then `SABScatterNG`. Direct kernels without Teff fall back to the free-gas extender.
- **DI_Sterile:** skipped.
- **DI_FreeGas:** `FreeGas` (not allowed together with ucnmode).

**inelas=freegas:** FreeGas for every composition element.

**inelas=vdosdebye:**
- requires knllux = -1 (legacy only, "fixme"), no ucnmode, a Debye temperature, and a crystalline material with atom info;
- builds `extractSABDataFromVDOSDebyeModel` per AtomInfo, then legacy `SABScatter`.

**Elastic gating related to dynamics:**
- Crystalline materials get ElIncScatter when `incoh_elas` is set.
- Non-crystalline materials with `@STATEOFMATTER solid` get incoherent-approximation elastic scattering using msd from dyninfo. Components without msd are reported in a warning.

---

## 5. NCMAT @DYNINFO (docs/ncmat_doc.md lines 142-450 and 740-1002)

**Types:**
- **sterile:** no scattering.
- **freegas:** free-gas model.
- **scatknl:** keys `alphagrid`/`betagrid`/`sab`, or `sab_scaled` (S*exp(-beta/2), optionally beta >= 0 only), or `qgrid`/`omegagrid`/`sqw` (the loader supports `KnlType::SQW`).
  - Grids must have 5 to 65534 points.
  - A `temperature` is mandatory, must be the same in all scatknl sections, and locks the material temperature.
  - Optional `egrid`: a full list (at least 10 points), or `Elow Eupper N` with 0 meaning automatic, or a single value meaning Eupper.
  - `<v>r<count>` repeat syntax is allowed.
- **vdos:** `vdos_egrid` (2-value shorthand allowed) and `vdos_density` (at least 5 points, auto-normalised). The grid is regularised. Optional `egrid`.
- **vdosdebye:** Debye VDOS. `debye_temp` allowed in-section from v5. Implied for all elements when v1-style crystals lack DYNINFO.

**NCMAT v8:**
- **`msd <val> [at_temperature T]`:** for vdosdebye and scatknl. It is converted to a Debye temperature via the isotropic Debye model, then back to an msd at the actual temperature (NCLoadNCMAT.cc ~line 570).
- **`debye_temp` in scatknl:** allowed only for clearly solid materials.
- **`effective_temperature`:** scatknl only, must be >= `temperature`. Used by the SCT extension. Automatic estimation applies to files of any version.
- **Automatic kernel msd:** only for v8+ files, so pre-v8 physics is unchanged.
- **Crystals with scatknl:** allowed. If msd cannot be resolved, loading fails. Concurrent analysis runs via `jobs_msdcalc` ("FIXME: Revisit for NCrystal 5.0.0 ... ABI breaking DynInfo_ScatKnlDirect update").
- **Other v8 rules:** @STATEOFMATTER mandatory for non-crystals; @DEBYETEMPERATURE removed; freegas/sterile banned in crystals.

---

## 6. cfg variables (inc/cfgutils/NCCfgVars.hh, NCCfgManip.hh)

| Variable | Default | Notes |
|---|---|---|
| `temp` | -1 = 293.15K, or the data's fixed T | Material temperature. |
| `inelas` | auto | auto / none / 0 / false / sterile / external / dyninfo / vdosdebye / freegas (see Sec. 4). |
| `vdoslux` | 3 (legacy) | User description is still the legacy table: 0..5 -> 100x50 .. 3200x1600 grid, Emax 0.5..12 eV, and the "Debye reduces lux by 3" note. Hidden 1000..1005 = legacy, 2000..2006 = next-gen. |
| `knllux` (new) | -1 | "expert-only". -6..-1 legacy (LegacySABAlgOpts); 0..6 next-gen; 100..106 next-gen with linear sigma(E); 200..206 next-gen with free-gas extender. |
| `ucnmode` | "" | refine / only / remove [:threshold, default 300 neV]; legacy kernels only. |

Related elastic switches: `incoh_elas`, `coh_elas`, `sans`.

---

## 7. Tests (/workspace/tests)

**C++ apps (tests/src):**
- `app_sctext`: validates SABSCTExtender (exact Teff=T limit, XS and moments against 2-D integration, comparison with the exact kernel).
- `app_debyemsd`: Debye msd formulas.
- `app_fmavdos`: FMA reproducibility of Gn eval/evalMany.
- `app_gnconvnoisefloor`: calibrates the FFT convolution noise floor against an O(n^2) reference.
- `app_gnerange`: `estimateGnErange` / `VDOSGn::eRange` on synthetic and real Gn spectra.
- `app_perfsab`: performance of next-gen SAB (vdoslux=2000).
- `app_perfvdos`: VDOS kernel construction timing.
- `app_sabsliceinteg`: log-lin slice integration (`createTailedBreakdown`).
- `app_sabutils`: SAB utilities.
- `app_vdos2knldiag`: VDOS -> kernel diagnostics at full precision.
- `app_vdosmergegrids`: `mergeGridsWithTol`.
- `app_vdosutils`: VDOS utils, including `findABExtentWithinKB` against mpmath.
- `app_gasmix`: gas mixtures and water vapour pressure.

**Python scripts (tests/scripts):**
- `sabxs_testdata`, `sabxs_trouble`, and `long_sabxs_{simple,light,heavy,debye,highdebye,lowtemp,maxlux,vdoslux6,knllux6}`: next-gen cross-section reference monitoring. Driver: tests/pypath/NCTestUtils/sabxs.py; references: tests/data/sabxs/*.txt; vdoslux 2000-2006 cycled with knllux 0-5.
- `sabsample`: kernel sampling.
- `sabcellsample`, `sabcelleval`: cell samplers and evaluation.
- `sabsurvey`: surveyor plots.
- `sabanalyse`: Teff/msd analyser against VDOSEval truth.
- `sabcache`: caching.
- `dkelas`: direct-kernel materials, msd-based elastic, state-of-matter gating.
- `ncmatv8msd`: v8 msd/debye_temp keywords and NCMATComposer.
- `multimsd`: msd handling with multiple atom sites.
- `vdosgn`: Gn against mpmath.
- `vdosexpand`, `vdoseval`, `vdos2knl`, `vdosignore` (vdos2sab exclude), `analysevdos`, `vdosregularise`, `vdosreg2`, `vdos2ncmat`.
- `hfgpy`: hydrogen-free-gas related.
- `kinutils`: kinematic limits against mpmath.
- `filtertable`, `filtertable_api`.
- `cfgstr`: cfg variables.

---

## 8. What is new or reworked on this branch

- **New files:**
  - phys_utils/NCSCTUtils; SABSCTExtender (in NCSABExtender); extd_utils/NCSABAnalyser; dyninfoutils/NCDynInfoWithTeffMSD;
  - sab/NCSABCfg, NCSABProcessor, NCSABExtended, NCSABSurveyor, NCSABCellInteg(+FMA), NCSABCellSample, NCSABKBCellSmpl, NCSABRefSampler, NCSABIdx; sabscatter/NCSABScatterNG;
  - vdos/NCVDOSExpand, NCVDOSKnlGrid, NCVDOSLux, NCVDOSUtils(+FMA), and FMA variants of Gn and ToScatKnl;
  - filter/NCFilterTable;
  - cfg var `knllux`.
- **Reworked:**
  - VDOSGn: noise floor, tapered truncation, smooth thinning, G1 thickening, 0-anchored grids, MaxLux/direct convolution.
  - VDOSToScatKnl: next-gen grids, Stirling series.
  - VDOSEval: Teff/msd and validation.
  - NCStdScatFact: next-gen branch, Teff.
  - NCLoadNCMAT: v8 msd/Teff.
  - ElIncScatter: solid non-crystalline scatknl.
  - KinUtils: Taylor alpha-minus.
  - `reducePtsByEquidistribution(Robust)` added to NCMath.
- **Open FIXMEs:**
  - vdosdebye only works with knllux=-1;
  - ucnmode only works with legacy;
  - MaxLux is not wired to a vdoslux level;
  - the Emax lower limit of 0.1 eV is a first guess;
  - `bcSamplingLargeSRatioThreshold` is unused;
  - NG/Legacy naming;
  - SABExtended requires SigmaBound=1;
  - isotropic outcome question in SABExtended;
  - vdoslux description still documents only the legacy levels.
