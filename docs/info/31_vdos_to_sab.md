# From a phonon density of states to $S(\alpha,\beta)$

## Summary

For solids, NCrystal generates the scattering kernel of each element at load
time from a one-dimensional phonon density of states (VDOS), for the
requested temperature. The method is Sjolander's multi-phonon expansion in
the incoherent and Gaussian (harmonic) approximations: the one-phonon
function $G_1$ follows from the VDOS, the $n$-phonon functions $G_n$ are
obtained by repeated self-convolution, and
$S(\alpha,\beta) = \sum_n e^{-x} x^n/n!\; kT\, G_n(\beta kT)$ with
$x \propto \alpha$ set by the mean-squared displacement. The implementation
decides how many orders are needed to cover the kinematic region up to a
target energy (5 eV), builds the $(\alpha,\beta)$ grid from the shapes of the
$G_n$, and evaluates the sum with care for very high orders (thousands are
common for light atoms at high temperature). An idealised Debye density of
states is used when only a Debye temperature is known. The expansion also
provides the mean-squared displacement and the effective temperature of the
atom.

Code: `ncrystal_core/src/vdos/` (`NCVDOSEval`, `NCVDOSGn`, `NCVDOSExpand`,
`NCVDOSToScatKnl`, `NCVDOSKnlGrid`), the FFT convolution in
`ncrystal_core/src/utils/NCFastConvolve.cc`, the Debye VDOS in
`ncrystal_core/src/dyninfoutils/NCDynInfoUtils.cc`.

## Physics

### The density of states

The input is a density of states $\rho(E)$ as a function of phonon energy
$E = \hbar\omega$, tabulated on an energy grid, from which NCrystal defines
the function it works with:

- linear interpolation between the tabulated points;
- below the first point $E_{min}$, $\rho(E) = \rho(E_{min})\,(E/E_{min})^2$,
  the Debye-like behaviour of acoustic phonons in a solid (so the first point
  should be placed where this is a fair description; it must be at least
  1e-5 eV);
- zero above the last point $E_{max}$;
- normalised to unit integral over $[0, E_{max}]$ (the overall scale of the
  input is irrelevant).

Internally the grid is equidistant and extends to exactly zero with a whole
number of bins; input not of that form is re-parameterised on loading.

### One-phonon function, $\gamma_0$, displacement and effective temperature

Following Sjolander (1958), three integrals of the VDOS set the scales of
the problem:

$$\gamma_0 = \int_0^{E_{max}} \frac{\rho(E)}{E}\coth\frac{E}{2kT}\,dE$$

(Sjolander eq. II.3 at $t=0$; NCrystal keeps the energies dimensional where
Sjolander normalises to $E_{max}$), from which the isotropic mean-squared
displacement is

$$\langle u^2\rangle = \frac{\hbar^2}{2M}\gamma_0$$

(the displacement along one axis, entering the Debye-Waller exponent as
$2W = q^2\langle u^2\rangle$), and the effective temperature, defined by
$kT_{eff}$ = mean energy per vibrational degree of freedom including
zero-point motion,

$$kT_{eff} = \int_0^{E_{max}} \rho(E)\,\frac{E}{2}\coth\frac{E}{2kT}\,dE \;\ge\; kT$$

(the ENDF File 7 effective temperature, used by the short-collision-time
model of chapter 33; $T \le T_{eff} \le T + E_{max}/2k$).

The one-phonon function is Sjolander's $\bar G_1$ (eq. II.29) made
asymmetric with the detailed-balance factor:

$$G_1^{sym}(E) = \frac{\rho(|E|)}{2\gamma_0\,|E|\,\sinh(|E|/2kT)}, \qquad
G_1(E) = G_1^{sym}(E)\,e^{-E/2kT} = \frac{\rho(|E|)}{\gamma_0\,E\,(e^{E/kT}-1)}.$$

Here positive $E$ is a neutron energy *gain* (suppressed by detailed
balance), negative $E$ an energy loss. $G_1$ has unit area; the parabolic
part below $E_{min}$ gives $G_1 \propto u/\sinh u$ with $u = E/2kT$.

### Multi-phonon functions (derivation)

Sjolander's eq. II.27 defines $G_n$ recursively as the convolution of
$G_{n-1}$ with $G_1$, so that $G_n$ is the $n$-fold convolution of $G_1$: the
energy distribution of a process exchanging $n$ phonons. All $G_n$ have unit
area and satisfy detailed balance. In the Gaussian approximation the
$\alpha$ dependence of the $n$-phonon term factorises as a Poisson weight:

$$S(\alpha,\beta) = \sum_{n=1}^{\infty} \frac{e^{-x} x^n}{n!}\; kT\,G_n(\beta kT), \qquad x = 2W = q^2\langle u^2\rangle = \alpha \cdot \frac{2 m_n kT}{\hbar^2}\langle u^2\rangle$$

(Sjolander eq. II.28 rewritten in $(\alpha,\beta)$ and without the factor
$\sigma kT/4E$ that NCrystal keeps outside $S$; the factor $kT$ converts the
energy density $G_n$ into a density in $\beta$). The $n = 0$ term is the
elastic line $e^{-x}\delta(\beta)$ and is deliberately left out: Bragg
diffraction and incoherent elastic scattering are separate processes with
the same Debye-Waller factor $e^{-x}$.

### The Debye model VDOS

When only a Debye temperature $T_D$ is known, the density of states is the
Debye spectrum $\rho(E) \propto E^2$ for $E \le kT_D$ and zero above, and it
goes through the same expansion as any VDOS. The Debye model is a crude
description of real phonon spectra, but the resulting cross sections and
angular distributions are reasonable even though the energy transfers are
not (BeO comparison in the 2.0.0 release notes); it is the fallback for
crystals specified with Debye temperatures only and for the quick `solid::`
materials.

## Inputs and configuration

- NCMAT @DYNINFO `type vdos` with `vdos_egrid` (two values = end points of
  an equidistant grid, else one value per density point) and `vdos_density`
  (at least 5 values, any normalisation); `type vdosdebye` with `debye_temp`
  (or, before NCMAT v8, the @DEBYETEMPERATURE section). The optional `egrid`
  field sets the energy grid of the processed kernel and, through its upper
  end, the target $E_{max}$ of the expansion (chapter 32).
- Material temperature `temp` (the expansion works from 0.1 K upwards).
- `vdoslux`: the luxury level of the expansion (table below).

## Implementation

### Expansion order and target energy

The expansion adds orders until the region of the $(\alpha,\beta)$ plane
where the last order contributes no longer intersects the kinematic region
$D(E_{target})$. The contributing region of order $n$ is the rectangle where
$G_n$ and the Poisson factor $e^{-x}x^n/n!$ both exceed a level-dependent
relative threshold of their maxima (1e-6 at the lowest level, 1e-9 at the
default, 1e-13 at the highest). The target energy is 5 eV (or the upper end
of an explicit `egrid`), so that the extension model of chapter 33 only
takes over where the free-gas-like behaviour is a fair description.
Resource limits (total and per-order spectrum sizes, growing with the
level) stop the expansion earlier for extreme inputs; the target energy is
then lowered until the available orders suffice, with a warning, and the
expansion fails with an error if the covered range would fall below 0.1 eV,
since extrapolating from a lower $E_{max}$ was found to give errors of
10-20%. The number of orders is typically tens to a few hundred for heavy
atoms at room temperature and thousands for hydrogen or hot materials; the
covered energy is recorded with the kernel as its `suggestedEmax`.

### Computing the $G_n$

- $G_1$ is tabulated on an equidistant grid from $-E_{max}$ to $+E_{max}$
  with a node at exactly zero. The grid is refined beyond the input grid so
  that the bin width is small compared with $kT$ (which the detailed-balance
  factor requires at low temperatures) and so that linear interpolation
  between the points reproduces $G_1$ (which matters for VDOS curves with
  structure narrow compared with both $kT$ and the input grid).
- Higher orders are obtained by FFT convolution of lower ones, pairing the
  largest available orders ($G_{2n} = G_n * G_n$, $G_{2n+1} = G_n * G_{n+1}$),
  so a given order needs only about $\log_2 n$ convolutions, which limits the
  accumulation of numerical error; the orders are computed concurrently in a
  thread pool, and every spectrum is renormalised to unit area.
- The spectra grow with the order and are kept compact in two ways: the
  far wings, below a fixed relative level of the maximum, are truncated,
  with the truncation edge located by a smooth fit rather than a hard
  threshold so that FFT round-off cannot move it; and from order 4 the bin
  width is doubled whenever the spectrum is smooth enough for linear
  interpolation on the coarser grid to reproduce it. Values that are pure
  FFT round-off are set to zero.
- All spectra live on lattices anchored at $E = 0$, so that truncation and
  thinning never shift a grid by a fraction of a bin; together with the
  explicit floating-point measures of chapter 70 this makes the expansion
  bit-identical across platforms.
- The legacy configuration reproduces the 2.x-4.x expansion (fixed point
  counts instead of smoothness-gated thinning, no $G_1$ refinement, no
  noise handling), and a direct $O(N^2)$ convolution exists as a reference
  for tests.

### Evaluating the sum

The Poisson factors $e^{-x}x^n/n!$ are evaluated recursively for small $n$
and, for large $n$, from a Stirling series for $n!$ with all cancellations
arranged inside the exponent, so that orders in the thousands and $x$ values
of hundreds are handled without overflow or loss of precision. For each
order the $G_n$ is interpolated linearly onto the $\beta$ grid and
accumulated into the table. An optional per-order weight lets plugins scale
or remove orders (e.g. the coherent part of the one-phonon term for models
that treat coherent one-phonon scattering themselves).

### The $(\alpha,\beta)$ grid

The grid dimensions are set by the level (table below). Points are placed
where the kernel has structure: $\beta$ points by equidistribution over a
weighted sum of the $G_n$ (weights falling roughly as $1/n^2$, so that the
sharp low-order features dominate), $\alpha$ points linearly spaced in
$\sqrt{x}$ (since the $n$-phonon terms tend to Gaussians of width
$\sqrt{n}$ in $x$, this gives uniform relative interpolation errors), and
extra points along the line $\alpha = \beta$ onto which the kinematic region
collapses as $E \to 0$. Both grids contain exactly zero and have exactly the
requested size. The legacy grids used fixed region schemes instead.

### Luxury levels

| `vdoslux` | grid $n_\alpha \times n_\beta$ | contribution threshold |
|---|---|---|
| 0 | 60 x 130 | 1e-6 |
| 1 | 110 x 240 | 1e-7 |
| 2 | 170 x 400 | 1e-8 |
| 3 (default) | 280 x 650 | 1e-9 |
| 4 | 470 x 1080 | 1e-11 |
| 5 | 770 x 1770 | 1e-12 |
| 6 | 2000 x 4000 | 1e-13 |

(On the development branch these are reached as `vdoslux` 2000..2006, with
0..5 still selecting the legacy expansion whose grids were
$100\cdot 2^{lux} \times 50\cdot 2^{lux}$ and target energies 0.5-12 eV.)
Initialisation takes of the order of 0.1 s per material at the default
level, and kernels are cached across materials and phases sharing the same
VDOS and element.

The Python API exposes the steps (`NCrystal.vdos.analyseVDOS`, `extractGn`,
`extractKnl`, the plotting helpers) and the JSON query `["vdos","expand"]`
does the same from the command line.

## Validation

- $G_n$ against high-precision (mpmath) references; the truncation-edge
  estimator and the FFT noise model against brute-force references
  (`tests/scripts/vdosgn.py`, `tests/src/app_gnerange`,
  `app_gnconvnoisefloor`, `app_vdosutils`).
- The `sabxs` reference tests monitor the cross sections of many stdlib and
  synthetic materials over all levels, including low temperatures, light and
  heavy atoms, high Debye temperatures and unusual VDOS curves collected from
  plugins (`tests/data/vdos_*.ncmat`).
- Physics comparison: BeO cross sections from VDOS at eight temperatures
  against the ENDF/B-VIII.0 kernels (2.0.0 release notes).

## Approximations and limitations

- Incoherent approximation and Gaussian (harmonic, isotropic) approximation:
  no anharmonicity, no dispersion, no anisotropic displacements; one VDOS
  per element (no partial VDOS per crystallographic site).
- The $E^2$ extrapolation below the first VDOS point is a solid-state
  assumption; materials with other low-energy behaviour must provide it
  explicitly or use direct kernels.
- The kernel describes the kinematic region up to its $E_{max}$ (5 eV by
  default); above it, the extension model applies.
- The expansion can fail for extreme combinations of VDOS, temperature and
  mass (the error suggests a lower temperature or higher `vdoslux`).

## New or changed in 5.0.0

Release 4.4.6 already generated kernels from VDOS curves at any temperature
with this expansion (since 2.0.0), with multithreaded convolutions and
caching across materials. The 5.0.0 expansion (the new default) adds:

- $G_1$ refinement that keeps detailed balance exact at low temperatures and
  resolves narrow structure (previously up to 10% cross-section errors in
  cold hydrogenous materials), and thinning gated on smoothness rather than
  point counts (previously up to tens of percent for VDOS curves with narrow
  features at low temperature).
- Noise-robust truncation and FFT noise handling, anchored grids, and the
  new $(\alpha,\beta)$ grid construction.
- Resource-based limits with a fixed 5 eV target, warnings whenever the
  range is reduced and an error below 0.1 eV, instead of silently reduced
  ranges; no reduced level for Debye-model kernels.
- The effective temperature from the VDOS, used by the SCT extension.
- Bit-identical results across platforms (chapter 70).

## References

- A. Sjolander, "Multi-phonon processes in slow neutron scattering by
  crystals", Arkiv for Fysik 14 (1958) 315 (eqs. II.3, II.27, II.28, II.29;
  NB: eq. II.1 is missing a factor $\hbar$, and Sjolander's energies are
  normalised to $E_{max}$).
- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, section 2.4.
- R. E. MacFarlane et al., "The NJOY Nuclear Data Processing System, Version
  2016", LA-UR-17-20093 (the LEAPR module uses the same expansion to produce
  ENDF File 7 kernels).
- NCrystal 2.0.0 release announcement (BeO comparison, `vdoslux`).
