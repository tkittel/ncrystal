# Cross sections and sampling from a tabulated $S(\alpha,\beta)$

## Summary

Given a kernel $S(\alpha,\beta)$ on an $(\alpha,\beta)$ grid (from a VDOS
expansion or from direct input), NCrystal provides the inelastic cross section
$\sigma(E)$ at any neutron energy and samples scattering outcomes $(E',\mu)$
with the correct joint distribution, directly from the table, at any energy,
without interpolating between pre-sampled distributions. Both rely on the fact
that the kinematically accessible region $D(E)$ grows monotonically with $E$:
the integral of $S$ over $D(E)$ gives $\sigma(E)$, and samples drawn from an
overlay distribution covering $D(E)$ and rejected outside it are exactly
distributed at $E$ (Cai et al. 2019). The 5.0.0 processing integrates $S$
exactly cell by cell with the standard log-linear interpolation in $\alpha$
and linear in $\beta$, tabulates the integral on an automatically determined
energy grid with monotone cubic interpolation in between, and samples by
picking a grid cell by its integral and a point within it, followed by the
kinematic rejection (with dedicated, tighter overlays for the cells crossed
by the kinematic boundary at very low energies). The 2.0.0-4.4.6 processing is kept as a compatibility
mode. Above the last energy of the table the extension models of chapter 33
take over.

Code: `ncrystal_core/src/sab/` (`NCSABProcessor`, `NCSABSurveyor`,
`NCSABCellSample`, `NCSABExtended`; legacy `NCSABIntegrator`, `NCSABSampler`),
`ncrystal_core/src/sabscatter/`, kinematics in
`ncrystal_core/include/NCrystal/internal/phys_utils/NCKinUtils.hh`.

## Physics

### Cross section as a phase-space integral

With the conventions of chapter 30,

$$\sigma(E) = \frac{C}{E}\, I(E), \qquad I(E) = \iint_{D(E)} S(\alpha,\beta)\,d\alpha\,d\beta, \qquad C = \frac{\sigma_b\,kT}{4},$$

where $D(E) = \{\beta \ge -E/kT,\ \alpha_-(\beta) \le \alpha \le \alpha_+(\beta)\}$
and

$$\alpha_\pm(\beta) = \frac{2E}{kT} + \beta \pm 2\sqrt{\frac{E}{kT}\Big(\frac{E}{kT}+\beta\Big)}$$

(sampling paper eqs. 3-4). Equivalently, a point $(\alpha,\beta)$ is
accessible at $E$ iff $(\alpha-\beta)^2 \le 4\alpha E/kT$, which is the
single test the samplers apply.

The joint density of outcomes at energy $E$ is
$p(\alpha,\beta|E) = C\,S(\alpha,\beta)\,\Theta_E(\alpha,\beta)/(\sigma(E)E)$
with $\Theta_E$ the indicator of $D(E)$. Energy transfer and scattering angle
follow from $E' = E + \beta kT$ and $\mu = (E + E' - \alpha kT)/(2\sqrt{EE'})$;
at the kinematic end point, where the accessible $\alpha$ interval
vanishes, the angle is isotropic.

### The rejection principle (derivation)

For $E_l < E_h$, $D(E_l) \subset D(E_h)$, hence

$$p(\alpha,\beta|E_l) = \frac{\sigma(E_h)E_h}{\sigma(E_l)E_l}\, p(\alpha,\beta|E_h)\,\Theta_{E_l}(\alpha,\beta).$$

A sample from the distribution at $E_h$ which happens to fall inside
$D(E_l)$ is therefore an exact sample at $E_l$; the acceptance rate is
$\sigma(E_l)E_l/(\sigma(E_h)E_h) = I(E_l)/I(E_h)$, which in the $1/v$ region
is $\sqrt{E_l/E_h}$ (sampling paper eqs. 9-11). Any overlay distribution
that is proportional to $S$ on $D(E_l)$ and known on a superset of it works
the same way; the 5.0.0 processing uses cell-wise overlays rather than a
full distribution at $E_h$.

### Interpolation of $S$ within a grid cell

Inside a cell $[\alpha_1,\alpha_2]\times[\beta_1,\beta_2]$ the table is
interpolated as ENDF prescribes: log-linear in $\alpha$ (linear in
$\ln S$), falling back to linear when either end value is zero, and linear
in $\beta$. The $\alpha$ integral of a log-linear segment is analytic,

$$\int_{\alpha_1}^{\alpha_2} S\,d\alpha = \frac{(\alpha_2-\alpha_1)(S_2-S_1)}{\ln(S_2/S_1)}$$

(sampling paper eq. 14), and so is its inverse CDF, which the $\alpha$
samplers use. The $\beta$ direction is integrated by Romberg quadrature
(5.0.0) or the trapezoid rule (legacy).

## Inputs and configuration

- The kernel: grids, values, temperature, $\sigma_b$, mass, and the energy
  up to which it is complete (`suggestedEmax`, set by the VDOS expansion;
  absent for direct kernels without `egrid`).
- `egrid` in the NCMAT @DYNINFO section: a complete list of energies (at
  least 10, used as given), or `Elow Eupper N` with 0 meaning automatic, or
  a single upper value (which for VDOS dynamics also sets the target
  $E_{max}$ of the expansion). Most data files leave it unset and get the
  automatic grid.
- `knllux` 0..6 (default 3): the luxury level of the processing. It sets the
  precision of the quadratures (adaptive Romberg with target precisions
  from 1e-3 at level 2 to 1e-12 at level 6; fixed low-order rules below),
  the nominal number of energy grid points (100, 140, 200, 300, 450, 600,
  2400; see below for the actual grid) and the size of the SCT lookup table
  (chapter 33). On the development branch
  the legacy processing is still reachable as negative values, and variants
  with linear instead of cubic interpolation of $I(E)$ (100..106) or with
  the free-gas instead of the SCT extender (200..206) exist for
  comparisons.
- `ucnmode` (legacy processing only; chapter 30).

## Implementation (5.0.0 processing)

### Surveying the grid

Every grid cell is classified by the two energies at which the growing
region $D(E)$ first touches it and first covers it completely, and cells are
sorted by these energies. At a given energy the touched cells form a prefix
of the list; cells touched but not covered are crossed by the kinematic
curves, and for those the cell is split into $\alpha$ regions in each of
which the integrand is smooth, so that the quadrature converges quickly.
Grids are limited to 65534 points per axis.

### Cell integrals and the energy grid

- The integrals of $S$ over whole cells (analytic in $\alpha$, linear in
  $\beta$) are computed once. $I(E)$ at a grid energy is the sum over the
  covered cells plus the quadrature over the crossed parts of the
  touched-but-not-covered cells; walking the grid upwards, the covered sum
  is carried over, and cells whose maximal possible contribution is
  negligible are skipped (through a smooth taper rather than a hard
  threshold, so that FFT noise in a VDOS-derived kernel cannot flip the
  decision).
- $E_{max}$ is the energy up to which the kernel is complete when known
  (from the VDOS expansion or `egrid`); otherwise it is taken, with a small
  safety margin, where the full integrals of the cells, accumulated in the
  order in which $D(E)$ touches them, have essentially reached their total,
  i.e. where the table has run out of content (capped at $10^4\,kT$ and
  1 keV).
- $E_{min}$ is the energy below which $I(E)/\sqrt{E}$ is constant to within a
  small fraction of its overall variation, i.e. where the $1/v$ law has set
  in, found by a root search to a level-dependent accuracy.
- The grid has 43 points geometrically spaced over the four decades below
  $E_{min}$ (a ratio of about 1.24 between neighbours, so that the
  $\sqrt{E}$ growth of $I$ keeps the sampling acceptance rate between
  neighbouring points near 0.9), followed by points geometrically spaced
  from $E_{min}$ to $E_{max}$, in total one fewer than the level's nominal
  number (an explicit `N` below 80 is raised to 80; FIXME: the "one
  fewer" is an off-by-one in the construction, to be fixed). Below the first
  point
  $I(E) \propto \sqrt{E}$, so $\sigma \propto 1/\sqrt{E}$; above the last
  point the extension of chapter 33 applies; in between $I(E)$ is
  interpolated with a monotone cubic (PCHIP) spline. The spline determines
  its own end slopes; the slope at the low end is checked against the
  $\sqrt{E}$ form (a kernel failing this check is rejected), while values
  above $E_{max}$ never use the spline.

### Sampling

Each energy grid point $E_h$ carries a sampling recipe, and a neutron at $E$
uses the recipe of the first grid point $E_h \ge E$ (below the lowest grid
point that of the lowest one; above $E_{max}$ the extension of chapter 33
takes over). A point is drawn from the interpolated $S$
within a cell exactly: one of the two $\beta$ edges is picked in proportion
to its $\alpha$ integral, $\alpha$ is drawn from the log-linear profile
along that edge, and $\beta$ from a linear density that is maximal at the
picked edge and zero at the other (the bilinear-in-$\beta$ interpolation is
the sum of these two terms). The final candidate is accepted iff
$(\alpha-\beta)^2 \le 4\alpha E/kT$, which makes the result exact at $E$
whatever the overlay. Which of two recipes a grid point uses depends on the
full-cell acceptance rate $I(E_h)$ divided by the summed full integrals of
all cells touched at $E_h$:

- *Full-cell sampling* (acceptance rate at least 0.15): a cell is picked
  among the touched cells in proportion to its full integral, a point is
  drawn in the whole cell, and it is kept if inside $D(E)$, otherwise the
  cell is picked anew. Covered cells dominate, so this is efficient over
  almost the whole grid. How far down it stays efficient depends on the
  size of the kernel's grid cells near $\alpha=\beta=0$, and so on the
  material, the temperature and, for VDOS-derived kernels, `vdoslux`: at
  $E = 10^{-3}\,kT$ the acceptance rate is about 0.2-0.3 at vdoslux 1,
  0.4-0.6 at vdoslux 3 and 0.65-0.85 at vdoslux 5 for hydrogen in
  polyethylene and for aluminium.
- *Bounded-cell sampling* (lower acceptance rates, i.e. at very low
  energies, where $D(E_h)$ is a thin sliver along $\alpha \approx \beta$
  and crosses cells much larger than the part of them inside it; for the
  same two materials below about $3\times10^{-4}\,kT$ at vdoslux 1,
  $5\times10^{-5}\,kT$ at vdoslux 3 and $5\times10^{-6}\,kT$ at vdoslux 5):
  a cell is picked in proportion to its integral within
  $D(E_h)$, and points are drawn in that cell until one falls inside
  $D(E_h)$, which makes the pair of steps an exact draw from $S$ on
  $D(E_h)$; the result is then kept if also inside $D(E)$, otherwise both
  steps are repeated. Covered cells, and crossed cells of which at least a
  fifth of the integral lies within $D(E_h)$, are drawn in the whole cell;
  the other crossed cells use an overlay restricted to the $\alpha$ range
  that $D(E_h)$ cuts out of the cell along each $\beta$ edge, precomputed
  and stored compactly per grid point. Cells whose possible contribution is
  negligible (the tapered shortcut above) are left out.

A slow reference implementation of both the cross section and the sampling
at any energy exists for tests and is exposed through the JSON query
`["sab","refeval",...]`.

### Process wrapping

The processed kernel (processor plus extender, chapter 33) is wrapped as an
isotropic scattering process scaled by $\sigma_b$ and the element fraction;
processes with the same kernel merge by adding their scales, and processed
kernels are cached on the kernel identity, energy grid, level and effective
temperature, so a kernel shared between materials or phases is processed
once.

## Implementation (legacy, 2.0.0-4.4.6)

The legacy processing is the one described in the sampling paper: at each
of 300 geometric energy grid points, the table is integrated analytically
along $\alpha$ within the accessible range of each $\beta$ row and by the
trapezoid rule along $\beta$ (with a $\sqrt{x}$ treatment of the first bin
above $-E/kT$, added in 3.1.0 after a study of low-energy artefacts), and
the per-row conditional distributions are stored; cross sections are
interpolated linearly in $E$. Sampling follows Algorithm 1 of the paper:
sample $\beta$ from the distribution at the next-higher grid energy, sample
$\alpha$ from the two bracketing row distributions with a common percentile
and interpolate, accept if within $[\alpha_-,\alpha_+]$. Above $E_{max}$ the
sampler mixes table and extender samples, correcting for the mismatch of
their integrals at $E_{max}$, and fails if the mismatch is too large.
`ucnmode` integrates the kernel exactly over the part of $D(E)$ leaving the
neutron below the UCN threshold, giving a dedicated UCN-production process
for biased simulations.

## Validation

- The processor against the reference implementation at every grid energy
  and in between, including sampled distributions
  (`tests/scripts/sabcelleval.py`, `sabcellsample.py`, `sabsample.py`).
- The `sabxs` reference cross sections over all levels for the stdlib and
  synthetic materials.
- The published validation of the sampling method: total cross sections
  within 0.05% (0.1 eV) and 0.2% (1 eV) of the analytic free gas,
  statistically identical energy-transfer distributions, water-sphere
  spectra within a fraction of a percent of MCNP6 (sampling paper); UCN
  production cross sections in superfluid helium (3.1.0 release notes).

## Approximations and limitations

- Results are exact for the *interpolated* table; the physics is as good as
  the table and the interpolation scheme, which is why the grid
  construction of chapter 31 matters.
- The interpolation of $I(E)$ between energy grid points is the only
  approximation in the cross section; it is monitored by the `sabxs` tests
  at the 1e-4 to 1e-6 level depending on the level.
- Sampling is less efficient at very low energies ($E \ll kT$), where
  $D(E)$ is thin, although the bounded-cell overlays limit the cost. Below
  the lowest grid point ($10^{-4} E_{min}$) the recipe of that point is
  used, so the acceptance rate there falls off as $\sqrt{E/E_{\rm lowest}}$.
- FIXME (development state, to be removed once fixed): the bookkeeping that
  selects bounded-cell sampling double-counts newly covered cells, so the
  selection condition can never be met and all grid points currently use
  full-cell sampling (exact, but with acceptance rates down to about 1% at
  the lowest grid points for hydrogen in polyethylene); an explicit `egrid`
  list currently needs more than 10 points.
- Direct kernels that end where $S$ is still large should carry an
  `egrid`, since the automatic $E_{max}$ then has to guess.

## New or changed in 5.0.0

Release 4.4.6 already provided exact-at-energy rejection sampling from
kernels (Algorithm 1 since 2.0.0, the low-energy refinements and UCN
production since 3.1.0, reference evaluators since 4.4.0). The 5.0.0
processing, now the default, replaces the row-based Algorithm 1 with exact
cell-wise phase-space integrals, automatic energy grids with monotone cubic
interpolation, and cell-based sampling, removing the interpolation and
mismatch artefacts of the legacy scheme; it is faster in the hot paths and
bit-reproducible across platforms (chapter 70). `ucnmode` remains available
only with the legacy processing and is scheduled for removal after 5.0.0.

## References

- X.-X. Cai, T. Kittelmann, E. Klinkby, J. I. Marquez Damian, J. Comput.
  Phys. 380 (2019) 400, doi:10.1016/j.jcp.2018.11.043.
- ENDF-6 Formats Manual (ENDF-102), File 7 interpolation rules for
  $S(\alpha,\beta)$.
- NCrystal 3.1.0 release announcement (first-bin treatment, UCN process).
- F. N. Fritsch, R. E. Carlson, "Monotone piecewise cubic interpolation",
  SIAM J. Numer. Anal. 17 (1980) 238 (PCHIP).
