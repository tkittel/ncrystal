# Beyond the kernel: free-gas and short-collision-time extension, and kernel analysis

## Summary

A tabulated kernel covers the kinematic region only up to some energy
$E_{max}$ (a few eV for VDOS-generated kernels, whatever the evaluator chose
for direct kernels). Above it NCrystal keeps the table for the part of phase
space it covers and models the rest with an analytical "extender": either
the short-collision-time (SCT) model, a free gas at the atom's effective
temperature $T_{eff}$ with detailed balance at the true temperature $T$
(the standard ENDF/NJOY recipe), or the plain free gas at $T$ (the
2.0.0-4.4.6 behaviour, exactly the SCT model with $T_{eff} = T$). The SCT
model needs $T_{eff}$, which comes from the VDOS for expanded kernels, from
the NCMAT v8 `effective_temperature` keyword for direct kernels, or is
estimated from the kernel itself through the moment sum rules of the
Gaussian approximation; the same analysis estimates the mean-squared
displacement that gives direct-kernel solids their elastic physics. Cross
sections above $E_{max}$ are $k/E + \sigma_{ext}(E)$ and sampling mixes
table and extender by their phase-space integrals, so the transition is
continuous and unbiased.

Code: `ncrystal_core/src/phys_utils/NCSCTUtils.cc`,
`ncrystal_core/src/sab/NCSABExtender.cc`, `NCSABExtended.cc`,
`ncrystal_core/src/extd_utils/NCSABAnalyser.cc`,
`ncrystal_core/src/dyninfoutils/NCDynInfoUtils.cc`.

## Physics

### Why a free-gas-like model at high energy

At high momentum transfer the collision is over so quickly that the struck
atom reacts as if free, moving with the velocity it had at that instant (the
impulse approximation, Lamb 1939). For a bound atom the velocity spread is
not Maxwellian at $T$: each vibrational mode of energy $\hbar\omega$ carries
a mean energy $(\hbar\omega/2)\coth(\hbar\omega/2kT)$, zero-point motion
included, half of it kinetic. Weighted with the density of states this
defines the effective temperature

$$kT_{eff} = \int \rho(\hbar\omega)\,\frac{\hbar\omega}{2}\coth\frac{\hbar\omega}{2kT}\,d(\hbar\omega) \;\ge\; kT,$$

so that the mean kinetic energy of the atom is $\frac32 kT_{eff}$ (the ENDF
File 7 effective temperature; chapter 31 gives the VDOS evaluation). At
large transfers a bound atom therefore scatters like a free gas at
$T_{eff}$, with the Doppler width of the recoil ridge set by $T_{eff}$. A
literal free gas at $T_{eff}$ would however obey detailed balance at
$T_{eff}$ instead of $T$ and strongly overestimate upscattering (Sutton et
al. 2009). The standard fix, due to the Gulf General Atomic team (Borgonovi
1970) and used in NJOY's THERMR, keeps the free-gas-at-$T_{eff}$ form on the
downscatter side, where the impulse approximation is right, and *defines*
the upscatter side by detailed balance at $T$.

### The SCT scattering law

In ENDF variables ($\alpha_E = \alpha/A$ with $A = M/m_n$, symmetric $S$),

$$S_{SCT}(\alpha_E, -|\beta|) = \frac{\exp\big(-\frac{(\alpha_E-|\beta|)^2\,T}{4\alpha_E T_{eff}} - \frac{|\beta|}{2}\big)}{\sqrt{4\pi\alpha_E T_{eff}/T}},$$

with $S_{SCT}(\alpha_E,+|\beta|)$ from the symmetry of the symmetric $S$,
i.e. detailed balance at $T$. Several published versions of this formula
contain misprints; the form above is the intended "GA" one (Sutton et al.
catalogue the variants). For $T_{eff} = T$ it is exactly the free-gas law,
so the free-gas extension is the special case of SCT. NCrystal works with
the asymmetric $S$ and its own $\alpha$ convention (chapter 30), and the
cross section is $\sigma = \sigma_b/(4c)\int S\,d\alpha\,d\beta$ over the
accessible region, $c = E/kT$.

### Sampling and cross section in the SCT model (derivation)

Write the transfers in units of $kT_{eff}$: with $r = T_{eff}/T$,
$(\alpha,\beta) = (r\alpha^*, r\beta^*)$. Substitution shows that on the
downscatter side the SCT law is identical, event for event, to the free-gas
law at $T_{eff}$ in its own units, and that the kinematic region maps onto
itself under this scaling. The upscatter side differs only by the
$\alpha$-independent factor $\exp(-\beta^*(r-1)) \le 1$. Hence:

- Sampling: run the free-gas sampler at $T_{eff}$, scale the result by $r$,
  and accept upscatter candidates with probability $e^{-\beta^*(r-1)}$
  (rejection; the acceptance rate is at least
  $\sigma_{SCT}/\sigma_{FG}(T_{eff})$, never small at the energies where the
  model is used).
- Cross section: by the same decomposition

$$\sigma_{SCT}(E) = \sigma_{FG,T_{eff}}(E) - \frac{\sigma_b A}{4c}\int_{\beta>0} I(\beta)\,\big(1-e^{-(r-1)\beta}\big)\,d\beta, \qquad c = E/kT_{eff},$$

  where $I(\beta)$ is the free-gas "$\beta$-kernel", the asymmetric
  free-gas $S$ integrated over the allowed $\alpha$ range, which with
  $u = \sqrt{\alpha_E}$ has the closed form

$$I(\beta) = \tfrac12\Big[\operatorname{erf}\big(y(u)\big) + e^{-\beta}\operatorname{erf}\big(z(u)\big)\Big]_{u_-}^{u_+}, \quad y(u) = \frac{u^2+\beta}{2u},\ z(u) = \frac{u^2-\beta}{2u},$$

  from the antiderivative $\int e^{-(u/2+\beta/2u)^2}du =
  \frac{\sqrt\pi}{2}[\operatorname{erf}(u/2+\beta/2u) +
  e^{-\beta}\operatorname{erf}(u/2-\beta/2u)]$. For $\min(cA, c/A)$ large
  the error functions saturate (corrections of order
  $\operatorname{erfc}\sqrt{cA} + \operatorname{erfc}\sqrt{c/A}$), leaving
  $I(\beta) = e^{-\beta}$ and the correction integral exactly $(r-1)/r$.
  Since $c \gtrsim 50$ whenever the model is used above a typical
  $E_{max}$ of several eV, this closed form covers most evaluations to
  double precision; the transition region at lower $c$ is covered by a
  spline table of the integral versus $\ln c$ built at initialisation.

The total cross section is only mildly affected by the difference between
SCT and a free gas at $T$ (an $O(kT_{eff}/E)$ effect); the important
improvement is in the secondary energy and angle distributions, whose
Doppler widths differ by the full factor $r$ (for hydrogen in room-temperature
solids $r$ is typically 3-5, and much larger at cryogenic temperatures;
FIXME: these $r$ values are an estimate, not a survey; a benchmark of the
SCT extension against evaluated data and across the stdlib is needed at
some point).

### Combining table and extender

Inside $D(E_{max})$ the table is used and outside it the extender; both are
proportional to $S$ on their domains, so with $I$ the phase-space integral
(chapter 32),

$$I(E) = I_{table}(E_{max}) - I_{ext}(E_{max}) + I_{ext}(E), \qquad
\sigma(E) = \frac{k}{E} + \sigma_{ext}(E), \quad k = \big(\sigma_{table}(E_{max}) - \sigma_{ext}(E_{max})\big)E_{max}.$$

Sampling at $E > E_{max}$ decides between "inside $D(E_{max})$" (table
sample at $E_{max}$, which is exact there since $D(E_{max}) \subset D(E)$)
and "outside" (extender sample at $E$, rejected if inside $D(E_{max})$)
with the probability $P_{inside} = I_{table}(E_{max})/I(E)$. Just above
$E_{max}$ the outside part is a thin shell which would be expensive to
reject-sample directly, so in that band the sampling is done at a slightly
higher energy and the result rejected if outside $D(E)$.

### Estimating $T_{eff}$ and the displacement from a kernel

In the incoherent/Gaussian approximation every $\alpha$ row of the
asymmetric kernel obeys exact moment sum rules, which make per-row
estimators possible:

- Wick's second-moment rule: $\mathrm{Var}_\beta(\alpha) = 2\alpha\,(m_n/M)\,(T_{eff}/T)$
  (in NCrystal's $\alpha$ convention), giving $T_{eff}$ from the variance
  of each row; the first-moment rule fixes the recoil centre at
  $-\alpha m_n/M$ and serves as a self-check of the row.
- The zeroth-moment (Debye-Waller) rule for kernels with the elastic line
  excluded: $N(\alpha) = \int S\,d\beta = 1 - e^{-c\alpha}$ with
  $c = (2m_n kT/\hbar^2)\langle u^2\rangle$ (the $x/\alpha$ ratio of the
  VDOS expansion), giving the displacement from the deficit $1 - N$. Kernels
  without such a deficit (liquids, $N \approx 1$) correctly yield "not
  inferable".

The analyser treats each tabulated $\alpha$ row as an independent
measurement: it computes the moments of the row, keeps only rows that are
resolved by the grid, covered by real data out to several standard
deviations and (for $T_{eff}$) consistent with the recoil-centre rule, which
excludes the small-$\alpha$ rows of stored kernels where finite precision
breaks the detailed-balance cancellation, as well as genuinely non-Gaussian
kernels such as cryogenic ortho-hydrogen. The final values are medians over
the admissible rows, with their spread as a robustness metric; an estimate
is accepted only if enough rows contribute and the spread is small, with a
stricter policy for the displacement (a refusal merely means no elastic
component, but a wrong value would add wrong physics) than for $T_{eff}$ (a
refusal falls back to the free gas at $T$, usually a far larger error than
any credible estimate). Against kernels expanded from VDOS curves of known
properties over a wide range of masses and temperatures, accepted estimates
are within about 1% for $T_{eff}$ and 5% for the displacement.

## Inputs and configuration

- The SCT extension is used by the 5.0.0 kernel processing whenever an
  effective temperature is available (a variant with the free-gas extender
  exists for comparisons); the legacy processing always extends with the
  free gas at $T$.
- NCMAT v8 `scatknl` keywords (`docs/ncmat_doc.md`): `effective_temperature`
  (used directly for SCT), `msd` or `debye_temp` (used directly for the
  displacement). Without them the estimates from the kernel are used, with
  one exception: displacement estimation is disabled for pre-v8 NCMAT files,
  so that the physics of existing files does not change (whereas $T_{eff}$
  estimation applies to them too, since it only changes the extension above
  the table, which was an approximation anyway).
- `vdos` and `vdosdebye` dynamics always provide $T_{eff}$ and the
  displacement from the density of states (chapter 31).
- The JSON query `["sab","analyse",...]` exposes the analyser's per-row
  moments and estimates.

## Implementation

The SCT cross section is the free-gas cross section at $T_{eff}$ minus the
upscatter correction, evaluated from the closed form $(r-1)/r$ where it
applies and from the spline table below; the sampler wraps the free-gas
sampler at $T_{eff}$ (chapter 34) with the $r$ scaling and the upscatter
acceptance. $T_{eff} = T$ reproduces the free gas exactly. The extender is
combined with the processed table as described above and a null extender
exists for kernels that must get nothing outside the table (e.g. coherent
single-phonon models from plugins). The resolution of $T_{eff}$ and
displacement per dynamics type is the table of chapter 30; when neither a
value nor a confident estimate exists, the extension falls back to the free
gas and the element gets no elastic component (with a warning for solids).

## Validation

- `tests/src/app_sctext`: $T_{eff} = T$ equivalence with the free gas, the
  closed-form $\beta$-kernel and the correction integral against
  quadrature, continuity between regimes, sampling distributions against the
  analytic law.
- `tests/scripts/sabanalyse.py`: estimates on VDOS-expanded kernels against
  the VDOS truth, on the stdlib direct kernels and on synthetic kernels; the
  acceptance criteria were tuned on a larger set of ENDF- and JENDL-converted
  kernels kept outside the repository.
- The `sabxs` tests with the free-gas and SCT extenders document the size of
  the effect per material.

## Approximations and limitations

- The SCT model is the impulse approximation plus the GA detailed-balance
  prescription; it is the ENDF/NJOY standard, not an exact theory, with
  percent-level accuracy in cross sections, better where $c$ is large.
- One effective temperature per element; molecular kernels describing
  several atoms get a single $T_{eff}$ and mass.
- The analyser needs a Gaussian-type kernel with resolved recoil ridge and
  real data in the wings; coarse or truncated evaluations yield refusals,
  which fall back to the free gas (for $T_{eff}$) or no elastic component
  (displacement).

## New or changed in 5.0.0

Release 4.4.6 extended all kernels with the free gas at the material
temperature, whereas NJOY-style evaluations use the SCT recipe. New in
5.0.0: the SCT model as the default extension, effective temperatures for
every VDOS-based material and for direct kernels (from the NCMAT v8 keyword
or the kernel analysis), displacements estimated from kernels and hence
elastic physics for direct-kernel solids, and the NCMAT v8 keywords and
state-of-matter rules that control all of this.

## References

- W. E. Lamb, Phys. Rev. 55 (1939) 190.
- G. M. Borgonovi, "Neutron scattering kernels calculations at epithermal
  energies", Gulf General Atomic report GA-9950 (1970).
- ENDF-6 Formats Manual (ENDF-102), section 7.4.
- R. E. MacFarlane, "New thermal neutron scattering files for ENDF/B-VI
  release 2", LA-12639-MS (1994); NJOY2016 manual, THERMR chapter.
- T. M. Sutton, T. H. Trumbull, C. R. Lubitz, "Comparison of some Monte
  Carlo models for bound hydrogen scattering", Proc. M&C 2009; and their
  CSEWG 2008 presentation "The short-collision-time approximation in ENDF".
- A. Sjolander, Arkiv for Fysik 14 (1958) 315.
- X.-X. Cai, T. Kittelmann, E. Klinkby, J. I. Marquez Damian, J. Comput.
  Phys. 380 (2019) 400.
