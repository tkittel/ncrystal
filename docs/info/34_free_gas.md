# The free-gas model

## Summary

The free-gas model describes neutrons scattering on unbound target atoms in
thermal equilibrium at temperature $T$: the internal structure of the
targets and their interactions are ignored, apart from maintaining the
Maxwellian velocity distribution. Its scattering function has a closed form,
and NCrystal implements it as a fully quantum-mechanical
$S(\alpha,\beta)$ model: the cross section follows from an analytic integral
over the kinematic region, and scattering outcomes are sampled from the exact
joint distribution at any energy with a dedicated rejection scheme, valid to
at least 1 GeV and running at about a million samples per second. It serves
three purposes: the inelastic model for gases and for atoms declared as
`freegas` (e.g. the oxygen of water), the extension model of tabulated
kernels above their range in the legacy processing and as the
$T_{eff} = T$ special case of the SCT model (chapter 33), and the building
block of the SCT sampler.

Code: `ncrystal_core/src/phys_utils/NCFreeGasUtils.cc`,
`ncrystal_core/src/freegas/NCFreeGas.cc`.

## Physics

### Scattering function

In NCrystal's variables (chapter 30), with $A = M/m_n$ the target-to-neutron
mass ratio,

$$S(\alpha,\beta) = \frac{1}{\sqrt{\alpha}}\exp\Big(-\frac{\beta^2}{4\alpha/A} - \frac{\alpha}{4A}\Big)$$

up to a constant (sampling paper eq. 19; in ENDF variables this is the
familiar $e^{-\beta/2}\exp(-(\alpha_E^2+\beta^2)/4\alpha_E)/\sqrt{4\pi\alpha_E}$;
see Squires chapter 4.5 for the $S(q,\omega)$ derivation). It satisfies
detailed balance at $T$, and its recoil ridge lies at $\beta = -\alpha/A$,
i.e. at $\beta = -\alpha_E$.

### Cross section

Integrating $S$ over the kinematic region $D(E)$ gives (sampling paper eq.
20, with the mass dependence made explicit)

$$\sigma(E) = \sigma_{free}\Big[\Big(1+\frac{1}{2a^2}\Big)\operatorname{erf}(a) + \frac{e^{-a^2}}{\sqrt{\pi}\,a}\Big], \qquad a^2 = \frac{A\,E}{kT}, \qquad \sigma_{free} = \sigma_b\Big(\frac{A}{1+A}\Big)^2.$$

The bracket is the usual free-gas Doppler factor, tending to 1 at high
energy and to the $1/v$ law at low energy; the prefactor relates the bound
cross section of the nuclide (chapter 43) to the free one.

### Energy-transfer distribution (derivation)

For a given $E$ the marginal distribution of $\beta$ follows by integrating
$S$ over $\alpha \in [\alpha_-(\beta),\alpha_+(\beta)]$ and involves
differences of error functions. For $A = 1$ it reduces to the textbook
result $f(E \to E') = e^{\beta}\operatorname{erf}(\sqrt{E/kT})/(E\sigma)$
for $E < E'$ and $\operatorname{erf}(\sqrt{E'/kT})/(E\sigma)$ for $E > E'$
(sampling paper eq. 21). Given $\beta$, the conditional distribution of
$\alpha$ is $S(\alpha,\beta)$ itself on $[\alpha_-,\alpha_+]$, which after
absorbing constants has the shape $f(x) = e^{-1/x - cx}/\sqrt{x}$: a single
peak falling off on both sides, with an essential singularity at $x \to 0$,
and no usable closed-form inverse.

## Inputs and configuration

- NCMAT @DYNINFO `type freegas` (no further fields; the mass and bound
  cross section come from the element, chapter 43); NCMAT v8 forbids it in
  crystalline materials since such atoms carry no displacement for the
  structure factors.
- `inelas=freegas` forces the model for every element of any material
  (`inelas=freegas;elas=0` gives a pure free gas).
- The quick-material syntax `freegas::<formula>/<density>` (e.g.
  `freegas::Ar/0.06perAa3`, with the temperature from `temp`), and the gas
  mixture syntax `gasmix::` (e.g. `gasmix::air/0.9relhumidity`,
  `gasmix::0.72xCO2+0.28xAr/massfractions/1.5atm/250K`) for mixtures given
  by formula or name, pressure, temperature and relative humidity (chapter
  42). The data library also contains free-gas files for the noble gases at
  STP (chapter 51).
- Material temperature `temp`; the model is valid at any temperature and
  energy of practical interest.

## Implementation

The cross section uses asymptotic expansions of the Doppler factor at small
and large $a$ and the full formula in between. Sampling proceeds in two
steps, both by rejection with overlays adapted on the fly, since every call
comes with a different energy and nothing can be precomputed: $\beta$ is
drawn from its marginal distribution with an overlay that is flat on the
downscatter side and exponential on the upscatter side, with the sampling
interval narrowed as rejected points reveal where the density is negligible,
and with cheap bounds on the density (from a lookup table of the error
function) deciding most candidates without an exact evaluation; $\alpha$ is
then drawn from $f(x)$ on $[\alpha_-,\alpha_+]$ with a flat overlay, an
$e^{-cx}/\sqrt{x}$ overlay for the tail, or a combination of the two,
depending on where the interval lies relative to the peak. At extremely high
energies the energy loss is drawn uniformly up to the kinematic maximum,
which is the correct limit. Near the kinematic end point the angle is
isotropic (chapter 32). The `FreeGas` process wraps provider and sampler
with the element's mass, temperature and $\sigma_{free}$; the SCT sampler of
chapter 33 reuses the same sampler at $T_{eff}$.

## Validation

- `tests/src/app_sctext` exercises the free-gas cross section and sampler
  directly (as the $T_{eff} = T$ limit of the SCT model): cross sections
  against the closed form and quadrature of $S$, sampled distributions
  against the analytic densities over a range of $A$ and $E/kT$; the
  `sabxs` tests cover the free gas as a kernel extender.
- The sampling paper used the free-gas kernel as the analytic benchmark of
  the kernel sampling (section 3): kernel-sampled results reproduce the
  analytic $\sigma(E)$ to 0.05% below 0.1 eV, and a free-gas water sphere
  reproduced MCNP6's classical free-gas treatment within a fraction of a
  percent between 1e-3 and 1 eV.

## Approximations and limitations

- Ignores molecular binding, rotations and vibrations: appropriate for
  monatomic gases and as a high-energy limit, a crude approximation for
  molecular gases at thermal energies (NCrystal has no dedicated molecular
  gas model; such kernels can be supplied as `scatknl`).
- Single target mass per element (the average isotope mass); no
  resonances.
- Classical transport codes sample the target velocity instead (Gelbard's
  scheme as implemented in OpenMC); whether to revisit NCrystal's
  $(\alpha,\beta)$-space sampler against those is an open item.

## New or changed in 5.0.0

The free-gas model has been in NCrystal since 2.0.0 and is unchanged in
physics. In 5.0.0 it became the special case of the SCT extension (chapter
33) and its sampler the engine of the SCT sampler, and it gained
bit-reproducibility across platforms with the rest of the physics utilities
(chapter 70).

## References

- X.-X. Cai, T. Kittelmann, E. Klinkby, J. I. Marquez Damian, J. Comput.
  Phys. 380 (2019) 400, eqs. 19-21 and section 3.
- G. L. Squires, "Introduction to the theory of thermal neutron
  scattering", chapter 4.5.
- NCrystal 2.0.0 release announcement (noble-gas data files).
