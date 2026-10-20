# Small-angle neutron scattering (SANS)

## Summary

Since 3.0.0 NCrystal has framework support for SANS physics as an additional
elastic component of multiphase materials: a scale factor
$C_{sans} = 2\pi(\phi/n_d)\Delta\rho^2$ built from the scattering-length
density contrast between two phases, their volume fractions and the number
density, and a per-atom cross section
$\sigma_{sans}(k) = C_{sans}(V_p/k^2)\int_0^{2k}Q\,P(Q)S(Q)\,dQ$ for
particles of volume $V_p$ with form factor $P(Q)$ and structure factor
$S(Q)$. The core ships one model, dilute hard spheres ($S = 1$, spherical
$P$), activated by an NCMAT `@CUSTOM_HARDSPHERESANS` section, mainly to
demonstrate and test the framework; a helper supports plugins that
describe SANS through tabulated $I(Q)$ curves (the SANSND nanodiamond
plugin). Native SANS input in NCMAT (without custom sections) is foreseen but
not implemented; `I(Q)` support in NCMAT is on the list for 5.0.0.

Code: `ncrystal_core/src/extd_utils/NCSANSUtils.cc`,
`ncrystal_core/src/sanshardsphere/`, `ncrystal_core/src/phys_utils/NCIofQHelper.cc`.

## Physics

### Cross section per atom (derivation)

The macroscopic SANS cross section of nanoparticles (number density
$N_p/V$, volume $V_p$, SLD contrast $\Delta\rho$) is
$\Sigma(Q) = (N_p/V)V_p^2\Delta\rho^2 P(Q)S(Q)$ (Hammouda, "The SANS
toolbox", eq. 1). Dividing by the atomic number density $n_d$ and using
$N_p/V = \phi/V_p$ with the particle volume fraction $\phi$ gives the
per-atom $d\sigma/d\Omega = (\phi/n_d)V_p\Delta\rho^2 P(Q)S(Q)$. For
elastic scattering $Q = 2k\sin(\theta/2)$, $d\Omega = 2\pi\,Q\,dQ/k^2$, and
integrating over $Q \in [0, 2k]$:

$$\sigma_{sans}(k) = \frac{2\pi\phi\Delta\rho^2}{n_d}\,\frac{V_p}{k^2}\int_0^{2k}Q\,P(Q)S(Q)\,dQ = C_{sans}\,\frac{V_p}{k^2}\int_0^{2k}Q\,P(Q)S(Q)\,dQ.$$

$C_{sans}$ is expressed in barn/Aa$^3$ = cm$^{-1}$. Sampling draws $Q$ with density $\propto Q P(Q) S(Q)$ on $[0, 2k]$
and sets the scattering angle from it; the energy is unchanged.

### The hard-sphere model

$S(Q) = 1$ (dilute) and $P(Q) = [3(\sin x - x\cos x)/x^3]^2$ with $x = QR$
for spheres of radius $R$. With $x = QR$ and $V_p = \frac43\pi R^3$ the
integral is closed-form:

$$\sigma_{sans}(k) = C_{sans}\,\frac{3\pi}{32}R^3\,(Rk)^{-6}\Big[32(Rk)^4 - 8(Rk)^2 + 4Rk\sin(4Rk) + \cos(4Rk) - 1\Big],$$

evaluated by a Taylor expansion at small $Rk$ for numerical stability;
sampling draws $x$ from the density $\propto[\sin x - x\cos x]^2/x^5$ on
$[0, 2kR]$. At long wavelengths ($Rk \to 0$, $P \to 1$) the cross section
tends to the constant $2C_{sans}V_p = 4\pi(\phi/n_d)\Delta\rho^2 V_p$
(isotropic scattering), and at short wavelengths it falls as
$3\pi C_{sans}R/k^2$.

### SLD contrast between phases

The scattering-length density of a phase is $\rho = n_d\sum_i f_i b_{coh,i}$.
For a `@CUSTOM_HARDSPHERESANS` section in phase A, the
contrast is taken with the phase following it in the phase list (phase B),
$\phi$ is the smaller of the two volume fractions, and more than one such
pair may occur (also in nested multiphase materials, searched recursively);
a B phase itself carrying the section is an error.

## Inputs and configuration

- NCMAT `@CUSTOM_HARDSPHERESANS` with one number, the sphere radius in Aa
  ($0 < R < 10^6$); only meaningful in a multiphase context, e.g. a nanoparticle phase in a matrix via
  `@OTHERPHASES` or `phases<...>`.
- `sans` (bool, default true; "primarily added to support future
  developments"); `elas=0` disables it with the other elastic parts;
  `comp=` accepts `sans`.
- Python: `NCMATComposer.add_hard_sphere_sans_model(sphere_radius)`.
- Plugins: the $I(Q)$ helper integrates and samples $Q I(Q)$ on $[0, 2k]$
  from tabulated points (flat $Q I(Q)$ below the first point, so the cross
  section is $1/v$ at low energies, zero above the last point; callers add
  points to change the extrapolation).

## Implementation

The SANS process is an isotropic elastic process with the scale factor and
radius; the factory wraps the standard scatter creation of the multiphase
material and appends one SANS process per detected phase pair. Native NCMAT
SANS input (not via `@CUSTOM_` sections) is foreseen at the same place.

## Validation

- `tests/scripts/mmc_sans.py` (a hard-sphere SANS material in MiniMC,
  chapter 61) and `ncmatphases.py` (multiphase NCMAT materials), plus the
  `nctool` dumps of multiphase test materials with hard-sphere sections.
- The SANSND plugin (out of tree) validated the $I(Q)$ route against
  nanodiamond measurements (Teshigawara et al. fit).

## Approximations and limitations

- Dilute approximation ($S = 1$) and monodisperse spheres in the core
  model; no polydispersity, no inter-particle structure factor, no
  non-spherical form factors.
- Contrast from coherent scattering lengths only (no magnetic SLD).
- SANS and the incoherent approximation for coherent elastic scattering of
  amorphous solids (chapter 24) overlap conceptually; the unofficial hack of
  chapter 24 lets a plugin take over the coherent part.
- NCMAT cannot yet carry a general $I(Q)$ or $S(Q)$ curve; this and the
  SANS handling in general are on the 5.0.0 review list.

## New or changed in 5.0.0

No changes to the SANS physics in 5.0.0; the framework and hard-sphere
model date from 3.0.0.

## References

- B. Hammouda, "The SANS toolbox", NIST (2012), eq. 1.
- NCrystal 3.0.0 release announcement (SANS framework, multiphase).
- SANSND plugin README (nanodiamond $I(Q)$ model).
