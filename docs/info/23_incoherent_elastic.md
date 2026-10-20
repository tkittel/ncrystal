# Incoherent elastic scattering

## Summary

Incoherent elastic scattering, the elastic diffuse scattering from the
spread of scattering lengths among isotopes and spin states (and, in the
incoherent approximation, from disorder), is modelled for every solid with
the isotropic harmonic Debye-Waller factor:
$d\sigma/d\Omega = \sum_j \sigma_{inc,j}/4\pi\,e^{-Q^2\langle u_j^2\rangle}$,
whose total cross section and angular distribution have closed forms. It
uses the same displacements as Bragg diffraction (chapter 35), so the
elastic physics of a material is self-consistent (since 2.0.0). For
amorphous solids the coherent part is added in the same form (chapter 24).
Cross sections and samples are essentially free compared with the other
processes.

Code: `ncrystal_core/src/elincscatter/NCElIncScatter.cc`,
`ncrystal_core/src/phys_utils/NCElIncXS.cc` (header
`ncrystal_core/include/NCrystal/internal/phys_utils/NCElIncXS.hh` carries
the derivation).

## Physics

### Cross section (derivation)

From the harmonic approximation (CPC 2020 eq. 23), the incoherent elastic
double-differential cross section of element $j$ is
$\sigma_{inc,j}/4\pi\,e^{-2W_j(Q)}\delta(\hbar\omega)$ with
$2W_j = Q^2\langle u_j^2\rangle$ for isotropic displacements (eq. 40).
With $Q^2 = 2k^2(1-\mu)$ for elastic scattering and integrating over the
azimuth,

$$\sigma_{inc,el}(E) = \frac{\sigma_{inc}}{2}\int_{-1}^{1}e^{-2k^2\langle u^2\rangle(1-\mu)}d\mu = \sigma_{inc}\,\frac{1-e^{-t}}{t}, \qquad t = 4k^2\langle u^2\rangle = \Big(\frac{4\pi\sqrt{\langle u^2\rangle}}{\lambda}\Big)^2$$

(elastic paper eqs. 56-57). At long wavelengths ($t \to 0$) the scattering
is isotropic with the full $\sigma_{inc}$; at short wavelengths it is
forward-peaked and the cross section falls as $\sigma_{inc}/t \propto
\lambda^2$. The scattering angle follows the density
$P(\mu) \propto e^{a\mu}$ with $a = t/2 = 2k^2\langle u^2\rangle$,
normalised by $N_t = t/(4\sinh(t/2))$.

For a polyatomic material the per-atom cross section is the fraction-
weighted sum over elements, each with its own displacement.

### What counts as incoherent

$\sigma_{inc}$ is the bound incoherent cross section of the element
(isotope mix and spin incoherence, chapter 43). In the incoherent
approximation the *coherent* bound cross section can be added to the same
formula, which is what NCrystal does for
non-crystalline solids where no Bragg diffraction describes the coherent
elastic part (chapter 24). The elastic paper shows (figure 35) that the
isotropic-displacement approximation is noticeably off at 1.8 Aa for many
materials, so anisotropic displacements remain a known improvement.

## Inputs and configuration

- `incoh_elas` (default on) enables the component for solids; `elas`,
  `comp` act on it too.
- The displacement of each element (chapter 35), from the structure for
  crystals or from the dynamics for non-crystalline solids; elements
  without a displacement (`freegas`, `sterile`, kernels that could not be
  analysed) are skipped with a warning, and if no element has one the
  component is left out.
- `temp` through the displacements; $\sigma_{inc}$ from the atom data
  (modifiable with `atomdb`, chapter 43).

## Implementation

The per-element terms are evaluated with a numerically safe form of
$(1-e^{-t})/t$ (series at small $t$, $1/t$ at large $t$), $\mu$ is sampled
for an element picked by its contribution, by rejection at small $t$ and
by inversion, $\mu = -1 + \frac{2}{t}\ln(1 + (e^t-1)R)$, otherwise. Both
cross section and sampling cost a handful of floating-point operations per
element.

## Validation

- Elastic paper section 6: total cross sections of vanadium and nickel
  against EXFOR data, and the discussion of the isotropic approximation.
- Repository tests: the `nctool` cross-section dumps of every stdlib
  material (`incoh_elas` contribution listed separately), `dkelas.py`
  (direct-kernel solids gaining the component), `multimsd.py`, and the
  amorphous-solid materials of chapter 24.

## Approximations and limitations

- Isotropic harmonic Debye-Waller factor; one displacement per element.
- Purely elastic: the quasi-elastic broadening of diffusive systems is
  not represented (liquids get no elastic component at all; the Liquids
  plugin adds a diffusion model, chapter 53).
- Magnetic diffuse scattering, nuclear-spin correlations and
  paramagnetic scattering are outside the model.

## New or changed in 5.0.0

The model is unchanged since 2.0.0 (when it was made consistent with the
Bragg Debye-Waller factors). In 5.0.0, solids described by direct kernels
get the component too (through the NCMAT v8 keywords or kernel analysis,
chapter 33), and elements without displacements are skipped with a warning
instead of an error.

## References

- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, section 2.3 (eqs. 23,
  40).
- T. Kittelmann, X.-X. Cai, CPC 267 (2021) 108082, section 6 (eqs. 56-57,
  sampling, EXFOR comparison, figure 35).
