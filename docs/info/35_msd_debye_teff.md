# Mean-squared displacements, Debye temperatures and effective temperatures

## Summary

Three scalar quantities per atom tie the dynamics of a material to the rest
of the physics: the isotropic mean-squared displacement $\langle u^2\rangle$
(the thermal vibration amplitude, which sets the Debye-Waller factors of
Bragg diffraction and incoherent elastic scattering and the Poisson weights
of the multi-phonon expansion), the Debye temperature $T_D$ (a one-parameter
description of the phonon spectrum, used to describe displacements at all
temperatures and as the fallback dynamics model), and the effective
temperature $T_{eff}$ (the mean kinetic energy per atom including zero-point
motion, used by the short-collision-time extension). NCrystal obtains them
from whichever dynamics information a material carries (VDOS, Debye
temperature, explicit NCMAT v8 values, or analysis of a tabulated kernel),
converts between displacement and Debye temperature through the isotropic
Debye model, and exposes them on the `Info` object.

Code: `ncrystal_core/src/vdos/NCVDOSEval.cc`,
`ncrystal_core/src/phys_utils/NCDebyeMSD.cc`,
`ncrystal_core/src/dyninfoutils/NCDynInfoUtils.cc`,
`ncrystal_core/src/ncmat/NCLoadNCMAT.cc`.

## Physics

### Definitions

For an atom vibrating harmonically and isotropically around its equilibrium
position, $\langle u^2\rangle$ denotes the mean-squared displacement *along
one axis* (so the total is $3\langle u^2\rangle$), which crystallographers
call $U_{iso}$ ($B_{iso} = 8\pi^2 U_{iso}$). The Debye-Waller factor of an
elastic process with momentum transfer $q$ is then $e^{-2W}$ with
$2W = q^2\langle u^2\rangle$; in the structure factors of chapter 12 each
atom's scattering length is multiplied by $e^{-W}$ with $q = 2\pi/d$, and the
same $\langle u^2\rangle$ sets the $\alpha$ scale of the phonon expansion
(chapter 31).

### From a density of states

With the unit-normalised VDOS $\rho(E)$ of chapter 31 (Sjolander 1958):

$$\langle u^2\rangle = \frac{\hbar^2}{2M}\int_0^{E_{max}}\frac{\rho(E)}{E}\coth\frac{E}{2kT}\,dE, \qquad
kT_{eff} = \int_0^{E_{max}}\rho(E)\,\frac{E}{2}\coth\frac{E}{2kT}\,dE.$$

Both include the zero-point contribution ($\coth \to 1$ at $T \to 0$), and
$T \le T_{eff} \le T + E_{max}/2k$.

### The isotropic Debye model

For a Debye spectrum with cut-off $kT_D$ the displacement integral has the
closed form (Glauber 1955, eqs. 11-12)

$$\langle u^2\rangle = \frac{3\hbar^2}{M k T_D}\Big[\frac14 + \Big(\frac{T}{T_D}\Big)^2\int_0^{T_D/T}\frac{x\,dx}{e^x-1}\Big],$$

valid down to $T = 0$ where it gives the zero-point value
$3\hbar^2/(4MkT_D)$, and growing linearly in $T$ at high temperature. The
inverse, $T_D$ from a displacement at a given temperature and mass, is
obtained numerically. The Debye model is also the fallback density of states
for the multi-phonon expansion (chapter 31).

### Debye temperatures as a transport of displacements between temperatures

NCrystal treats a Debye temperature as a temperature-independent material
parameter and a displacement as a quantity valid at one temperature. When a
displacement is given at temperature $T_1$ (NCMAT v8 `msd X at_temperature
T1`) and the material is used at $T_2$, the value is converted to $T_D$ at
$T_1$ and back to a displacement at $T_2$ through the Debye model. The
temperature attached to the displacement therefore restricts nothing; it
only pins the value. A VDOS, in contrast, gives the displacement directly at
any temperature.

## Inputs and configuration

| Dynamics type | Displacement | $T_{eff}$ |
|---|---|---|
| `vdos` | from the VDOS at the material temperature (a @DEBYETEMPERATURE section in v4-v7 files takes precedence, with a warning) | from the VDOS |
| `vdosdebye` | Debye model at the material temperature, from `debye_temp`, from `msd ... at_temperature` (v8), or from the @DEBYETEMPERATURE section (v1-v7) | not needed (legacy processing only) |
| `scatknl` | v8 `msd` (transported to the kernel temperature if given at another) or `debye_temp`; otherwise estimated from the kernel for v8 files and non-NCMAT sources (chapter 33); none for pre-v8 files | v8 `effective_temperature`, else estimated from the kernel |
| `freegas`, `sterile` | none (v8 forbids them in crystals; v1-v7 crystals took the @DEBYETEMPERATURE value) | none |

NCMAT v8 (`docs/ncmat_doc.md`) removes @DEBYETEMPERATURE, adds
`msd ... at_temperature`, allows `msd`/`debye_temp` on `scatknl` sections of
solids, adds `effective_temperature` on `scatknl` sections (at least the
kernel temperature), requires @STATEOFMATTER for non-crystalline materials,
and demands a usable displacement for every atom of a crystal at load time
(explicit keyword or confident analysis), with an error otherwise. The
material temperature `temp` enters all displacement evaluations; direct
kernels fix it to the kernel temperature.

## Implementation

The displacements and Debye temperatures of a crystal are resolved when the
`Info` object is built (concurrently for the elements of a material, since
the VDOS integrals and kernel analyses take milliseconds each), and stored
per atom together with the structure; a refused estimate is a load-time
error for crystals, naming the elements concerned. For non-crystalline
solids and for the SCT extension the values are resolved from the dynamics
when the processes are built. Python exposes them as
`AtomInfo.meanSquaredDisplacement` and `.debyeTemperature`, through
`NCrystal.vdos.analyseVDOS` (displacement, Debye temperature, $T_{eff}$ and
the integrals of a VDOS), `debyeIsotropicMSD` and
`debyeTempFromIsotropicMSD`; the `NCMATComposer` writes `msd` or Debye
temperatures, and `nctool --dump` prints both per atom.

## Validation

- `tests/scripts/vdoseval.py` and `analysevdos.py` (the VDOS integrals
  against references), `multimsd.py` (materials mixing several displacement
  sources), `ncmatv8msd.py` (the v8 keyword rules and the Debye-model round
  trip to machine precision), `sabanalyse.py` (estimates from kernels
  against VDOS truth).
- Historically, the Debye-temperature-based displacements of 1.0.0 were
  compared with measured Bragg edges and total cross sections (CPC 2020
  paper, section 4).

## Approximations and limitations

- Isotropic and harmonic displacements only: no anisotropic displacement
  parameters ($U_{ij}$), one value per element rather than per
  crystallographic site, no anharmonic temperature dependence beyond what a
  VDOS or the Debye model gives.
- The Debye model is a crude phonon spectrum; it is adequate for
  displacements but not for energy transfers (chapter 31).
- The transport of a displacement to another temperature goes through the
  Debye model and inherits its shape assumptions.
- Displacements estimated from kernels are provisional pending comparisons
  with evaluated reference data (chapter 33).

## New or changed in 5.0.0

Release 4.4.6 already derived displacements from VDOS curves (since 2.0.0)
and from Debye temperatures, per element. New in 5.0.0: the effective
temperature as a first-class quantity (from the VDOS, the v8 keyword, or
kernel analysis), displacements for direct kernels (v8 keywords or kernel
analysis, and hence crystals and solids with `scatknl` dynamics), the v8
`msd ... at_temperature` keyword replacing @DEBYETEMPERATURE, and the
`NCMATComposer` writing `msd` directly. Anisotropic displacements and site
occupancies are on the ABI discussion list for 5.0.0.

## References

- R. J. Glauber, "Time-dependent displacement correlations and inelastic
  scattering by crystals", Phys. Rev. 98 (1955) 1692, eqs. 11-12.
- A. Sjolander, Arkiv for Fysik 14 (1958) 315.
- V. F. Sears, Acta Cryst. A53 (1997) 35 (the Debye-Waller exponent).
- ENDF-6 Formats Manual (ENDF-102), section 7.4 (effective temperatures).
- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, sections 2.2-2.4.
