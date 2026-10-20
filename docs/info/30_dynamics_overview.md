# Atomic dynamics and inelastic scattering: overview

## Summary

NCrystal models the inelastic scattering of thermal neutrons on the atoms of
a material through scattering kernels, i.e. tabulated scattering functions
$S(\alpha,\beta)$, one per element (more precisely one per @DYNINFO section
of the material), within the incoherent approximation and for isotropic
materials. The kernels are either given directly in the input data, or, for
solids, generated at load time from a phonon density of states (VDOS) with a
multi-phonon expansion, or, as a fallback, generated from an idealised Debye
density of states. Each kernel is turned into a scattering process that gives
cross sections and samples scattering events at any neutron energy: inside
the tabulated range by integration and sampling of the table, and above it by
an analytical extension model. The same dynamics data also provide the
atomic mean-squared displacements that the elastic models need.

This chapter gives the whole picture; the following chapters describe each
step in detail:

- `31_vdos_to_sab.md`: VDOS to $S(\alpha,\beta)$ (the Sjolander expansion).
- `32_sab_cross_sections_and_sampling.md`: cross sections and sampling from
  a tabulated kernel.
- `33_kernel_extension_sct.md`: what happens above the tabulated range (free
  gas and short-collision-time models), and the estimation of effective
  temperatures and displacements from a kernel.
- `34_free_gas.md`: the free-gas model.
- `35_msd_debye_teff.md`: mean-squared displacements, Debye temperatures and
  effective temperatures.

## Physics

### The scattering function and the incoherent approximation

For an isotropic material the double-differential cross section of a neutron
of kinetic energy $E$ scattering into energy $E'$ and direction cosine $\mu$
is written in terms of the scattering function as

$$\frac{d^2\sigma}{dE'\,d\Omega} = \frac{\sigma_b}{4\pi\,kT}\sqrt{\frac{E'}{E}}\,S(\alpha,\beta)$$

with the dimensionless momentum and energy transfers

$$\alpha = \frac{E + E' - 2\mu\sqrt{EE'}}{kT}, \qquad \beta = \frac{E'-E}{kT}.$$

$\sigma_b$ is the bound scattering cross section of the atom. Note the
NCrystal convention for $\alpha$: it is $\hbar^2 q^2/(2 m_n kT)$ with the
neutron mass, i.e. it does *not* divide by the target-to-neutron mass ratio
$A = M/m_n$ as the ENDF convention does ($\alpha_{ENDF} = \alpha/A$). This
choice was made in the sampling paper (Cai et al. 2019) because a target mass
is ambiguous for kernels describing several isotopes or a whole molecule.
Negative $\beta$ means the neutron loses energy ("downscatter").

The function $S$ used throughout NCrystal is the *asymmetric* one, satisfying
detailed balance as $S(\alpha,\beta) = e^{-\beta} S(\alpha,-\beta)$. The
symmetric form $S_{sym} = e^{\beta/2} S$ is accepted on input and converted.

Integrating over the kinematically accessible region gives the cross section

$$\sigma(E) = \frac{\sigma_b\,kT}{4E}\iint_{D(E)} S(\alpha,\beta)\,d\alpha\,d\beta,$$

where $D(E)$ is the region $\beta \ge -E/kT$,
$\alpha_-(\beta) \le \alpha \le \alpha_+(\beta)$ with
$\alpha_\pm = (2E + kT\beta \pm 2\sqrt{E(E+kT\beta)})/kT$. $D(E)$ grows
monotonically with $E$, which is what the sampling algorithms exploit.

Coherent inelastic scattering is treated in the *incoherent approximation*:
the scattering function of the coherent part is approximated by that of the
incoherent part, so that one kernel per atom species, with the full bound
cross section $\sigma_b = \sigma_{coh} + \sigma_{inc}$, describes all
inelastic scattering (Placzek and Van Hove 1955; Sjolander 1958; CPC 2020
paper section 2.4). Coherent one-phonon effects, which this approximation
misses, can be supplied by plugins (chapter 53); the kernel machinery
supports removing the coherent part of individual phonon orders for that
purpose.

The elastic part of the scattering function (the $\beta = 0$ line of Bragg
diffraction and incoherent elastic scattering) is *not* part of the
kernels: it is handled by the dedicated elastic processes (chapters 20-24),
which take their Debye-Waller factors from the same dynamics information.

### Material temperature

A kernel is valid at one temperature. Kernels given directly in NCMAT data
lock the material temperature to the one declared with them. Kernels
generated from a VDOS are generated for the temperature actually requested,
so VDOS-based materials work at any temperature.

## Inputs and configuration

### The kinds of dynamics information

The input data describe the dynamics of each element with one of five
models. In NCMAT files these are the `type` values of @DYNINFO sections (see
`docs/ncmat_doc.md`); on the loaded `Info` object they become the
corresponding `DynamicInfo` classes, which any data source, not only NCMAT,
can provide:

| NCMAT `type` | Meaning |
|---|---|
| `scatknl` | A tabulated kernel on $(\alpha,\beta)$ grids (asymmetric, or symmetric and scaled), at a fixed temperature, optionally with its own energy grid and, from NCMAT v8, explicit displacement and effective-temperature data |
| `vdos` | A phonon density of states |
| `vdosdebye` | An idealised Debye density of states from a Debye temperature |
| `freegas` | The free-gas model |
| `sterile` | No scattering at all on this element |

The `fraction` of each section weights the element's process in the material
and must be consistent with the unit cell content for crystals. A material
may mix the types freely (e.g. a kernel for the hydrogen and free gas for
the oxygen of water). Crystals specified with Debye temperatures only get
`vdosdebye` dynamics for every element.

### Choosing the inelastic model: the `inelas` parameter

The cfg-string parameter `inelas` (default `auto`) selects how the inelastic
component of a material is built: `auto` uses the dynamics information of
the data when present, else the Debye model for solids with Debye
temperatures, else nothing for crystals and the free gas otherwise;
`vdosdebye` and `freegas` force that model for every element (so
`inelas=freegas;elas=0` turns any material into a pure free gas);
`external` selects an externally provided cross-section curve (only `.nxs`
files through the nxslib plugin); `none`/`0` disables inelastic scattering.
The pseudo-parameter `comp=` selects which components (`elas`, `coh_elas`,
`incoh_elas`, `inelas`, `sans`) to keep, and `elas=0` removes all elastic
components.

### Luxury levels

Two integer parameters trade precision for resources: `vdoslux` controls
the VDOS expansion and the size of the generated kernel (chapter 31),
`knllux` the processing of any kernel into cross sections and samplers
(chapter 32); both default to 3 on a scale of 0 to 6. Both exist in two
generations: the algorithms of releases 2.0.0 to 4.4.6 ("legacy") and the
5.0.0 algorithms, which are the 5.0.0 defaults and, together with the
space-group symmetry work of chapters 10-12, one of the two headline
features of the release (the SCT extension of chapter 33 and NCMAT v8 being
the other major additions). At the time of writing the development branch still
reaches the 5.0.0 algorithms through transitional values (`vdoslux`
2000..2006, `knllux` 0..6 with the legacy ones at negative values); the
chapters describe the 5.0.0 algorithms as the model and the legacy ones as
the compatibility mode.

### Ultra-cold neutrons

The `ucnmode` parameter (`refine`, `only`, `remove`, with an optional
threshold, default 300 neV) replaces the kernel treatment of events that
leave the neutron below the threshold with a dedicated UCN-production
process, which integrates the kernel exactly over the UCN part of phase
space and supports biased simulations of UCN production (release 3.1.0). It
is only available with the legacy kernel processing and is going to be
removed in a release after 5.0.0.

## Implementation

### The pipeline

```
NCMAT @DYNINFO (or any data source)
   |
   v
Info: dynamics per element (VDOS, Debye VDOS, direct kernel, free gas, sterile)
   |
   |     VDOS / Debye VDOS -> phonon expansion (chapter 31) -> kernel
   |     direct kernel     -> conversion to the standard asymmetric form
   v
kernel (asymmetric S on an alpha,beta grid; T, sigma_b, M, Emax)
   |
   |     processing (chapter 32) + extension model (chapter 33)
   v
scattering process of the element, scaled by sigma_b and the element
fraction, combined with the elastic processes into the material's process.
```

The free-gas type bypasses the kernel machinery and gives a free-gas process
directly (chapter 34). Kernels and their processed forms are cached on the
data rather than the material, so identical dynamics in different materials
or phases are expanded and processed once.

### Quantities derived from the dynamics

Two scalar quantities are extracted from the dynamics of each element:

| Quantity | VDOS | Debye VDOS | direct kernel |
|---|---|---|---|
| Mean-squared displacement | from the VDOS | Debye model at the material temperature | the explicit NCMAT v8 `msd`/`debye_temp` value, otherwise estimated from the kernel (NCMAT v8 and non-NCMAT sources only) |
| Effective temperature | from the VDOS | not needed | the explicit NCMAT v8 `effective_temperature`, otherwise estimated from the kernel |

The displacement feeds the Debye-Waller factors of the elastic processes
(Bragg structure factors and incoherent elastic scattering), the effective
temperature feeds the short-collision-time extension of the kernel. For
crystals the displacement must be known for every atom at load time, since
the structure factors need it. Chapter 35 gives the formulas, chapter 33
the kernel analysis.

### Which elements get elastic physics from their dynamics

- Crystalline materials: the structure factors and the incoherent elastic
  cross sections use the displacement of every atom, whatever its dynamics
  type (NCMAT v8 bans `freegas` and `sterile` from crystals for this reason).
- Non-crystalline solids (`@STATEOFMATTER solid`): each element with a
  displacement gets incoherent elastic scattering in the incoherent
  approximation of the coherent part (chapter 24); elements without one are
  left out, with a warning.
- Liquids and gases: no elastic scattering from the dynamics.

## Validation

- The VDOS expansion against high-precision references, the kernel
  processing against a slow reference implementation, and the cross
  sections of many materials over all luxury levels against stored
  references (chapters 31-32).
- The free gas, the SCT model and the kernel analyser have dedicated tests
  (chapters 33-34).
- Historically, kernels generated by NCrystal from VDOS curves were compared
  with the ENDF/B-VIII.0 kernels for BeO at eight temperatures (release 2.0.0
  notes), and the kernel sampling method was benchmarked against the
  analytic free gas and against MCNP6 and OpenMC (Cai et al. 2019).

## Approximations and limitations

- The incoherent approximation: no coherent inelastic (one-phonon,
  dispersion) effects in the core; isotropic scattering functions only.
- The kernels of VDOS origin assume harmonic lattice dynamics (the Gaussian
  approximation), isotropic displacements, and a VDOS that behaves as $E^2$
  below its first point (chapter 31).
- One kernel per element: molecules or unit cells are not described as a
  whole, and a tabulated kernel is only valid at its own temperature.
- $S(q,\omega)$ input is foreseen in the NCMAT specification but not
  implemented.
- Nuclear resonances and magnetic scattering are not modelled (magnetic
  scattering exists as a plugin for one specific system).

## New or changed in 5.0.0

Release 4.4.6 already modelled inelastic scattering through kernels expanded
from VDOS curves at any temperature, with exact rejection-based sampling,
free-gas extension and dedicated UCN production, which few transport codes
offer. The 5.0.0 development adds:

- A new generation of kernel processing (exact cell-wise phase-space
  integrals, automatic energy grids, cell-based sampling; chapter 32) and of
  VDOS expansion (exact detailed balance at low temperatures, robust
  handling of narrow spectral structure and of numerical noise, resource
  rather than order-count limits; chapter 31), both now the defaults.
- The short-collision-time extension with an effective temperature above
  the kernel range for all kernels, replacing the free gas at the material
  temperature (chapter 33).
- Effective temperatures and displacements estimated directly from tabulated
  kernels, giving elastic physics to solids described by direct kernels
  (chapter 33).
- NCMAT v8 keywords for explicit displacement and effective-temperature
  data, kernels in crystals, and the state-of-matter rules that decide
  elastic physics (chapter 50).
- Bit-reproducible results across platforms and faster hot paths (chapter
  70).

## References

- X.-X. Cai, T. Kittelmann, "NCrystal: A library for thermal neutron
  transport", Comput. Phys. Commun. 246 (2020) 106851,
  doi:10.1016/j.cpc.2019.07.015 (section 2.4 for the inelastic framework).
- X.-X. Cai, T. Kittelmann, E. Klinkby, J. I. Marquez Damian,
  "Rejection-based sampling of inelastic neutron scattering", J. Comput.
  Phys. 380 (2019) 400, doi:10.1016/j.jcp.2018.11.043 (the $(\alpha,\beta)$
  convention, kinematic limits, the sampling algorithm).
- A. Sjolander, "Multi-phonon processes in slow neutron scattering by
  crystals", Arkiv for Fysik 14 (1958) 315.
- G. Placzek, L. Van Hove, "Interference effects in the total neutron
  scattering cross-section of crystals", Nuovo Cimento 1 (1955) 233.
- NCrystal release announcements 2.0.0 and 3.1.0,
  https://github.com/mctools/ncrystal/wiki.
