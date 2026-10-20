# Multiphase materials, composition, density and temperature

## Summary

A material in NCrystal is either a single phase (one `Info` with a
composition, density, temperature, state of matter and optionally a crystal
structure and dynamics) or a multiphase material: an ordered list of phases
with volume fractions, composed either in the cfg-string
(`phases<0.3*A&0.7*B>;common`) or inside NCMAT data (`@OTHERPHASES`), nested
to any depth. Physics objects of a multiphase material are the fraction- and
number-density-weighted sums of the phases' processes (per-atom cross
sections), and the combined composition drives absorption and macroscopic
conversions. Density, temperature and the state of matter are validated and
completed by the `InfoBuilder`, and can be overridden or scaled from the
cfg-string.

Code: `ncrystal_core/src/factories/NCMatCfg.cc`,
`ncrystal_core/src/infobld/NCInfoBldComposition.cc`,
`ncrystal_core/src/stdscatfactory/NCStdMPScatFact.cc`.

## Concepts

- Composition: a list of (fraction, `AtomData`) by atom count, derived
  from the unit cell when there is one, else from the `@DYNINFO` fractions,
  else given directly (non-NCMAT sources).
- Density: mass density and number density are linked by the average atom
  mass; for crystals both follow from the cell, for others `@DENSITY`
  provides one and the other is computed. The `density` cfg parameter
  overrides (`2.3gcm3`, `1e22perAa3`, `1200kgm3`) or scales (`0.9x`,
  cumulative) it, which changes the macroscopic cross sections but not the
  per-atom physics.
- Temperature: `temp` (default -1 = 293.15 K unless the data lock it; units
  K, C, F; range 0.001 K to 1e6 K); NCMAT v7 `@TEMPERATURE` sets a default
  or locks the value, and `scatknl` sections lock it to the kernel
  temperature. It enters the Debye-Waller factors, the VDOS expansion and
  the free-gas model.
- State of matter: `Unknown`, `Solid`, `Liquid`, `Gas`; solid is implied by
  crystal structure or VDOS-type dynamics, otherwise declared
  (`@STATEOFMATTER`, mandatory in v8 for non-crystals); it decides whether a
  non-crystalline material gets elastic scattering (chapter 24).
- Multiphase: phases with volume fractions $\phi_i$; the per-atom cross
  section of the mixture is $\sum_i w_i\sigma_i$ with
  $w_i = \phi_i n_{d,i}/\sum_j\phi_j n_{d,j}$, i.e. weights by atom count,
  which is why the multiphase factory needs the number densities and not
  only the fractions. Multiphase `Info` objects expose the phase list, the
  combined composition, the combined density and a state of matter only if
  all phases agree.

## Inputs and configuration

- Cfg syntax `phases<FRAC1*CFG1&...&FRACN*CFGN>[;COMMONCFG]`: fractions by
  volume, each phase's final cfg is `CFGi;COMMONCFG`, except that
  `phasechoice` and `density` in the common part apply to the top level
  only. A single phase collapses to a normal single-phase configuration.
  Examples: `phases<0.5*Al_sg225.ncmat&0.5*Mg_sg194.ncmat>;temp=20K`,
  `phases<0.99*Al_sg225.ncmat&0.01*void.ncmat>` (porosity via the
  `void.ncmat` stdlib file: a sterile zero-cross-section material of
  density 1e-200), `phases<0.01*Al2O3_sg167_Corundum.ncmat&0.99*LiquidWaterH2O_T293.6K.ncmat>`.
- NCMAT v6 `@OTHERPHASES`: one `fraction cfgstr` line per secondary phase
  (any loadable cfg-string, even nested `phases<...>`); the file's own
  material becomes phase 0 with the remaining volume fraction.
- `phasechoice=N` selects one phase of a multiphase material for inspection
  (top-level pseudo-parameter).
- `@DENSITY` (v2) in `atoms_per_aa3`, `kg_per_m3`, `g_per_cm3` for
  non-crystals; `@TEMPERATURE` (v7); `@STATEOFMATTER` (v5, v8 rules).
- Gas mixtures: `gasmix::<spec>`: molar (default) or `massfractions` mixtures such as `0.7xCO2+0.3xAr`, the
  `air` alias, pressure (Pa, mbar, hPa, kPa, bar, atm; default 1 atm) or a
  density, temperature in K/C/F, `relhumidity` using the saturated vapour
  pressure of doi:10.1175/JAMC-D-17-0334.1, the ideal gas law for the
  density, and `atomdb` lines; emits a locked `@TEMPERATURE` and free-gas
  dynamics. Example: `gasmix::BF3/2atm/25C/B_is_0.95_B10_0.05_B11`.

## Implementation

When an `Info` object is built, the composition is derived or validated,
temperatures are checked for consistency (including kernel temperatures),
the average atom mass and the densities are completed (unit cell first,
then given values, with compatibility checks), and the state of matter is
completed. Cfg-level overrides such as `density` create variant `Info`
objects sharing the heavy data (hkl lists, dynamics) with the original.
The scatter object of an Info-level multiphase material is built from the
phases' scatter objects (recursively for nested phases) merged with the
atom-count weights; cfg-level multiphase materials are handled by the
configuration object itself, which forwards parameters to all phases.
`Info`, `Scatter` and `Absorption` objects are cached per phase, so a phase
shared by many materials is loaded once.

Python: `NCrystal.createInfo/createScatter/createAbsorption/load(cfg)`,
`Info.isMultiPhase()`, `Info.getPhases()`, `Info.getComposition()`,
`Info.getDensity()`, `Info.getNumberDensity()`, `Info.getTemperature()`,
`Info.stateOfMatter()`, `Info.getSLD()`;
`NCMATComposer.add_secondary_phase`, `set_density`, `set_state_of_matter`,
`lock_temperature`, `set_default_temperature`; `nctool --phases`.

## Validation

- `tests/src/app_cfgmp`, `app_cfg`, `app_cfgparse`, `app_extracfg`: cfg
  parsing incl. `phases<>`, `density`, `temp`; `tests/scripts/cfgstr.py`,
  `ncmatphases.py` (`@OTHERPHASES`), `mmcphases.py` (multiphase in MiniMC),
  `quickfact.py` (`freegas::`, `solid::`, `gasmix::`), and the `nctool`
  dumps of multiphase test materials such as
  `tests/data/Al_sg225_with_voids.ncmat`.
- NCrystal 3.0.0 release announcement (multiphase introduction, Al with
  voids example).

## Approximations and limitations

- Phases are mixed per atom without any geometric structure: no grain or
  particle shapes, no interfaces (apart from SANS, chapter 41), no
  correlations between phases.
- The state of matter is a single label per phase; mixed states (e.g.
  slurries) are just multiphase.
- `density` overrides change macroscopic quantities only; they do not
  rescale the lattice (use `strain` for that).
- Temperature overrides are refused for materials with direct kernels at a
  fixed temperature.

## New or changed in 5.0.0

NCMAT v8 state-of-matter rules (mandatory for non-crystals; crystals
always solid) and the browse tooling listing phases and physics
properties; multiphase itself is unchanged since 3.0.0.

## References

- `ncrystal_core/include/NCrystal/factories/NCMatCfg.hh` (the multiphase
  cfg-string notes), `docs/ncmat_doc.md` (@DENSITY, @OTHERPHASES,
  @TEMPERATURE, @STATEOFMATTER).
- NCrystal 3.0.0 release announcement.
