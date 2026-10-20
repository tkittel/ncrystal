# Absorption

## Summary

NCrystal is a scattering library, but applications such as McStas or the
virtual-API users need the full interaction, so a simple absorption process
is provided: the $1/v$ law, $\sigma_{abs}(E) = \sigma_{2200}\sqrt{E_{2200}/E}$,
with the thermal (2200 m/s) capture cross section of each element taken from
the atom database and averaged over the composition. No resonances, no
tabulated energy dependence. Absorption factories are pluggable, so a plugin
can replace the model.

Code: `ncrystal_core/src/absoov/NCAbsOOV.cc`,
`ncrystal_core/src/absfact/NCStdAbsFact.cc`.

## Physics

For a neutron of speed $v$ and a $1/v$ absorber,
$\sigma_{abs}(v) = \sigma_{2200}\,(2200\ \mathrm{m/s})/v$, i.e.
$\sigma_{abs}(E) = c/\sqrt{E}$ with $c = \sigma_{2200}\sqrt{E_{2200}}$,
$E_{2200} = 0.0253$ eV. This is exact for pure $s$-wave capture far from
resonances and is the standard thermal treatment (Geant4's and McStas'
thermal absorption are identical to it, CPC 2020 paper section 5); it is
wrong near resonances (e.g. the 1.3 eV resonance of Sn visible in the CPC
2020 EXFOR comparison, or Cd, Gd and In above their resonance energies).

For a material, $\sigma_{2200}$ is the composition-weighted average
$\sum_i f_i\sigma_{2200,i}$ over the elements and isotopes (chapter 43), with
the fractions by atom count; multiphase materials carry the combined
composition on their `Info` object, so the same formula applies.

## Inputs and configuration

- The capture cross sections of the atom database, overridable per file (`@ATOMDB`) or per cfg-string (`atomdb=`, chapter
  43), e.g. `atomdb=B:is:0.95:B10:0.05:B11` for enriched boron.
- `absnfactory` selects another registered absorption factory (plugins).
- `comp`/`elas`/`inelas` do not affect absorption; it is a separate
  `Absorption` process, created by `createAbsorption(cfg)` and combined with
  scattering by the application (MiniMC applies it as weight attenuation
  with its `absorption=0|1` engine option).

## Implementation

The process stores $c$ and returns $c/\sqrt{E}$, merges additively with
other instances (phase mixing), and is tabulated exactly by the filter
tables (chapter 61). Multiphase materials need no special handling since
the combined composition is on the `Info` object.

## Validation

- CPC 2020 paper section 5: Geant4 plugin absorption identical to Geant4's
  own; EXFOR total cross sections (Mo, Ni, Zn, UO$_2$, Sn) with absorption
  included.
- Repository: the `nctool` cross-section dumps and `.log` tests for all
  stdlib materials include the absorption column; `tests/scripts/atomdbpy.py`
  and the `@ATOMDB` cases of `ncmatparser.py` cover the composition
  averaging.

## Approximations and limitations

- $1/v$ only: no resonances, no $(n,\alpha)$ or $(n,p)$ versus $(n,\gamma)$
  distinction, no energy-dependent tabulation; applications needing
  resonance absorption must combine NCrystal scattering with their own
  absorption data (as OpenMC and Geant4 do, using NCrystal only for
  scattering below a cross-over energy).

## New or changed in 5.0.0

Unchanged since 1.0.0 (section 2.6 of the CPC 2020 paper). New in 5.0.0:
the capture cross sections of the atom database are browsable
(`ncrystal browse --atomdb -w "absxs > 1000"`), and absorption appears in
the filter tables.

## References

- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, section 2.6 and section
  5 (Geant4 comparison, EXFOR).
- V. F. Sears, "Neutron scattering lengths and cross sections", Neutron
  News 3 (1992) 26 (the source of the thermal capture cross sections via
  the ATI Vienna table).
