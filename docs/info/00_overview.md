# NCrystal: scope, process model and material classes

## Summary

NCrystal is a C++11 library (with C and Python interfaces and command-line
tools, no third-party dependencies, Apache 2.0) that provides the
thermal-neutron physics of materials to Monte Carlo transport and
ray-tracing codes: for a material described by a short configuration
string, it returns cross sections and samples scattering outcomes for
coherent and incoherent, elastic and inelastic scattering, plus $1/v$
absorption, in crystal powders, mosaic and layered single crystals,
amorphous solids, liquids and gases, including multiphase and isotopically
modified materials and (as a framework) SANS. It is integrated in OpenMC,
Geant4 and McStas, usable directly from C++, C and Python, and ships a
validated data library. Release 5.0.0 keeps this scope and ABI-breaks to
make room for the improvements described in chapter 90.

## The process model

Everything a transport code needs is expressed through two kinds of
immutable, thread-safe *processes* created from a configuration string
(`createScatter(cfg)`, `createAbsorption(cfg)`), each with

- `crossSection(E, direction)` (per atom, in barn; `crossSectionIsotropic(E)`
  for isotropic materials),
- `sampleScatter(rng, E, direction)` giving the outgoing energy and
  direction (or `(E', mu)` for isotropic materials),
- a `domain()` (energy range with possibly non-zero cross section) and
  flags (isotropic material or not, scattering or absorption, pure elastic
  or not),
- `createMerged` for combining compatible processes with scale factors
  (phases, components),
- a JSON description of the model and its parameters.

A material's scatter process is a composition of components: coherent
elastic (`coh_elas`: Bragg diffraction for crystals, incoherent
approximation for amorphous solids), incoherent elastic (`incoh_elas`),
inelastic (`inelas`: kernels, free gas), SANS (`sans`), and plugin-added
processes; each can be switched off (`elas=0`, `comp=...`) to study it in
isolation, and the same cfg-string gives the same physics in every
interface. "Isotropic material" is a property of the material (cross
sections independent of the neutron direction), not of the scattering: a
powder is isotropic although Bragg scattering is not; every single crystal
is anisotropic. Caches of neutron-state-dependent data live in per-thread
`CachePtr` objects, and random numbers come from an `RNG` interface that the
application may replace (`ncrystal_samplescatter_rs` in the C API lets the
caller pass its own stream per particle). The macroscopic conversion uses
the material's number density (`Info`). Conversions to `Info` objects
(`createInfo(cfg)`) expose the structural and dynamic data behind the
processes (chapters 10, 30, 42, 43).

## Material classes and what is modelled

| Class | Defined by | Elastic | Inelastic | Chapters |
|---|---|---|---|---|
| Crystal powder / polycrystal | unit cell, space group, positions or sites, displacements | powder Bragg + incoherent elastic | kernels from VDOS (or Debye fallback), direct kernels | 10-12, 20, 23, 30-35 |
| Mosaic single crystal | same + `dir1`, `dir2`, `mos` | Gaussian mosaic Bragg (+ powder below `sccutoff`) + incoherent elastic | same | 21 |
| Layered single crystal (PG) | same + `lcaxis` | rotation-averaged mosaic Bragg | same | 22 |
| Amorphous solid | composition, density, `solid`, VDOS/Debye/kernel dynamics | incoherent-approximation elastic | kernels | 24, 30-35 |
| Liquid | composition, density, kernels or free gas (Liquids plugin for diffusion) | none (plugin: diffusion, $S(Q)$) | kernels, SCT-extended | 30-35, 53 |
| Gas | composition, density (or `gasmix::`) | none | free gas | 34, 42 |
| Multiphase | phases with volume fractions | sum of the phases + SANS between phases | sum | 41, 42 |

Absorption is $1/v$ for everything (chapter 40); atomic data and mixtures
in chapter 43; the input formats in chapters 50-53; tools in 60-61;
numerics in 70.

## Interfaces

- C++: `NCrystal.hh` (managed processes with RNG), `NCProcImpl.hh` and
  `NCPluginMgmt.hh` (unmanaged, for frameworks and plugins), `NCInfo.hh`,
  the `MatCfg` configuration class.
- C: `ncrystal.h`, the stable ABI used by the Python module and the
  Geant4/OpenMC/McStas integrations.
- Python: `NCrystal` (`createScatter`, `createInfo`, `load`, `NCMATComposer`,
  `vdos`, `plot`, `minimc`, `browse`, `cli`), command-line tools `nctool`,
  `ncrystal <tool>` / `ncrystal_<tool>`, `ncrystal-config`.
- Virtual API (`NCVirtAPI_Type1_v1.hh`) for ABI-stable embedding (OpenMC).
- Integrations: OpenMC (`openmc.Material.from_ncrystal`), Geant4
  (`ncrystal-geant4`), McStas (`NCrystal_sample`, Union process,
  `.laz`/`.lau` and filter-table exports), ENDF export for MCNP/PHITS-style
  pipelines, and the NJOY+NCrystal predecessor.

## Configuration strings

`"<data>;par1=val1;par2=val2"`: the data name (a stdlib file, a path,
`source::name`, a quick material, or `phases<...>`), then parameters:
`temp`, `density`, `dcutoff`, `dcutoffup`, `mos`, `dir1`, `dir2`, `dirtol`,
`mosprec`, `lcaxis`, `lcmode`, `sccutoff`, `strain`, `coh_elas`,
`incoh_elas`, `inelas`, `sans`, `vdoslux`, `knllux`, `atomdb`, `ucnmode`,
`infofactory`, `scatfactory`, `absnfactory`, and the pseudo-parameters
`bragg`, `elas`, `comp`, `phasechoice`, `bkgd` (obsolete). The reference
is `nctool --doc` / `NCrystal.generateCfgStrDoc()`; files can embed
defaults with `NCRYSTALMATCFG[...]`.

## References

- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851 (the framework paper).
- `README.md`, the NCrystal wiki (Home, Using-NCrystal, Included-physics).
