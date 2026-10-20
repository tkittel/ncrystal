# Tools for analysis, simulation and export

## Summary

Around the physics models NCrystal provides tools to inspect materials
(`nctool`, the JSON `query` interface), to simulate neutrons in simple
geometries with its own fast Monte Carlo (`MiniMC`), to export models to
other codes (`ncmat2endf` for ENDF-6 thermal scattering files,
`ncmat2hkl` for McStas reflection lists, `mcstasunion` for McStas-Union
instrument code, `ncmat2cpp` for embedding data in C++), and, new on the
5.0.0 branch, to tabulate total cross sections as piecewise-linear
"filter tables" for attenuation-only components. The `ncrystal browse` tool
(chapter 51) finds data; it embodies no physics itself.

Code: `ncrystal_python/src/NCrystal/` (`_cli_*.py`, `minimc.py`,
`ncmat2endf.py`, `mcstasutils.py`, `filter.py`), `ncrystal_core/src/minimc/`,
`ncrystal_core/src/filter/`, `ncrystal_core/src/query/`.

## nctool

The main investigation tool: `--dump` prints the `Info` object (structure,
composition, dynamics, hkl list, densities), the default mode plots cross sections per component
(`--coh_elas`, `--incoh_elas`, `--inelas`, `--sans`, `--absorption`
toggles) and scattering distributions, `--phases` lists phases, `--doc`
prints the cfg-string reference, `--bench` times processes, `--mc` runs a
MiniMC scenario, `--extract` shows the text data of a material (e.g. the
generated NCMAT of `solid::`), `--plugins` lists plugins, `--browse`
lists data (superseded by `ncrystal browse`). Everything it shows is
produced by the same library calls an application makes.

## JSON queries (`ncrystal query`)

A low-level interface returning JSON for developer and test use: unit
conversions and numerical diagnostics, data and plugin listings, the atom
database, symmetry operations and equivalent reflections, site/position
conversions, the kernel-processing internals (reference evaluation and
sampling, cell surveys, kernel analysis), VDOS expansions, filter tables,
and per-object queries on `Info` (structure, sites). The physics chapters
name the queries that expose their internals.

## MiniMC

A Monte Carlo transport engine for a single convex volume of one material,
public since 4.3.0: geometries `sphere`, `box`, `slab` (infinite), `cyl`
(about the y axis, finite or infinite); sources `constant` (pencil),
`circular`, `isotropic` (point or sphere), with fixed, log-normal, uniform or
Maxwellian energy or wavelength distributions; neutrons are split at every
step into a transmitted part (tallied) and a forced collision, absorption is
applied as weight attenuation, Russian roulette and a scattering-order limit
control the depth; tallies of energy, energy transfer, wavelength, $q$,
scattering angle, scattering order, weight, with breakdowns by order, scored
as neutrons leave the volume; JSON persistence; multithreaded and vectorised
(above 10 MHz tally rates). Scenario strings such as
`"1.8Aa pencil on 1mfp sphere 1e6 times"` (energy in Aa, eV, meV, keV, MeV,
GeV, kT, or BT = Bragg threshold; thickness as a length or in mean free
paths) give quick interactive use; `geomcfg`/`srccfg`/`enginecfg` strings or
the Python objects give full control; `nctool --mc` and `ncrystal minimc`
front it. It is validated against a reference implementation in the test
suite and used there to check multiple scattering of single-crystal Ge and
PG, SANS materials, phases and absorption.

## ENDF export (`ncmat2endf`)

Since 4.2.0, `ncrystal_ncmat2endf` / `NCrystal.ncmat2endf.ncmat2endf(cfg)`
writes ENDF-6 thermal scattering (File 7) files from any single-phase
material whose dynamics are VDOS- or kernel-based, using endf-parserpy to
format and check the files: the $S(\alpha,\beta)$ kernel expanded by
NCrystal at the requested temperature (default $E_{max}$ 5 eV, optional
extra temperatures with the caveat that ENDF imposes one grid on all of
them), the elastic part in one
of three modes (`greater`: only the major component; `scaled`, the default:
the major component scaled to the total elastic cross section, with the
smallest-incoherent element carrying the coherent part in polyatomic
systems; `mixed`: both components, the modern but not universally supported
form), optional generalised information file sections and isotopic
expansion, and the metadata fields through `EndfMetaData`. This replaces
the earlier NJOY+NCrystal route and gives MCNP/PHITS users access to the
NCrystal models, with the documented caveats that NJOY processing follows,
that temperatures are fixed, and that the direct integrations (OpenMC,
Geant4, McStas) are preferable.

## McStas exports (`ncmat2hkl`, `mcstasunion`)

`ncrystal_ncmat2hkl` writes `.laz`/`.lau` reflection lists (with the
`spacegroup_setting` header since 5.0.0) for the McStas `PowderN` and
`Single_crystal` components from any crystalline material (chapter 52);
`ncrystal_mcstasunion` and `NCrystal.mcstasutils` generate McStas-Union
`NCrystal_process` instrument code and detect which components a material
needs (`cfgstr_detect_components`). McStas also has the direct
`NCrystal_sample` component and the Union process using the library itself.

## Filter tables (new)

A new facility (C++, C API, Python `NCrystal.filter.NCrystalFilter`)
produces a piecewise-linear table of the total cross section (scattering
plus absorption, per atom, with the number density for the macroscopic
conversion) of an isotropic material as a function of wavelength from 0 to
500 Aa, accurate to a relative tolerance of 1e-3: the cross section is
sampled densely (including points just below and above every known
discontinuity, i.e. Bragg edges and process energy-domain boundaries,
which are represented as two table points at the same wavelength), the
sample is reduced to the fewest points reproducing it within tolerance,
and the result is refined and verified until the tolerance holds
everywhere. Features narrower than the sampling cannot be found by
sampling, so a warning is issued for processes not known to be free of them
(none of the current ones). The intended consumer is GPU-friendly
attenuation-only filter components in McStas, hence the name.

## Other

`ncrystal_ncmat2cpp` embeds NCMAT files as C++ source (how the standard
library is compiled into the library); `nctool --bench` and developer
benchmarks track performance; the virtual API (chapter 52) serves OpenMC;
the Geant4 interface is a separate package (`ncrystal-geant4`).

## New or changed in 5.0.0

Filter tables (whole feature), the greatly extended JSON query system and
the `ncrystal query` tool, `ncmat2hkl` with settings, and `ncmat2endf`
exporting the 5.0.0 kernels. MiniMC is unchanged (public since 4.3.0).

## References

- NCrystal wiki pages "ncmat2endf" and "minimc" (with `minimc_geom`,
  `minimc_src`, `minimc_engine`, `minimc_scenario`).
- G. Schnabel, D. L. Aldama, R. Capote, "How to explain ENDF-6 to
  computers", arXiv:2312.08249 (endf-parserpy); ENDF-6 Formats Manual
  (ENDF-102).
- McStas `PowderN`, `Single_crystal`, `NCrystal_sample` and Union
  documentation.
