# Tools for creating materials

## Summary

NCrystal's Python package turns the physics models of the previous chapters
into a material-creation workflow: `NCMATComposer` assembles and validates
NCMAT data programmatically from crystal structures (CIF files, online
databases, ASE objects, existing materials), compositions, densities and
dynamics; `PhononDOSAnalyser` prepares phonon densities of states from
first-principles outputs; `ncrystal_cif2ncmat`, `ncrystal_hfg2ncmat`,
`ncrystal_vdos2ncmat` and `ncrystal_endf2ncmat` are the command-line fronts
for the common cases (crystals from CIF, hydrogen-rich amorphous materials
from functional groups, VDOS curves from text, ENDF thermal scattering
libraries). All of them end in NCMAT text that the core validates fully, so
the physics statements of chapters 10-43 apply to user-made materials
exactly as to the standard library.

Code: `ncrystal_python/src/NCrystal/` (`ncmat.py`, `cifutils.py`,
`vdos.py`, `hfg2ncmat.py`, and the `_cli_*.py` tools).

## NCMATComposer

- Sources: `from_cif` (local files, `codid::<id>` and `mpid::<id>` online
  entries from the Crystallography Open Database and the Materials Project,
  with `mp_apikey`), `from_ase`, `from_info` (any loaded material),
  `from_cfgstr`, `from_ncmat` (text), `from_hfg`.
- Structure: `set_cellsg` / `set_cellsg_cubic`, `set_atompos`,
  `refine_crystal_structure` and `verify_crystal_structure` (spglib-based
  standardisation and checks; the load-time check of chapter 10 is the
  authority).
- Composition and density: `set_composition`, `set_fraction`,
  `remap_atom`, `update_atomdb` (isotopes, enrichment, doping, chapter 43;
  partial occupancy emulated with a sterile pseudo-isotope, chapter 10),
  `set_density`.
- Dynamics (chapters 30-35): `set_dyninfo_vdos`, `set_dyninfo_vdosdebye`,
  `set_dyninfo_debyetemp`, `set_dyninfo_msd(temperature=...)` and
  `set_dyninfo_uiso` (writing the v8 `msd` keyword, hence v8 output),
  `set_dyninfo_scatknl`, `set_dyninfo_freegas`, `set_dyninfo_sterile`,
  `set_dyninfo_from_object`, `transfer_dyninfo_objects`,
  `allow_fallback_dyninfo`.
- Other: `set_state_of_matter`, `lock_temperature` /
  `set_default_temperature`, `add_secondary_phase`,
  `add_hard_sphere_sans_model`, custom sections (`add_raw_content`,
  `set_custom_section_data`), `set_plotlabel`.
- Output and use: `create_ncmat()` (v8 for crystals with a space group,
  with `@ATOMSITES` and the declared setting; the lowest version otherwise,
  full-precision cell parameters), `write(file)`, `register_as(name)`
  (in-memory, chapter 52), `load(cfg_params)`, `inspect()`, `plot_xsect()`.
- Rules enforced: positions must have the symmetry of the space group in
  its standard origin (an error at `create_ncmat`, with
  `refine_crystal_structure` as the remedy), `msd` data imply a solid,
  fractions and labels must be consistent.

## CIF import (`cif2ncmat`, `cifutils`)

`CIFLoader` reads CIF text, files or online entries with gemmi, uses spglib
to standardise and verify the structure, merges symmetry-equivalent atoms,
handles isotropic displacement parameters (`Uiso`, converted to the v8 `msd`
keyword through `set_dyninfo_uiso` with `--uisotemp`), and hands the result
to the composer. `ncrystal_cif2ncmat` adds dynamics options (a Debye
temperature, default 300 K, or VDOS data to attach) and validation plots
of the generated cross sections. 5.0.0 made the import robust across all
230 space groups and many CIF styles (tested with synthetic CIF data).

## Phonon densities of states (`PhononDOSAnalyser`, `vdos2ncmat`)

`PhononDOSAnalyser` is an immutable object reading DOS curves from raw
`(label, egrid, density)` lists, `AnyVDOS`/`DI_VDOS` objects, or Quantum
ESPRESSO `matdyn.dos` files (total and partial DOS columns), with methods to
inspect (`plot`, `plot_gn` for the Sjolander $G_n$ of chapter 31), `merge`,
`drop`, `update_label`, apply a low-energy cut-off (`apply_cutoff`, judged
with `plot_cutoff_effects` and `plot_cutoff_effects_on_xsects`, which show
the effect on the $G_n$ and on cross sections, since the $E^2$
extrapolation below the first point replaces whatever the calculation
produced there), regularise the grid (`apply_regularisation`, so that the
grid is equidistant with 0 as an implied point, as the loader will
otherwise do at load time), and finally `apply_to(composer)` mapping labels
to the composer's atoms. `ncrystal_vdos2ncmat` is the older command-line
route for a single curve from a text file ("somewhat outdated").
`NCrystal.vdos` also exposes `analyseVDOS`, `extractGn`, `extractKnl`,
`createVDOSDebye`, `debyeIsotropicMSD`, `debyeTempFromIsotropicMSD` for
scripting the physics of chapters 31 and 35 directly.

## Hydrogen-rich amorphous materials (`hfg2ncmat`)

`ncrystal_hfg2ncmat` builds an NCMAT file for a hydrogen-rich amorphous
material from a chemical formula, a density and the count of hydrogen
functional groups (CHali, CHaro, CH2, CH3, NH, NH2, NH3, OH, SH), using the
measured group VDOS curves of Romanelli et al. (J. Phys.: Condens. Matter
2021, doi:10.1088/1361-648X/abfc13) for hydrogen and a Debye model (default
400 K) for the other atoms; `@STATEOFMATTER solid` and `@DENSITY` complete
the file. Most of the 16 amorphous stdlib materials were made this way
(chapter 51). Also available as `NCMATComposer.from_hfg`.

## ENDF import (`endf2ncmat`)

`ncrystal_endf2ncmat` converts ENDF-6 thermal scattering (File 7) data to
`scatknl` @DYNINFO sections through PyNE, one file per temperature, for
expert users who have evaluated kernels (the 2.0.0 release converted the
whole ENDF/B-VIII.0 thermal sublibrary this way into `ncrystal-extra`).
Converted kernels are direct kernels at fixed temperatures (chapter 30); with
NCMAT v8 they can carry `effective_temperature` and displacement data, or
have them estimated (chapter 33).

## Verification helpers

`ncrystal_verifyatompos` checks that atom positions in an NCMAT file are
consistent with its space group using ASE (largely superseded by the
load-time check of chapter 10); the composer's `verify_crystal_structure`
and the `inspect`/`plot_xsect` methods, `nctool --dump` and the
`produce_validation_plot` functions serve the same purpose from the physics
side.

## Validation

- `tests/scripts/composer.py`, `composer2.py`, `cif2ncmat.py`,
  `cif2ncmatgen.py`, `long_cif2ncmatgen.py`, `verifyatompos.py`, `hfgpy.py`,
  `vdos2ncmat.py`, `analysevdos.py`, and the ncrystal-notebooks repository's
  material-creation notebooks (the most complete user documentation).

## Approximations and limitations

- The composer cannot express what NCMAT cannot (chapter 50): occupancy
  (emulated), anisotropic displacements (reduced to isotropic),
  non-standard settings (standardised with spglib).
- `cif2ncmat` depends on the quality of the CIF (disorder, partial sites
  and mixed sites need manual decisions); `hfg2ncmat` assumes the group
  VDOS curves are transferable between molecules.
- `endf2ncmat` carries the ENDF evaluation's limitations (fixed
  temperatures, grids) into NCMAT.

## New or changed in 5.0.0

The composer writes NCMAT v8 for crystals (sites, settings) and for `msd`
data, with full-precision cells; structure verification and refinement are
built on the core's symmetry machinery; `cif2ncmat` is more robust. Planned:
proper site occupancy support, and always emitting v8.

## References

- G. Romanelli et al., "Thermal neutron cross sections of amorphous
  hydrogenous materials from vibrational spectra of functional groups",
  J. Phys.: Condens. Matter 33 (2021) 285901, doi:10.1088/1361-648X/abfc13.
- spglib (A. Togo et al.), gemmi (M. Wojdyr), ASE; the Crystallography
  Open Database and the Materials Project; PyNE (ENDF reading).
- The ncrystal-notebooks repository (material creation tutorials).
