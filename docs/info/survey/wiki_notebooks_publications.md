<!-- Unedited AI-generated survey inventory (2026-10-20), kept as raw source notes for docs/info. Not reviewed. To be removed from history before merging to develop. -->

# Inventory of public NCrystal documentation sources (surveyed 2026-10-20)

Sources: shallow clones of `https://github.com/mctools/ncrystal.wiki.git` (-> `/home/ubuntu/.claude/jobs/cd3edb19/tmp/wiki`) and `https://github.com/mctools/ncrystal-notebooks.git` (-> `/home/ubuntu/.claude/jobs/cd3edb19/tmp/notebooks`), plus WebFetch/WebSearch/Crossref/arXiv API queries. Nothing else downloaded.

## 1. GitHub wiki (https://github.com/mctools/ncrystal/wiki)

### 1.1 Complete page list (from `_Sidebar.md` + repo files)

Sidebar order: Home, Get-NCrystal, Using-NCrystal, Data-library, CfgRefDoc ("Cfg string parameters"), Included-physics, NCMAT-format, minimc, ncmat2endf, Plugins, CuratedPlugins, PluginsDevelopment, Announcement-Release2.0.0, -2.1.0, -2.5.0, -2.7.0, -3.0.0, -3.1.0, -3.2.0, About, Contact.
Not in sidebar but present: minimc_engine, minimc_geom, minimc_scenario, minimc_src, VDOSAnimations, CfgStrings (stub: "under construction, see Using-NCrystal"), For-developers (stale travis-ci badges), plus asset dirs `datalib/` (per-material .ncmat + inspect/VDOS/onlinedb-validation PDF/PNG plots + dump txt), `vdosanim/`, and announcement image dirs.
Raw markdown: `https://raw.githubusercontent.com/wiki/mctools/ncrystal/<Page>.md`. Sizes: NCMAT-format 897 lines, CfgRefDoc 262, Announcement-3.0.0 (largest), Announcement-2.5.0 252, Data-library 164, ncmat2endf 139, Using-NCrystal 113, minimc_* 30-112 each, Included-physics 39, Home 22.

### 1.2 Physics-relevant pages, summarised

- **Home** (https://github.com/mctools/ncrystal/wiki/Home): one-paragraph overview (coherent/incoherent, elastic/inelastic, powders, mosaic & layered single crystals, amorphous, liquids, gases, multiphase, enriched, SANS); integrations (OpenMC, Geant4, McStas); the 3 citation papers; a release-by-release feature chronology v1.0 -> v4.4 (v2.0 inelastic/kernels, v2.1 atomdb, v2.2 plugins, v2.5 C++/thread-safety, v2.6 VDOS everywhere/Debye temps obsoleted, v2.7 amorphous + big datalib growth, v3.0 multiphase+SANS, v3.1 UCN/VCN, v3.2 gasmix, v3.6 Python API/material creation utils, v3.7 PyPI, v3.9/4.3 MiniMC, v4.0 Windows, v4.1 virtual API for OpenMC runtime plugin, v4.2 ncmat2endf, v4.4 `strain`).
- **Included-physics** (https://github.com/mctools/ncrystal/wiki/Included-physics): short page. Elastic: Bragg diffraction from explicit hkl lists with structure factors from isotropic Debye-Waller factors (MSD from VDOS, fallback Debye temperature); powder, Gaussian-mosaic single crystal, layered crystal (PG); back-scattering / rocking-curve asymmetry handled consistently. Incoherent elastic: same DW model, eq. 23 of CPC 2020 paper. SANS: framework since v3.0 with only a hard-sphere model. Inelastic: S(alpha,beta) kernels in incoherent approximation, kernels either embedded in NCMAT or generated on the fly from 1D phonon DOS (points to Announcement-2.0.0); notes "we have yet to publish a detailed publication on the inelastic models", only the JCP 2019 sampling paper. Multiphase/enriched: refs 3.0.0 and 2.1.0 announcements. Absorption: 1/v model only.
- **Using-NCrystal** (https://github.com/mctools/ncrystal/wiki/Using-NCrystal): setup; the uniform cfg-string syntax with worked examples (temp units, `density=0.45x`/`gcm3`, `gasmix::...`, mosaic Ge monochromator with `mos`/`dir1`/`dir2`, `bragg=0` for single-crystal filters, multiphase `phases<...>`); links to CfgRefDoc.
- **CfgRefDoc** (https://github.com/mctools/ncrystal/wiki/CfgRefDoc): auto-generated reference of all cfg parameters, grouped: Base (`atomdb, dcutoff, dcutoffup, infofactory, strain, temp`), Basic scattering (`coh_elas, incoh_elas, inelas, sans, scatfactory, vdoslux`, pseudo-params `bkgd, bragg, comp, elas`), Advanced single-crystal (`dir1, dir2, dirtol, lcaxis, lcmode, mos, mosprec, sccutoff, ucnmode`), Absorption (`absnfactory`), Special (`density, phasechoice`). Each with type, units, default, description.
- **Data-library** (https://github.com/mctools/ncrystal/wiki/Data-library): static table for release v4.2.10 of 134 NCMAT files (110 crystalline, 16 amorphous solids, 5 gases, 2 liquids, 1 "unknown") with state, space group, density, per-element inelastic model (VDOS / Debye / kernel) and links to inspect plots, VDOS plots and online-database validation plots (`datalib/`); also `nctool --browse` usage, plugin data, gasmix syntax, cif2ncmat/hfg2ncmat pointers, link to dgcode material cookbook.
- **NCMAT-format** (https://github.com/mctools/ncrystal/wiki/NCMAT-format): the full format spec v1..v7: @CELL, @ATOMPOSITIONS, @SPACEGROUP, @DEBYETEMPERATURE (v1); @DYNINFO with dynamic models Sterile / Free gas / Scattering kernel / VDOS / VDOSDebye, @DENSITY, "effects of @DYNINFO on inelastic modelling" (v2); @ATOMDB, @CUSTOM_* (v3); cell/Debye changes (v4); @DYNINFO changes, @STATEOFMATTER and its physics effects (v5); @OTHERPHASES (v6); @TEMPERATURE (v7); emacs mode. (Same content as repo `docs/ncmat_doc.md`.)
- **Announcement-Release2.0.0** ("Inelastic physics gets some love"): the de-facto inelastic physics write-up: S(alpha,beta)/S(q,w) kernels, kinematic region, extrapolation beyond tabulated range, on-the-fly Sjolander phonon expansion from VDOS (`vdoslux`, `vdos2ncmat`), Debye-model fallback VDOS, BeO comparison of direct ENDF kernel vs VDOS vs Debye; NCMAT v2 and cfg-string updates.
- **Announcement-Release2.1.0** ("It's time to get dirty"): flexible atomic definitions (`atomdb`: isotopes, enrichment, overriding scattering lengths/xs), in-memory NCMAT, physics-extension hooks.
- **Announcement-Release2.5.0**: technical (thread safety, RNG streams, TextData/data-source factories, modern C++); little physics.
- **Announcement-Release2.7.0**: amorphous solids support, `ncrystal_hfg2ncmat` (hydrogen-rich functional-group approximation), 16 new amorphous + 64 new crystalline materials, NCMAT creation/verification tools, init speedup.
- **Announcement-Release3.0.0**: multiphase materials (cfg `phases<..>` and NCMAT @OTHERPHASES), SANS framework with hard-sphere model (`@CUSTOM_HARDSPHERESANS` example of diamond in D2O), density overrides, quick unstructured materials (`solid::AgBr/6.48gcm3`), new `comp`/`sans` params, hkl list improvements, McStas .laz/.lau export/import, license cleanup.
- **Announcement-Release3.1.0** ("UCN/VCN"): improved kernel sampling near the kinematic end-point (superfluid He4 example), dedicated biased UCN-production process (`ucnmode`), exact S(alpha,beta) evaluation/integration utilities.
- **Announcement-Release3.2.0**: `gasmix::` on-the-fly gas mixtures (fractions, pressure, temperature, humidity).
- **VDOSAnimations**: webm animations of how phonon orders in the Sjolander expansion build up S(alpha,beta) (Al, PE, graphite, V).
- **minimc** + minimc_engine/geom/scenario/src: MiniMC built-in stepping Monte Carlo (box/sphere/cylinder/slab, pencil/circular/isotropic sources, fixed/log-normal/uniform/Maxwell spectra, tallies incl. energy, q, angle, energy loss, scatter order; >10 MHz vectorised/multithreaded); cfg-string references for engine (parameters, tallies), geometry, source and shorthand "scenario".
- **ncmat2endf** (https://github.com/mctools/ncrystal/wiki/ncmat2endf): v4.2 export of NCMAT materials to ENDF-6 (via endf-parserpy), caveats (NJOY processing, only VDOS-based single-phase incoherent-approx materials, fixed temperatures), elastic modes `greater`/`scaled`/`mixed`, Python/CLI usage, references (ND2025 Madrid proceedings, NJOY+NCrystal NIMA 2022, endf-parserpy arXiv:2312.08249).
- **Plugins / PluginsDevelopment**: mechanism (since v2.2, pip-installable since v4.0); developing data-only vs physics plugins.
- **CuratedPlugins** (https://github.com/mctools/ncrystal/wiki/CuratedPlugins), each with install/PyPI/source URL: GraphiticData (graphite oxides, C60), UraniumOxideData (UO2 VDOS per temperature, Noguere 2020), WaterData (H2O/D2O wide T range from ENDF/B-VIII.0 CAB kernels), BzScope (coherent inelastic kernels beyond incoherent approx, IHEP), ClathrateHydrateData, CrysExtn (Sabine and Becker-Coppens extinction for polycrystals, github.com/XuShuqi7/ncplugin-CrysExtn), CrysText (texture), LiquidHeliumData (He-4 near 1 K, NIMA 2023 doi:10.1016/j.nima.2023.168284), Liquids (on-the-fly liquids, ICANS XXV 2026), MagScat (paramagnetic scattering on O2 in clathrates, Zimmer PRC 93 035503), SANSND (nanodiamond SANS, doi:10.1080/00295639.2023.2196926), nxslib (.nxs loader), DummyDataPlugin, DummyPlugin.
- **About**: citation requests (3 papers), funding (BrightnESS 676548, HighNESS 951782). **Get-NCrystal**: conda/pip/source/other. **Contact**: discussions, issues, ncrystal-developers@cern.ch.

Gaps worth noting: no wiki page on absorption beyond "1/v"; no "Validation" page (validation only as per-material plots in datalib and the BeO plots in the 2.0.0 announcement); no page on strain beyond CfgRefDoc; the inelastic model is documented only via the 2.0.0/3.1.0 announcements.

## 2. ncrystal-notebooks (https://github.com/mctools/ncrystal-notebooks)

README groups them as basic / material creation / misc / contributed; all have Colab links. Titles from the first markdown cell:

Basic:
- `ncrystal1_basic_01_Introduction_and_Python_API` -- "Introduction and basic usage of the NCrystal Python API": Info/Scatter/Absorption objects, cross sections, sampling.
- `ncrystal1_basic_02_Data_Infrastructure_and_StdDataLib` -- data library, data sources, browsing/dumping NCMAT files.
- `ncrystal1_basic_03_Scatter_patterns_with_the_builtin_MiniMC_framework` -- "The NCrystal MiniMC simulation engine": geometry + multiple scattering in a single convex volume.
- `ncrystal1_basic_04_NCrystal_as_backend_for_third_party_apps` -- using NCrystal inside McStas/OpenMC/Geant4.

Material creation:
- `ncrystal2_advanced_01_Creating_materials_and_the_NCMATComposer` -- building NCMAT with `NCMATComposer`.
- `ncrystal2_advanced_02_Import_crystal_structure_from_CIF_or_databases` -- CIF / COD / Materials Project import (v3.6+).
- `ncrystal2_advanced_03_Add_phonon_info_with_PhononDOSAnalyser_with_QuantumEspresso_example` -- curating VDOS curves (Quantum Espresso) with `PhononDOSAnalyser`; stresses VDOS matters for elastic (DW) as well as inelastic.

Misc:
- `ncrystal2_advanced_04_VDOS2KNL_Connecting_phonons_to_inelastic_scattering` -- interactive VDOS -> S(alpha,beta) expansion (Sjolander method) widgets.
- `ncrystal2_advanced_05_export_ENDF` -- ncmat2endf ENDF export example.
- `misc/ncrystal_sapphire_filter` -- sapphire single-crystal filter transmission; `bragg=0` quick model vs proper single-crystal modelling.

Contributed, NEUWAVE-12 tutorial (`notebooks/contributed/NEUWAVE-12/`):
- `Neutron_filters_exercise` -- xs and transmission of polycrystalline and single-crystal filters vs experimental data (Freund).
- `Sample_container_transmission_exercise` -- xs/transmission of metals used for sample holders.
- `Water_Ice` -- wavelength dependence of total xs of H2O and ice vs temperature (neutron imaging application, Siegwart).
- `Extinction_correction_exercise` -- computing coherent elastic xs from hkl info and applying Sabine extinction.
- `Installing_Plugins_Texture_exercise` -- CrysText texture plugin, modification of coherent elastic xs.
- `In_scattering` -- computing in-scattering.
- `Polymers_with_Functional_Group_Approximation` -- PLA via the functional-group approximation (hfg2ncmat).
- `Transmission_with_NCrystal_and_McStas` -- conda install + simple McStas transmission instrument.

Contributed, ND2025 (`notebooks/contributed/ND2025/`):
- `ND2025_ncmat2endf` -- export to ENDF and use in OpenMC.
- `ND2025_OpenMC_NCrystal` -- NCrystal materials directly in OpenMC.
- `ND2025_Extinction_ENDF_library` -- CrysExtn extinction plugin + ncmat2endf to produce an extinction-corrected ENDF library.

## 3. Website https://mctools.github.io/ncrystal/

Exists but is a bare landing page: single heading "NCrystal" and one link "Click here to access the NCrystal wiki" (https://github.com/mctools/ncrystal/wiki/). No sections, no API/Sphinx docs there.

## 4. Publications

Core (requested citations, from About page):
- X.-X. Cai, T. Kittelmann, "NCrystal: A library for thermal neutron transport", Comput. Phys. Commun. 246 (2020) 106851, https://doi.org/10.1016/j.cpc.2019.07.015 (arXiv:1901.08890).
- T. Kittelmann, X.-X. Cai, "Elastic neutron scattering models for NCrystal", Comput. Phys. Commun. 267 (2021) 108082, https://doi.org/10.1016/j.cpc.2021.108082 (arXiv:2012.04294).
- X.-X. Cai, T. Kittelmann et al., "Rejection-based sampling of inelastic neutron scattering", J. Comput. Phys. 380 (2019) 400-407, https://doi.org/10.1016/j.jcp.2018.11.043.

Related / tooling / plugins:
- K. Ramic, J.I. Marquez Damian, T. Kittelmann, D.D. Di Julio, D. Campi, M. Bernasconi, G. Gorini, V. Santoro, "NJOY+NCrystal: An open-source tool for creating thermal neutron scattering libraries with mixed elastic support", NIM A 1027 (2022) 166227, https://doi.org/10.1016/j.nima.2021.166227 (arXiv:2108.11737).
- J.I. Marquez Damian, D.D. DiJulio, G. Muhrer, "Nuclear data development at the European Spallation Source", UCANS-Web 2020, arXiv:2103.06133.
- K. Ramic et al., "Advances in Nuclear Data and Software Development for the HighNESS Project", Nucl. Sci. Eng. 198(1) (2024) 74-82, https://doi.org/10.1080/00295639.2023.2184196.
- N. Rizzi et al., nanodiamond SANS models/benchmark (basis of SANSND plugin), Nucl. Sci. Eng. 2023, https://doi.org/10.1080/00295639.2023.2196926.
- Liquid helium kernel (basis of LiquidHeliumData plugin), NIM A 2023, https://doi.org/10.1016/j.nima.2023.168284.
- L. Zanini et al., "Moderator developments in HighNESS...", EPJ Web Conf. 298 (2024) 03001, https://doi.org/10.1051/epjconf/202429803001 (uses NCrystal-based data).

2024-2026 (extinction, ENDF export, third-party integrations):
- S. Xu et al., "Impact of extinction effects on neutron transmission in solid beryllium metal", J. Appl. Cryst. (2025), https://doi.org/10.1107/S1600576725007939.
- T. Kittelmann, D. DiJulio, S. Xu, J.I. Marquez Damian, "Revisiting Becker-Coppens (1974): updated recipes for estimating extinction factors in spherical crystallites", Acta Cryst. A 82 (2026) 163-178, https://doi.org/10.1107/S2053273326001245.
- J.I. Marquez Damian, T. Kittelmann, D. Di Julio, S. Xu, D. Roubtsov, G. Muhrer, "Open source pipelines for thermal scattering data generation: from NCrystal to ENDF-6", ND2025 (16th Int. Conf. Nuclear Data for Science and Technology), Madrid, June 2025 (proceedings DOI not yet found on Crossref).
- Grammer et al., "An NCrystal extension for thermal neutron scattering in MCNP", NIM A (2026), https://doi.org/10.1016/j.nima.2025.171079.
- M. Tang, Z.-Y. Pan, N. Yang, T. Kittelmann, X.-X. Cai, "BzScope: an absolute cross section calculator for neutron-phonon scattering", arXiv:2511.01178 (Nov 2025; NCrystal plugin for coherent 1- and 2-phonon scattering).
- Jalil et al., "Validation of ENDF/B-VIII.0 and NCrystal thermal scattering laws for graphite in criticality safety analysis using OpenMC", Nuclear Analysis (2026), https://doi.org/10.1016/j.nucana.2026.100219 (user-side validation).
- D.L. Christensen et al., "Simulating neutron diffraction from deformed mosaic crystals in McStas", arXiv:2511.19081 (2025; mentions NCrystal).
- Conference talks, ICANS XXV (Malmo, 13-17 Apr 2026): T. Kittelmann, NCrystal status/roadmap (https://indico.ess.eu/event/3811/contributions/24451/); S. Xu et al., "Quantitative data analysis for neutron scattering experiments using NCrystal" (clathrate magnetic scattering, Be extinction at HIPPO, hydrogen imaging at IMAT/RADEN) (https://indico.ess.eu/event/3811/contributions/24505/); D. Di Julio et al., "A liquids plugin for NCrystal" (on-the-fly liquids: water, LH2, LD2; proceedings submitted) (https://indico.ess.eu/event/3811/contributions/24554/).
- NEUWAVE-12 tutorial slides (Sept 2024): https://indico.ess.eu/event/3439/sessions/7330/attachments/15695/30086/2024-09-05-neuwave-workshop.pdf.

Not found: any dedicated published paper on NCrystal's inelastic/VDOS-expansion model, on MiniMC, on multiphase/SANS, or on texture (CrysText); the wiki itself states the inelastic models are not yet the subject of a dedicated publication.
