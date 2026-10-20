# Physics provided by plugins

## Summary

Since 2.2.0 NCrystal can be extended at run time by plugins: shared
libraries registering text-data, info, scatter or absorption factories,
typically activated by an `@CUSTOM_<NAME>` section in an NCMAT file and
adding a process on top of the standard physics. Since 4.0.0 plugins are
installed with `pip` and discovered automatically through
`ncrystal-pluginmanager`. Several physics capabilities that are *not* in the
core exist as plugins maintained by the NCrystal developers and
collaborators: extinction in polycrystals (CrysExtn), texture (CrysText),
coherent one-phonon kernels beyond the incoherent approximation (BzScope),
liquid dynamics with diffusion and $S(Q)$ (Liquids), paramagnetic scattering
(MagScat), nanodiamond SANS (SANSND), legacy `.nxs` input (nxslib), and
data-only plugins with extra materials. This chapter records what each
provides and how it attaches to the core, so that the capabilities of the
ecosystem can be stated separately from those of the library itself.

Code: `ncrystal_core/include/NCrystal/plugins/NCPluginMgmt.hh`,
`ncrystal_pypluginmgr/`, templates in `examples/plugin` and
`examples/plugin_dataonly`; the curated list in `devel/plugin_database.yml`
and the wiki page "CuratedPlugins".

## Mechanism

A plugin registers factories with a priority; for scatter factories the
usual pattern is to create the standard physics of the material and add
the plugin's process when the custom section is present (SANS plugins
instead declare multiphase capability and use the phase-pairing helpers of
chapter 41). Plugins can also replace models or provide only data files
(`plugins::<Name>/...`). Discovery is through `ncrystal-pluginmanager`
(plugins installed as Python packages), or explicit environment variables;
`nctool --plugins` and `ncrystal browse --plugins` list them. The custom
section's lines reach the plugin as raw strings through the `Info` object,
so the core never interprets them.

## The curated plugins

The list below follows the curated plugin page of the NCrystal wiki (which
carries installation instructions and source links); each plugin's own
repository is the reference for its status and usage.

- **CrysExtn**: extinction models for polycrystalline materials (Sabine's
  and Becker-Coppens' models among others), applied as per-plane
  extinction factors on the Bragg contributions; activated by a
  `@CUSTOM_CRYSEXTN` section with the model and its crystallite, grain and
  mosaic parameters. Used in the beryllium transmission study of Xu et al.
  (2025); the Becker-Coppens recipes were revisited in Kittelmann et al.
  (2026).
- **CrysText**: texture (preferred orientation) models for polycrystals,
  replacing the isotropic orientation average of the powder model.
- **BzScope**: coherent inelastic scattering kernels for isotropic
  materials computed with the BzScope tool, i.e. beyond the incoherent
  approximation; uses the core's ability to remove the coherent one-phonon
  part from the VDOS expansion (chapter 31) and the null extender (chapter
  33). Upstream at CINEMA/IHEP, mirrored at mctools.
- **Liquids**: on-the-fly modelling of liquids, combining a VDOS for the
  solid-like part with a diffusion model and an optional tabulated $S(Q)$
  for the coherent part, and a rotational model for liquid H$_2$/D$_2$;
  activated by `@CUSTOM_LIQUIDS` (DiJulio et al., ICANS XXV, 2026).
- **MagScat**: paramagnetic scattering on O$_2$ in clathrate hydrates and
  O$_2$-C$_{60}$ after Zimmer (2016); activated by `@CUSTOM_MAGSCAT`.
- **SANSND**: small-angle scattering of nanodiamond powders from a measured
  structure factor combined with a bulk diamond kernel, through the $I(Q)$
  helper of chapter 41 (Rizzi et al., 2023).
- **nxslib**: loading of `.nxs` files with their empirical non-Bragg
  cross-section curves (`inelas=external`); separate from the core since
  3.0.0 for licence reasons.
- **Data plugins**: LiquidHeliumData (liquid helium near 1 K, Granada et
  al. 2023), WaterData (light and heavy water kernels at many
  temperatures), UraniumOxideData (temperature-optimised UO$_2$ VDOS),
  GraphiticData (graphite oxides, fullerite, intercalation compounds),
  ClathrateHydrateData (clathrate hydrates).

The `ncrystal-extra` repository additionally holds large or experimental
data (e.g. the ENDF/B-VIII.0 kernels converted to NCMAT) and utilities.

## What this means for the capability statement

- Core NCrystal: kinematic Bragg diffraction (powder, mosaic single
  crystal, layered), incoherent elastic, kernel-based inelastic from VDOS or
  kernels, free gas, $1/v$ absorption, hard-sphere SANS, multiphase.
- With plugins: extinction, texture, coherent one-phonon inelastic, liquid
  diffusion models with $S(Q)$, paramagnetic scattering, nanodiamond SANS,
  `.nxs` input, and more materials.
- Plugins depend on the core's extension points, which the 5.0.0 work
  touched: the resolution of displacements and effective temperatures for
  plugin-built kernels (estimation allowed), the null extender, the
  per-order weights of the VDOS expansion, the switch handing coherent
  elastic scattering of amorphous solids to a plugin, and the
  state-of-matter rules (a liquid with a `vdos` section must remain
  expressible for the Liquids plugin: an open discussion point).

## Validation

- Each plugin carries its own self-test, the curated plugins are tested
  against new NCrystal releases by a dedicated CI, and the
  `ncrystal-verify` package runs the core's Python tests in installed
  environments.
- Physics validation is in the plugin publications.

## New or changed in 5.0.0

No change to the plugin mechanism itself since 4.0.0 (pip installation,
plugin manager). The 5.0.0 physics changes affect plugins through the hkl
lists (settings), the symmetry check of structures, the
displacement/effective-temperature resolution for plugin-built kernels, the
state-of-matter rules, and the new inelastic processing; plugin and data
updates will follow the release.

## References

- NCrystal wiki pages "Plugins", "PluginsDevelopment", "CuratedPlugins";
  `devel/plugin_database.yml`.
- S. Xu et al., J. Appl. Cryst. 58 (2025) 1957 (CrysExtn); T. Kittelmann et
  al., Acta Cryst. A82 (2026) (Becker-Coppens recipes); N. Rizzi et al.,
  Nucl. Sci. Eng. (2023) doi:10.1080/00295639.2023.2196926 (SANSND); O.
  Zimmer, Phys. Rev. C 93 (2016) 035503 (MagScat); J. R. Granada et al.,
  NIM A 1053 (2023) 168284 (liquid helium); M. Boin, J. Appl. Cryst. 45
  (2012) 603 (nxslib); D. D. DiJulio et al., ICANS XXV (2026) (Liquids).
