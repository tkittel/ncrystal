# Glossary

- **alpha, beta**: dimensionless momentum and energy transfer of a
  scattering event in units of $kT$; NCrystal's $\alpha$ omits the mass
  ratio $A$ of the ENDF convention (chapter 01).
- **alpha2x**: the constant converting $\alpha$ to $x = 2W$ in the phonon
  expansion, $2m_n kT\langle u^2\rangle/\hbar^2$ (chapter 31).
- **Analyser (SABAnalyser)**: estimator of effective temperature and
  displacement from a tabulated kernel (chapter 33).
- **Bounded-cell / full-cell sampling**: the two sampling recipes of the
  next-gen kernel processing (chapter 32). Full-cell sampling picks a
  touched grid cell by its full integral and draws in the whole cell; it is
  used wherever its acceptance rate is at least 0.15. Bounded-cell sampling,
  for very low energies, picks cells by their integral within the
  kinematically accessible region and uses overlays restricted to that
  region for the cells it crosses.
- **Bragg threshold / cut-off**: the longest wavelength $2d_{max}$ at which
  Bragg diffraction occurs (chapter 20).
- **cfg-string**: `"<data>;par=val;..."` material configuration (chapter
  00).
- **Choice code / setting**: the ITVB code identifying one of the 530
  space-group settings, e.g. `227:2` (chapter 11).
- **Demi-normal**: one of each pair of opposite plane normals (chapter 12).
- **Direct kernel**: an $S(\alpha,\beta)$ table given in the data
  (`scatknl`), as opposed to one expanded from a VDOS (chapter 30).
- **dcutoff / dcutoffup**: d-spacing limits of the hkl list (chapter 12).
- **Debye-Waller factor**: $e^{-2W}$, $2W = Q^2\langle u^2\rangle$
  (chapter 35).
- **DynamicInfo (DI_*)**: the per-element dynamics description on an `Info`
  object: sterile, free gas, direct kernel, VDOS, VDOS-Debye (chapter 30).
- **Effective temperature ($T_{eff}$)**: mean kinetic energy per atom in
  temperature units, incl. zero-point motion; drives the SCT model
  (chapters 33, 35).
- **egrid**: the neutron energy grid of a processed kernel, or the NCMAT
  field setting it (chapter 32).
- **Extender**: the analytical model (free gas or SCT) used above a
  kernel's energy range (chapter 33).
- **Filter table**: piecewise-linear total cross-section table for
  attenuation-only components (chapter 61).
- **$G_n$**: Sjolander's $n$-phonon energy-transfer functions (chapter 31).
- **Incoherent approximation**: using the incoherent (self) scattering
  function for coherent inelastic scattering too (chapter 30), and, for
  amorphous solids, for coherent elastic scattering (chapter 24).
- **Info**: the immutable object holding structure, composition, dynamics,
  density, temperature and hkl list of a material (chapter 00).
- **knllux / vdoslux / sablux**: luxury levels of the kernel processing,
  the VDOS expansion and the SAB configuration (chapters 31, 32).
- **Kinematic region $D(E)$**: the accessible $(\alpha,\beta)$ area at
  energy $E$ (chapter 32).
- **Legacy / next-gen**: the 2.0.0-4.4.6 inelastic algorithms versus the
  5.0.0 ones (chapters 30-32).
- **MiniMC**: NCrystal's Monte Carlo engine for single volumes (chapter
  61).
- **msd**: mean-squared displacement along one axis, $U_{iso}$ (chapter
  35).
- **Multiplicity**: number of symmetry-equivalent planes in an hkl family,
  Friedel pairs included (chapter 12).
- **NCMAT**: NCrystal's material text format, versions v1-v8 (chapter 50).
- **Orbit / site**: the set of positions generated from one site by the
  space group; the site is its representative (chapter 10).
- **Overlay**: a distribution sampled instead of the target one and
  corrected by rejection (chapters 32, 34).
- **Plugin**: a dynamically loaded factory extending NCrystal (chapter 53).
- **Process**: a `Scatter` or `Absorption` object with cross sections and
  sampling (chapter 00).
- **Quick material**: `freegas::`, `solid::`, `gasmix::` generated data
  (chapter 52).
- **SABData / ScatKnlData**: standard and input forms of a kernel table
  (chapter 30).
- **SCT**: short-collision-time approximation (chapter 33).
- **sccutoff**: d-spacing below which single-crystal planes are treated as
  powder (chapter 21).
- **SLD**: scattering-length density, for SANS contrast (chapter 41).
- **State of matter**: solid, liquid, gas or unknown; decides elastic
  physics of non-crystals (chapters 24, 42).
- **suggestedEmax**: the energy up to which an expanded kernel is complete
  (chapter 31).
- **Teff**: see effective temperature.
- **ucnmode**: legacy UCN-production split of the inelastic process
  (chapter 30).
- **VDOS**: vibrational (phonon) density of states (chapter 31).
- **Virtual file**: in-memory data registered under a file name (chapter
  52).
