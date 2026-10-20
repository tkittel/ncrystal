# Bragg diffraction in powders and polycrystals

## Summary

For a crystalline material without orientation information NCrystal models
coherent elastic scattering as ideal powder (non-textured polycrystal)
Bragg diffraction: every reflection family $hkl$ with $d_{hkl} \ge \lambda/2$
contributes $\lambda^2 d_{hkl}|F_{hkl}|^2/(2V_{uc}n_a)$ to the cross section
per atom, and a scattering event picks a family with probability
proportional to its contribution and scatters the neutron elastically onto
the Debye-Scherrer cone of that family. The hkl list with its structure
factors is prepared once from the crystal structure (chapter 12). Cross
sections and samples cost one binary search each, at tens of millions per
second. The model is the one documented in the 2020 CPC paper (section 2.2)
and the 2021 elastic paper (section 3); no extinction, texture, grain-size
or strain-broadening effects are included in the core (plugins exist, see
chapter 53).

Code: `ncrystal_core/include/NCrystal/internal/powderbragg/NCPowderBragg.hh`,
`ncrystal_core/src/powderbragg/NCPowderBragg.cc`; assembled by
`ncrystal_core/src/stdscatfactory/NCStdScatFact.cc`.

## Physics

### Cross section (derivation)

Coherent elastic scattering of a crystal is, per unit cell,

$$S^{el}_{coh}(\vec Q,\omega) = \frac{(2\pi)^3}{V_{uc}}\delta(\hbar\omega)\sum_{hkl}\delta(\vec Q - \vec\tau_{hkl})\,|F(\vec\tau_{hkl})|^2,$$

with the unit-cell structure factor $F(\vec Q) = \sum_i b_i e^{-W_i(Q)}e^{i\vec Q\cdot\vec p_i}$
(CPC 2020 eqs. 25, 31). A reflection on the plane family $hkl$ requires
$\lambda \le 2d_{hkl}$ and occurs at the Bragg angle $\sin\theta_B = \lambda/2d_{hkl}$.
Averaging the orientation of the crystallites over the sphere, each
reciprocal-lattice vector $\vec\tau_{hkl}$ of length $2\pi/d_{hkl}$
contributes $\lambda^2 d_{hkl}|F_{hkl}|^2/(16\pi^3)$ per unit cell per
plane (CPC 2020 eq. 34), and summing over the plane normals in the family
(multiplicity $m_{hkl}$, i.e. all symmetry-equivalent $\vec\tau$, both signs
counted) gives the per-atom cross section

$$\sigma_{powder}(\lambda) = \frac{\lambda^2}{2V_{uc}\,n_a}\sum_{hkl:\ \lambda\le 2d_{hkl}} d_{hkl}\,|F_{hkl}|^2$$

(CPC 2020 eq. 35, elastic paper eq. 13), where the sum runs over all
planes, equivalently over families weighted by multiplicity, and $n_a$ is
the number of atoms per unit cell. In the general mosaic formulation of the
elastic paper (section 2) this is $\sigma_{hkl} = q_{hkl}\,g_{hkl}$ with
$q_{hkl} = \lambda^3|F|^2/(V_{uc}n_a\sin 2\theta_B)$ and the geometric
factor $g = \sin\alpha/2$ for an isotropic orientation density
$W = 1/4\pi$, $\alpha = \pi/2 - \theta_B$.

The cross section is a sum of step functions of $\lambda$ (the Bragg
edges at $\lambda = 2d_{hkl}$) times $\lambda^2$, vanishing above the
largest $2d$ (the Bragg cut-off, below which only incoherent elastic and
inelastic scattering remain).

### Scattering outcome

The neutron keeps its energy and scatters by $2\theta_B$ with uniform
azimuth around its direction: $\mu = \cos 2\theta_B = 1 - 2\sin^2\theta_B
= 1 - \lambda^2/2d^2$. Since all planes of a family and all crystallite
orientations contribute equally in a powder, nothing further needs to be
sampled.

### Debye-Waller factors

$W_i(Q) = \frac12 Q^2\langle u_i^2\rangle$ with the isotropic displacement
of chapter 35, so $e^{-W_i(\tau_{hkl})} = e^{-2\pi^2\langle u_i^2\rangle/d_{hkl}^2}$.
The temperature dependence of the Bragg cross sections enters only through
these factors (and through nothing else: the lattice is not expanded with
temperature unless the user applies `strain`).

## Inputs and configuration

- The hkl list of the `Info` object (`hklList()`: $d$, $|F|^2$,
  multiplicity, and the list type; chapter 12) and `StructureInfo`
  ($V_{uc}$, $n_a$). Data sources: NCMAT crystal structures (structure
  factors computed), `.laz`/`.lau` files (structure factors read), or
  plugin-provided lists.
- `dcutoff` (default automatic: 0.1 Aa, or 0.2 Aa for cells with more
  than 40 atoms): planes with $d$ below it are not generated; `dcutoffup`
  removes planes above it (chapter 12).
- `coh_elas` / `bragg` (default on) enables the component; `elas=0` or
  `comp=...` remove it with the other elastic parts.
- `strain` (range [-0.5, 0.5]) scales the cell lengths by $1 + $`strain`
  (isotropic), shifting all Bragg edges.
- `temp` through the displacements.
- Orientation parameters (`mos`, `dir1`, `dir2`, `lcaxis`) switch to the
  single-crystal models of chapters 21-22; in a single crystal, planes with
  $d <$ `sccutoff` (default 0.4 Aa) are still treated with this powder
  model (next section).

## Implementation

Planes are sorted by decreasing $d$ (equal $d$ merged), and two arrays are
built: the Bragg-edge energies $E_{edge} = h^2/(8m_n d^2)$ and the
cumulative sums of $d|F|^2 m$. A cross section is then a binary search for
the last edge below $E$ followed by one multiplication (the $\lambda^2$
prefactor becomes $1/E$), and sampling a second binary search in the
cumulative sums to pick the plane, with the outgoing direction built from
$\sin^2\theta_B = E_{edge}/E$ and an isotropic azimuth without
trigonometric calls (elastic paper eq. 15). The process declares the Bragg
threshold as the lower bound of its energy domain, so callers can skip it
for colder neutrons, and instances merge by merging their edge lists
(phases, or the powder part of a single crystal). Cross sections cost about
50 MHz and samples about 10 MHz on a 2020-era core (elastic paper section
7); 5.0.0 adds branchless searches and multi-energy evaluation for the
filter tables (chapter 61).

## Validation

- CPC 2020 paper section 4: hkl lists and positions cross-checked against
  nxslib/SgInfo for all data-library files; McStas-simulated powder
  patterns refined with GSAS-II (four materials, R below 1.9%) and
  idealised patterns refined with FullProf for 17 materials (R below
  2.7%), recovering lattice parameters, positions and displacements;
  measured $|F|^2$ for monoclinic SiLu$_2$O$_5$ reproduced; EXFOR total
  cross sections (Mo, Ni, Zn, UO$_2$, Sn) agree in the Bragg-edge region.
- Elastic paper: an isotropically illuminated single crystal reproduces
  the powder cross section to better than 1e-3 (1e-4 at `mosprec` 1e-5),
  tying the two models together.
- Repository tests: `tests/src/app_fillhkl` (Bragg threshold, partial
  lists), `tests/scripts/fillhkl.py`, `hklobj.py`, `nosg.py`, `strain.py`,
  `dkelas.py`, the many `.log`-checked `nctool` dumps and cross-section
  curves of the stdlib materials, and the symmetry self-check of structure
  factors at load time (chapter 12).

## Approximations and limitations

- Ideally imperfect crystallites in the Born approximation: no extinction
  (primary or secondary), no multiple scattering within grains.
- No preferred orientation (texture), no finite grain size or strain
  broadening of the Bragg edges, no lattice-parameter distributions.
- Isotropic harmonic Debye-Waller factors, one displacement per element
  (chapter 35).
- Extinction and texture for powders are available through the CrysExtn
  and CrysText plugins (chapter 53), which scale or replace the per-plane
  contributions; the Becker-Coppens recipes used there have been revisited
  by the NCrystal authors (Kittelmann et al. 2026).

## New or changed in 5.0.0

The powder model itself is unchanged since 1.0.0. What changed is its
input, the hkl list (chapter 12): exact bit-identical d-spacings, the
Laue-class fix for seven trigonal space groups, all species included in
$|F|^2$, the structure-factor symmetry self-check, correct lattice matrices
for triclinic and rhombohedral cells (issue #383), and space-group settings
and @ATOMSITES input (chapters 10-11). The filter tables of chapter 61
tabulate the powder Bragg edges exactly.

## References

- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, sections 2.2, 3.2 and
  4.
- T. Kittelmann, X.-X. Cai, "Elastic neutron scattering models for
  NCrystal", CPC 267 (2021) 108082, sections 2-3 and 7.
- T. Kittelmann, D. D. DiJulio, S. Xu, J. I. Marquez Damian, "Revisiting
  Becker-Coppens (1974)", Acta Cryst. A82 (2026) (extinction recipes, with
  NCrystal-computed Be cross sections).
