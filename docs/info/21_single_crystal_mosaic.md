# Bragg diffraction in single crystals with Gaussian mosaicity

## Summary

Given an orientation (`dir1`, `dir2`) and a mosaic spread (`mos`), NCrystal
models a single crystal as an ensemble of crystallites whose plane normals
follow a truncated, properly normalised two-dimensional Gaussian on the unit
sphere around the nominal normal. The cross section of each reflection is
the integral of that density along the "Bragg circle" of normals satisfying
the Bragg condition for the incoming direction, evaluated per plane by
Gauss-Legendre quadrature or, where accurate, by a closed-form approximation
derived for NCrystal; scattering picks a plane and then an actual crystallite
normal on the circle by rejection, and reflects the neutron specularly. The
model is exact for its assumptions (no interpolation grids in energy or
angle), anisotropic in the sense that cross sections depend on the neutron
direction, and fast enough for ray tracing (1-10 MHz above 1.5 Aa). Weak
short-d planes are delegated to the powder model (`sccutoff`). This is the
model published in the 2021 elastic paper (section 4).

Code: `ncrystal_core/src/scbragg/NCSCBragg.cc`,
`ncrystal_core/src/phys_utils/NCGaussMos.cc`, `NCGaussOnSphere.cc`,
`ncrystal_core/include/NCrystal/internal/extd_utils/NCOrientUtils.hh`,
`NCPlaneProvider.hh`; assembled in `NCStdScatFact.cc`.

## Physics

### Mosaic model

The orientation density of crystallite normals around the nominal normal
$\hat n$ is

$$W(\delta) = N\,e^{-\delta^2/2\sigma^2}\,\Theta(\tau - \delta),$$

where $\delta$ is the angular distance on the sphere, $\sigma$ the mosaic
spread (the user gives the FWHM, `mos` $= 2\sqrt{2\ln 2}\,\sigma$), $\tau$
the truncation angle and $N$ the normalisation over the whole sphere,
$N^{-1} = 2\pi\int_0^\tau e^{-\theta^2/2\sigma^2}\sin\theta\,d\theta$
(evaluated numerically to 1e-12). The truncation is set from the precision
parameter $\epsilon$ (`mosprec`) as $\tau = \max(3, 1.1\sqrt{-2\ln\epsilon})\,\sigma$
(4.09 $\sigma$ at the default $\epsilon = 10^{-3}$), and must satisfy
$\tau < \pi/2$, which limits the mosaicity to about 14 degrees at the
tightest precision.

### Cross section (derivation)

For a plane with spacing $d$ and the neutron direction $\hat k_i$ at
wavelength $\lambda \le 2d$, define the Bragg angle $\sin\theta_B = \lambda/2d$
and $\alpha = \pi/2 - \theta_B$, the angle a normal must make with $-\hat k_i$
for reflection. The reflecting normals form the Bragg circle of angular
radius $\alpha$ around $-\hat k_i$. With $\gamma$ the angle between the
nominal normal and $-\hat k_i$, the cross section per atom is
(elastic paper eqs. 10-12, 24)

$$\sigma_{hkl} = q_{hkl}\,g_{hkl}, \qquad q_{hkl} = \frac{\lambda^3|F_{hkl}|^2}{V_{uc}n_a\sin 2\theta_B}, \qquad
g_{hkl} = 2N\sin\alpha\int_0^{t'}\exp\Big(-\frac{\arccos^2(\sin\alpha\sin\gamma\cos t + \cos\alpha\cos\gamma)}{2\sigma^2}\Big)dt,$$

with $t' = \arccos\max(-1, (\cos\tau - \cos\alpha\cos\gamma)/(\sin\alpha\sin\gamma))$
the half-length of the arc of the circle inside the truncation cap
($g = 0$ when the circle and the cap do not meet, the branch-free test
$(1-\cos^2\alpha)(1-\cos^2\gamma) > \max(0,\cos\tau - \cos\alpha\cos\gamma)^2$).
Both the normal and the anti-normal of a plane are considered. The
powder limit is recovered with $W = 1/4\pi$.

When the circle is long compared with $\sigma$ and far from the
truncation edge, the integral has the closed-form approximation derived in
the elastic paper (eq. 35),

$$g \approx \frac{N}{1/(2\pi\sigma^2)}\,\frac{e^{-\delta_0^2/2\sigma^2}}{\sqrt{2\pi\sigma^2}}\,\operatorname{erf}\sqrt{\frac{\tau^2-\delta_0^2}{2\sigma^2}}\,\sqrt{\frac{\sin\alpha}{\sin\gamma}}, \qquad \delta_0 = |\alpha-\gamma|,$$

a Gaussian rocking curve in $\delta_0$ times a truncation factor and a
curvature factor $\sqrt{\sin\alpha/\sin\gamma}$ that makes the rocking
curves asymmetric (to lowest order the $1 - \eta\cot\alpha/2$ of Wuttke).
It is implemented as $g = S(\cos(\alpha-\gamma))\sqrt{\sin\alpha/\sin\gamma}$
with a spline $S$ (eq. 36), used when $\tau < \psi$ and
$\cos\psi\sin\gamma\sin\alpha + \cos\alpha\cos\gamma < \cos\tau$ (with
$\psi \approx 22.6$ degrees at $\epsilon = 10^{-3}$), and numerical
integration otherwise. The latter uses Gauss-Legendre quadrature of the exact
density along the circle, with the number of points (8-24, on segments of
at most $10\sigma$) chosen from a table of measured maximum errors, so that
these are at least 10 times below the requested accuracy. When the closed
form is never valid (large mosaicities or small `mosprec`), a 2D table of
such integrals is set up instead, covering $\alpha$ from 1.2 times the
truncation angle to $\pi/2$. It is indexed by $\sin(\gamma-\alpha)$ and
$\tan(\alpha/2)$ (no transcendental functions needed), and holds
$\ln(g/(\sqrt{1-z^2}\sqrt{\sin\alpha/\sin\gamma}))$ with
$z=\sin(\gamma-\alpha)/\sin\tau$, which is smooth. It is interpolated with
4x4-point Lagrange polynomials, and is validated and refined at setup until
its maximum error is 10 times below the requested accuracy (for values not
negligible), or else not used.

### Scattering

A plane (and sign of normal) is chosen with probability proportional to
its cross section, a point $t \in [-t', t']$ on its Bragg circle is sampled
by rejection from a constant overlay (acceptance 15-35%), which selects a
crystallite, and the neutron is reflected specularly in that crystallite:
$\vec k_f = \vec k_i + \vec\tau\hat n$ (energy unchanged, exact
double-angle formulas, no trigonometric calls, then rotation back to the
lab frame).

## Inputs and configuration

- `dir1=@crys:c1,c2,c3@lab:l1,l2,l3` and `dir2=...`, with `dirtol`
  (default 1e-4 rad) bounding the allowed mismatch between the opening
  angles of the two directions in the two frames. The crystal direction is
  a Cartesian vector in the crystal frame of chapter 01 (`@crys:`), the
  normal of a plane (`@crys_hkl:h,k,l`), or a lattice direction
  $u\vec a+v\vec b+w\vec c$ (`@crys_uvw:u,v,w`).
- `dir2=arbitrary` when the rotation of the crystal around `dir1` is of no
  importance (an unspecified but reproducible rotation is then used).
- `rot=a1,a2,a3:angle`: a right-handed rotation of the crystal around an
  axis in the lab frame, applied after `dir1` and `dir2` (e.g.
  `rot=0,1,0:30deg`).
- `mos` (FWHM, in (0, pi/2]), `mosprec` (default 1e-3, range [1e-7, 0.1]).
- `sccutoff` (default 0.4 Aa; 0 disables): planes with $d$ below it are
  treated as having infinite mosaicity, i.e. as powder (chapter 20), with
  a one to two orders of magnitude speed-up at short wavelengths.
- `dcutoff`, `dcutoffup`, `coh_elas`, `temp`, `strain` as in chapter 20.
- `lcaxis` selects the layered model of chapter 22 instead.

## Implementation

The model is set up in the frame of the crystal and shared between all
materials which only differ in the orientation of the crystal (`dir1`,
`dir2`, `rot`): each of these merely wraps it in a `RotatedScatter` process,
which rotates neutron directions to and from the crystal frame. Many
differently oriented crystals of the same material thus cost little memory
and initialisation time.

The planes of the material are grouped into families of equal $(d, |F|^2)$
with their normals, sorted by $d$, and
processed only while $\lambda \le 2d$. For a given neutron state (cached
per thread), the few planes able to contribute are first preselected: a
normal $\hat n$ only contributes if $||\hat n\cdot\hat k_i|-\lambda/2d|$
is at most $2\sin(\theta_t/2)$, with $\theta_t$ the truncation angle,
which is tested in single precision (with a safety margin) for up to 64
planes at a time, in a loop written to be vectorised by compilers. For the
surviving planes, an exact test on normal and anti-normal together decides
whether the Bragg circle meets the truncated Gaussian at all, the geometric
factor is then obtained from the closed-form approximation where it is
accurate and from numerical integration otherwise,
and the per-plane cross sections are accumulated for later plane
selection. The Gaussian density, its circle integrals and the closed form
are served from splined lookup tables sized by `mosprec`, so that no
transcendental function is evaluated per plane; points on circles are
sampled exactly by rejection. Throughput is 1-10 MHz above 1.5 Aa and
0.1-1 MHz below, where many planes contribute (elastic paper section 7),
with `sccutoff` making the difference at short wavelengths. A d-spacing
spread parameter exists in the interfaces but is not implemented.

## Validation

- Elastic paper: rocking curves of eq. 24 against 100-digit `mpmath`
  references for $\sigma$ of 1 degree and 1 arcminute including
  back-scattering; sampled azimuthal distributions against the density;
  a 1e10-step "zig-zag" walk for consistency between cross sections and
  sampling; isotropic illumination reproducing the powder cross section to
  1e-3; reflectivity of Ge-511 slabs against Sears' Darwin-Hamilton
  solution (in plane) and Wuttke's Monte Carlo (out of plane), where the
  curvature factor matters.
- Repository tests: `tests/src/app_backscat` (back-scattering geometry),
  `tests/scripts/azimuth.py` (azimuthal distributions), `mmc_scge.py`
  (single-crystal Ge in MiniMC, chapter 61), the `nctool` cross-section
  dumps of oriented materials, the orientation parsing tests in
  `app_cfg`/`app_cfgparse`, `tests/src/app_rotscat` and
  `tests/scripts/scshare.py` (models in the crystal frame, shared between
  orientations), `tests/src/app_orientutils` and `tests/scripts/scorient.py`
  (equivalent ways of specifying an orientation).

## Approximations and limitations

- Isotropic Gaussian mosaicity only (no anisotropic mosaic, no
  Lorentzian or other shapes); no d-spacing distribution; $\tau < \pi/2$
  limits the usable mosaicity.
- Ideally imperfect crystal: Born approximation, no primary or secondary
  extinction, no dynamical diffraction (the Darwin-Hamilton transfer
  problem is left to the transport code, which handles multiple
  scattering).
- Spline artefacts of order `mosprec` near the truncation edges.

## New or changed in 5.0.0

The model is unchanged and as published for 2.0.0-4.4.6, but is now set up
in the crystal frame and shared between orientations (see above). New ways
of specifying orientations: `@crys_uvw:`, `dir2=arbitrary` and `rot`.
Bug fix: planes with a normal exactly along the neutron were ignored near
back-scattering (`GaussMos::calcCrossSections`).
Numerical integrations of circles use Gauss-Legendre quadrature of the exact
density rather than Romberg integration of a splined one, and are mostly
replaced by a table when the closed form is never valid: 5-9 times faster
in such cases, and more precise. Cross sections are evaluated faster due to the
preselection of planes (results are unchanged): typically 2-3 times for a mosaicity of 0.5 degrees
(7 times for crystals of low symmetry), and more for smaller mosaicities.
Indirect changes: the hkl input (chapter 12), orientation handling for
non-default space-group settings (chapter 11), and `.lau` files with setting
information (chapter 52).

## References

- T. Kittelmann, X.-X. Cai, CPC 267 (2021) 108082, section 4 (model,
  eqs. 24-36, validation).
- The Sears (Darwin-Hamilton in-plane reflectivity) and Wuttke
  (out-of-plane Monte Carlo, asymmetric rocking curves) references as
  cited in section 4 of the elastic paper.
- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, section 3 (cfg
  parameters `mos`, `dir1`, `dir2`, `dirtol`).
