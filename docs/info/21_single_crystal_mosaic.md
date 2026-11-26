# Bragg diffraction in single crystals with Gaussian mosaicity

## Summary

Given an orientation (`dir1`, `dir2`) and a mosaic spread (`mos`), NCrystal
models a single crystal as an ensemble of crystallites whose plane normals
follow a truncated, properly normalised two-dimensional Gaussian on the unit
sphere around the nominal normal. The cross section of each reflection is
the integral of that density along the "Bragg circle" of normals satisfying
the Bragg condition for the incoming direction, evaluated per plane by a
closed-form approximation derived for NCrystal where that is accurate, and
otherwise by Gauss-Legendre quadrature or a validated table of such
integrals. Scattering picks a plane and then an actual crystallite
normal on the circle by rejection, and reflects the neutron specularly. The
model is exact for its assumptions (no interpolation grids in energy or
angle), anisotropic in the sense that cross sections depend on the neutron
direction, and fast enough for ray tracing (1-10 MHz above 1.5 Aa). Weak
short-d planes are delegated to the powder model (`sccutoff`). This is the
model published in the 2021 elastic paper (section 4). The model is set up
in the frame of the crystal and shared between all orientations of the same
crystal.

Code: `ncrystal_core/src/scbragg/NCSCBragg.cc`,
`ncrystal_core/src/phys_utils/NCGaussMos.cc`, `NCGaussOnSphere.cc`,
`ncrystal_core/include/NCrystal/internal/extd_utils/NCOrientUtils.hh`,
`NCPlaneProvider.hh`, `ncrystal_core/src/utils/NCLatticeUtils.cc`
(crystal directions and orientations) and `NCGaussLegendre.hh`, with cfg
parsing in `ncrystal_core/src/cfgutils/`, all assembled in
`NCStdScatFact.cc`.

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
$\tau < \pi/2$, which limits $\sigma$ to about 14 degrees (a FWHM of
about 34 degrees) at the tightest precision.

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
It is implemented as
$g = \sqrt{u\sin\alpha/\sin\gamma}\,F(\cos\delta_0)\,(1+u\,\mathrm{E}[v^2]/(4\sin\alpha\sin\gamma))$,
with $u = \cos\delta_0-\cos\tau$ and a spline $F$ (rather than the spline
$S = \sqrt{u}F$ of eq. 36, which has a square root edge at $u=0$). Unlike
eq. 35, $F$ uses the exact density, and with
$s = \sqrt{2\sin\alpha\sin\gamma}\,\sin(t/2) = \sqrt{u}\,v$ along the
circle the only approximation is
$1/\cos(t/2) \approx 1+s^2/(4\sin\alpha\sin\gamma)$, with
$\mathrm{E}[v^2]$ (the mean of $v^2$ with weight $W$) from a rational
approximation. This removes relative errors of order $\sigma^2$ (in radians),
and an underestimate of up to about $\epsilon/2$. The closed
form is used when $\tau < \psi$ and
$\cos\psi\sin\gamma\sin\alpha + \cos\alpha\cos\gamma < \cos\tau$ (eqs.
37-38, with $\psi \approx 18.7$ degrees at $\epsilon = 10^{-3}$), and when
in addition a first-order estimate of its relative error is below
$\epsilon/2$. The closed form is never used where rounding of its input
$\cos\delta_0$ (relative errors of about $4\varepsilon_{mach}/\sigma^2$)
would exceed $\epsilon/2$, i.e. for $\sigma$ below 0.3 arcseconds at
$\epsilon = 10^{-3}$ or 30 arcseconds at $\epsilon = 10^{-7}$.

Elsewhere, circle integrals are evaluated numerically to an accuracy of
$\epsilon/10$ (clamped to $[10^{-7}, 10^{-4}]$), by Gauss-Legendre
quadrature of the exact density along the circle, with the number of
points (8-24, on segments of at most $10\sigma$) chosen from a table of
measured maximum errors, so that these are at least 10 times below that
accuracy. Where the closed form fails at the peak of the rocking curve for
$\alpha = 1.2\tau$ (large mosaicities or small `mosprec`), a 2D table of
such integrals is set up and used for $\alpha$ from $1.2\tau$ to
$\pi/2$. The table is indexed by $\sin(\gamma-\alpha)$ and
$\tan(\alpha/2)$ (no transcendental functions needed), and holds
$\ln(g/(\sqrt{1-z^2}\sqrt{\sin\alpha/\sin\gamma}))$ with
$z=\sin(\gamma-\alpha)/\sin\tau$, which is smooth. Values are
interpolated with 4x4-point Lagrange polynomials. At setup the table
(64, 128 or 256 points in $\alpha$, four times as many in $z$) is
validated against numerical integration and refined until its maximum error
is 10 times below the numerical accuracy (for values not negligible), or
else not used.

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
  normal of a plane (`@crys_hkl:h,k,l`, via the reciprocal lattice), or a
  lattice direction $u\vec a+v\vec b+w\vec c$ (`@crys_uvw:u,v,w`, via the
  direct lattice). The components of `dir2` parallel to `dir1` are ignored.
- `dir2=arbitrary` when the rotation of the crystal around `dir1` is of no
  importance: a deterministic `dir2` orthogonal to `dir1` (in both frames)
  is then chosen, which needs the lattice. The `dir2` parameter can also be
  left out for layered crystals with `lcaxis` along `dir1` (chapter 22).
- `rot=a1,a2,a3:angle` (angle in rad, deg, arcmin or arcsec): a right-handed
  rotation of the material around an axis in the lab frame, applied after
  the orientation from `dir1` and `dir2` (e.g. `rot=0,1,0:30deg`). The `rot`
  parameter has no effect on materials which are not single crystals, so it
  can be added to any cfg-string.
- A single crystal needs `mos` and `dir1`, and `dir2` unless `lcaxis` is
  along `dir1`. Otherwise the error is
  `Single crystals need mos and dir1, and also dir2 unless lcaxis is along dir1`.
- `mos` (FWHM, in (0, pi/2]), `mosprec` (default 1e-3, range [1e-7, 0.1]).
- `sccutoff` (default 0.4 Aa; 0 disables): planes with $d$ below it are
  treated as having infinite mosaicity, i.e. as powder (chapter 20), with
  a one to two orders of magnitude speed-up at short wavelengths.
- `dcutoff`, `dcutoffup`, `coh_elas`, `temp`, `strain` as in chapter 20.
- `lcaxis` selects the layered model of chapter 22 instead.
- C++: `MatCfg::get_dir1()`, `get_dir2()` and `get_rot()` return what the
  cfg-string specifies (`NullOpt` if not set, null vectors for
  `dir2=arbitrary`) and never need the lattice. `ScatterRequest` resolves
  the orientation, which needs a material with a crystal structure
  (`MissingInfo` otherwise): `get_dir2()` resolves an arbitrary or left out
  `dir2` (before `rot`), and `createSCOrientation()` includes `rot`.
  `ScatterRequest::isSingleCrystal()` is only true for single phase
  crystalline materials, so orientation parameters are ignored for
  non-crystalline phases of multiphase materials.
- JSON: `dir1` and `dir2` are
  `{"basis":"cartesian"|"hkl"|"uvw","crystal":[..],"lab":[..]}` (or
  `"arbitrary"`), and `rot` is `{"axis":[..],"angle":rad}`.

## Implementation

The model is set up in the frame of the crystal (`SCBragg::createModel`)
and shared between all materials which only differ in the orientation of
the crystal (`dir1`, `dir2`, `rot`): each `SCBragg` process holds the shared
model and its own crystal-to-lab rotation (which must be a proper rotation),
with which neutron directions are rotated to and from the crystal frame.
The factory caches models (and the powder process for planes below
`sccutoff`) by material, `mos`, `mosprec` and `sccutoff`. Many differently
oriented crystals of the same material thus cost little memory and
initialisation time. The JSON description of an `SCBragg` process holds
`model_uid` (the same for processes sharing a model) and `crystal2lab`.

The planes of the material are grouped into families of equal $(d, |F|^2)$
with their normals, sorted by $d$, and
processed only while $\lambda \le 2d$ (planes for which rounding would
give $\sin\theta_B > 1$ just below the threshold are excluded). For a
given neutron state (cached per thread), the few planes able to contribute
are first preselected: a normal $\hat n$ only contributes if
$||\hat n\cdot\hat k_i|-\lambda/2d|$ is at most $2\sin(\theta_t/2)$, with
$\theta_t$ the truncation angle, which is tested in single precision (with
a safety margin, so results are unaffected) in chunks of 64 planes (16 for
the tail), in a loop written to be vectorised by compilers. For the
surviving planes, an exact test on normal and anti-normal together decides
whether the Bragg circle meets the truncated Gaussian at all, the geometric
factor is then obtained from the closed-form approximation where it is
accurate and from the table or numerical integration otherwise,
and the per-plane cross sections are accumulated for later plane
selection. The cache key is the neutron energy and the direction
normalised and snapped to a grid of $2^{-47}$, compared bitwise, so results
depend only on the input and not on earlier calls. A neutron which is not
deflected keeps its exact input direction. The Gaussian density and the
closed form are served from splined lookup tables sized by `mosprec`, so
that no transcendental function is evaluated per plane where the closed
form or the table of circle integrals applies. Points on circles are
sampled exactly by rejection. Throughput was 1-10 MHz above 1.5 Aa and
0.1-1 MHz below, where many planes contribute (elastic paper section 7,
before the plane preselection), with `sccutoff` making the difference at
short wavelengths. A d-spacing
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
- Repository tests (C++ in `tests/src/`, Python in `tests/scripts/`):
  - `app_backscat`: back-scattering geometry, including neutrons exactly
    along plane normals.
  - `app_gaussmosplanes`: the plane preselection against a loop over all
    planes (also for non-default `mosprec` and at the Bragg threshold).
  - `app_scbraggcache`: results independent of the call history, and
    undeflected neutrons keeping their input direction.
  - `app_goscircleint`: closed form, table and numerical circle integrals
    against `mpmath` references, and the expected table sizes.
  - `app_math`: the Gauss-Legendre rules of `NCGaussLegendre.hh`.
  - `app_rotscat`, `scshare.py`: models shared between orientations.
  - `app_orientutils`: lattice directions, the zone law, `dir2=arbitrary`,
    parallel crystal directions, `SCOrientation` output.
  - `scorient.py`: equivalent ways of specifying an orientation,
    including `rot` and `dir2=arbitrary`.
  - `app_cfg`, `app_cfgparse`, `app_cfgmp`: parsing, normalisation,
    getters, setters and error messages of the cfg parameters.
  - `azimuth.py` (azimuthal distributions), `mmc_scge.py` (single-crystal
    Ge in MiniMC, chapter 61), and the `nctool` cross-section dumps of
    oriented materials.
- The `mpmath` reference values (independent of NCrystal) are generated by
  `devel/misc/gen_mosaicrefs.py` into `tests/data/mosaicrefs/`. The
  developer benchmark `devel/simplebuild/static/NCDev/app_benchscbragg`
  measures `SCBragg` speed.

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

The model is as published for 2.0.0-4.4.6, but is now set up in the
crystal frame and shared between orientations (see above). New ways of
specifying orientations: `@crys_uvw:`, `dir2=arbitrary` and `rot`.
Bug fixes: planes with a normal exactly along the neutron were ignored near
back-scattering (`GaussMos::calcCrossSections`), and the JSON description
of `SCBragg` gave the mosaicity under the key `dmax` instead of `mos`.
Results no longer depend on the call history of a cache.
Numerical integrations of circles use Gauss-Legendre quadrature of the exact
density rather than Romberg integration of a splined one (also without the
reduced accuracy formerly used for very small mosaicities), and are mostly
replaced by a table where the closed form is often invalid: 5-9 times
faster in such cases, and more precise. The closed form uses the exact
density, a first order correction and an error estimate, which changes some
cross sections at large mosaicities by up to about $\sigma^2$. Cross sections are evaluated faster
due to the preselection of planes (results are unchanged): typically 2-3
times for a mosaicity of 0.5 degrees (7 times for crystals of low
symmetry), and more for smaller mosaicities.
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
