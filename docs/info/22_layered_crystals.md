# Bragg diffraction in layered crystals (pyrolytic graphite)

## Summary

Pyrolytic graphite and similar layered crystals have crystallites whose
$c$ axes are aligned (with a mosaic spread) but which are rotated at random
around that axis. NCrystal's `LCBragg` model takes the Gaussian mosaic model
of chapter 21 and averages it analytically over the rotation angle $\phi$
about the layer axis (`lcaxis`): plane families are grouped by $(d, \theta_n)$
with $\theta_n$ the angle between the normal and the axis, the $\phi$
intervals where a family can reflect are located geometrically, and the
cross section is integrated over them, mostly via a series expansion with
exact moments of the mosaicity distribution, otherwise numerically. Sampling
draws $\phi$ from a piecewise-constant overlay and then proceeds as in a single crystal. Two
slow reference models (averages over explicit rotations) are available
through `lcmode` for validation. The same model describes fibre textures
(`dir1avg`) and single crystals on a spinning sample holder (`lcaxis` in
the lab frame). As for chapter 21, the model is set up in the frame of the
crystal and shared between orientations. Documented in the 2021 elastic
paper, section 5.

Code: `ncrystal_core/src/lcbragg/NCLCBragg.cc`,
`ncrystal_core/src/extd_utils/NCLCUtils.cc` (`LCHelper`, `LCPlaneSet`,
`LCROIFinder`, `LCStdFrame`), `NCLCRefModels.cc` (`LCBraggRef`,
`LCBraggRndmRot`).

## Physics

With $L$ the layer axis, $\theta_k$ the angle between $-\hat k_i$ and $L$
and $\theta_n$ the angle between a plane normal and $L$, a crystallite
rotated by $\phi$ about $L$ presents its normal at an angle $\gamma(\phi)$
to $-\hat k_i$ given by

$$\cos\gamma(\phi) = \sin\theta_k\sin\theta_n\cos\phi + \cos\theta_n\cos\theta_k,$$

and the geometric factor of the layered crystal is the average of the
single-crystal one (chapter 21) over $\phi$:

$$g_{LC}(\alpha,\hat k_i) = \frac{1}{2\pi}\int_{-\pi}^{\pi} g_{SC}\big(\alpha,\gamma(\phi)\big)\,d\phi,$$

with $\sigma_{hkl} = q_{hkl}\,g_{LC}$ as before (elastic paper eqs.
43-50). Planes whose normals are on the axis ($\theta_n = 0$) are
unaffected by the rotation and treated as in a single crystal
(the $00l$ planes of graphite, the ones a monochromator uses); planes
with $\theta_n$ and $\pi - \theta_n$ are processed jointly. The powder
limit is not reached: a layered crystal is anisotropic in $\theta_k$ even
with an isotropic in-plane average, and the structure of the cross section
(plateaus $\propto \lambda^2/\sqrt{1-(\lambda/2d)^2}$, triangular regions,
elastic paper figures 19-22) is distinctive.

The integrand is non-zero only where the Bragg circle meets the truncated
Gaussian, i.e. for $\phi$ in the interval(s) where

$$\frac{\cos(\alpha+\tau) - \cos\theta_k\cos\theta_n}{\sin\theta_k\sin\theta_n} \le \cos\phi \le \frac{\cos\max(0,\alpha-\tau) - \cos\theta_k\cos\theta_n}{\sin\theta_k\sin\theta_n}$$

(eqs. 49-50), after a cheap pre-check that
$\cos\theta_n \in [\cos(\theta_k+\alpha), \cos\max(0,\theta_k-\alpha)]$
allows any contribution.

## Inputs and configuration

- `lcaxis=c1,c2,c3` (a Cartesian vector in the crystal frame of chapter 01,
  also written `@crys:c1,c2,c3`, e.g. `0,0,1` for graphite, which the stdlib
  `C_sg194_pyrolytic_graphite.ncmat` embeds as `NCRYSTALMATCFG[lcaxis=0,0,1]`
  so that the axis is set whenever an orientation is given), together with
  `dir1`, `dir2` and `mos` as in chapter 21, and `mosprec`, `sccutoff`,
  `dirtol`, `rot`. The axis can also be the normal of a plane
  (`lcaxis=@crys_hkl:h,k,l`), a lattice direction (`@crys_uvw:u,v,w`), an
  axis in the lab frame (`@lab:l1,l2,l3`, e.g. for a single crystal on a
  spinning sample holder), or the crystal direction of `dir1` (`@dir1`,
  which needs `dir1`). A lab axis is converted to the crystal frame with the
  orientation before `rot`, so `rot` rotates the whole material including
  that axis.
- `dir2` can be left out when `lcaxis` is a direction in the crystal which
  is parallel to that of `dir1` (cos^2 of their angle within 1e-12 of 1,
  which can need the lattice), or `@dir1`, since the rotation around `dir1`
  is then of no importance. With the stdlib PG file,
  `mos=1deg;dir1=@crys_hkl:0,0,1@lab:0,1,0` is thus enough.
- Pseudo-parameter `dir1avg=X` (alias `fibre=X`): exactly
  `dir1=X;lcaxis=@dir1`, i.e. a crystal averaged over all rotations around
  a primary direction (a fibre texture).
- `lcmode` (default 0): 0 selects `LCHelper`; a positive $N$ selects
  `LCBraggRef`, the average of $N$ evenly rotated single crystals
  (accurate only for very large $N$: $2\pi/N$ must be about 5% of the
  FWHM, so $N \sim 7\times10^3$ at 1 degree and $2.6\times10^6$ at 10
  arcseconds); a negative $-N$ selects `LCBraggRndmRot`, which draws $N$
  random rotations at every call (not thread-safe, for reference only; the
  McStas `Single_crystal` "PG mode" corresponds to $N = 1$ and is shown
  inadequate in the elastic paper).
- Without an orientation, a layered polycrystal is treated as an ordinary
  powder (chapter 20).
- C++: `MatCfg::get_lcaxis()` returns the `LCAxisState` given (form and
  vector, `NullOpt` if not set), while `ScatterRequest::get_lcaxis()` always
  gives a vector in the crystal frame (and needs a crystal structure for
  all but the Cartesian form). JSON:
  `{"form":"cartesian"|"hkl"|"uvw"|"lab"|"primary","axis":[..]}` (no axis
  for `primary`).

## Implementation

Normals with the same $d$ and $\theta_n$ are merged into plane sets with
their $|F|^2$ summed and sorted by $d$. All evaluations take place in a
frame with the layer axis along $z$, where a crystallite rotation is a
rotation of the normal about $z$ by $\phi$; for each plane set the $\phi$
intervals with non-zero cross section are located from the inequalities
above (with the degenerate cases of neutron or normal on the axis), the
single-crystal cross section is integrated over each interval. With
$u = \cos\gamma = A\cos\phi + B$, the integral becomes
$\int g(u)\,(A^2-(u-B)^2)^{-1/2}\,du$, where $g$ is only non-zero in a narrow
window around $\lambda/2d$. When that window is well inside $[B-|A|, B+|A|]$
(the usual case), the square-root factor is expanded to 6th order, which
needs the moments of $g$. These are exact and simple functions of
$\lambda/2d$ and six averages over the mosaicity distribution, computed once.
The method is used only if a strict bound on the neglected terms is below
`mosprec`/10. Otherwise the trapezoidal rule is used, starting with 17
points and halving the step until two estimates agree to within `mosprec`
(clamped to $[10^{-7}, 10^{-2}]$), with the circle integrals of chapter 21
inside, evaluated in batches of 64. Convergence can be like $h^{1.5}$ (with
$h$ the step), so unless the last two differences show fast convergence,
0.55 times the last one is added (a Richardson step). Cumulative
cross sections and per-interval overlays are kept for sampling.
Sampling picks an interval by cross section, draws $\phi$ by rejection
from a piecewise-constant overlay (acceptance 30-40%), and then samples the
single-crystal outcome at that rotation. The model is about ten times
slower than the single-crystal model (elastic paper section 7).

As for `SCBragg`, an `LCBragg` process holds a model set up in the crystal
frame (`LCBragg::createModel`) and its own crystal-to-lab rotation. The
factory caches models by material, `mos`, `mosprec`, `sccutoff`, `lcaxis`
(as a unit vector in the crystal frame) and `lcmode`. The JSON description
holds `model_uid` and `crystal2lab`, and the plane set fields
(`nplanesets`, `dmin`, `dmax`) are null for `lcmode` other than 0.

## Validation

- Elastic paper section 5: `LCHelper` against `LCBraggRef` with very large
  $N$ for PG at several wavelengths and orientations; the structure of the
  cross sections analysed analytically.
- Repository tests (C++ in `tests/src/`, Python in `tests/scripts/`):
  `app_lc` (the three modes against each other), `app_lcinteg` (series
  and numerical integration against `mpmath` references from
  `devel/misc/gen_mosaicrefs.py`), `app_rotscat` and `scshare.py` (models
  shared between orientations, all `lcaxis` forms, `lcmode`, JSON),
  `scorient.py` (`lcaxis` forms, `dir2` left out, `dir1avg`/`fibre`),
  `app_cfg`/`app_cfgparse` (parsing and errors), `pgbe.py` (issue #171, PG
  and Be geometries), `mmc_scpg.py` (PG monochromator in MiniMC, chapter
  61).

## Approximations and limitations

- The same assumptions as chapter 21 (Gaussian mosaic, no extinction, no
  d-spacing spread), plus a uniform rotation distribution about the axis
  (no in-plane texture) and a single layer axis.
- The reference modes are slow and `lcmode<0` is not thread-safe.

## New or changed in 5.0.0

No physics changes, but as for chapter 21 the model is now set up in the
crystal frame and shared between orientations, and `LCBragg` processes
have a JSON description and summary string. New forms of `lcaxis`
(`@crys:`, `@crys_hkl:`, `@crys_uvw:`, `@lab:`, `@dir1`), the
pseudo-parameters `dir1avg` and `fibre`, and `dir2` is optional for an
`lcaxis` along `dir1`. The description of the plain form of `lcaxis` was
corrected (a Cartesian vector, not direct-lattice coordinates). Indirect
changes through the hkl list and settings, and through the circle
integrals, as for chapter 21. The integration over $\phi$ now mostly uses
the expansion described above, and otherwise the trapezoidal rule (rather
than Romberg integration) with batched evaluations and a final
extrapolation step. Cross sections are
evaluated 3-15 times faster, and change slightly (well within `mosprec`,
typically in the 5th significant digit).

## References

- T. Kittelmann, X.-X. Cai, CPC 267 (2021) 108082, section 5 (eqs. 43-50,
  figures 19-22, reference models).
