# Bragg diffraction in layered crystals (pyrolytic graphite)

## Summary

Pyrolytic graphite and similar layered crystals have crystallites whose
$c$ axes are aligned (with a mosaic spread) but which are rotated at random
around that axis. NCrystal's `LCBragg` model takes the Gaussian mosaic model
of chapter 21 and averages it analytically over the rotation angle $\phi$
about the layer axis (`lcaxis`): plane families are grouped by $(d, \theta_n)$
with $\theta_n$ the angle between the normal and the axis, the $\phi$
intervals where a family can reflect are located geometrically, and the
cross section is integrated over them with an adaptive Simpson scheme using
the single-crystal closed form inside; sampling draws $\phi$ from a
piecewise-constant overlay and then proceeds as in a single crystal. Two
slow reference models (averages over explicit rotations) are available
through `lcmode` for validation. Documented in the 2021 elastic paper,
section 5.

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

- `lcaxis=c1,c2,c3` (direct-lattice coordinates, e.g. `0,0,1` for
  graphite; the stdlib `C_sg194_pyrolytic_graphite.ncmat` embeds
  `NCRYSTALMATCFG[lcaxis=0,0,1]` so that the axis is set whenever an
  orientation is given), together with `dir1`, `dir2` and `mos` as in
  chapter 21, and `mosprec`, `sccutoff`, `dirtol`.
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

## Implementation

Normals with the same $d$ and $\theta_n$ are merged into plane sets with
their $|F|^2$ summed and sorted by $d$. All evaluations take place in a
frame with the layer axis along $z$, where a crystallite rotation is a
rotation of the normal about $z$ by $\phi$; for each plane set the $\phi$
intervals with non-zero cross section are located from the inequalities
above (with the degenerate cases of neutron or normal on the axis), the
single-crystal cross section is integrated over each interval with an
adaptive Simpson scheme using the closed form of chapter 21 inside, and
cumulative cross sections and per-interval overlays are kept for sampling.
Sampling picks an interval by cross section, draws $\phi$ by rejection
from a piecewise-constant overlay (acceptance 30-40%), and then samples the
single-crystal outcome at that rotation. The model is about ten times
slower than the single-crystal model (elastic paper section 7).

## Validation

- Elastic paper section 5: `LCHelper` against `LCBraggRef` with very large
  $N$ for PG at several wavelengths and orientations; the structure of the
  cross sections analysed analytically.
- Repository tests: `tests/src/app_lc` (the three modes against each
  other), `tests/scripts/pgbe.py` (issue #171, PG and Be geometries),
  `mmc_scpg.py` (PG monochromator in MiniMC, chapter 61).

## Approximations and limitations

- The same assumptions as chapter 21 (Gaussian mosaic, no extinction, no
  d-spacing spread), plus a uniform rotation distribution about the axis
  (no in-plane texture) and a single layer axis.
- The reference modes are slow and `lcmode<0` is not thread-safe.

## New or changed in 5.0.0

No physics changes; indirect changes through the hkl list and settings as
for chapter 21.

## References

- T. Kittelmann, X.-X. Cai, CPC 267 (2021) 108082, section 5 (eqs. 43-50,
  figures 19-22, reference models).
