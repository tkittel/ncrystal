# Reflection lists, d-spacings and structure factors

## Summary

From the unit cell, the symmetry and the atomic displacements NCrystal
generates, at load time, the list of reflection families $hkl$ with
$d_{hkl}$ in $[d_{cut}, d_{cut,up}]$, each with its multiplicity and
squared structure factor $|F_{hkl}|^2$ (Debye-Waller damped), which is the
input of every coherent elastic model (chapters 20-22). With a space group
the families are exactly the symmetry-equivalent groups; without one,
planes are grouped by equal $(d, |F|^2)$. d-spacings are computed so that
planes equal by symmetry get bit-identical values on every platform, the
lists are sorted strictly by $(d, \text{multiplicity}, hkl)$, and a
symmetry self-check on low-index planes guards against wrong space groups.
Lists are generated lazily and can be requested in part (e.g. only the
Bragg threshold).

Code: `ncrystal_core/src/extd_utils/NCFillHKL.cc`, `KSqCalc` in
`ncrystal_core/include/NCrystal/internal/utils/NCLatticeUtils.hh`,
`NCEqRefl` (chapter 11), the lazy hkl machinery of `Info`.

## Physics

### Structure factors

For the reciprocal-lattice vector $\vec\tau_{hkl}$ of length
$|k| = 2\pi/d_{hkl}$,

$$|F_{hkl}|^2 = \Big|\sum_{\text{species } s} b_{coh,s}\,e^{-W_s}\sum_{\text{positions } p\in s} e^{2\pi i(hx_p+ky_p+lz_p)}\Big|^2, \qquad W_s = \tfrac12|k|^2\langle u_s^2\rangle$$

(CPC 2020 eq. 25; Sears, Acta Cryst. A53 (1997) 35 for the exponent),
with the coherent scattering length $b_{coh}$ of the species (an element,
isotope or mixture, chapter 43) and the isotropic displacement of chapter
35. The multiplicity $m_{hkl}$ counts all equivalent planes including the
Friedel pairs, so that $\sum_{hkl} d|F|^2 m$ over families is the sum over
all reflections; the hkl list stores one representative per family and
expands it on request.

### d-spacings

$1/d^2$ is the textbook quadratic form in $h, k, l$ with the reciprocal
metric of the cell. It is evaluated with the integer parts formed exactly
and exact trigonometry for 90 and 120 degree angles, so that only the final
floating-point operations round and planes with equal $d$ by the cell
symmetry (cubic 333 and 511, hexagonal 700 and 530) get bit-identical
values on all platforms.

### Cut-offs

- `dcutoff` (default automatic: 0.1 Aa, or 0.2 Aa for cells with more than
  40 atoms): planes below it are not generated; the number of planes grows
  as $d_{cut}^{-3}$, but the initialisation time tends to a constant thanks
  to an early rejection of weak planes (below). `dcutoffup` removes long-d
  planes (useful to study the lower edges in isolation).
- Reflections with $|F|^2$ below 1e-5 barn (forbidden by symmetry, or
  negligible) are dropped.

## Inputs and configuration

- The validated structure (chapter 10) with a displacement for every atom.
- `dcutoff`, `dcutoffup`.
- Alternative sources of lists: `.laz`/`.lau` files (structure factors and
  multiplicities read and validated against the space group, chapter 52),
  and `Info` objects built by plugins with explicit lists or normals.

## Implementation

- The loop runs over half of the hkl box bounded by $|h| \le a/d_{cut}$
  etc. (since $hkl$ and $-h-k-l$ always share a family), with a cheap
  $|k|^2$ preselection before any expensive work.
- With a space group, each point is mapped to its family representative
  (chapter 11) and a bit table sized to the hkl range marks representatives
  already handled; $|F|^2$ is computed once per family. Before the loop, a
  self-check computes $|F|^2$ for all members of the low-index families and
  rejects the structure if members differ by more than 1%: a structure
  whose positions do not have the symmetry of the declared space group is
  thus caught instead of silently getting wrong Bragg physics (the
  position-based check of chapter 10 is the first line of defence; this one
  is independent of it).
- Without a space group, accepted points are sorted by $d$ and merged into
  families of equal $(d, |F|^2)$ within a small relative tolerance set by
  the rounding noise between truly equivalent planes; families take the
  average values of their members.
- A cheap upper bound on $|F|^2$ (all contributions in phase) rejects weak
  planes before any trigonometry, which keeps the load time bounded as
  $d_{cut} \to 0$ (CPC 2020 section 3.1). All species are always included
  in $|F|^2$ (skipping heavily damped ones, as before 5.0.0, biased weak
  planes by up to a few percent).
- Debye-Waller factors, phases and $|k|^2$ are evaluated with the
  platform-independent arithmetic of chapter 70, so lists are bit-identical
  across platforms and compilers.
- On the `Info` object the list is generated lazily and partially: a
  d-spacing range can be requested, and the Bragg threshold (the largest
  d-spacing) is available without building the full list. Entries are
  sorted by $d$ (descending), then multiplicity, then $hkl$. Any list type
  (symmetry families, explicit index lists, explicit normals) is turned
  into the stream of plane normals that the single-crystal models need.
- Python: `Info.hklList()` (`all_indices=True` expands families),
  `Info.hklObjects()`, `Info.braggthreshold`, `Info.dspacingFromHKL(h,k,l)`;
  `nctool --dump` prints the list; `ncmat2hkl` exports it (chapter 61).
- Lists are generated in milliseconds to tens of milliseconds per stdlib
  material at the default cut-off.

## Validation

- `tests/src/app_fillhkl` (thresholds, partial lists), `tests/scripts/fillhkl.py`,
  `hklobj.py`, `nosg.py` (no space group), `strain.py`; the `nctool` dumps of
  the stdlib checked in `.log` files reproduce every family.
- `app_symchecks` and `app_eqrefltables` for the symmetry side (chapter
  11); `latticecells.py` for d-spacings against metric-tensor calculations
  and spglib.
- CPC 2020 paper: lists cross-checked against nxslib/SgInfo for all files;
  powder-pattern refinements (GSAS-II, FullProf) recovering structures;
  measured $|F|^2$ of a monoclinic structure reproduced.

## Approximations and limitations

- Nuclear coherent scattering lengths only: no anomalous scattering, no
  magnetic structure factors, no occupancies or anisotropic displacements
  in $F$ (chapter 10).
- Isotropic Debye-Waller factor per species.
- The no-space-group path groups by numerical coincidence of $(d, |F|^2)$,
  which may merge accidentally degenerate families (harmless for powders).
- Reflection conditions are found numerically (the $|F|^2$ cut), not from
  the symmetry.

## New or changed in 5.0.0

NCrystal 4.4.6 already generated symmetry-correct families for all 230
groups (default settings), with early rejection of weak planes, and the
lists were validated against independent codes and refinements. New in
5.0.0: exact, bit-identical d-spacings and lists across platforms,
equivalent reflections derived from the operations of the actual setting
(fixing Laue class $\bar31m$ and making non-default settings correct), the
structure-factor symmetry self-check, all species included in $|F|^2$, a
more robust family grouping without space group, lower memory use, and the
corrected lattice matrix of triclinic and rhombohedral-axes cells (issue
#383).

## References

- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, sections 2.2 and 3.1
  (structure factors, cut-offs, early rejection, validation).
- V. F. Sears, Acta Cryst. A53 (1997) 35 (the Debye-Waller exponent).
- International Tables Vol. B Table A1.4.2.7 (settings; chapter 11).
