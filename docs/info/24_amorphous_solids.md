# Elastic scattering in amorphous solids

## Summary

Non-crystalline solids (glasses, polymers, powders described without a
unit cell) have no Bragg diffraction, but their atoms still vibrate around
fixed positions, so elastic scattering exists. NCrystal models it, since
2.7.0, with the incoherent-elastic machinery of chapter 23 applied to the
*total* bound cross section $\sigma_{coh} + \sigma_{inc}$ of each element
(the incoherent approximation for the coherent elastic part), using
displacements derived from the dynamics sections. The static structure
factor $S(Q)$ of the amorphous structure is not modelled in the core; a
plugin can take over the coherent part. Which materials are solids follows
from the NCMAT @STATEOFMATTER rules, mandatory in NCMAT v8 for
non-crystalline files.

Code: the amorphous-solid branch of
`ncrystal_core/src/stdscatfactory/NCStdScatFact.cc` and `ElIncScatter`
(chapter 23); state-of-matter resolution in
`ncrystal_core/src/infobld/NCInfoBldComposition.cc`.

## Physics

In the incoherent approximation the coherent elastic scattering of a
disordered solid is replaced by its angular average, i.e. by the same
Debye-Waller-damped isotropic form as incoherent scattering:

$$\frac{d\sigma}{d\Omega} = \sum_j \frac{w_j(\sigma_{coh,j}+\sigma_{inc,j})}{4\pi}\,e^{-Q^2\langle u_j^2\rangle}, \qquad
\sigma_{el}(E) = \sum_j w_j(\sigma_{coh,j}+\sigma_{inc,j})\,\frac{1-e^{-t_j}}{t_j},$$

with $t_j = 4k^2\langle u_j^2\rangle$ (chapter 23). This reproduces the
correct total elastic cross section at long wavelengths (where $S(Q) \to$
its average over the accessible $Q$ range) and the correct Debye-Waller
suppression at short wavelengths, but it ignores the oscillations of
$S(Q)$ (the "amorphous halo") and therefore the angular structure and the
small-$Q$ behaviour of the coherent part. The inelastic part of such
materials is modelled by the usual kernels (chapters 30-35), for which the
incoherent approximation is the standard choice anyway.

## Inputs and configuration

- A non-crystalline NCMAT file (no @CELL/@ATOMPOSITIONS) whose state is
  `solid`: automatically so when any @DYNINFO section is of type `vdos` or
  `vdosdebye`; otherwise an explicit `@STATEOFMATTER solid` (NCMAT v7
  optional, v8 mandatory for all non-crystalline files, which removes the
  "unknown" state). Liquids and gases get no elastic component.
- Displacements per element from the dynamics (chapter 35):
  `vdos` (from the curve), `vdosdebye` (`debye_temp` or v8 `msd`), `scatknl`
  with v8 `msd`/`debye_temp` or confident kernel analysis (v8 files and
  plugin data only); `freegas`/`sterile` elements and kernels that could
  not be analysed contribute nothing, with a warning naming them.
- `coh_elas` (default on) adds the coherent part, `incoh_elas` the
  incoherent part; both off removes the component.
- The quick material `solid::<formula>/<density>[/TDebye<T>K]` (e.g.
  `solid::Al2O3/3.9gcm3/TDebye500K`) builds an amorphous solid with
  Debye-model dynamics (default Debye temperature 300 K), giving a quick
  elastic plus inelastic model of any composition and density.
- A plugin can take over the coherent part with its own $S(Q)$-based model
  (an unofficial switch disables the incoherent approximation for it).

## Implementation

Non-crystalline materials whose state is solid get the incoherent elastic
process of chapter 23 built with the total bound cross section of each
element (coherent plus incoherent, each part switchable), provided at least
one element has a displacement. The state of matter is implied solid by a
crystal structure or VDOS-type dynamics and must be declared otherwise;
multiphase materials keep a state only if all phases agree.

## Validation

- Stdlib amorphous materials (polyethylene, acrylic glass, epoxy, the
  nylons, polystyrene, polycarbonate, Kapton, rubber, ...; 16 of them were
  added with 2.7.0) have their elastic contributions in the `nctool`
  cross-section dumps checked by `.log` tests; `dkelas.py` exercises
  direct-kernel solids; `ncmatmisc.py` and `ncmatparser.py` cover the
  @STATEOFMATTER rules; `quickfact.py` covers `solid::`.
- There is no benchmark of the incoherent approximation against measured
  $S(Q)$ of glasses in the repository; the approximation is the standard
  one of thermal scattering libraries and its limitations are known
  (below).

## Approximations and limitations

- No static structure factor: coherent elastic scattering is isotropic
  (Debye-Waller-damped), so diffraction halos and the strong small-$Q$
  suppression of coherent scattering in dense amorphous solids are absent.
  For hydrogenous polymers the incoherent part dominates and the error is
  small; for e.g. amorphous SiO$_2$ or metallic glasses it is not.
- The same displacement assumptions as chapters 23 and 35.
- Materials of unknown state (NCMAT v7 and earlier without the section,
  non-NCMAT sources) get no elastic scattering, silently.
- Planned: $I(Q)$ or $S(Q)$ curves in NCMAT for SANS and coherent elastic
  models (reminder recorded for 5.0.0), which would replace the incoherent
  approximation here.

## New or changed in 5.0.0

Elastic scattering for amorphous solids described by direct kernels (via
the NCMAT v8 keywords and kernel analysis, chapters 33 and 35), the v8
requirement of an explicit state of matter for every non-crystalline file,
and partial components with warnings when some elements lack displacements.
Otherwise unchanged since 2.7.0, when the state-of-matter concept and the
amorphous-solid elastic component were introduced.

## References

- NCMAT format documentation, sections "The @STATEOFMATTER section" and
  "Effect of state of matter on NCrystal physics modelling"
  (`docs/ncmat_doc.md`).
- T. Kittelmann, X.-X. Cai, CPC 267 (2021) 108082, section 6 (the elastic
  formula reused here).
