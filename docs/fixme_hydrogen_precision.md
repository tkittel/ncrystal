# Precision of VDOSDebye phonon expansions (most visible for hydrogen)

Temporary handover report for continuing the investigation. Delete once
resolved.

## Summary

Inelastic cross sections from NCrystal's phonon expansions of VDOSDebye
kernels (the normal path: `extractSABDataFromVDOSDebyeModel` -> `SABScatter`)
deviate from the exact Debye-model cross sections by up to ~2e-3, and these
deviations do **not** shrink with increasing vdoslux (they converge to an
offset, and at the lowest energies even grow with vdoslux). Two patterns:

1. **Hydrogen (and to lesser degree deuterium), thermal and above**: expansions
   are low by 1.4e-3 (0.0253 eV) to 2.2e-3 (2 eV) at vdoslux=6.
2. **All elements, lowest energies (E=1e-5 eV, 1/v region)**: expansions are
   high by ~2e-3 at vdoslux=6 (O, C, Al, Pb), growing from ~1.2e-3 at
   vdoslux=3.

The exact values are from `VDOS::VDOSDebyeXS` (exact one-phonon part plus a
fitted multi-phonon part), which agrees with the independent reference
implementation `DebyeRef` in
`devel/simplebuild/static/NCDev/scripts/fit_vdosdebye_xs` to within 4e-4 in all
cases checked (mostly ~1e-4), while the expansions deviate from that same
reference by up to 2.6e-3. So the expansions (or the VDOS they expand) are
the likely culprit, not `VDOSDebyeXS`.

## Measurements

T=293.15K, single-element materials (one vdosdebye `@DYNINFO`, `comp=inelas`).
Values are (expansion/exact - 1) for vdoslux=3,4,5,6, where "exact" is
`VDOSDebyeXS` (via `NCRYSTAL_MULTIVDOSDEBYE=always`):

```
el TDebye  E[eV]  lux3      lux4      lux5      lux6
H    600   1e-05: -8.79e-04 -7.64e-04 -7.10e-04 -6.87e-04
H    600  0.0253: -1.51e-03 -1.45e-03 -1.42e-03 -1.40e-03
H    600     0.1: -1.94e-03 -1.83e-03 -1.78e-03 -1.76e-03
H    600     0.5: -3.06e-03 -2.40e-03 -2.19e-03 -2.10e-03
H    600       2: -3.94e-03 -2.78e-03 -2.40e-03 -2.21e-03
H    300   1e-05: -1.13e-03 -9.30e-04 -8.96e-04 -8.79e-04
H    300  0.0253: -1.53e-03 -1.56e-03 -1.56e-03 -1.57e-03
H    300     0.1: -2.43e-03 -2.16e-03 -2.09e-03 -2.06e-03
H    300     0.5: -3.60e-03 -2.71e-03 -2.45e-03 -2.33e-03
H    300       2: -4.30e-03 -2.99e-03 -2.57e-03 -2.32e-03
D    600   1e-05: +3.17e-04 +4.44e-04 +4.99e-04 +5.36e-04
D    600  0.0253: -5.25e-04 -4.32e-04 -3.92e-04 -3.68e-04
D    600     0.1: -7.11e-04 -6.02e-04 -5.60e-04 -5.39e-04
D    600     0.5: -2.21e-03 -1.65e-03 -1.46e-03 -1.38e-03
D    600       2: -2.85e-03 -1.96e-03 -1.66e-03 -1.52e-03
Li   400   1e-05: +9.79e-04 +1.33e-03 +1.49e-03 +1.56e-03
Li   400  0.0253: +1.19e-04 +2.43e-04 +3.07e-04 +3.39e-04
Li   400     0.1: -2.84e-04 -1.79e-04 -1.40e-04 -1.21e-04
Li   400     0.5: -9.23e-04 -6.36e-04 -5.28e-04 -4.90e-04
Li   400       2: -1.07e-03 -6.62e-04 -5.30e-04 -4.91e-04
C    400   1e-05: +1.15e-03 +1.57e-03 +1.75e-03 +1.86e-03
C    400  0.0253: +2.23e-04 +3.98e-04 +4.76e-04 +5.28e-04
C    400     0.1: -1.68e-05 +8.21e-05 +1.21e-04 +1.43e-04
C    400     0.5: -3.87e-04 -2.38e-04 -1.72e-04 -1.52e-04
C    400       2: -6.23e-04 -3.77e-04 -2.99e-04 -3.03e-04
O    400   1e-05: +1.15e-03 +1.71e-03 +1.92e-03 +2.05e-03
O    400  0.0253: +1.08e-04 +3.76e-04 +4.78e-04 +5.48e-04
O    400     0.1: -5.79e-05 +6.15e-05 +1.07e-04 +1.34e-04
O    400     0.5: -4.68e-04 -3.60e-04 -3.10e-04 -2.96e-04
O    400       2: -5.59e-04 -3.73e-04 -3.14e-04 -3.23e-04
Al   400   1e-05: +1.17e-03 +1.80e-03 +2.03e-03 +2.17e-03
Al   400  0.0253: +1.38e-04 +3.90e-04 +5.15e-04 +5.89e-04
Al   400     0.1: +2.40e-04 +3.99e-04 +4.64e-04 +5.10e-04
Al   400     0.5: -3.32e-04 -2.83e-04 -2.53e-04 -2.41e-04
Al   400       2: -3.14e-04 -2.01e-04 -1.65e-04 -1.87e-04
Pb   100   1e-05: +1.36e-03 +2.12e-03 +2.41e-03 +2.54e-03
Pb   100  0.0253: -1.87e-04 +9.68e-05 +2.33e-04 +3.15e-04
Pb   100     0.1: +3.02e-04 +4.04e-04 +4.46e-04 +4.73e-04
Pb   100     0.5: +3.47e-04 +3.53e-04 +3.66e-04 +3.68e-04
Pb   100       2: +4.93e-05 -2.69e-05 -2.29e-05 -3.99e-05
```

Against the reference `DebyeRef` (sigma/sigma_bound; E = 1e-5, 0.0253, 0.1,
0.5, 2 eV):

```
H  600K exact/ref-1: +5.0e-05 -9.2e-05 -1.2e-04 -4.2e-05 +2.5e-05 | lux6/ref-1: -6.4e-04 -1.5e-03 -1.9e-03 -2.1e-03 -2.2e-03
O  400K exact/ref-1: -5.4e-05 -2.2e-04 -1.9e-06 +1.0e-04 +7.9e-05 | lux6/ref-1: +2.0e-03 +3.2e-04 +1.3e-04 -2.0e-04 -2.4e-04
Pb 100K exact/ref-1: +1.1e-04 +3.9e-04 -1.8e-04 -3.3e-04 +3.4e-05 | lux6/ref-1: +2.6e-03 +7.0e-04 +2.9e-04 +3.4e-05 -6.3e-06
C  400K exact/ref-1: +3.2e-05 -3.3e-04 -1.9e-04 -1.4e-04 -3.3e-05 | lux6/ref-1: +1.9e-03 +2.0e-04 -4.3e-05 -3.0e-04 -3.4e-04
```

(Also at 1e-3 eV for H: exact/ref-1=-8.4e-05, lux6/ref-1=-6.9e-04.)

## Consequence for vdoslux=6

Materials with VDOSDebye entries use exact cross sections (`VDOSDebyeXS`, via
the grouped VDOSDebye treatment) at vdoslux 0-5, but vdoslux=6 does not use the
grouping and therefore keeps the cross sections of the expansions. Until the
deviations above are fixed, the highest vdoslux thus has the least accurate
VDOSDebye cross sections (by up to ~2e-3). A possible remedy is to use the
grouped treatment at vdoslux=6 as well, but with no merging of entries.

Once the expansions are fixed, vdoslux=6 (and perhaps lower values) should
ideally be more accurate than the `VDOSDebyeXS` parameterisation (1e-3
target, 1e-2 for T/TDebye in [50,1000]). Then revisit: (a) which cross
sections the grouped treatment uses at high vdoslux, and (b) the `CalcError`
for T/TDebye>1000, whose message currently suggests vdoslux=100+value (no
grouping, i.e. a separate phonon expansion for each entry) as an alternative
"which has no such hard limit, but can be slower to initialise".

## Hypotheses (unverified)

- **Low-energy VDOS representation.** `extractFromDIVDOSDebyeNoCache` (in
  `NCDynInfoUtils.cc`) builds the Debye VDOS with points only in the upper 50%
  of [0,k*TDebye], relying on the quadratic extrapolation of `VDOSEval` below
  the first point. The E->0 (1/v) cross section is governed by the low-energy
  VDOS (one-phonon, |beta| -> 0) and the Debye-Waller factor. An offset there
  could explain pattern 2. Check by comparing the one-phonon G1 and the
  Debye-Waller integral from the VDOS used in the expansion against the
  analytic Debye forms (`VDOSDebyeXS` and `DebyeRef` have them).
- **Gn convolution resolution.** `VDOSGn` uses min_nbins=400, which earlier
  notes suggested limits VDOS cross sections to ~1e-3. A resolution limit
  would give an offset that does not shrink with vdoslux, as seen.
- **Light atoms, high energies.** For A~1 the multi-phonon terms matter up to
  high orders already at thermal energies, and the SCT extension above Emax
  and the order truncation may both contribute (pattern 1). Compare the
  expansion's sigma_1 and sum of sigma_n>=2 separately against
  `VDOSDebyeXS::sigmaParts` to see which part is off.
- **SAB grid / integration.** The kernel is integrated on (alpha,beta) grids;
  for H the kinematic region at fixed E covers large alpha ranges. Check by
  integrating the expanded S(alpha,beta) directly with a fine quadrature
  (bypassing SABProcessor) to separate kernel errors from integration errors.

## Suggested next steps

1. Split the deviation into one-phonon and multi-phonon parts: get sigma_1 from
   the expansion (e.g. vdos2sab exclude flags, cf. `vdos2sabExcludeFlag` in
   `NCDynInfoUtils.hh`, to keep only order 1) and compare with
   `VDOSDebyeXS::sigmaParts`.
2. Compare the Debye-Waller factor (lambda) and G1 used by the expansion with
   the analytic ones.
3. Redo with a VDOS grid covering all of [0,k*TDebye] to test the first
   hypothesis.
4. Check whether raising VDOSGn min_nbins removes the offsets.

## Reproduction

Exact values via `NCRYSTAL_MULTIVDOSDEBYE=always` (MultiVDOSDebyeScatter with
one entry uses `VDOSDebyeXS` for cross sections); expansions via `=never`.
Run with `ncdevtool sbrun python3 script.py`:

```python
import NCrystalDev as NC
from NCTestUtils.env import ncsetenv
el, td = 'H', 600.0
NC.registerInMemoryFileData( 'chk.ncmat', '\n'.join( [
    'NCMAT v7', '@STATEOFMATTER', '  solid', '@DENSITY', '  1.0 g_per_cm3',
    '@DYNINFO', f'  element {el}', '  fraction 1', '  type vdosdebye',
    f'  debye_temp {td}' ] ) )
def xs( mode, lux, e ):
    ncsetenv( 'MULTIVDOSDEBYE', mode )
    NC.clearCaches()
    sc = NC.createScatter( f'chk.ncmat;vdoslux={lux};comp=inelas' )
    return sc.xsect( wl = NC.ekin2wl( e ) )
for e in ( 1e-5, 0.0253, 0.1, 0.5, 2.0 ):
    ex = xs( 'always', 3, e )
    print( e, [ xs( 'never', lux, e ) / ex - 1 for lux in ( 3, 4, 5, 6 ) ] )
```

Reference values (sigma/sigma_bound) from the fit script (plain python3 with
numpy; takes a few minutes; A is the mass in neutron masses, us are
E/(k*TDebye)):

```python
import importlib.machinery, importlib.util
p = 'devel/simplebuild/static/NCDev/scripts/fit_vdosdebye_xs'
ldr = importlib.machinery.SourceFileLoader( 'fitmod', p )
m = importlib.util.module_from_spec( importlib.util.spec_from_loader( 'fitmod', ldr ) )
ldr.exec_module( m )
k, T, td, mass = 8.617333262e-5, 293.15, 600.0, 1.007975973752
es = [ 1e-5, 0.0253, 0.1, 0.5, 2.0 ]
s1, s2 = m.DebyeRef( td / T ).parts( [ e / ( k * td ) for e in es ],
                                     mass / 1.00866491606 )
print( [ a + b for a, b in zip( s1, s2 ) ] )
```
