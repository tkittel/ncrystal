# Atomic data: scattering lengths, cross sections, masses, mixtures

## Summary

Every element or isotope in a material is represented by an `AtomData`
object carrying its average mass, coherent scattering length $b_{coh}$,
bound incoherent cross section $\sigma_{inc}$ and thermal capture cross
section $\sigma_{abs}$, from which the coherent, total and free cross
sections follow. The built-in database has 348 entries (natural elements and
isotopes, Sears-based scattering data and AME2016 masses); files and
cfg-strings can override values, define unknown isotopes, and build mixtures
of elements and isotopes (enrichment, doping, substitutional disorder) whose
effective constants, including the disorder contribution to the incoherent
cross section, are computed exactly. The data can be browsed from the
command line and Python.

Code: `ncrystal_core/include/NCrystal/interfaces/NCAtomData.hh`,
`ncrystal_core/src/atomdb/`.

## Physics

For a mixture of components with fractions $f_i$ (by count):

- mass $\bar m = \sum_i f_i m_i$ and coherent scattering length
  $\bar b = \sum_i f_i b_i$ (averages);
- coherent cross section $\sigma_{coh} = 4\pi\bar b^2$;
- incoherent cross section
  $\sigma_{inc} = \sum_i f_i\sigma_{inc,i} + 4\pi\big(\sum_i f_i b_i^2 - \bar b^2\big)$,
  i.e. the average of the components' spin/isotope incoherence plus the
  variance of the scattering length over the mixture, which is exactly the
  Laue diffuse scattering of random substitution (CPC 2020 eqs. 13-16);
- total bound scattering $\sigma_b = \sigma_{coh} + \sigma_{inc}$, and free
  cross sections $\sigma_{free} = \sigma_b (A/(1+A))^2$ with $A = \bar m/m_n$
  (used by the free-gas model, chapter 34);
- capture $\sigma_{abs}$ at 2200 m/s, averaged (chapter 40).

Scattering lengths are treated as real constants (absorption is handled as
a separate process), valid to $O(10^{-4})$ for thermal neutrons (CPC 2020
section 2.1). Internally lengths are stored in $\sqrt{\mathrm{barn}} = 10$ fm
so that $\sigma = 4\pi b^2$ without conversion factors.

## The database

- 348 entries keyed by $(Z, A)$ with $A = 0$ for the natural element:
  mass, $b_{coh}$, $\sigma_{inc}$, $\sigma_{abs}$.
- Scattering data from the ATI Vienna "Table of coherent scattering lengths
  and cross sections" (the recommended values of Landolt-Bornstein I/16A
  and Rauch/Waschkowski, following Sears 1992), parsed from its HTML table;
  masses and natural abundances from AME2016 (Huang et al., Chinese Physics
  C 41 (2017) 030002), natural-element masses computed from the isotopes
  (slightly different from CIAAW values).
- Manual fixes documented at the end of the generated table: four legacy
  element entries kept from before 2.1.0; Xe $\sigma_{inc}$ set to 1.16133 b
  (was 0 in the source); Ca $b_{coh}$ 4.7 -> 5.795 fm and Ca-40 5.9096 fm
  per NIM A 2019 (doi:10.1016/j.nima.2019.02.072).
- Natural abundances are present in the generator but not used; isotopic
  compositions other than natural are specified explicitly as mixtures.
- Labels: element symbols, `El<A>` for isotopes (`B10`, `Si29`), `D` and
  `T` as aliases of `H2` and `H3`, and the generic labels `X`, `X1`..`X99`
  which may only denote mixtures (never carry constants directly).

## Inputs and configuration

- NCMAT `@ATOMDB` (v3 and later; `docs/ncmat_doc.md`): lines
  `Sym mass u b_coh fm sigma_inc b sigma_abs b` to define or override an
  element or isotope, `X is f1 A f2 B ...` to define mixtures (fractions
  summing to 1; `X is Al` as an alias; isotope labels cannot be redefined
  as mixtures), evaluated in order so that later lines can build on earlier
  ones (`B is 0.9 B10 0.1 B11` followed by `B is 0.999 B 0.001 C`), and
  `nodefaults` on the first line to disable the built-in database (for
  validation). Examples: enriched boron, chromium doping of aluminium
  positions, a hydrogen-oxygen pseudo-atom `X`.
- The `atomdb=` cfg parameter with the same syntax, `:` for spaces and `@`
  between lines (`atomdb=B:is:0.95:B10:0.05:B11`), appended after the file's
  own `@ATOMDB`; and the `_`-separated form in the quick factories
  (`solid::B4C/2.5gcm3/B_is_0.95_B10_0.05_B11`).

## Implementation

`AtomData` objects are immutable and classified as a single isotope, a
natural element, or a composite with components; identical definitions in
different files are shared. Python: `NCrystal.atomdata` (`atomDB(Z, A)`,
`iterateAtomDB`), `AtomData` objects from `Info.getComposition()` and
`AtomInfo.atomData` (`.coherentScatLen`, `.incoherentXS`, `.captureXS`,
`.averageMassAMU`, `.isNaturalElement`, ...), `ncrystal browse --atomdb`
with physics selections such as `-w "absxs > 1000 and natural"`, and
`NCMATComposer.update_atomdb` / `remap_atom`.

## Validation

- `tests/scripts/atomdbpy.py` (the Python database access) and
  `ncmatparser.py` (syntax and semantics of `@ATOMDB`), the `nctool` dumps
  showing per-atom constants, and the composition tests of multiphase
  materials; `tests/data/Al_multirole.ncmat` exercises several roles of one
  element.
- The sources of the constants are reference tables; NCrystal adds no
  evaluation of its own beyond the three documented fixes.

## Approximations and limitations

- Real, energy-independent scattering lengths; no resonance or
  energy-dependent $b$, no magnetic form factors, no nuclear polarisation.
- Natural abundances are not applied automatically to make isotopic
  mixtures (natural elements use the tabulated natural-abundance values of
  $b$ and $\sigma$ directly, which is what the tables provide).

## New or changed in 5.0.0

The database and mixture machinery date from 2.1.0 (user-modifiable atom
data) and are unchanged in physics. New in 5.0.0: the database is
browsable (`ncrystal browse --atomdb`, `NCrystal.browse.AtomDBBrowser`).
Partial site occupancy, currently emulated with sterile pseudo-isotope
mixtures by the composer, is to become a proper per-site property (chapter
10).

## References

- V. F. Sears, Neutron News 3 (1992) 26; Landolt-Bornstein New Series
  I/16A (2000), chapter 6; the ATI Vienna scattering-length table.
- W. J. Huang et al., "The AME2016 atomic mass evaluation (I)", Chinese
  Physics C 41 (2017) 030002 (as cited in the database generator; the mass
  table `mass16.txt`).
- The calcium scattering-length re-evaluation, NIM A (2019),
  doi:10.1016/j.nima.2019.02.072.
- X.-X. Cai, T. Kittelmann, CPC 246 (2020) 106851, section 2.1.
