# Numerics, performance and cross-platform reproducibility

## Summary

NCrystal's physics is implemented so that results are the same on every
platform, compiler and build configuration: identical to the last bit for
the whole VDOS-to-kernel pipeline, structure factors and the mosaic model
since 5.0.0, which replaced every platform-dependent ingredient (system
`exp`/`erf`, compiler-chosen fused multiply-add contraction, noise-driven
discrete decisions) with portable implementations and explicit fused
operations. The same work made the hot paths faster (runtime FMA dispatch
for portably built binaries, branchless searches). This chapter states the
reproducibility policy, the numerical toolbox behind the physics chapters,
and the performance characteristics.

Code: `ncrystal_core/include/NCrystal/internal/utils/` and
`ncrystal_core/src/utils/`; the FMA dispatch rules in
`docs/devel_fma_attribute.md`.

## Reproducibility policy

A difference between platforms is treated as a defect when it is a
genuinely different result (a different grid point chosen, a family split
differently, a cut placed in another bin), not when it is last-digit noise.
Three classes of measures achieve this:

- Portable elementary functions: exponentials, error functions and the
  trigonometric kernels of the structure-factor calculation are evaluated by
  NCrystal's own implementations (correctly rounded to within one or two
  units in the last place, and faster than typical system libraries), so
  results do not depend on the platform's mathematical library.
- Explicit fused multiply-adds wherever a compiler could contract an
  expression in more than one way (quadrature updates, interpolation,
  kinematic limits, sampling utilities, the $|k|^2$ and $|F|^2$ kernels),
  and integer-first evaluation of exactly representable quantities such as
  $1/d^2$.
- Noise-robust discrete decisions: numerical noise is snapped to zero where
  it is known to be noise (FFT round-off), thresholds are replaced by
  smooth tapers or fits where a hard threshold could flip under noise
  (truncation edges, negligible cells, grid point selection), and
  catastrophic cancellation is avoided in all evaluations of the kinematic
  limits.

The policy includes the build: compiler contraction settings, forced FMA
rebuilds, 32-bit x87 versus SSE2 targets and WebAssembly are tested to give
the same reference cross sections, and `long double` is never used.

## Toolbox

Compensated summation (used for cumulative cross sections, kernel
integrals and composition averages) and double-double arithmetic where
extra precision matters; fixed and adaptive Romberg and Simpson quadrature;
monotone cubic (PCHIP) and natural cubic spline interpolation; branchless
binary searches on sorted grids (3-5 times faster than the standard
library, up to 24% faster Bragg sampling and 15% faster kernel sampling);
FFT convolutions with calibrated noise floors; noise-tolerant point
reduction for grid construction; the xoroshiro128+ random number generator
with per-thread streams (applications may supply their own); and strong
types for physical quantities with validation.

## FMA dispatch

Portably built binaries (wheels, conda packages) cannot assume hardware
fused multiply-add. NCrystal therefore compiles its hot loops in several
variants and selects the hardware one at run time (through the compiler's
function multi-versioning on GCC and Clang, and through a dedicated AVX2
build of the same functions on MSVC, where software FMA emulation would
otherwise dominate the VDOS expansion time). The decorated functions are
audited so that every floating-point expression is either an explicit
fused operation or safe under contraction, which is what makes the
dispatch invisible in the results.

## Performance characteristics

- Elastic models (2020-era core, elastic paper section 7): empty-call
  baseline 500 MHz, powder Bragg and incoherent elastic cross sections about
  50 MHz, powder sampling about 10 MHz, mosaic single crystal 1-10 MHz
  above 1.5 Aa (0.1-1 MHz below), layered crystals ten times slower.
- Kernels: a cross section is one spline evaluation after a branchless
  search; a sample costs a few random numbers and a cell lookup at thermal
  energies; the free gas gives about 1 MHz for a full $(\alpha,\beta)$
  sample.
- Initialisation: hkl lists in milliseconds to tens of milliseconds, VDOS
  expansions of the order of 0.1 s at the default level (multithreaded and
  cached across materials), kernel processing similar; kernels of order
  10 MB at the default level. (FIXME: these expansion figures are from the
  4.4.6 `vdoslux` documentation; measure time and memory of the 5.0.0
  expansion and processing per level before quoting them.)
- Threading: all physics objects are immutable after construction and
  shared; per-thread caches and random streams; a thread pool parallelises
  the orders of the VDOS expansion and the loading of materials.

## Validation

Unit tests of every numerical utility (elementary functions, quadrature,
interpolation, searches, FFT, point reduction, random numbers), the
kinematic utilities and the $G_n$ functions against high-precision
references, the reference cross sections of many materials across the CI
platforms, and forced-contraction, 32-bit and WebAssembly builds in
developer setups.

## Approximations and limitations

- Reproducibility is guaranteed for IEEE double with round-to-nearest and
  relies on correctly rounded fused multiply-add (available on all
  supported platforms).
- Last-digit differences can still occur between platforms in places not
  yet audited; only those that change discrete outcomes are fixed.

## New or changed in 5.0.0

Everything under the reproducibility policy and the FMA dispatch is new in
5.0.0, together with the branchless searches, the portable elementary
functions, the per-thread C API error state, and the performance pass over
the kernel and hkl hot paths. Release 4.4.6 was already multithread-safe
with caching and factory threads.

## References

- `docs/devel_fma_attribute.md`.
- W. J. Cody, W. Waite, "Software manual for the elementary functions"
  (1980) (range reduction); the xoroshiro128+ generator of Blackman and
  Vigna.
- T. Kittelmann, X.-X. Cai, CPC 267 (2021) 108082, section 7 (performance
  figures).
