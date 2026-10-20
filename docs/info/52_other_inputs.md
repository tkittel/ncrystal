# Other inputs: data sources, .laz/.lau files, quick materials, in-memory data

## Summary

Besides NCMAT files from the standard library, NCrystal materials can come
from files on disk, from in-memory text registered by the application, from
on-the-fly generators (`freegas::`, `solid::`, `gasmix::`), from McStas
`.laz`/`.lau` reflection lists, from plugin-provided data and factories
(chapter 53), and from `Info` objects built directly in C++. All go through
the same factory infrastructure: text-data factories resolve a name to text,
info factories turn text into an `Info`, and scatter/absorption factories
turn `Info` objects into processes. The cfg-string names the data
(`name::file` to pin a source) and carries the parameters.

Code: `ncrystal_core/src/factories/`, `ncrystal_core/src/lazlau/`,
`ncrystal_core/src/quickfact/`, `ncrystal_core/src/gasmix/`;
public customisation API in `NCDataSources.hh`.

## Data sources

| Source name | What it serves | Notes |
|---|---|---|
| `stdlib` | the standard library (chapter 51) | embedded or at `NCRYSTAL_DATADIR`; can be disabled |
| `stdpath` | directories in `NCRYSTAL_DATA_PATH` (plus a compile-time list) | |
| `plugins` | data directories of plugins, `plugins::<Name>/file.ncmat` | only on explicit request |
| `relpath`, `abspath` | files relative to the working directory, absolute paths | can be disabled for sandboxing |
| `virtual` | in-memory files registered by the application | names should end in `.ncmat` |
| `freegas`, `solid`, `gasmix` | on-the-fly NCMAT text from the request string | chapters 34, 24, 42 |

Priorities decide which source wins for a bare name; `name::file` selects
one source. Custom search directories and in-memory files can be added at
run time from C++ and Python (`NCrystal.registerInMemoryFileData`), and
`ncrystal browse` lists what every source provides.

## .laz / .lau reflection lists

NCrystal reads the files used by the McStas `PowderN` (`.laz`) and
`Single_crystal` (`.lau`) components: a header with the unit cell, space
group and (since 5.0.0) optionally the setting, and per-plane lines with
$hkl$, $d$, $|F|^2$ and multiplicity. The hkl list is validated and
normalised against the space group (consistent $|F|^2$ within a family,
listed versus calculated $d$, multiplicities), the cell is completed and
checked, and the dynamics are Debye-model with a Debye temperature
(default 300 K) at the given temperature (default 293.15 K). Files without
a setting are assumed in the default setting; rhombohedral axes are not
supported. The reverse direction, `ncrystal_ncmat2hkl`, writes such files
from any NCrystal material for McStas users (chapter 61).

## Quick materials

- `freegas::<formula>/<density>[/<atomdb>]`: free-gas dynamics for each
  element, `@STATEOFMATTER gas`.
- `solid::<formula>/<density>[/TDebye<T>K][/<atomdb>]`: amorphous solid
  with Debye-model dynamics (default Debye temperature 300 K), e.g.
  `solid::B4C/2.52gcm3`.
- `gasmix::<spec>`: real gas mixtures from formula, pressure, temperature,
  humidity (chapter 42).
- Density units `perAa3`, `kgm3`, `gcm3`; `atomdb` lines with `_`
  separators (`B_is_0.95_B10_0.05_B11`). The generated NCMAT text can be
  inspected with `nctool --extract`.

## Direct construction and the C++/Python API

C++ code (plugins, tests) can build `Info` objects directly from unit cell,
dynamics, composition and density, with the same validation as the NCMAT
loader, and processes directly from their implementations. The Python
`NCMATComposer` (chapter 60) is the usual way to create materials
programmatically and register them in memory. The virtual API (used by
OpenMC) is an ABI-stable facade over the same factories.

## Validation

- `tests/scripts/ncmat2hkl.py`, `ncmat2hkl2.py` (laz/lau round trips incl.
  settings), `tests/src/app_inmemfiles`, `app_textdata`, `app_mtfileread`
  and `tests/scripts/removedatasrc.py` (registration, search paths,
  multithreaded reading), `quickfact.py`, `mcstasunion.py`, and
  `app_vapit1v1` for the virtual API.

## Approximations and limitations

- `.laz`/`.lau` files carry no dynamics beyond a Debye temperature; the
  reader is not a full McStas-format parser.
- Quick materials use idealised dynamics (free gas, Debye model) and
  literature-free densities given by the user.
- `.nxs` files need the nxslib plugin (chapter 53); CIF files are converted
  by the Python tools, not read by the core.

## New or changed in 5.0.0

The space-group setting in `.laz`/`.lau` headers (export and import),
exact d-spacings and cell checks when reading them, and `ncrystal browse`
listing all sources; otherwise the data-source machinery is as in 4.x.

## References

- McStas `PowderN` and `Single_crystal` component documentation (the
  `.laz`/`.lau` formats).
