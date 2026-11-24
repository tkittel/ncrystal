# Changelog input: CIF parsing and Python messages

To be curated into the CHANGELOG at release time.

## CIF data

* New JSON query `["cif","structure",...]`: the structure and provenance
  data items of CIF data, read by a new C++ CIF 1.1 parser, and the space
  group setting chosen for it (now in C++, with the same rules and
  messages), with any warnings returned in the result.
* C interface: error messages are no longer cut off after 511 characters
  (now 4095).
* CIF data loads faster: sites are collected and merged in C++ (e.g. a large
  cell taking 85 s before now loads in about 1 s). The
  `["cif","structure",...]` query also returns the cell snapped to the
  setting and the atoms. Messages about CIF sites name the data and line.
* Without spglib refinement (and without `override_spacegroup`), CIF data
  whose cell and atoms have a higher symmetry in the same coordinates (within
  `tol_symprec`, e.g. P1 data of a P-1 structure) is loaded in that setting,
  with a warning. Displacements contradicting it are an error, or keep the
  given setting in permissive mode.
* `override_spacegroup` (and cif2ncmat `--spacegroup`) also accepts Hall
  symbols, and the prefixes "HM:" and "Hall:". Symbols which mean different
  settings as Hermann-Mauguin and as Hall symbols (e.g. "P 2") are an error
  without prefix.
* New `cifutils.spacegroup_settings()` and `ncrystal_cif2ncmat
  --list-spacegroups` list the space group settings with their
  Hermann-Mauguin and Hall symbols.
* Warnings in the comments of NCMAT data created from CIF data
  (cif2ncmat) are wrapped, at most 10 lines each. They depend only on the
  input, not on quiet flags, logging or the warnings already shown, and no
  longer include warnings from other Python packages (gemmi, spglib, ...).
* `NCMATComposer.refine_crystal_structure(quiet=True)` also hides the
  warnings about dropped anisotropic displacements.
* CIF data is read with the new C++ CIF 1.1 parser instead of gemmi, which is
  no longer needed. Error messages name the data, block and line. In strict
  mode, anisotropic displacements with missing tags, with labels not
  belonging to a site, or with "?" components are an error (permissive mode
  ignores them with a warning), and sites with "?" coordinates are an error
  even if unoccupied. Space group symbols which are "?" are skipped in
  favour of alternative tags. Anisotropic displacements given as B are
  converted to U (B = 8 pi^2 U), like Biso.
* `CIFLoader.raw_cif_dict` is removed (incompatible change for code using
  it).
* Fixed author names in the comments of NCMAT data created from CIF data with
  a single `_publ_author_name` value (e.g. "F, et al." instead of the name).
* The DOI in the comments of NCMAT data created from CIF data can also come
  from the primary citation (`_citation_doi`).
* Data downloaded from materialsproject.org gets the data item
  `_[local]_ncrystal_mp_database_code`.
* Repeated identical warnings while loading CIF data are shown once (the
  NCMAT comments list them with their count).

## Messages and warnings (Python)

* Notices, warnings and raw output from NCrystal (from both the Python and
  C++ code) now go through one channel. Warnings are printed as
  "NCrystal WARNING: ..." (and notices as "NCrystal: ...") via the
  redirectable print function, instead of being emitted via Python's
  warnings module, and they are no longer deduplicated. The new
  `NCrystal.use_logging()` sends notices and warnings (optionally also raw
  output) to the logger "NCrystal" instead (GitHub issues #114 and #361).
* New module `NCrystal.msg` with `use_logging` and the `MsgCtrl` context
  manager, which sets prefixes or a logger, and collects, suppresses,
  deduplicates or prints despite logging notices, warnings (incl. Python
  warnings, as the warnings module shows them, without changing its
  filters) and raw output.
* Status lines of library functions and command-line tools (e.g.
  "Wrote <file>") are now notices rather than raw output.
* Deprecated functionality now emits `NCrystalFutureWarning` (a
  `FutureWarning`, via the warnings module). The `NCrystalUserWarning`
  class is removed (incompatible change for code using it).

## Command-line tools

* Command-line tools print notices without prefix and warnings as
  "WARNING: ...". `NCrystal.cli.run` gets `msg_style='api'` (default) or
  `'cli'` (output exactly as in the terminal).
* New `--quiet` option for ncrystal_atomdb and ncrystal_browse, which
  silences all output except the requested data.
* The "-o stdout" output of ncrystal_cif2ncmat, ncrystal_ncmat2hkl and
  ncrystal_minimc (and cif2ncmat `--showcif`) now goes through the
  redirectable print function, e.g. when using `NCrystal.cli.run`.
* The obsolete "python -mNCrystal.mcstasutils" usage now fails; use
  ncrystal_mcstasunion or ncrystal_ncmat2hkl instead.
