# Distribution packaging

## Debian

The `debian/` subdirectory contains the template for building Debian packages of
NCrystal (also for derivatives like Ubuntu). It needs Debian 13 (trixie) or
Ubuntu 26.04, or later. For Ubuntu 24.04, see below.

To build the packages, run the following in the root directory of an NCrystal
source tree. Nothing needs to be installed beforehand: `apt build-dep` installs
everything needed for the build (dpkg-buildpackage, debhelper, cmake, compilers,
Python build tools, ...), as listed in `debian/control`.

```
devel/packaging/prepare-debian
sudo apt build-dep ./
dpkg-buildpackage -b -us -uc
```

`prepare-debian` copies the template to `debian/`, names the library package
after the major version in the `VERSION` file (`libncrystal<major>`), and writes
`debian/changelog` with that version.

The packages (`.deb` files) end up in the parent directory. To check them with
`lintian *.changes`, or to run the package tests in `debian/tests` on the
installed packages with `autopkgtest`, first `sudo apt install lintian
autopkgtest`.

For Ubuntu 24.04, whose Python build tools are older, use
`devel/packaging/prepare-debian --ubuntu2404` instead, which also applies
`ubuntu2404.patch`.

Binary packages (`apt install ncrystal` installs them all):

- `libncrystal<major>`: the shared library (named after its SONAME).
- `libncrystal-dev`: headers and CMake files.
- `python3-ncrystal`: the NCrystal Python module.
- `python3-ncrystal-verify`: the `ncrystal-verify` command, also used by the
  package tests.
- `ncrystal`: all command-line tools (incl. `ncrystal-config`).

Notes:

- Command-line tools and their man pages come from `[project.scripts]` in
  `ncrystal_python/pyproject.toml`, so new tools need no changes here.
- Packages are built with `-O3 -DNDEBUG`, as in upstream's default builds.
  `DEB_BUILD_OPTIONS=noopt` gives a debug build (`-O0`, assertions enabled).
  Debug symbols go into the `-dbgsym` packages (or stay in the binaries with
  `DEB_BUILD_OPTIONS=nostrip`).
