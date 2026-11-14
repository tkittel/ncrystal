# Distribution packaging

## Debian

The `debian/` subdirectory contains the files for building Debian packages of
NCrystal (also for derivatives like Ubuntu). It needs Debian 13 (trixie) or
Ubuntu 26.04, or later.

To build the packages, run the following in the root directory of an NCrystal
source tree. Nothing needs to be installed beforehand: `apt build-dep` installs
everything needed for the build (dpkg-buildpackage, debhelper, cmake, compilers,
Python build tools, ...), as listed in `debian/control`.

```
cp -r devel/packaging/debian debian
sudo apt build-dep ./
dpkg-buildpackage -b -us -uc
```

The packages (`.deb` files) end up in the parent directory. To check them with
`lintian *.changes`, or to run the package tests in `debian/tests` on the
installed packages with `autopkgtest`, first `sudo apt install lintian
autopkgtest`.

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
- A new SONAME (major version) means renaming the `libncrystal<major>` package
  in `control` and its `.install` file.
- `debian/changelog` only holds a template entry, which is replaced with the
  actual version when building packages for distribution.
