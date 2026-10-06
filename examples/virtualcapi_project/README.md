Virtual C API example
=====================

Small example showing how a C++ application can use NCrystal without any
build-time dependency on NCrystal, via type 2 version 1 of NCrystal's virtual C
API (`NCrystal/virtualapi/ncvirtapi.h`). This API provides material information
(composition, density, temperature, ...), scattering (cross sections and
sampling), and absorption cross sections, and is intended for applications
like FLUKA and Geant4.

Instead of working with the plain C API directly, this example uses the
`NCrystalVAPIType2V1` class, defined in `src/NCrystalVAPIType2V1.hh` and
`src/NCrystalVAPIType2V1.cc`. Applications can simply copy these two files
(which are in the public domain) into their own projects, and compile them with
their other code. The class takes care of locating NCrystal at runtime (via
the `ncrystal-config` command), loading its shared library, and accessing the
API, and it gives a C++ interface with standard library types and exceptions:

```
auto nc = NCrystalVAPIType2V1::load();//throws if NCrystal is not available
auto info = nc->createInfo( "stdlib::Al_sg225.ncmat;temp=200K" );
auto scatter = nc->createScatter( "stdlib::Al_sg225.ncmat;temp=200K" );
double density = info.density();
NCrystalVAPIType2V1::Neutron n{ 0.025, 0.0, 0.0, 1.0 };//ekin [eV], direction
double xs = scatter.crossSection( n );
```

The application can thus have optional support for NCrystal without
complicated build-time flags and dependency management, and distributed
packages of the application work with future releases of NCrystal as they come
out.

To build and run the example (NCrystal is only needed for running it):

```
cmake -S . -B build
cmake --build build
./build/ncrystal-virtcapi-example
```

Applications which prefer to use the plain C API can instead include a copy of
`NCrystal/virtualapi/ncvirtapi.h`, and use `NCrystalVAPIType2V1.cc` as an
example of how to load and access it. For an example using the older C++
virtual API (type 1 version 1), see `examples/virtualapi_project`.
