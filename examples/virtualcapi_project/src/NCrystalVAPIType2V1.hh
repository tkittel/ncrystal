
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// As an addition to the usual NCrystal license pasted above, note that THIS  //
// PARTICULAR FILE is also placed into the Public Domain, so that it might    //
// be easily adopted into the code-base of any project wishing to use it:     //
//                                                                            //
// This is free and unencumbered software released into the public domain.    //
//                                                                            //
// Anyone is free to copy, modify, publish, use, compile, sell, or            //
// distribute this software, either in source code form or as a compiled      //
// binary, for any purpose, commercial or non-commercial, and by any          //
// means.                                                                     //
//                                                                            //
// In jurisdictions that recognize copyright laws, the author or authors      //
// of this software dedicate any and all copyright interest in the            //
// software to the public domain. We make this dedication for the benefit     //
// of the public at large and to the detriment of our heirs and               //
// successors. We intend this dedication to be an overt act of                //
// relinquishment in perpetuity of all present and future rights to this      //
// software under copyright law.                                              //
//                                                                            //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,            //
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF         //
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.     //
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR          //
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,      //
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR      //
// OTHER DEALINGS IN THE SOFTWARE.                                            //
//                                                                            //
// For more information, please refer to <https://unlicense.org/>             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef NCrystalVAPIType2V1_hh
#define NCrystalVAPIType2V1_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// A C++ class for using NCrystal through type 2 version 1 of its virtual C   //
// API (see NCrystal/virtualapi/ncvirtapi.h), without any build-time          //
// dependency on NCrystal. Copy this file and NCrystalVAPIType2V1.cc into     //
// your project, and compile them with your other code (C++11 or later; on    //
// Linux and macOS, link with the dl library, e.g. ${CMAKE_DL_LIBS} in CMake).//
//                                                                            //
// At runtime, NCrystalVAPIType2V1::load() locates the NCrystal shared        //
// library via the ncrystal-config command (or, as in NCrystal's Python API,  //
// via the NCRYSTAL_LIB environment variable, with the namespace in           //
// NCRYSTAL_LIB_NAMESPACE_PROTECTION if needed), loads it, and accesses the   //
// interface. All errors are thrown as exceptions derived from                //
// std::runtime_error (NCrystalVAPIType2V1::Error for errors reported by      //
// NCrystal).                                                                 //
//                                                                            //
// Units: eV for neutron energies, barn per atom for cross sections, g/cm^3   //
// for densities, atoms/Aa^3 for number densities, and kelvin for             //
// temperatures.                                                              //
//                                                                            //
// Threads: Info objects can be used by any number of threads at the same     //
// time. Scatter and Absorption objects own caches, and must therefore only   //
// be used by one thread at a time: give each thread its own clone (clones    //
// share the physics, but have separate caches).                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class NCrystalVAPIType2V1 final {
public:

  //Load NCrystal and access the interface (the first call loads it, and later
  //calls return the same object). Throws std::runtime_error if NCrystal can
  //not be found or loaded, or does not provide the interface:
  static std::shared_ptr<const NCrystalVAPIType2V1> load();

  //Errors reported by NCrystal, with the NCrystal error type (e.g.
  //"FileNotFound" or "BadInput"):
  class Error final : public std::runtime_error {
  public:
    Error( const std::string& type, const std::string& message );
    const std::string& type() const noexcept { return m_type; }
  private:
    std::string m_type;
  };

  //Neutron state: kinetic energy (eV) and direction (unit vector):
  struct Neutron {
    double ekin;
    double ux, uy, uz;
  };

  //An entry of the composition of a material (A=0 for natural elements, and
  //fractions by number of atoms):
  struct Component {
    unsigned long Z;
    unsigned long A;
    double fraction;
  };

  //Natural isotope abundances provided by the application: for a given Z, the
  //isotopes as (A,fraction) pairs, or an empty list if not known:
  using NaturalAbundances = std::function<
    std::vector<std::pair<unsigned long,double>>( unsigned long Z ) >;

  //Material information:
  class Info final {
  public:
    //An id which is the same for objects of the same material (e.g. created
    //from the same cfg-string), useful to avoid duplicate materials:
    std::uint64_t uniqueID() const;
    double density() const;//g/cm^3
    double numberDensity() const;//atoms/Aa^3
    double temperature() const;//kelvin (throws if the material has none)

    //The composition. By default (preferNaturalElements=true), natural
    //elements are returned as such (with A=0), and natural abundances are
    //not needed, so most applications can simply call composition() without
    //arguments. Natural abundances are only needed (and an error is thrown
    //if they are not provided) for breaking down natural elements into
    //isotopes. This happens for all elements if preferNaturalElements=false,
    //or otherwise only for the (unusual) materials which contain the same
    //element both as the natural element and as specific isotopes (e.g. a
    //mixture of natural and enriched boron):
    std::vector<Component> composition( bool preferNaturalElements = true,
                                        const NaturalAbundances& = nullptr )
                                        const;

    Info( Info&& ) noexcept;
    Info& operator=( Info&& ) noexcept;
    ~Info();
  private:
    friend class NCrystalVAPIType2V1;
    Info( std::shared_ptr<const NCrystalVAPIType2V1>, void * );
    std::shared_ptr<const NCrystalVAPIType2V1> m_api;
    void * m_h;
  };

  //Scattering:
  class Scatter final {
  public:
    Scatter clone() const;//For usage in another thread.

    //Whether the scattering depends on the neutron direction (e.g. for single
    //crystals), in which case directions are in the frame of the material:
    bool isOriented() const;

    double crossSection( const Neutron& ) const;//barn per atom

    //Sample a scattering, updating the neutron. Random numbers must be
    //uniformly distributed in [0,1) or (0,1]. Exceptions thrown by the
    //random number generator are propagated:
    void sampleScatter( const std::function<double()>& rng, Neutron& ) const;
    void sampleScatter( double (*rng)( void * ), void * rng_state,
                        Neutron& ) const;

    Scatter( Scatter&& ) noexcept;
    Scatter& operator=( Scatter&& ) noexcept;
    ~Scatter();
  private:
    friend class NCrystalVAPIType2V1;
    Scatter( std::shared_ptr<const NCrystalVAPIType2V1>, void * );
    std::shared_ptr<const NCrystalVAPIType2V1> m_api;
    void * m_h;
  };

  //Absorption:
  class Absorption final {
  public:
    Absorption clone() const;//For usage in another thread.
    double crossSection( const Neutron& ) const;//barn per atom

    Absorption( Absorption&& ) noexcept;
    Absorption& operator=( Absorption&& ) noexcept;
    ~Absorption();
  private:
    friend class NCrystalVAPIType2V1;
    Absorption( std::shared_ptr<const NCrystalVAPIType2V1>, void * );
    std::shared_ptr<const NCrystalVAPIType2V1> m_api;
    void * m_h;
  };

  //Create objects from NCrystal cfg-strings (e.g. "Al_sg225.ncmat;temp=200K"):
  Info createInfo( const std::string& cfgstr ) const;
  Scatter createScatter( const std::string& cfgstr ) const;
  Absorption createAbsorption( const std::string& cfgstr ) const;

  //Clear NCrystal's internal caches (e.g. of loaded materials), to release
  //memory. Existing objects stay valid:
  void clearCaches() const;

  NCrystalVAPIType2V1( const NCrystalVAPIType2V1& ) = delete;
  NCrystalVAPIType2V1& operator=( const NCrystalVAPIType2V1& ) = delete;
  NCrystalVAPIType2V1( NCrystalVAPIType2V1&& ) = delete;
  NCrystalVAPIType2V1& operator=( NCrystalVAPIType2V1&& ) = delete;
  ~NCrystalVAPIType2V1() = default;

private:
  explicit NCrystalVAPIType2V1( const void * raw_api ) : m_raw( raw_api ) {}
  const void * m_raw;//The struct of the C API.
};

#endif
