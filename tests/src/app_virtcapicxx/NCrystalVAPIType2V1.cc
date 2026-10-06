
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

#include "NCrystalVAPIType2V1.hh"
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <mutex>

#if !defined(NCVAPI_WINDOWS) && ( defined (_WIN32) || defined (WIN32) )
#  define NCVAPI_WINDOWS
#endif
#ifdef NCVAPI_WINDOWS
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

// The following is an unmodified copy of the declarations of type 2 version 1
// from NCrystal/virtualapi/ncvirtapi.h (only whitespace and comments may
// differ), and must NOT be changed:

// BEGIN COPY FROM ncvirtapi.h
#include <stddef.h>
extern "C" {

  /* Error information. The strings are always null-terminated (and possibly */
  /* truncated):                                                             */
  typedef struct {
    int code;            /* Non-zero when an error occurred.                */
    char type[64];       /* Error type (e.g. "BadInput" or "FileNotFound"). */
    char message[1024];  /* Error message.                                  */
  } ncrystal_vapi_error_t;

  /****************************************************************************/
  /* Type 2, version 1 (interface id 2001, available from NCrystal 4.5.0).    */
  /*                                                                          */
  /* Material information, scattering and absorption, as needed by e.g.       */
  /* FLUKA and Geant4. As in NCrystal itself, there are three kinds of        */
  /* handles, created independently from cfg-strings:                         */
  /*                                                                          */
  /* * Info handles: material information (composition, density,              */
  /*   temperature, ...). They have no caches, and can therefore be used by   */
  /*   any number of threads at the same time.                                */
  /* * Scatter and absorption handles: cross sections (and for scattering,    */
  /*   sampling). They own caches, exactly as the scatter handles of type 1   */
  /*   version 2 (see above), so each material in each thread needs its own   */
  /*   handle, which multi-threaded applications obtain by cloning.           */
  /*                                                                          */
  /* The conventions are as in type 1 version 2: barn per atom for cross      */
  /* sections, eV for neutron energies, a "neutron" parameter is an array     */
  /* (ekin,ux,uy,uz), and functions returning int return 0 on success and     */
  /* non-zero on failure (in which case output parameters are not modified),  */
  /* while functions returning handles return NULL on failure. Functions      */
  /* without an error parameter can not fail, but must be given valid         */
  /* (non-NULL) handles.                                                      */
  /****************************************************************************/

  typedef struct {
    unsigned long interface_id;/* Always 2001.                              */
    size_t struct_size;/* sizeof(ncrystal_vapi_type2_v1_t) in NCrystal.     */

    /* Info handles (deallocating a NULL handle does nothing):                */
    struct ncrystal_vapi_t2v1_info *
      (*create_info)( const char * cfgstr, ncrystal_vapi_error_t * );
    void (*deallocate_info)( struct ncrystal_vapi_t2v1_info * );

    /* A unique id of the material (an integral value in [0,2^53]). Handles   */
    /* with the same id have the same material, and handles created from the  */
    /* same cfg-string get the same id (as long as NCrystal's caches are not  */
    /* cleared), so it can be used to avoid creating duplicate materials in   */
    /* the application:                                                       */
    double (*info_unique_id)( const struct ncrystal_vapi_t2v1_info * );

    double (*info_density)( const struct ncrystal_vapi_t2v1_info * );
                                                         /* g/cm^3        */
    double (*info_number_density)( const struct ncrystal_vapi_t2v1_info * );
                                                         /* atoms/Aa^3    */

    /* The temperature in kelvin (fails if the material does not have one):   */
    int (*info_temperature)( const struct ncrystal_vapi_t2v1_info *,
                             double * temperature,
                             ncrystal_vapi_error_t * );

    /* The composition as (Z,A,fraction) entries, with A=0 for natural        */
    /* elements and fractions by number of atoms (summing to 1). It writes    */
    /* at most "capacity" entries to the arrays Z, A and fraction (which may  */
    /* be NULL if capacity is 0), and returns the total number of entries (so */
    /* a first call with capacity 0 gives the needed capacity), or 0 on       */
    /* failure. With prefer_natural_elements=1, elements which are only       */
    /* present as the natural element are returned with A=0, and all other   */
    /* elements (and with prefer_natural_elements=0 all elements) are broken  */
    /* down into isotopes, which needs natural isotope abundances from the    */
    /* application: natabund must then fill in at most "capacity" (at least   */
    /* 128) isotopes of element Z, as mass numbers and fractions by number of */
    /* atoms, and return their number (or 0 if the abundances of the element  */
    /* are not known). It is passed natabund_state unchanged, and may only be */
    /* NULL if no break down is needed.                                       */
    size_t (*info_composition)( const struct ncrystal_vapi_t2v1_info *,
                                int prefer_natural_elements,
                                size_t (*natabund)( void * natabund_state,
                                                    unsigned long Z,
                                                    unsigned long * A,
                                                    double * fraction,
                                                    size_t capacity ),
                                void * natabund_state,
                                unsigned long * Z,
                                unsigned long * A,
                                double * fraction,
                                size_t capacity,
                                ncrystal_vapi_error_t * );

    /* Scatter handles (as in type 1 version 2):                              */
    struct ncrystal_vapi_t2v1_scatter *
      (*create_scatter)( const char * cfgstr, ncrystal_vapi_error_t * );
    struct ncrystal_vapi_t2v1_scatter *
      (*clone_scatter)( const struct ncrystal_vapi_t2v1_scatter *,
                        ncrystal_vapi_error_t * );
    void (*deallocate_scatter)( struct ncrystal_vapi_t2v1_scatter * );

    /* Whether the scattering depends on the neutron direction (e.g. for      */
    /* single crystals), in which case directions are in the frame of the     */
    /* material. Returns 1 if so, and 0 for isotropic materials:              */
    int (*scatter_is_oriented)( const struct ncrystal_vapi_t2v1_scatter * );

    int (*scatter_cross_section)( struct ncrystal_vapi_t2v1_scatter *,
                                  const double * neutron,
                                  double * xsect,
                                  ncrystal_vapi_error_t * );

    int (*sample_scatter)( struct ncrystal_vapi_t2v1_scatter *,
                           double (*rng)( void * rng_state ),
                           void * rng_state,
                           double * neutron,
                           ncrystal_vapi_error_t * );

    /* Absorption handles (with caches, like scatter handles):                */
    struct ncrystal_vapi_t2v1_absorption *
      (*create_absorption)( const char * cfgstr, ncrystal_vapi_error_t * );
    struct ncrystal_vapi_t2v1_absorption *
      (*clone_absorption)( const struct ncrystal_vapi_t2v1_absorption *,
                           ncrystal_vapi_error_t * );
    void (*deallocate_absorption)( struct ncrystal_vapi_t2v1_absorption * );

    int (*absorption_cross_section)( struct ncrystal_vapi_t2v1_absorption *,
                                     const double * neutron,
                                     double * xsect,
                                     ncrystal_vapi_error_t * );

    /* Clear NCrystal's internal caches (e.g. of loaded materials), to        */
    /* release memory. Existing handles stay valid:                           */
    void (*clear_caches)( void );
  } ncrystal_vapi_type2_v1_t;

}
// END COPY FROM ncvirtapi.h

namespace {

  using API = ncrystal_vapi_type2_v1_t;
  using Error = NCrystalVAPIType2V1::Error;

  //Collects errors from the C API, and throws them as exceptions:
  class ErrorCollector {
  public:
    ErrorCollector() { std::memset( &m_err, 0, sizeof(m_err) ); }
    ncrystal_vapi_error_t * ptr() noexcept { return &m_err; }
    void throwIf( bool failed )
    {
      if ( failed || m_err.code != 0 )
        throw Error( m_err.code ? m_err.type : "UnknownError",
                     m_err.code ? m_err.message : "Unknown error" );
    }
  private:
    ncrystal_vapi_error_t m_err;
  };

  void toArray( const NCrystalVAPIType2V1::Neutron& n, double * a )
  {
    a[0] = n.ekin; a[1] = n.ux; a[2] = n.uy; a[3] = n.uz;
  }

  //Context of the trampolines below, which capture exceptions, since they
  //must not propagate through the C API:
  struct RNGContext {
    const std::function<double()> * fct;
    std::exception_ptr exception;
  };

  struct NatAbundContext {
    const NCrystalVAPIType2V1::NaturalAbundances * fct;
    std::exception_ptr exception;
  };

  ////////////////////////////////////////////////////////////////////////////
  // Locating and loading the NCrystal library.

  std::string getEnv( const char * name )
  {
    const char * v = std::getenv( name );
    return v ? std::string( v ) : std::string();
  }

  void runNCrystalConfig( std::string& shlibpath, std::string& ns )
  {
#ifdef NCVAPI_WINDOWS
    FILE* pipe = _popen( "ncrystal-config --show shlibpath namespace", "r" );
#else
    FILE* pipe = popen( "ncrystal-config --show shlibpath namespace"
                        " 2>/dev/null", "r" );
#endif
    if ( !pipe )
      throw std::runtime_error( "Could not run the ncrystal-config command" );
    auto readLine = [pipe]( std::string& tgt ) -> bool
    {
      //Read line, and discard trailing whitespace (including newlines):
      char buffer[4096];
      if ( std::fgets( buffer, sizeof(buffer), pipe ) == nullptr )
        return false;
      tgt = buffer;
      while ( !tgt.empty()
              && std::isspace( static_cast<unsigned char>( tgt.back() ) ) )
        tgt.pop_back();
      return true;
    };
    const bool ok = readLine( shlibpath ) && !shlibpath.empty()
      && readLine( ns );
#ifdef NCVAPI_WINDOWS
    const int ec = _pclose( pipe );
#else
    const int ec = pclose( pipe );
#endif
    if ( !ok || ec != 0 )
      throw std::runtime_error( "Could not locate NCrystal: the command"
                                " \"ncrystal-config --show shlibpath"
                                " namespace\" failed (is NCrystal"
                                " installed?)" );
  }

  const API * loadRawAPI()
  {
    //The library and its symbol namespace, from NCRYSTAL_LIB (as in
    //NCrystal's Python API) or else from ncrystal-config:
    std::string shlibpath = getEnv( "NCRYSTAL_LIB" );
    std::string ns;
    if ( !shlibpath.empty() )
      ns = getEnv( "NCRYSTAL_LIB_NAMESPACE_PROTECTION" );
    else
      runNCrystalConfig( shlibpath, ns );

#ifdef NCVAPI_WINDOWS
    HMODULE handle = LoadLibraryA( shlibpath.c_str() );
    if ( !handle )
      throw std::runtime_error( "Could not load the NCrystal library: "
                                + shlibpath );
#else
    dlerror();//clear previous errors
    void * handle = dlopen( shlibpath.c_str(), RTLD_LOCAL | RTLD_LAZY );
    if ( !handle ) {
      const char * e = dlerror();
      throw std::runtime_error( "Could not load the NCrystal library: "
                                + shlibpath + ( e ? std::string(" (")
                                                + e + ")" : "" ) );
    }
#endif

    const std::string symbol = "ncrystal" + ns + "_access_virtual_c_api";
#ifdef NCVAPI_WINDOWS
    void * addr = reinterpret_cast<void*>( GetProcAddress( handle,
                                                           symbol.c_str() ) );
#else
    void * addr = dlsym( handle, symbol.c_str() );
#endif
    if ( !addr )
      throw std::runtime_error( "The NCrystal library " + shlibpath
                                + " does not provide the virtual C API (the"
                                " function " + symbol + " was not found)."
                                " A newer version of NCrystal is needed." );

    using AccessFct = const void * (*)( unsigned long );
    auto access = reinterpret_cast<AccessFct>( addr );
    auto api = static_cast<const API*>( access( 2001 ) );
    if ( !api )
      throw std::runtime_error( "The NCrystal library " + shlibpath
                                + " does not provide type 2 version 1 of the"
                                " virtual C API (interface id 2001). A newer"
                                " version of NCrystal is needed." );
    if ( api->interface_id != 2001 || api->struct_size != sizeof( API ) )
      throw std::runtime_error( "Inconsistent definition of type 2 version 1"
                                " of NCrystal's virtual C API (was the copy"
                                " in NCrystalVAPIType2V1.cc modified?)" );
    return api;
  }

}

//The trampolines get C language linkage, as the function pointer types of
//the C API:
extern "C" {
  static double ncvapi_rng_trampoline( void * p )
  {
    auto& ctx = *static_cast<RNGContext*>( p );
    try {
      return ( *ctx.fct )();
    } catch ( ... ) {
      //Make NCrystal fail, and rethrow the exception afterwards:
      ctx.exception = std::current_exception();
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  static size_t ncvapi_natabund_trampoline( void * p, unsigned long Z,
                                            unsigned long * A,
                                            double * fraction,
                                            size_t capacity )
  {
    auto& ctx = *static_cast<NatAbundContext*>( p );
    try {
      auto isotopes = ( *ctx.fct )( Z );
      if ( isotopes.size() > capacity )
        throw std::runtime_error( "Too many isotopes returned by the"
                                  " natural abundance function" );
      for ( std::size_t i = 0; i < isotopes.size(); ++i ) {
        A[i] = isotopes[i].first;
        fraction[i] = isotopes[i].second;
      }
      return isotopes.size();
    } catch ( ... ) {
      ctx.exception = std::current_exception();
      return 0;
    }
  }
}

//Access to the C API struct of an NCrystalVAPIType2V1 object:
struct NCrystalVAPIType2V1::Access {
  static const API& api( const NCrystalVAPIType2V1& o )
  {
    return *static_cast<const API*>( o.m_raw );
  }
};

NCrystalVAPIType2V1::Error::Error( const std::string& type,
                                   const std::string& message )
  : std::runtime_error( "NCrystal error (" + type + "): " + message ),
    m_type( type )
{
}

std::shared_ptr<const NCrystalVAPIType2V1> NCrystalVAPIType2V1::load()
{
  static std::mutex mtx;
  static std::shared_ptr<const NCrystalVAPIType2V1> instance;
  std::lock_guard<std::mutex> lock( mtx );
  if ( !instance )
    instance = std::shared_ptr<const NCrystalVAPIType2V1>
      ( new NCrystalVAPIType2V1( loadRawAPI() ) );
  return instance;
}

NCrystalVAPIType2V1::Info
NCrystalVAPIType2V1::createInfo( const std::string& cfgstr ) const
{
  ErrorCollector err;
  auto h = Access::api( *this ).create_info( cfgstr.c_str(), err.ptr() );
  err.throwIf( h == nullptr );
  return Info( shared_from_this(), h );
}

NCrystalVAPIType2V1::Scatter
NCrystalVAPIType2V1::createScatter( const std::string& cfgstr ) const
{
  ErrorCollector err;
  auto h = Access::api( *this ).create_scatter( cfgstr.c_str(), err.ptr() );
  err.throwIf( h == nullptr );
  return Scatter( shared_from_this(), h );
}

NCrystalVAPIType2V1::Absorption
NCrystalVAPIType2V1::createAbsorption( const std::string& cfgstr ) const
{
  ErrorCollector err;
  auto h = Access::api( *this ).create_absorption( cfgstr.c_str(),
                                                   err.ptr() );
  err.throwIf( h == nullptr );
  return Absorption( shared_from_this(), h );
}

void NCrystalVAPIType2V1::clearCaches() const
{
  Access::api( *this ).clear_caches();
}

namespace {
  //The handles are stored as void pointers in the classes:
  struct ncrystal_vapi_t2v1_info * infoH( void * h )
  {
    return static_cast<struct ncrystal_vapi_t2v1_info *>( h );
  }
  struct ncrystal_vapi_t2v1_scatter * scatterH( void * h )
  {
    return static_cast<struct ncrystal_vapi_t2v1_scatter *>( h );
  }
  struct ncrystal_vapi_t2v1_absorption * absorptionH( void * h )
  {
    return static_cast<struct ncrystal_vapi_t2v1_absorption *>( h );
  }
}

////////////////////////////////////////////////////////////////////////////////
// Info

NCrystalVAPIType2V1::Info::Info( std::shared_ptr<const NCrystalVAPIType2V1> api,
                                 void * h )
  : m_api( std::move( api ) ), m_h( h )
{
}

NCrystalVAPIType2V1::Info::Info( Info&& o ) noexcept
  : m_api( std::move( o.m_api ) ), m_h( o.m_h )
{
  o.m_h = nullptr;
}

NCrystalVAPIType2V1::Info&
NCrystalVAPIType2V1::Info::operator=( Info&& o ) noexcept
{
  if ( this != &o ) {
    if ( m_h )
      Access::api( *m_api ).deallocate_info( infoH( m_h ) );
    m_api = std::move( o.m_api );
    m_h = o.m_h;
    o.m_h = nullptr;
  }
  return *this;
}

NCrystalVAPIType2V1::Info::~Info()
{
  if ( m_h )
    Access::api( *m_api ).deallocate_info( infoH( m_h ) );
}

std::uint64_t NCrystalVAPIType2V1::Info::uniqueID() const
{
  return static_cast<std::uint64_t>
    ( Access::api( *m_api ).info_unique_id( infoH( m_h ) ) );
}

double NCrystalVAPIType2V1::Info::density() const
{
  return Access::api( *m_api ).info_density( infoH( m_h ) );
}

double NCrystalVAPIType2V1::Info::numberDensity() const
{
  return Access::api( *m_api ).info_number_density( infoH( m_h ) );
}

double NCrystalVAPIType2V1::Info::temperature() const
{
  ErrorCollector err;
  double t = 0.0;
  err.throwIf( Access::api( *m_api ).info_temperature( infoH( m_h ), &t,
                                                       err.ptr() ) != 0 );
  return t;
}

std::vector<NCrystalVAPIType2V1::Component>
NCrystalVAPIType2V1::Info::composition( bool preferNaturalElements,
                                        const NaturalAbundances& natab ) const
{
  const API& api = Access::api( *m_api );
  NatAbundContext ctx{ &natab, nullptr };
  auto call = [&]( unsigned long * Z, unsigned long * A, double * fraction,
                   std::size_t capacity )
  {
    ErrorCollector err;
    const std::size_t n
      = api.info_composition( infoH( m_h ), preferNaturalElements ? 1 : 0,
                              natab ? ncvapi_natabund_trampoline : nullptr,
                              &ctx, Z, A, fraction, capacity, err.ptr() );
    if ( ctx.exception )
      std::rethrow_exception( ctx.exception );
    err.throwIf( n == 0 );
    return n;
  };
  //Usually a single call is enough, otherwise call again with the needed
  //capacity:
  std::vector<unsigned long> Z( 64 ), A( 64 );
  std::vector<double> fraction( 64 );
  std::size_t n = call( Z.data(), A.data(), fraction.data(), Z.size() );
  if ( n > Z.size() ) {
    Z.resize( n );
    A.resize( n );
    fraction.resize( n );
    if ( call( Z.data(), A.data(), fraction.data(), n ) != n )
      throw std::runtime_error( "Inconsistent composition from NCrystal" );
  }
  std::vector<Component> res;
  res.reserve( n );
  for ( std::size_t i = 0; i < n; ++i )
    res.push_back( Component{ Z[i], A[i], fraction[i] } );
  return res;
}

////////////////////////////////////////////////////////////////////////////////
// Scatter

NCrystalVAPIType2V1::Scatter::Scatter
( std::shared_ptr<const NCrystalVAPIType2V1> api, void * h )
  : m_api( std::move( api ) ), m_h( h )
{
}

NCrystalVAPIType2V1::Scatter::Scatter( Scatter&& o ) noexcept
  : m_api( std::move( o.m_api ) ), m_h( o.m_h )
{
  o.m_h = nullptr;
}

NCrystalVAPIType2V1::Scatter&
NCrystalVAPIType2V1::Scatter::operator=( Scatter&& o ) noexcept
{
  if ( this != &o ) {
    if ( m_h )
      Access::api( *m_api ).deallocate_scatter( scatterH( m_h ) );
    m_api = std::move( o.m_api );
    m_h = o.m_h;
    o.m_h = nullptr;
  }
  return *this;
}

NCrystalVAPIType2V1::Scatter::~Scatter()
{
  if ( m_h )
    Access::api( *m_api ).deallocate_scatter( scatterH( m_h ) );
}

NCrystalVAPIType2V1::Scatter NCrystalVAPIType2V1::Scatter::clone() const
{
  ErrorCollector err;
  auto h = Access::api( *m_api ).clone_scatter( scatterH( m_h ), err.ptr() );
  err.throwIf( h == nullptr );
  return Scatter( m_api, h );
}

bool NCrystalVAPIType2V1::Scatter::isOriented() const
{
  return Access::api( *m_api ).scatter_is_oriented( scatterH( m_h ) ) != 0;
}

double NCrystalVAPIType2V1::Scatter::crossSection( const Neutron& n ) const
{
  ErrorCollector err;
  double a[4];
  toArray( n, a );
  double xs = 0.0;
  err.throwIf( Access::api( *m_api ).scatter_cross_section
               ( scatterH( m_h ), a, &xs, err.ptr() ) != 0 );
  return xs;
}

void NCrystalVAPIType2V1::Scatter::sampleScatter( double (*rng)( void * ),
                                                  void * rng_state,
                                                  Neutron& n ) const
{
  ErrorCollector err;
  double a[4];
  toArray( n, a );
  err.throwIf( Access::api( *m_api ).sample_scatter
               ( scatterH( m_h ), rng, rng_state, a, err.ptr() ) != 0 );
  n.ekin = a[0]; n.ux = a[1]; n.uy = a[2]; n.uz = a[3];
}

void
NCrystalVAPIType2V1::Scatter::sampleScatter( const std::function<double()>& rng,
                                             Neutron& n ) const
{
  RNGContext ctx{ &rng, nullptr };
  ErrorCollector err;
  double a[4];
  toArray( n, a );
  const int ret = Access::api( *m_api ).sample_scatter
    ( scatterH( m_h ), ncvapi_rng_trampoline, &ctx, a, err.ptr() );
  if ( ctx.exception )
    std::rethrow_exception( ctx.exception );
  err.throwIf( ret != 0 );
  n.ekin = a[0]; n.ux = a[1]; n.uy = a[2]; n.uz = a[3];
}

////////////////////////////////////////////////////////////////////////////////
// Absorption

NCrystalVAPIType2V1::Absorption::Absorption
( std::shared_ptr<const NCrystalVAPIType2V1> api, void * h )
  : m_api( std::move( api ) ), m_h( h )
{
}

NCrystalVAPIType2V1::Absorption::Absorption( Absorption&& o ) noexcept
  : m_api( std::move( o.m_api ) ), m_h( o.m_h )
{
  o.m_h = nullptr;
}

NCrystalVAPIType2V1::Absorption&
NCrystalVAPIType2V1::Absorption::operator=( Absorption&& o ) noexcept
{
  if ( this != &o ) {
    if ( m_h )
      Access::api( *m_api ).deallocate_absorption( absorptionH( m_h ) );
    m_api = std::move( o.m_api );
    m_h = o.m_h;
    o.m_h = nullptr;
  }
  return *this;
}

NCrystalVAPIType2V1::Absorption::~Absorption()
{
  if ( m_h )
    Access::api( *m_api ).deallocate_absorption( absorptionH( m_h ) );
}

NCrystalVAPIType2V1::Absorption NCrystalVAPIType2V1::Absorption::clone() const
{
  ErrorCollector err;
  auto h = Access::api( *m_api ).clone_absorption( absorptionH( m_h ),
                                                   err.ptr() );
  err.throwIf( h == nullptr );
  return Absorption( m_api, h );
}

double NCrystalVAPIType2V1::Absorption::crossSection( const Neutron& n ) const
{
  ErrorCollector err;
  double a[4];
  toArray( n, a );
  double xs = 0.0;
  err.throwIf( Access::api( *m_api ).absorption_cross_section
               ( absorptionH( m_h ), a, &xs, err.ptr() ) != 0 );
  return xs;
}
