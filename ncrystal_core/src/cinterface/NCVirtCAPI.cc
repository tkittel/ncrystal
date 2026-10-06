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

//Implementation of the virtual C APIs defined in NCrystal/virtualapi/ncvirtapi.h

#include "NCrystal/cinterface/ncrystal.h"
#include "NCrystal/virtualapi/ncvirtapi.h"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include "NCrystal/misc/NCCompositionUtils.hh"
#include "NCrystal/core/NCMem.hh"
#include <cstdio>
#include <cstdlib>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace VirtCAPI {
    namespace {

      //Scatter and absorption handles (the C types are only declared, and
      //handles are reinterpret_cast'ed to and from these structs):
      struct ProcHandle {
        ProcImpl::ProcPtr proc;
        CachePtr cache;
      };

      //Info handles:
      struct InfoHandle {
        InfoPtr info;
      };

      template<class TCHandle>
      ProcHandle& extract( TCHandle * h )
      {
        return *reinterpret_cast<ProcHandle*>( h );
      }

      template<class TCHandle>
      const ProcHandle& extract( const TCHandle * h )
      {
        return *reinterpret_cast<const ProcHandle*>( h );
      }

      const InfoHandle& extractInfo( const struct ncrystal_vapi_t2v1_info * h )
      {
        nc_assert_always( h != nullptr );
        return *reinterpret_cast<const InfoHandle*>( h );
      }

      template<class TCHandle>
      TCHandle * wrap( ProcHandle * sp )
      {
        return reinterpret_cast<TCHandle*>( sp );
      }

      void copyStr( char * dest, std::size_t destsize, const char * src )
      {
        nc_assert( destsize > 0 );
        std::size_t n = std::min( std::strlen( src ), destsize - 1 );
        std::memcpy( dest, src, n );
        dest[n] = '\0';
      }

      void setError( ncrystal_vapi_error_t * err,
                     const char * type, const char * msg ) noexcept
      {
        if ( !err ) {
          std::fprintf( stderr, "NCrystal ERROR [%s]: %s\n", type, msg );
          std::fflush( stderr );
          std::exit( 1 );
        }
        err->code = 1;
        copyStr( err->type, sizeof(err->type), type );
        copyStr( err->message, sizeof(err->message), msg );
      }

      //Run the function, converting any exception to an error (returns 0 on
      //success and 1 on failure):
      template<class TFct>
      int guarded( ncrystal_vapi_error_t * err, TFct&& fct ) noexcept
      {
        try {
          fct();
          return 0;
        } catch ( Error::Exception& e ) {
          setError( err, e.getTypeName(), e.what() );
        } catch ( std::exception& e ) {
          setError( err, "std::exception", e.what() );
        } catch ( ... ) {
          setError( err, "UnknownError", "Unknown error" );
        }
        return 1;
      }

      class RNGFromFct final : public RNGStream {
      public:
        RNGFromFct( double (*fct)( void * ), void * state )
          : m_fct( fct ), m_state( state ) {}
      protected:
        double actualGenerate() override
        {
          //NCrystal needs numbers in (0,1]:
          double v = m_fct( m_state );
          if ( v > 0.0 && v <= 1.0 )
            return v;
          if ( v == 0.0 )
            return std::numeric_limits<double>::min();
          NCRYSTAL_THROW2( BadInput, "Random number function returned"
                           " invalid value: " << v );
        }
      private:
        double (*m_fct)( void * );
        void * m_state;
      };

      void checkNeutron( const double * n )
      {
        if ( !( n[0] >= 0.0 && std::isfinite( n[0] ) ) )
          NCRYSTAL_THROW2( BadInput, "Invalid neutron energy: " << n[0] );
        const double u2 = n[1]*n[1] + n[2]*n[2] + n[3]*n[3];
        if ( !( u2 > 0.0 && std::isfinite( u2 ) ) )
          NCRYSTAL_THROW2( BadInput, "Invalid neutron direction: ("
                           << n[1] << ", " << n[2] << ", " << n[3] << ")" );
      }

      //Implementations shared by the scatter and absorption handles of the
      //different APIs:

      template<class TCHandle>
      TCHandle * createProc( const char * cfgstr, bool scatter,
                             ncrystal_vapi_error_t * err )
      {
        ProcHandle * res = nullptr;
        guarded( err, [&res,cfgstr,scatter]()
        {
          if ( !cfgstr )
            NCRYSTAL_THROW( BadInput, "Invalid (null) cfgstr" );
          auto proc = ( scatter
                        ? FactImpl::createScatter( cfgstr )
                        : FactImpl::createAbsorption( cfgstr ) );
          res = new ProcHandle{ std::move( proc ), nullptr };
        } );
        return wrap<TCHandle>( res );
      }

      template<class TCHandle>
      TCHandle * cloneProc( const TCHandle * h, bool scatter,
                            ncrystal_vapi_error_t * err )
      {
        ProcHandle * res = nullptr;
        guarded( err, [&res,h,scatter]()
        {
          if ( !h )
            NCRYSTAL_THROW2( BadInput, "Invalid (null) "
                             << ( scatter ? "scatter" : "absorption" )
                             << " handle" );
          res = new ProcHandle{ extract( h ).proc, nullptr };
        } );
        return wrap<TCHandle>( res );
      }

      template<class TCHandle>
      void deallocateProc( TCHandle * h )
      {
        delete reinterpret_cast<ProcHandle*>( h );
      }

      template<class TCHandle>
      int crossSection( TCHandle * h, const double * n, double * xsect,
                        ncrystal_vapi_error_t * err )
      {
        return guarded( err, [h,n,xsect]()
        {
          if ( !h || !n || !xsect )
            NCRYSTAL_THROW( BadInput, "Invalid (null) argument" );
          checkNeutron( n );
          auto& ph = extract( h );
          *xsect = ph.proc->crossSection( ph.cache,
                                          NeutronEnergy{ n[0] },
                                          NeutronDirection( n[1], n[2], n[3] )
                                          ).dbl();
        } );
      }

      template<class TCHandle>
      int sampleScatter( TCHandle * h, double (*rng)( void * ),
                         void * rng_state, double * n,
                         ncrystal_vapi_error_t * err )
      {
        return guarded( err, [h,rng,rng_state,n]()
        {
          if ( !h || !rng || !n )
            NCRYSTAL_THROW( BadInput, "Invalid (null) argument" );
          checkNeutron( n );
          auto& ph = extract( h );
          RNGFromFct rngwrapper( rng, rng_state );
          auto out = ph.proc->sampleScatter( ph.cache, rngwrapper,
                                             NeutronEnergy{ n[0] },
                                             NeutronDirection( n[1],
                                                               n[2],
                                                               n[3] ) );
          n[0] = out.ekin.dbl();
          n[1] = out.direction[0];
          n[2] = out.direction[1];
          n[3] = out.direction[2];
        } );
      }

      //The type of natabund is a template parameter, since the function
      //pointer type of the C API has C language linkage:
      template<class TNatAbund>
      std::size_t composition( const struct ncrystal_vapi_t2v1_info * h,
                               int prefer_natural_elements,
                               TNatAbund natabund,
                               void * natabund_state,
                               unsigned long * resZ,
                               unsigned long * resA,
                               double * resFrac,
                               std::size_t capacity,
                               ncrystal_vapi_error_t * err )
      {
        std::size_t res = 0;
        guarded( err, [&]()
        {
          if ( !h )
            NCRYSTAL_THROW( BadInput, "Invalid (null) info handle" );
          if ( capacity > 0 && !( resZ && resA && resFrac ) )
            NCRYSTAL_THROW( BadInput, "Invalid (null) output array" );
          CompositionUtils::NaturalAbundanceProvider natprov{ nullptr };
          if ( natabund ) {
            natprov = [natabund,natabund_state]( unsigned Z )
            {
              constexpr std::size_t bufsize = 128;
              unsigned long bufA[bufsize];
              double bufFrac[bufsize];
              const std::size_t niso = natabund( natabund_state, Z,
                                                 bufA, bufFrac, bufsize );
              if ( niso > bufsize )
                NCRYSTAL_THROW2( BadInput, "Natural abundance function"
                                 " returned too many isotopes for Z="<<Z );
              //NB: niso==0 means that the abundances are not known, which
              //translates to an empty result:
              std::vector<std::pair<unsigned,double>> result;
              for ( std::size_t i = 0; i < niso; ++i ) {
                if ( bufFrac[i] == 0.0 )
                  continue;
                if ( bufA[i] < Z || bufA[i] > 999 )
                  NCRYSTAL_THROW2( BadInput, "Invalid (Z,A) value returned"
                                   " from natural abundance function: Z="
                                   <<Z<<", A="<<bufA[i] );
                if ( !( bufFrac[i] > 0.0 && bufFrac[i] <= 1.0 ) )
                  NCRYSTAL_THROW2( BadInput, "Invalid fraction returned from"
                                   " natural abundance function: "
                                   <<bufFrac[i] );
                result.emplace_back( static_cast<unsigned>( bufA[i] ),
                                     bufFrac[i] );
              }
              return result;
            };
          }
          auto bd = CompositionUtils::createFullBreakdown
            ( extractInfo( h ).info->getComposition(), natprov,
              ( prefer_natural_elements
                ? CompositionUtils::PreferNaturalElements
                : CompositionUtils::ForceIsotopes ) );
          //Count first, so nothing is written on failure:
          std::size_t n = 0;
          for ( auto& e : bd )
            n += e.second.size();
          std::size_t i = 0;
          for ( auto& e : bd ) {
            for ( auto& af : e.second ) {
              if ( i < capacity ) {
                resZ[i] = e.first;
                resA[i] = af.first;
                resFrac[i] = af.second;
              }
              ++i;
            }
          }
          res = n;
        } );
        return res;
      }

    }
  }
}

namespace NCVC = NCrystal::VirtCAPI;

extern "C" {

  static struct ncrystal_vapi_t1v2_scatter *
  ncrystal_vapi_t1v2_create_scatter( const char * cfgstr,
                                     ncrystal_vapi_error_t * err )
  {
    return NCVC::createProc<struct ncrystal_vapi_t1v2_scatter>( cfgstr, true,
                                                                err );
  }

  static struct ncrystal_vapi_t1v2_scatter *
  ncrystal_vapi_t1v2_clone_scatter( const struct ncrystal_vapi_t1v2_scatter * h,
                                    ncrystal_vapi_error_t * err )
  {
    return NCVC::cloneProc( h, true, err );
  }

  static void
  ncrystal_vapi_t1v2_deallocate_scatter( struct ncrystal_vapi_t1v2_scatter * h )
  {
    NCVC::deallocateProc( h );
  }

  static int
  ncrystal_vapi_t1v2_cross_section( struct ncrystal_vapi_t1v2_scatter * h,
                                    const double * n,
                                    double * xsect,
                                    ncrystal_vapi_error_t * err )
  {
    return NCVC::crossSection( h, n, xsect, err );
  }

  static int
  ncrystal_vapi_t1v2_sample_scatter( struct ncrystal_vapi_t1v2_scatter * h,
                                     double (*rng)( void * ),
                                     void * rng_state,
                                     double * n,
                                     ncrystal_vapi_error_t * err )
  {
    return NCVC::sampleScatter( h, rng, rng_state, n, err );
  }

  static struct ncrystal_vapi_t2v1_info *
  ncrystal_vapi_t2v1_create_info( const char * cfgstr,
                                  ncrystal_vapi_error_t * err )
  {
    NCVC::InfoHandle * res = nullptr;
    NCVC::guarded( err, [&res,cfgstr]()
    {
      if ( !cfgstr )
        NCRYSTAL_THROW( BadInput, "Invalid (null) cfgstr" );
      res = new NCVC::InfoHandle{ NC::FactImpl::createInfo( cfgstr ) };
    } );
    return reinterpret_cast<struct ncrystal_vapi_t2v1_info*>( res );
  }

  static void
  ncrystal_vapi_t2v1_deallocate_info( struct ncrystal_vapi_t2v1_info * h )
  {
    delete reinterpret_cast<NCVC::InfoHandle*>( h );
  }

  static double
  ncrystal_vapi_t2v1_info_unique_id( const struct ncrystal_vapi_t2v1_info * h )
  {
    const auto uid = NCVC::extractInfo( h ).info->getUniqueID().value;
    //The ids are counters, and can not realistically exceed 2^53:
    nc_assert_always( uid <= ( std::uint64_t(1) << 53 ) );
    return static_cast<double>( uid );
  }

  static double
  ncrystal_vapi_t2v1_info_density( const struct ncrystal_vapi_t2v1_info * h )
  {
    return NCVC::extractInfo( h ).info->getDensity().dbl();
  }

  static double
  ncrystal_vapi_t2v1_info_number_density
  ( const struct ncrystal_vapi_t2v1_info * h )
  {
    return NCVC::extractInfo( h ).info->getNumberDensity().dbl();
  }

  static int
  ncrystal_vapi_t2v1_info_temperature( const struct ncrystal_vapi_t2v1_info * h,
                                       double * temperature,
                                       ncrystal_vapi_error_t * err )
  {
    return NCVC::guarded( err, [h,temperature]()
    {
      if ( !h || !temperature )
        NCRYSTAL_THROW( BadInput, "Invalid (null) argument" );
      auto& info = *NCVC::extractInfo( h ).info;
      if ( !info.hasTemperature() )
        NCRYSTAL_THROW( MissingInfo, "The material does not have a"
                        " temperature" );
      *temperature = info.getTemperature().dbl();
    } );
  }

  static size_t
  ncrystal_vapi_t2v1_info_composition
  ( const struct ncrystal_vapi_t2v1_info * h,
    int prefer_natural_elements,
    size_t (*natabund)( void *, unsigned long, unsigned long *, double *,
                        size_t ),
    void * natabund_state,
    unsigned long * Z,
    unsigned long * A,
    double * fraction,
    size_t capacity,
    ncrystal_vapi_error_t * err )
  {
    return NCVC::composition( h, prefer_natural_elements, natabund,
                              natabund_state, Z, A, fraction, capacity, err );
  }

  static int
  ncrystal_vapi_t2v1_scatter_is_oriented
  ( const struct ncrystal_vapi_t2v1_scatter * h )
  {
    nc_assert_always( h != nullptr );
    return NCVC::extract( h ).proc->isOriented() ? 1 : 0;
  }

  static struct ncrystal_vapi_t2v1_scatter *
  ncrystal_vapi_t2v1_create_scatter( const char * cfgstr,
                                     ncrystal_vapi_error_t * err )
  {
    return NCVC::createProc<struct ncrystal_vapi_t2v1_scatter>( cfgstr, true,
                                                                err );
  }

  static struct ncrystal_vapi_t2v1_scatter *
  ncrystal_vapi_t2v1_clone_scatter( const struct ncrystal_vapi_t2v1_scatter * h,
                                    ncrystal_vapi_error_t * err )
  {
    return NCVC::cloneProc( h, true, err );
  }

  static void
  ncrystal_vapi_t2v1_deallocate_scatter( struct ncrystal_vapi_t2v1_scatter * h )
  {
    NCVC::deallocateProc( h );
  }

  static int
  ncrystal_vapi_t2v1_scatter_cross_section
  ( struct ncrystal_vapi_t2v1_scatter * h,
    const double * n,
    double * xsect,
    ncrystal_vapi_error_t * err )
  {
    return NCVC::crossSection( h, n, xsect, err );
  }

  static int
  ncrystal_vapi_t2v1_sample_scatter( struct ncrystal_vapi_t2v1_scatter * h,
                                     double (*rng)( void * ),
                                     void * rng_state,
                                     double * n,
                                     ncrystal_vapi_error_t * err )
  {
    return NCVC::sampleScatter( h, rng, rng_state, n, err );
  }

  static struct ncrystal_vapi_t2v1_absorption *
  ncrystal_vapi_t2v1_create_absorption( const char * cfgstr,
                                        ncrystal_vapi_error_t * err )
  {
    return NCVC::createProc<struct ncrystal_vapi_t2v1_absorption>( cfgstr,
                                                                   false,
                                                                   err );
  }

  static struct ncrystal_vapi_t2v1_absorption *
  ncrystal_vapi_t2v1_clone_absorption
  ( const struct ncrystal_vapi_t2v1_absorption * h,
    ncrystal_vapi_error_t * err )
  {
    return NCVC::cloneProc( h, false, err );
  }

  static void
  ncrystal_vapi_t2v1_deallocate_absorption
  ( struct ncrystal_vapi_t2v1_absorption * h )
  {
    NCVC::deallocateProc( h );
  }

  static int
  ncrystal_vapi_t2v1_absorption_cross_section
  ( struct ncrystal_vapi_t2v1_absorption * h, const double * n, double * xsect,
    ncrystal_vapi_error_t * err )
  {
    return NCVC::crossSection( h, n, xsect, err );
  }

  static void ncrystal_vapi_t2v1_clear_caches()
  {
    NC::clearCaches();
  }

}

namespace {
  const ncrystal_vapi_type1_v2_t s_vapi_type1_v2 = {
    1002,
    sizeof( ncrystal_vapi_type1_v2_t ),
    ncrystal_vapi_t1v2_create_scatter,
    ncrystal_vapi_t1v2_clone_scatter,
    ncrystal_vapi_t1v2_deallocate_scatter,
    ncrystal_vapi_t1v2_cross_section,
    ncrystal_vapi_t1v2_sample_scatter
  };

  const ncrystal_vapi_type2_v1_t s_vapi_type2_v1 = {
    2001,
    sizeof( ncrystal_vapi_type2_v1_t ),
    ncrystal_vapi_t2v1_create_info,
    ncrystal_vapi_t2v1_deallocate_info,
    ncrystal_vapi_t2v1_info_unique_id,
    ncrystal_vapi_t2v1_info_density,
    ncrystal_vapi_t2v1_info_number_density,
    ncrystal_vapi_t2v1_info_temperature,
    ncrystal_vapi_t2v1_info_composition,
    ncrystal_vapi_t2v1_create_scatter,
    ncrystal_vapi_t2v1_clone_scatter,
    ncrystal_vapi_t2v1_deallocate_scatter,
    ncrystal_vapi_t2v1_scatter_is_oriented,
    ncrystal_vapi_t2v1_scatter_cross_section,
    ncrystal_vapi_t2v1_sample_scatter,
    ncrystal_vapi_t2v1_create_absorption,
    ncrystal_vapi_t2v1_clone_absorption,
    ncrystal_vapi_t2v1_deallocate_absorption,
    ncrystal_vapi_t2v1_absorption_cross_section,
    ncrystal_vapi_t2v1_clear_caches
  };
}

const void * ncrystal_access_virtual_c_api( unsigned long interface_id )
{
  if ( interface_id == 1002 )
    return &s_vapi_type1_v2;
  if ( interface_id == 2001 )
    return &s_vapi_type2_v1;
  return nullptr;
}
