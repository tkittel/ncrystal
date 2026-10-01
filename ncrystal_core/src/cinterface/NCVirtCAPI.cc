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
#include <cstdio>
#include <cstdlib>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace VirtCAPI {
    namespace {

      //Handle of the type1_v2 API (the C type is only declared, and handles
      //are reinterpret_cast'ed to and from this struct):
      struct T1V2Scatter {
        ProcImpl::ProcPtr proc;
        CachePtr cache;
      };

      T1V2Scatter& extract( ncrystal_vapi_t1v2_scatter_t * h )
      {
        return *reinterpret_cast<T1V2Scatter*>( h );
      }

      const T1V2Scatter& extract( const ncrystal_vapi_t1v2_scatter_t * h )
      {
        return *reinterpret_cast<const T1V2Scatter*>( h );
      }

      ncrystal_vapi_t1v2_scatter_t * wrap( T1V2Scatter * sp )
      {
        return reinterpret_cast<ncrystal_vapi_t1v2_scatter_t*>( sp );
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

    }
  }
}

namespace NCVC = NCrystal::VirtCAPI;

extern "C" {

  static ncrystal_vapi_t1v2_scatter_t *
  ncrystal_vapi_t1v2_create_scatter( const char * cfgstr,
                                     ncrystal_vapi_error_t * err )
  {
    NCVC::T1V2Scatter * res = nullptr;
    NCVC::guarded( err, [&res,cfgstr]()
    {
      if ( !cfgstr )
        NCRYSTAL_THROW( BadInput, "Invalid (null) cfgstr" );
      auto proc = NC::FactImpl::createScatter( cfgstr );
      res = new NCVC::T1V2Scatter{ std::move( proc ), nullptr };
    } );
    return NCVC::wrap( res );
  }

  static ncrystal_vapi_t1v2_scatter_t *
  ncrystal_vapi_t1v2_clone_scatter( const ncrystal_vapi_t1v2_scatter_t * h,
                                    ncrystal_vapi_error_t * err )
  {
    NCVC::T1V2Scatter * res = nullptr;
    NCVC::guarded( err, [&res,h]()
    {
      if ( !h )
        NCRYSTAL_THROW( BadInput, "Invalid (null) scatter handle" );
      res = new NCVC::T1V2Scatter{ NCVC::extract( h ).proc, nullptr };
    } );
    return NCVC::wrap( res );
  }

  static void
  ncrystal_vapi_t1v2_deallocate_scatter( ncrystal_vapi_t1v2_scatter_t * h )
  {
    delete reinterpret_cast<NCVC::T1V2Scatter*>( h );
  }

  static int
  ncrystal_vapi_t1v2_cross_section( ncrystal_vapi_t1v2_scatter_t * h,
                                    const double * n,
                                    double * xsect,
                                    ncrystal_vapi_error_t * err )
  {
    return NCVC::guarded( err, [h,n,xsect]()
    {
      if ( !h || !n || !xsect )
        NCRYSTAL_THROW( BadInput, "Invalid (null) argument" );
      NCVC::checkNeutron( n );
      auto& sc = NCVC::extract( h );
      *xsect = sc.proc->crossSection( sc.cache,
                                      NC::NeutronEnergy{ n[0] },
                                      NC::NeutronDirection( n[1], n[2], n[3] )
                                      ).dbl();
    } );
  }

  static int
  ncrystal_vapi_t1v2_sample_scatter( ncrystal_vapi_t1v2_scatter_t * h,
                                     double (*rng)( void * ),
                                     void * rng_state,
                                     double * n,
                                     ncrystal_vapi_error_t * err )
  {
    return NCVC::guarded( err, [h,rng,rng_state,n]()
    {
      if ( !h || !rng || !n )
        NCRYSTAL_THROW( BadInput, "Invalid (null) argument" );
      NCVC::checkNeutron( n );
      auto& sc = NCVC::extract( h );
      NCVC::RNGFromFct rngwrapper( rng, rng_state );
      auto out = sc.proc->sampleScatter( sc.cache, rngwrapper,
                                         NC::NeutronEnergy{ n[0] },
                                         NC::NeutronDirection( n[1],
                                                               n[2],
                                                               n[3] ) );
      n[0] = out.ekin.dbl();
      n[1] = out.direction[0];
      n[2] = out.direction[1];
      n[3] = out.direction[2];
    } );
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
}

const void * ncrystal_access_virtual_c_api( unsigned interface_id )
{
  if ( interface_id == 1002 )
    return &s_vapi_type1_v2;
  return nullptr;
}
