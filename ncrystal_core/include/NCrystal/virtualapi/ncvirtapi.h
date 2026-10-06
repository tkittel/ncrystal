#ifndef NCrystal_virtapi_h
#define NCrystal_virtapi_h

/******************************************************************************/
/*                                                                            */
/*  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   */
/*                                                                            */
/*  Copyright 2015-2026 NCrystal developers                                   */
/*                                                                            */
/*  Licensed under the Apache License, Version 2.0 (the "License");           */
/*  you may not use this file except in compliance with the License.          */
/*  You may obtain a copy of the License at                                   */
/*                                                                            */
/*      http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                            */
/*  Unless required by applicable law or agreed to in writing, software       */
/*  distributed under the License is distributed on an "AS IS" BASIS,         */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/*  See the License for the specific language governing permissions and       */
/*  limitations under the License.                                            */
/*                                                                            */
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/* As an addition to the usual NCrystal license pasted above, note that THIS  */
/* PARTICULAR FILE is also placed into the Public Domain, so that it might    */
/* be easily adopted into the code-base of any project wishing to use it:     */
/*                                                                            */
/* This is free and unencumbered software released into the public domain.    */
/*                                                                            */
/* Anyone is free to copy, modify, publish, use, compile, sell, or            */
/* distribute this software, either in source code form or as a compiled      */
/* binary, for any purpose, commercial or non-commercial, and by any          */
/* means.                                                                     */
/*                                                                            */
/* In jurisdictions that recognize copyright laws, the author or authors      */
/* of this software dedicate any and all copyright interest in the            */
/* software to the public domain. We make this dedication for the benefit     */
/* of the public at large and to the detriment of our heirs and               */
/* successors. We intend this dedication to be an overt act of                */
/* relinquishment in perpetuity of all present and future rights to this      */
/* software under copyright law.                                              */
/*                                                                            */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,            */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF         */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.     */
/* IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR          */
/* OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,      */
/* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR      */
/* OTHER DEALINGS IN THE SOFTWARE.                                            */
/*                                                                            */
/* For more information, please refer to <https://unlicense.org/>             */
/*                                                                            */
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/* Virtual C APIs of NCrystal.                                                */
/*                                                                            */
/* These allow applications to use NCrystal without any build-time dependency */
/* on it: the application includes a copy of this file, and at runtime loads  */
/* the NCrystal shared library (e.g. located via "ncrystal-config --show      */
/* shlibpath namespace") and looks up the function:                           */
/*                                                                            */
/*   const void * ncrystal<ns>_access_virtual_c_api( unsigned long id );      */
/*                                                                            */
/* where <ns> is the NCrystal symbol namespace (usually empty). It returns a  */
/* pointer to a struct of function pointers (one struct type per interface    */
/* id, defined below), or NULL if the interface is not available in the       */
/* loaded NCrystal library. The returned struct lives until the program ends. */
/*                                                                            */
/* The interfaces only use plain C types, and are thus independent of the C++ */
/* standard library and compiler used to build NCrystal. No exceptions        */
/* propagate through the functions. Instead, functions which can fail take a  */
/* pointer to a caller-provided ncrystal_vapi_error_t, which is filled in     */
/* only when an error occurs. If that pointer is NULL, NCrystal instead       */
/* prints the error message and terminates the program.                       */
/*                                                                            */
/* Once released, an interface struct never changes (except comments and     */
/* formatting). Changed or added functionality will result in new interface   */
/* structs with new interface ids.                                            */
/*                                                                            */
/* Objects created through an interface (e.g. materials) are referred to by   */
/* handles, which are pointers to incomplete struct types (e.g. "struct       */
/* ncrystal_vapi_t2v1_info *"). These struct types are only declared by their */
/* first use, as the return type of the function creating the objects, so     */
/* each interface is defined by a single struct.                              */
/*                                                                            */
/* To be compatible with C90, the following conventions are used for numbers: */
/*                                                                            */
/* * Counts and array lengths are size_t.                                     */
/* * Unique ids are returned as double, and are guaranteed to be integral     */
/*   values between 0 and 2^53 (so they are represented exactly).             */
/* * Status codes are int, and other integers are unsigned long.              */
/*                                                                            */
/******************************************************************************/

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

  /* Error information. The strings are always null-terminated (and possibly */
  /* truncated):                                                             */
  typedef struct {
    int code;            /* Non-zero when an error occurred.                */
    char type[64];       /* Error type (e.g. "BadInput" or "FileNotFound"). */
    char message[1024];  /* Error message.                                  */
  } ncrystal_vapi_error_t;

  /****************************************************************************/
  /* Type 1, version 2 (interface id 1002, available from NCrystal 4.5.0).    */
  /*                                                                          */
  /* Scattering cross sections and sampling, as needed by e.g. OpenMC (which  */
  /* gets the material composition, density and temperature elsewhere, and    */
  /* handles absorption itself). It supersedes the C++ VirtAPI_Type1_v1.      */
  /*                                                                          */
  /* * Units: barn per atom for cross sections, eV for neutron energies.      */
  /* * The "neutron" parameter is an array (ekin,ux,uy,uz), with (ux,uy,uz)   */
  /*   the unit direction vector. When sampling, it is updated in place.      */
  /* * Each scatter handle owns a cache, which makes repeated calls fast. A   */
  /*   handle must therefore not be used by more than one thread at a time,   */
  /*   and multi-threaded applications should give each thread its own clone */
  /*   (clones share the underlying physics, but have separate caches).       */
  /*   Since caches are specific to a material, each material in each thread  */
  /*   needs its own handle. Creating, cloning and deallocating handles is    */
  /*   thread-safe.                                                           */
  /* * The random number function must return numbers uniformly distributed */
  /*   in [0,1) or (0,1]. The rng_state pointer is passed on to it, and it is */
  /*   only called during the sample_scatter call.                            */
  /* * Functions returning int return 0 on success and non-zero on failure    */
  /*   (in which case the neutron array is not modified). Functions returning */
  /*   handles return NULL on failure.                                        */
  /****************************************************************************/

  typedef struct {
    unsigned long interface_id;/* Always 1002.                              */
    size_t struct_size;/* sizeof(ncrystal_vapi_type1_v2_t) in NCrystal.     */

    struct ncrystal_vapi_t1v2_scatter *
      (*create_scatter)( const char * cfgstr, ncrystal_vapi_error_t * );
    struct ncrystal_vapi_t1v2_scatter *
      (*clone_scatter)( const struct ncrystal_vapi_t1v2_scatter *,
                        ncrystal_vapi_error_t * );
    /* Deallocating a NULL handle is allowed (and does nothing): */
    void (*deallocate_scatter)( struct ncrystal_vapi_t1v2_scatter * );

    int (*cross_section)( struct ncrystal_vapi_t1v2_scatter *,
                          const double * neutron,
                          double * xsect,
                          ncrystal_vapi_error_t * );

    int (*sample_scatter)( struct ncrystal_vapi_t1v2_scatter *,
                           double (*rng)( void * rng_state ),
                           void * rng_state,
                           double * neutron,
                           ncrystal_vapi_error_t * );
  } ncrystal_vapi_type1_v2_t;

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

#ifdef __cplusplus
}
#endif

#endif
