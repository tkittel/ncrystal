#!/usr/bin/env python3

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

# Test the virtual C API (NCrystal/virtualapi/ncvirtapi.h) like a client
# would use it: loading the NCrystal shared library at runtime (here with
# ctypes), looking up ncrystal<ns>_access_virtual_c_api, and calling the
# functions via the returned struct of function pointers. The ctypes
# structures below are written from the description in ncvirtapi.h.

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC
from NCrystalDev._locatelib import get_libpath_and_namespace
import ctypes
import math
import random

class ErrorInfo(ctypes.Structure):
    _fields_ = [ ('code', ctypes.c_int),
                 ('type', ctypes.c_char * 64),
                 ('message', ctypes.c_char * 1024) ]

handle_t = ctypes.c_void_p
err_p = ctypes.POINTER(ErrorInfo)
dbl_p = ctypes.POINTER(ctypes.c_double)
rngfct_t = ctypes.CFUNCTYPE( ctypes.c_double, ctypes.c_void_p )

class APIType1V2(ctypes.Structure):
    _fields_ = [
        ('interface_id', ctypes.c_ulong),
        ('struct_size', ctypes.c_size_t),
        ('create_scatter', ctypes.CFUNCTYPE( handle_t, ctypes.c_char_p, err_p )),
        ('clone_scatter', ctypes.CFUNCTYPE( handle_t, handle_t, err_p )),
        ('deallocate_scatter', ctypes.CFUNCTYPE( None, handle_t )),
        ('cross_section', ctypes.CFUNCTYPE( ctypes.c_int, handle_t, dbl_p,
                                            dbl_p, err_p )),
        ('sample_scatter', ctypes.CFUNCTYPE( ctypes.c_int, handle_t, rngfct_t,
                                             ctypes.c_void_p, dbl_p, err_p )),
    ]

class Component(ctypes.Structure):
    _fields_ = [ ('Z', ctypes.c_ulong),
                 ('A', ctypes.c_ulong),
                 ('fraction', ctypes.c_double) ]

natabund_t = ctypes.CFUNCTYPE( ctypes.c_size_t, ctypes.c_void_p, ctypes.c_ulong,
                               ctypes.POINTER(ctypes.c_ulong), dbl_p,
                               ctypes.c_size_t )
info_t = ctypes.c_void_p

class APIType2V1(ctypes.Structure):
    _fields_ = [
        ('interface_id', ctypes.c_ulong),
        ('struct_size', ctypes.c_size_t),
        ('create_info', ctypes.CFUNCTYPE( info_t, ctypes.c_char_p, err_p )),
        ('deallocate_info', ctypes.CFUNCTYPE( None, info_t )),
        ('info_unique_id', ctypes.CFUNCTYPE( ctypes.c_double, info_t )),
        ('info_density', ctypes.CFUNCTYPE( ctypes.c_double, info_t )),
        ('info_number_density', ctypes.CFUNCTYPE( ctypes.c_double, info_t )),
        ('info_temperature', ctypes.CFUNCTYPE( ctypes.c_int, info_t, dbl_p,
                                               err_p )),
        ('info_composition', ctypes.CFUNCTYPE( ctypes.c_size_t, info_t,
                                               ctypes.c_int, natabund_t,
                                               ctypes.c_void_p,
                                               ctypes.POINTER(Component),
                                               ctypes.c_size_t, err_p )),
        ('create_scatter', ctypes.CFUNCTYPE( handle_t, ctypes.c_char_p, err_p )),
        ('clone_scatter', ctypes.CFUNCTYPE( handle_t, handle_t, err_p )),
        ('deallocate_scatter', ctypes.CFUNCTYPE( None, handle_t )),
        ('scatter_is_oriented', ctypes.CFUNCTYPE( ctypes.c_int, handle_t )),
        ('scatter_cross_section', ctypes.CFUNCTYPE( ctypes.c_int, handle_t,
                                                    dbl_p, dbl_p, err_p )),
        ('sample_scatter', ctypes.CFUNCTYPE( ctypes.c_int, handle_t, rngfct_t,
                                             ctypes.c_void_p, dbl_p, err_p )),
        ('create_absorption', ctypes.CFUNCTYPE( handle_t, ctypes.c_char_p,
                                                err_p )),
        ('clone_absorption', ctypes.CFUNCTYPE( handle_t, handle_t, err_p )),
        ('deallocate_absorption', ctypes.CFUNCTYPE( None, handle_t )),
        ('absorption_cross_section', ctypes.CFUNCTYPE( ctypes.c_int, handle_t,
                                                       dbl_p, dbl_p, err_p )),
        ('clear_caches', ctypes.CFUNCTYPE( None )),
    ]

def load_api( interface_id ):
    libpath, ns = get_libpath_and_namespace()
    lib = ctypes.CDLL( str(libpath) )
    fct = getattr( lib, f'ncrystal{ns}_access_virtual_c_api' )
    fct.restype = ctypes.c_void_p
    fct.argtypes = [ ctypes.c_ulong ]
    return fct( interface_id )

class CountingRNG:
    def __init__( self, seed ):
        self.rng = random.Random( seed )
        self.ncalls = 0
    def __call__( self, _state ):
        self.ncalls += 1
        return self.rng.random()

def neutron( wl, direction ):
    return ( ctypes.c_double * 4 )( NC.wl2ekin(wl), *direction )

def main():
    for bad_id in ( 0, 1001, 1003 ):
        assert load_api( bad_id ) is None
    addr = load_api( 1002 )
    assert addr is not None
    api = APIType1V2.from_address( addr )
    assert api.interface_id == 1002
    assert api.struct_size == ctypes.sizeof( APIType1V2 )
    print('Interface 1002 found, with struct_size matching ctypes: OK')

    err = ErrorInfo()
    for cfg in ( 'stdlib::Al_sg225.ncmat',
                 'stdlib::Polyethylene_CH2.ncmat;temp=250K' ):
        h = api.create_scatter( cfg.encode(), ctypes.byref(err) )
        assert h
        c = api.clone_scatter( h, ctypes.byref(err) )
        assert c
        sc = NC.createScatter( cfg )
        for wl in ( 0.5, 1.8, 4.0, 10.0 ):
            for hh in ( h, c ):
                xs = ctypes.c_double( -1.0 )
                n = neutron( wl, (0.0, 0.0, 1.0) )
                assert api.cross_section( hh, n, ctypes.byref(xs),
                                          ctypes.byref(err) ) == 0
                assert xs.value == sc.xsect( wl = wl ), (xs.value, sc.xsect( wl = wl ))
        print(f'{cfg}: xs(1.8Aa)={sc.xsect(wl=1.8):.6g} barn, same as the'
              ' Python API: OK')

        #Sampling with a Python random number generator:
        rng = CountingRNG( 123 )
        c_rng = rngfct_t( rng )
        for _ in range( 100 ):
            n = neutron( 1.8, (0.0, 0.0, 1.0) )
            assert api.sample_scatter( h, c_rng, None, n,
                                       ctypes.byref(err) ) == 0
            assert n[0] >= 0.0
            assert abs( math.sqrt( n[1]**2 + n[2]**2 + n[3]**2 ) - 1 ) < 1e-12
        assert rng.ncalls > 0
        print(f'{cfg}: 100 samples with a Python RNG: OK')
        api.deallocate_scatter( h )
        api.deallocate_scatter( c )

    #Errors:
    h = api.create_scatter( b'stdlib::doesnotexist.ncmat', ctypes.byref(err) )
    assert not h
    assert err.code != 0 and err.type == b'FileNotFound'
    print(f'Expected error: [{err.type.decode()}] {err.message.decode()}')
    h = api.create_scatter( b'stdlib::Al_sg225.ncmat', ctypes.byref(err) )
    err = ErrorInfo()
    xs = ctypes.c_double( -1.0 )
    n = neutron( 1.8, (0.0, 0.0, 0.0) )
    assert api.cross_section( h, n, ctypes.byref(xs), ctypes.byref(err) ) != 0
    assert err.type == b'BadInput' and xs.value == -1.0
    print(f'Expected error: [{err.type.decode()}] {err.message.decode()}')
    api.deallocate_scatter( h )
    api.deallocate_scatter( None )
    test_type2()
    print('All OK')

def natabund_py( _state, Z, A, fraction, capacity ):
    #Natural abundances from a Python dict:
    data = { 1 : [ (1, 0.99985), (2, 0.00015) ],
             6 : [ (12, 0.99), (13, 0.01) ],
             32 : [ (70, 0.2052), (72, 0.2745), (73, 0.0776), (74, 0.3652),
                    (76, 0.0775) ] }
    assert capacity >= 128
    for i, (a, f) in enumerate( data.get( Z, [] ) ):
        A[i] = a
        fraction[i] = f
    return len( data.get( Z, [] ) )

def test_type2():
    addr = load_api( 2001 )
    assert addr is not None
    api = APIType2V1.from_address( addr )
    assert api.interface_id == 2001
    assert api.struct_size == ctypes.sizeof( APIType2V1 )
    print('Interface 2001 found, with struct_size matching ctypes: OK')
    err = ErrorInfo()
    c_natabund = natabund_t( natabund_py )
    for cfg in ( 'stdlib::Polyethylene_CH2.ncmat;temp=250K',
                 ( 'stdlib::Ge_sg227.ncmat;mos=40arcmin'
                   ';dir1=@crys_hkl:5,1,1@lab:0,0,1'
                   ';dir2=@crys_hkl:0,-1,1@lab:0,1,0' ) ):
        #Info, compared with the Python API:
        h = api.create_info( cfg.encode(), ctypes.byref(err) )
        assert h
        info = NC.createInfo( cfg )
        assert api.info_density( h ) == info.density
        assert api.info_number_density( h ) == info.numberdensity
        t = ctypes.c_double( -1.0 )
        assert api.info_temperature( h, ctypes.byref(t), ctypes.byref(err) ) == 0
        assert t.value == info.getTemperature()
        n = api.info_composition( h, 0, c_natabund, None, None, 0,
                                  ctypes.byref(err) )
        assert n > 0, err.message
        cmps = ( Component * n )()
        assert api.info_composition( h, 0, c_natabund, None, cmps, n,
                                     ctypes.byref(err) ) == n
        compos = [ ( c.Z, c.A, round( c.fraction, 8 ) ) for c in cmps ]
        print(f'{cfg}: T={t.value}K, density={api.info_density(h):.6g}g/cm3,'
              f' composition (all isotopes): {compos}')
        api.deallocate_info( h )

        #Scatter and absorption, compared with the Python API:
        sc = api.create_scatter( cfg.encode(), ctypes.byref(err) )
        ab = api.create_absorption( cfg.encode(), ctypes.byref(err) )
        assert sc and ab
        pysc, pyab = NC.createScatter( cfg ), NC.createAbsorption( cfg )
        assert api.scatter_is_oriented( sc ) == ( 1 if pysc.isOriented() else 0 )
        for wl in ( 0.5, 1.8, 4.0 ):
            xs_sc, xs_ab = ctypes.c_double( -1.0 ), ctypes.c_double( -1.0 )
            nn = neutron( wl, (0.0, 0.0, 1.0) )
            assert api.scatter_cross_section( sc, nn, ctypes.byref(xs_sc),
                                              ctypes.byref(err) ) == 0
            assert api.absorption_cross_section( ab, nn, ctypes.byref(xs_ab),
                                                 ctypes.byref(err) ) == 0
            assert xs_sc.value == pysc.xsect( wl = wl, direction = (0,0,1) )
            assert xs_ab.value == pyab.xsect( wl = wl, direction = (0,0,1) )
        print(f'  is_oriented={api.scatter_is_oriented(sc)}, scatter and'
              ' absorption cross sections same as the Python API: OK')
        api.deallocate_scatter( sc )
        api.deallocate_absorption( ab )

    #A material without temperature:
    h = api.create_info( b'phases<0.5*stdlib::Al_sg225.ncmat;temp=200K'
                         b'&0.5*stdlib::Cu_sg225.ncmat;temp=300K>',
                         ctypes.byref(err) )
    assert h
    t = ctypes.c_double( -1.0 )
    assert api.info_temperature( h, ctypes.byref(t), ctypes.byref(err) ) != 0
    assert err.type == b'MissingInfo' and t.value == -1.0
    print(f'Expected error: [{err.type.decode()}] {err.message.decode()}')
    api.deallocate_info( h )
    api.clear_caches()

if __name__ == '__main__':
    main()
