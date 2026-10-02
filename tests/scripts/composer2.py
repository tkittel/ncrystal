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

# NEEDS: spglib
import NCrystalDev as NC
from NCTestUtils.common import ensure_error


def main():
    print('\n\n  ================> Al example 1\n\n')

    c_Al = NC.NCMATComposer()
    c_Al.set_cellsg_cubic( 4.05 )
    c_Al.set_atompos( [ ('Al',0,0,0),
                        ('Al',0,1/2,1/2),
                        ('Al',1/2,0,1/2),
                       ('Al',1/2,1/2,0)])
    c_Al.allow_fallback_dyninfo()
    with ensure_error(NC.NCBadInput,
                      "Must provide a space group number (or invoke"
                      " .refine_crystal_structure()) before it is"
                      " possible to verify a crystal structure"):
        c_Al.create_ncmat()

    c_Al.refine_crystal_structure()

    print( c_Al() )
    c_Al.load().dump()

    print('\n\n  ================> Al example 2\n\n')
    c_Al2 = NC.NCMATComposer()
    c_Al2.set_cellsg_cubic( 4.05, spacegroup=225 )
    c_Al2.set_atompos( [('Al',0,0,0),('Al',0,1/2,1/2),('Al',1/2,0,1/2),('Al',1/2,1/2,0)])
    c_Al2.allow_fallback_dyninfo()
    c_Al2.set_composition('Al','0.99 Al 0.01 Cr')
    print(c_Al2())

    print('\n\n  ================> Al example 3\n\n')

    c_Al3 = NC.NCMATComposer()
    c_Al3.set_cellsg_cubic( 4.05 )
    c_Al3.set_atompos( [('tight_atom',0,0,0),('loose_atom',0,1/2,1/2),('loose_atom',1/2,0,1/2),('loose_atom',1/2,1/2,0)])
    c_Al3.set_dyninfo_msd('tight_atom',msd=0.005, temperature=200)
    c_Al3.set_dyninfo_msd('loose_atom',msd=0.02, temperature=200)
    c_Al3.set_composition('tight_atom','Al')
    c_Al3.set_composition('loose_atom','Al')
    c_Al3.refine_crystal_structure()#Detect spacegroup
    print(c_Al3())

    print('\n\n  ================> Al SANS example\n\n')

    c=NC.NCMATComposer('Al_sg225.ncmat') #<--- Can init from cfg-string.
    c.add_secondary_phase(0.01,'void.ncmat')
    c.add_hard_sphere_sans_model(50)
    print(c())
    c.load().dump()



def test_v8_scatknl_keywords():
    #The effective_temperature/msd/debye_temp parameters of
    #set_dyninfo_scatknl emit NCMAT v8 keywords (switching the header
    #to v8), flowing to the extension/elastic machinery on load. A
    #non-crystalline v8 material without a declared state of matter is
    #refused already at create_ncmat time:
    import NCrystalDev as NC
    import numpy as np
    from NCrystalDev.misc import evaluate_query as q
    def mkc( **kwargs ):
        c = NC.NCMATComposer()
        c.set_composition('Al','Al')
        c.set_density( 2.7, 'g/cm3' )
        a = np.linspace(1e-6,100.0,10)
        b = np.linspace(0.0,100.0,10)#symmetric-scaled form
        c.set_dyninfo_scatknl( 'Al', alphagrid=a, betagrid=b,
                               temperature=300.0,
                               sab_scaled=np.ones(100), **kwargs )
        return c
    try:
        mkc( effective_temperature = 450.0 ).create_ncmat()
    except NC.NCBadInput as e:
        print(f'=> correctly refused stateless v8: {e}')
    c = mkc( effective_temperature = 450.0, msd = 0.05 )
    c.set_state_of_matter('solid')
    txt = c.create_ncmat()
    assert txt.startswith('NCMAT v8')
    assert '\n  effective_temperature 450\n' in txt
    assert '\n  msd .05\n' in txt
    r = q( ['sab','analyse',c.register_as('comp_v8.ncmat',''),''] )
    r = r['sabanalyse']
    assert r['teff_extender'] == 450.0
    print('=> v8 keywords roundtrip (teff_extender=450): OK')
    #no keywords -> still v7:
    c2 = mkc()
    assert c2.create_ncmat().startswith('NCMAT v7')
    print('=> plain scatknl composer output remains NCMAT v7: OK')

if __name__ == '__main__':
    main()
    test_verify_origin_and_triclinic()
    test_refine_corrections()
    test_v8_scatknl_keywords()
