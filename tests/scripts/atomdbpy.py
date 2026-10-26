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

# NEEDS: numpy

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.stabilise_ncpprint # noqa F401
import NCrystalDev as NC
import NCrystalDev.atomdata as ad
from NCrystalDev.misc import evaluate_query
from NCTestUtils.common import ensure_error

def test_atomdbbrowser():
    import json

    import NCrystalDev._common as nc_common
    import NCrystalDev.misc as nc_misc
    #The static database is only queried once (and cached):
    nqueries = [0]
    orig_q = nc_misc.evaluate_query
    def counting_q( q ):
        nqueries[0] += ( q == ['util','atomdb'] )
        return orig_q( q )
    nc_misc.evaluate_query = counting_q
    try:
        ad.AtomDBBrowser('He')
        ad.AtomDBBrowser().where('absxs > 1')
        assert nqueries[0] <= 1
    finally:
        nc_misc.evaluate_query = orig_q
    #Iteration gives copies (modifying them does not affect the cache):
    a = ad.AtomDBBrowser()
    d = list(a)
    d[0]['mass'] = -1.0
    assert next(iter(ad.AtomDBBrowser()))['mass'] > 0.0
    print('atomdb fields:',[ n for n,_ in ad.atomdb_fields_doc() ])
    assert len(a) == len(d) and repr(a) == f'AtomDBBrowser({len(d)} entries)'
    print('He + b10 (case-insensitive):',ad.AtomDBBrowser('He','b10').names())
    print('Li:',a.match('Li').names())
    print('where:',a.where('absxs > 1000','A=0').names())
    sel = ad.AtomDBBrowser('He','B10')
    print('==> text table:')
    print(sel.table(),end='')
    print('==> json:')
    nc_common.ncpprint( json.loads( sel.to_json() ) )
    assert list(sel)[1]['label'] == 'He3'
    print('==> atomdb_lines():')
    for ll in ad.AtomDBBrowser('He','N15','In113','Gd157').atomdb_lines():
        print(ll)
    print('==> atomdb_lines(align=False):')
    for ll in ad.AtomDBBrowser('He','N15').atomdb_lines(align=False):
        print(ll)
    #The lines reproduce the database values exactly when used in NCMAT:
    for d, ll in zip( a, a.atomdb_lines() ):
        lbl = d['label']
        ncmat = ( 'NCMAT v7\n@DENSITY\n  1 g_per_cm3\n@DYNINFO\n'
                  f'  element {lbl}\n  fraction 1\n  type freegas\n'
                  f'@ATOMDB\n  nodefaults\n  {ll}\n' )
        adata = NC.directLoad( ncmat ).info.composition[0][1]
        ref = NC.atomDB( d['Z'], d['A'] ) if d['A'] else NC.atomDB( lbl )
        assert adata.averageMassAMU() == ref.averageMassAMU(), lbl
        assert adata.coherentScatLenFM() == ref.coherentScatLenFM(), lbl
        assert adata.incoherentXS() == ref.incoherentXS(), lbl
        assert adata.captureXS() == ref.captureXS(), lbl
    print('atomdb_lines roundtrip OK for',len(a),'entries')
    def bad( msg, fct ):
        with ensure_error(NC.NCBadInput,msg):
            fct()
    names = ', '.join( n for n,_ in ad.atomdb_fields_doc( conditions = True ) )
    bad('Invalid property "foo" in condition "foo > 1" (must be one'
        ' of: '+names+')',lambda : a.where('foo > 1'))
    bad('Invalid property "COHSL" in condition "COHSL<0" (names are'
        ' case-sensitive, did you mean "cohsl"?)',lambda : a.where('COHSL<0'))
    bad('Invalid property "a" in condition "a=0" (names are case-sensitive,'
        ' did you mean "A"?)',lambda : a.where('a=0'))
    bad('Invalid condition "label=He" ("label" can not be used in'
        ' conditions)',lambda : a.where('label=He'))
    assert 'label' not in [ n for n,_ in ad.atomdb_fields_doc(conditions=True) ]
    print('AtomDB OK')

for a in NC.iterateAtomDB():
    print(a)
    print(repr(a.to_atomdb_str()))
for a in NC.iterateAtomDB(objects=False):
    print(a)

#Same formatting as the @ATOMDB lines from the C++ code:
for d in evaluate_query(['util','atomdb']):
    a = NC.atomDB( d['Z'], d['A'] ) if d['A'] else NC.atomDB( d['label'] )
    assert d['atomdbline'] == d['label'] + ' ' + a.to_atomdb_str(), d['label']
print('to_atomdb_str() identical to the C++ formatting for all entries')

test_atomdbbrowser()
