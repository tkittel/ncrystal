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
import NCrystalDev.cfgstr as nc_cfgstr
from NCrystalDev.exceptions import NCBadInput, NCFileNotFound
from NCTestUtils.common import ensure_error

def main():
    import pprint
    for c in [ 'stdlib::Al_sg225.ncmat',
               'Al_sg225.ncmat;temp=293.15K',
               'Al_sg225.ncmat;temp=300K',
               'Al_sg225.ncmat;temp=20C',
               'Al_sg225.ncmat;temp=293.15K; vdoslux =3;dcutoff=   0.1 nm',
               'Al_sg225.ncmat;temp=293.15K; vdoslux =3;dcutoff=   0.001 mu',
               'Al_sg225.ncmat;temp=300K',
               'C_sg194_pyrolytic_graphite.ncmat',#embedded cfg
               ('phases<0.7*C_sg227_Diamond.ncmat&0.28*solid::V/6.12gcm3'
                '&0.02*Polyethylene_CH2.ncmat;comp=inelas>;temp=380K'),
               ('phases<0.5*solid::B4C/2.52gcm3/B_is_0.95_B10_0.05_B11'
                '&0.5*gasmix::BF3/2bar>'),
               'phases<1.0*C_sg227_Diamond.ncmat>;temp=380K',
               ('phases<0.4*C_sg227_Diamond.ncmat'
                '&0.6*C_sg227_Diamond.ncmat>;temp=380K'),
              ]:
        print(f"Decoding cfg {c!r}:")
        pprint.pprint(nc_cfgstr.decodeCfg(c),
                      indent=4,
                      sort_dicts=False)

    with ensure_error(NCFileNotFound,
                      'Could not find data: "foobarnonexistent.ncmat"'):
        nc_cfgstr.decodeCfg('foobarnonexistent.ncmat;temp=300K')

    with ensure_error(NCBadInput,
                      'Syntax error - invalid value'
                      ' "warm" provided for parameter "temp"'):
        nc_cfgstr.decodeCfg('stdlib::Al_sg225.ncmat;temp=warm')

def test_print_redirection():
    #generateCfgStrDoc('print') must use the (redirectable) NCrystal print:
    import contextlib
    import io

    from NCrystalDev._common import capture_print_ctxmgr
    leaked = io.StringIO()
    with contextlib.redirect_stdout(leaked), capture_print_ctxmgr() as cap:
        nc_cfgstr.generateCfgStrDoc('print')
    assert not leaked.getvalue(), 'output bypassed NCrystal print'
    assert ''.join(cap.data) == nc_cfgstr.generateCfgStrDoc('txt_full') + '\n'
    print('generateCfgStrDoc("print") output redirected OK')

def test_type_erasure():
    #decodeCfg and decodecfg_vdoslux accept anything MaterialSource accepts,
    #with optional cfg_params (incl. any embedded NCRYSTALMATCFG in the data):
    import NCrystalDev as NC
    from NCrystalDev.exceptions import NCCalcError
    from NCrystalDev.misc import MaterialSource
    txt = NC.createTextData('stdlib::Al_sg225.ncmat').rawData
    txt_emb = txt + '# NCRYSTALMATCFG[vdoslux=1]\n'
    def show( lbl, material, cfg_params = None ):
        d = nc_cfgstr.decodeCfg( material, cfg_params )
        v = nc_cfgstr.decodecfg_vdoslux( material, cfg_params )
        d2 = nc_cfgstr.decodeCfg( material, cfg_params = cfg_params )
        #raw data gets a new TextData (and uid) per decoding:
        for e in (d,d2):
            del e['textdata_uid']
        assert d == d2
        print(f'{lbl}: data_name={d["data_name"]!r}'
              f' textdata_name={d["textdata_name"]!r} pars={d["pars"]}'
              f' vdoslux={v}')
    show( 'cfgstr', 'stdlib::Al_sg225.ncmat' )
    show( 'cfgstr+params', 'stdlib::Al_sg225.ncmat', 'temp=20K;vdoslux=5' )
    show( 'ncmat', txt )
    show( 'ncmat+embedded', txt_emb )
    show( 'ncmat+embedded+params', txt_emb, 'temp=50K;vdoslux=4' )
    show( 'TextData', NC.createTextData('stdlib::Al_sg225.ncmat'), 'vdoslux=2' )
    import pathlib
    p = pathlib.Path('tmp_cfgstr_weird #name&.ncmat')
    p.write_text(txt)
    show( 'file path', p, 'temp=30K' )
    p.unlink()
    #(non-crystalline, so no spglib needed):
    show( 'NCMATComposer', NC.NCMATComposer.from_cfgstr('gasmix::He'),
          'vdoslux=6' )
    show( 'MaterialSource', MaterialSource(txt_emb,cfg_params='temp=10K') )
    show( 'MaterialSource+params',
          MaterialSource(MaterialSource(txt_emb,cfg_params='temp=10K'),
                         cfg_params='vdoslux=5') )
    js = nc_cfgstr.decodeCfg( txt_emb, 'temp=50K', asJSONStr = True )
    import json
    assert ( json.loads(js)['allpars']
             == nc_cfgstr.decodeCfg( txt_emb, 'temp=50K' )['allpars'] )
    print('multiphase vdoslux:', nc_cfgstr.decodecfg_vdoslux(
        'phases<0.5*Al_sg225.ncmat&0.5*Cu_sg225.ncmat>;vdoslux=2' ))
    with ensure_error(NCCalcError,
                      'Could not determine unique value of parameter'
                      ' "vdoslux" on multiphase material (different values'
                      ' found in different phases).'):
        nc_cfgstr.decodecfg_vdoslux('phases<0.5*Al_sg225.ncmat;vdoslux=1'
                                    '&0.5*Cu_sg225.ncmat;vdoslux=4>')
    with ensure_error(NCBadInput,
                      'Can not decode cfg of a MaterialSource which is'
                      ' preloaded.'):
        nc_cfgstr.decodeCfg( NC.createInfo('stdlib::Al_sg225.ncmat') )

def test_doc_modes():
    import json
    doc = nc_cfgstr.generateCfgStrDoc('python')
    assert doc == json.loads(nc_cfgstr.generateCfgStrDoc('json'))
    short = nc_cfgstr.generateCfgStrDoc('txt_short')
    full = nc_cfgstr.generateCfgStrDoc('txt_full')
    assert isinstance(short,str) and isinstance(full,str)
    assert len(full) > len(short)
    print('generateCfgStrDoc modes OK')

if __name__ == '__main__':
    main()
    test_print_redirection()
    test_type_erasure()
    test_doc_modes()
