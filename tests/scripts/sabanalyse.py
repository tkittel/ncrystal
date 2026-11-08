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

#Exercises the ["sab","analyse",...] query (NCSABAnalyser): Teff/msd
#moment-estimators on S(alpha,beta) kernels, validated against the
#VDOSEval-derived truth for VDOS-based materials. NB: printed values are
#deliberately low-precision (tolerance-asserted instead), keeping the
#reference log robust against last-bit platform differences.

import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query as q

def analyse( cfgstr, lbl='', diag=False ):
    query = ['sab','analyse',cfgstr,lbl]
    if diag:
        query += ['diag']
    return q(query)['sabanalyse']

def check_vdos_case( cfgstr, lbl, teff_tol, msd_tol ):
    r = analyse( cfgstr, lbl )
    te, tv = r['teff'], r['teff_vdos']
    ms, mv = r['msd'], r['msd_vdos']
    assert te is not None and tv is not None
    assert ms is not None and mv is not None
    rel_te = abs( te/tv - 1.0 )
    rel_ms = abs( ms/mv - 1.0 )
    ok_te = rel_te < teff_tol
    ok_ms = rel_ms < msd_tol
    print(f'==> {cfgstr}//{lbl or "mono"}:')
    print(f'    teff ~ {te:.3g} K (vdos {tv:.3g} K), relerr<{teff_tol:g}:'
          f' {"OK" if ok_te else "FAIL"}'
          f' [nrows={r["teff_nrows"]}]')
    print(f'    msd  ~ {ms:.2g} Aa^2 (vdos {mv:.2g} Aa^2),'
          f' relerr<{msd_tol:g}: {"OK" if ok_ms else "FAIL"}'
          f' [nrows={r["msd_nrows"]}]')
    assert ok_te and ok_ms, f'tolerances exceeded ({rel_te=}, {rel_ms=})'
    #These well-behaved kernels must pass the turn-key acceptance
    #policy, so "auto" simply repeats the (default-options) raw values:
    assert r['auto'] == { 'teff': te, 'msd': ms }

def main():
    check_vdos_case('stdlib::Al_sg225.ncmat','', 0.02, 0.02)
    check_vdos_case('stdlib::CaH2_sg62_CalciumHydride.ncmat;temp=20','H',
                    0.03, 0.08)
    check_vdos_case('stdlib::CaH2_sg62_CalciumHydride.ncmat;temp=20','Ca',
                    0.02, 0.02)
    check_vdos_case('stdlib::Fe_sg229_Iron-alpha.ncmat;temp=800','',
                    0.02, 0.02)
    check_vdos_case('stdlib::Polyethylene_CH2.ncmat','H', 0.02, 0.03)

    #A genuine direct (scatknl) kernel: a liquid. Teff must come out
    #sane (no vdos reference exists), while the msd estimate must be
    #self-flagged as unreliable (no Debye-Waller deficit in a liquid):
    r = analyse('stdlib::LiquidWaterH2O_T293.6K.ncmat','H')
    assert r['teff_vdos'] is None and r['msd_vdos'] is None
    T = r['temperature']
    ok_teff = ( r['teff'] is not None and 2.0*T < r['teff'] < 10.0*T
                and r['teff_relspread'] < 0.5 )
    ok_msd = r['msd'] is None or r['msd_relspread'] > 0.5
    #Turn-key policy: teff accepted, unreliable msd rejected:
    ok_auto = ( r['auto']['teff'] == r['teff']
                and r['auto']['msd'] is None )
    print('==> liquid water (direct kernel):')
    print(f'    teff in (2T,10T) with modest spread:'
          f' {"OK" if ok_teff else "FAIL"}')
    print(f'    msd self-flagged unreliable: {"OK" if ok_msd else "FAIL"}')
    print(f'    auto: teff accepted, msd rejected:'
          f' {"OK" if ok_auto else "FAIL"}')
    assert ok_teff and ok_msd and ok_auto

    #NCMAT v8 explicit effective_temperature takes precedence over the
    #kernel-based estimate (query teff is the extender-level value):
    txt = NC.createTextData('benzene_solid_100K_sabsmall.ncmat').rawData
    parts = txt.replace('NCMAT v7','NCMAT v8',1).split('@DYNINFO')
    for i, p in enumerate(parts):
        if 'element     C' in p:
            parts[i] = p + '  effective_temperature 600\n'
    NC.registerInMemoryFileData( 'v8efft_benzene.ncmat',
                                 '@DYNINFO'.join(parts) )
    r = analyse('virtual::v8efft_benzene.ncmat','C')
    ok_efft = r['teff_extender'] == 600.0 and abs(r['teff']/578.1-1)<0.01
    print(f'==> explicit effective_temperature honoured:'
          f' {"OK" if ok_efft else "FAIL"}')
    assert ok_efft

    #Trimmed JENDL-5 fixtures (tests/data, ENDF truths in the file
    #headers): teff agrees with the ENDF values, msd is accepted for C in
    #benzene@100K and C-from-benzene@20K, but refused for H in benzene@100K:
    def teff_ok( r, endf ):
        t = r['auto']['teff']
        return t is not None and abs( t/endf - 1.0 ) < 0.03
    rC = analyse('benzene_solid_100K_sabsmall.ncmat','C')
    rH = analyse('benzene_solid_100K_sabsmall.ncmat','H')
    ok_C = teff_ok( rC, 572.9425 ) and rC['auto']['msd'] is not None
    ok_H = teff_ok( rH, 1184.981 ) and rH['auto']['msd'] is None
    r20 = analyse('C_from_benzene_solid_20K_sabsmall.ncmat')
    ok_20 = ( teff_ok( r20, 556.8517 )
              and r20['auto']['msd'] is not None
              and 0.012 < r20['auto']['msd'] < 0.016 )
    print(f'==> JENDL benzene fixtures: 100K C teff+msd:'
          f' {"OK" if ok_C else "FAIL"}, 100K H teff, no msd:'
          f' {"OK" if ok_H else "FAIL"}, 20K C teff+msd:'
          f' {"OK" if ok_20 else "FAIL"}')
    assert ok_C and ok_H and ok_20

    #Diagnostics mode: aligned per-row arrays and valid status codes:
    d = analyse('stdlib::Al_sg225.ncmat','',diag=True)['diagnostics']
    n = len(d['alpha'])
    assert n > 100
    for k in ('m0','mean','variance','teff_row','msd_row',
              'teff_row_status','msd_row_status'):
        assert len(d[k]) == n, k
    assert all( 0 <= int(s) <= 6 for s in d["teff_row_status"] )
    assert all( 0 <= int(s) <= 6 for s in d["msd_row_status"] )
    assert d['floor_value'] > 0.0
    nok = sum( 1 for s in d['teff_row_status'] if int(s)==0 )
    assert nok > 100
    print(f'==> diagnostics: {n} rows, aligned arrays and'
          ' valid status codes: OK')

    #Consumer-level check: the auto-detected Teff must reach the factory
    #layer, i.e. next-gen direct-kernel materials get the SCT extension
    #by default (free-gas only when explicitly selected with the
    #["sab","xs",...] query; NB the total xs above Emax is
    #continuity-anchored and hence insensitive to the extender model --
    #only sampling differs -- so we check the chosen model instead of xs
    #values):
    cfgstr = ( 'stdlib::LiquidWaterH2O_T293.6K.ncmat;comp=inelas;vdoslux=1'
               ';knllux=1' )
    s = NC.load( cfgstr ).scatter.getSummary()
    #NB: only H has a kernel (the O component is plain free gas):
    c = [ c for _,c in s['specific']['components']
          if c['name'] == 'SABScatterNG' ]
    assert len(c) == 1
    m_sct = c[0]['specific']['extension_method']
    m_fg = q( [ 'sab', 'xs', cfgstr, 'H', 'freegas', '1.0' ] )['extension_method']
    print(f'==> direct-kernel extension methods: knllux=1 -> {m_sct},'
          f' freegas query -> {m_fg}')
    assert ( m_sct, m_fg ) == ( 'sct', 'freegas' )

if __name__ == '__main__':
    main()
