#!/usr/bin/env python3

#Temporary CI diagnostic: where does the Windows NaN in sabsample test 20 come from?
import NCTestUtils.enable_fpe # noqa F401
import NCTestUtils.enable_testdatapath # noqa F401
import math
import NCrystalDev as NC
from NCrystalDev.misc import evaluate_query as q
cfg = 'benzene_solid_100K_sabsmall.ncmat;vdoslux=2001'
def tryq(lbl, query):
    try:
        r = q(query, huge_arrays=True)
        print(f'{lbl}: OK', flush=True)
        return r
    except Exception as e:
        print(f'{lbl}: {type(e).__name__}: {e}', flush=True)
tryq('proc E=7.0 nsample=0', ['sab','proc',cfg,'H','0','1','7.0'])
r = tryq('proc E=7 nsample=0 (again)', ['sab','proc',cfg,'H','0','1','7.0'])
if r:
    for k,v in r['sabproc'].items():
        if isinstance(v,list) and v and isinstance(v[0],float):
            print('  ',k,len(v),min(v),max(v))
for seed in (1,2,3):
    for n in (1,10,100,1000,10000,100000):
        if not tryq(f'proc E=7 seed={seed} nsample={n}',
                    ['sab','proc',cfg,'H',str(n),str(seed),'7.0']):
            break
tryq('refsample E=7', ['sab','refsample',cfg,'H','100000','1','7.0','ref'])
sc = NC.createScatter(cfg)
for ekin in (7.0, 5.0, 10.0):
    nbad = 0
    for i in range(100000):
        e,mu = sc.sampleScatterIsotropic(ekin)
        if not (math.isfinite(e) and math.isfinite(mu)):
            nbad += 1
            if nbad <= 3:
                print(f'  createScatter E={ekin}: sample {i}: E={e} mu={mu}')
    print(f'createScatter E={ekin}: {nbad} non-finite of 100000', flush=True)
raise SystemExit(1)#always fail, so ctest shows the output
