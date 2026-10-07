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

# Test 64bit RNG seeds and stream indices (passed to the C API as two 32bit
# halves, so they work also where unsigned long is only 32bit).

import NCTestUtils.enable_fpe # noqa F401
import NCrystalDev as NC

cfg = 'stdlib::Al_sg225.ncmat;dcutoff=0.5'

def samples( sc ):
    return [ sc.sampleScatterIsotropic(0.025)[1] for i in range(5) ]

def main():
    s7 = samples(NC.createScatterIndependentRNG(cfg,seed=7))
    s7b = samples(NC.createScatterIndependentRNG(cfg,seed=7))
    s7hi = samples(NC.createScatterIndependentRNG(cfg,seed=7+2**32))
    smax = samples(NC.createScatterIndependentRNG(cfg,seed=2**64-1))
    assert s7 == s7b, "same seed must give same stream"
    assert s7 != s7hi, "seeds differing above bit 32 must give different streams"
    assert smax != s7
    for bad in (-1, 2**64, 1.5):
        try:
            NC.createScatterIndependentRNG(cfg,seed=bad)
        except NC.NCBadInput:
            pass
        else:
            raise SystemExit(f'seed={bad!r} should have been rejected')
    sc = NC.createScatter(cfg)
    c1 = samples(sc.clone(rng_stream_index=3))
    c2 = samples(sc.clone(rng_stream_index=3+2**32))
    assert c1 != c2, "stream indices differing above bit 32 must differ"
    try:
        sc.clone(rng_stream_index=2**64)
    except NC.NCBadInput:
        pass
    else:
        raise SystemExit('rng_stream_index=2**64 should have been rejected')
    print('all ok')

if __name__ == '__main__':
    main()
