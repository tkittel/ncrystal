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

def main():
    #Checking https://github.com/mctools/ncrystal/issues/262
    import NCrystalDev as NC
    import NCrystalDev.misc as nc_misc
    #At 5000K the cheaper inelas probe (vdoslux<=2) fails, so the retry with
    #the material's own vdoslux is also tested:
    for cfg in [ "CaO2H2_sg164_CalciumHydroxide.ncmat;temp=3000K",
                 "CaO2H2_sg164_CalciumHydroxide.ncmat;temp=5000K",
                 "CaO2H2_sg164_CalciumHydroxide.ncmat;temp=2000K;vdoslux=1",
                 "Al_sg225.ncmat;comp=inelas",
                 "phases<0.5*Al_sg225.ncmat;vdoslux=1&0.5*Cu_sg225.ncmat>",
                 "gasmix::He",
                ]:
        print(cfg)
        print('  ',nc_misc.detect_scattering_components(cfg))
        for c,lbl in nc_misc.cfgstr_detect_components(cfg):
            print(f'   {lbl}: {c}')
    #Raw data (with cfg_params) via MaterialSource:
    txt = NC.createTextData('stdlib::Al_sg225.ncmat').rawData
    ms = nc_misc.MaterialSource(txt,cfg_params='temp=5000K;vdoslux=1')
    print('raw data:',nc_misc.detect_scattering_components(ms))

if __name__ == '__main__':
    main()
