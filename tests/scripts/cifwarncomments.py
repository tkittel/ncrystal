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

# Test of the formatting of warnings in the comments of NCMAT data created
# from CIF data (wrapped, with repeats merged and long warnings truncated).

import NCTestUtils.enable_fpe # noqa F401
from NCrystalDev._cifcomposerimpl import _fmt_warnings_for_comments as fmt

def show( title, warnings, nwmax = 10 ):
    print( f'==> {title}:' )
    for line in fmt( warnings, 'in this test', nwmax = nwmax ):
        print( f'# {line}'.rstrip() )

def main():
    show( 'none', [] )
    show( 'short ones, with a repeat', [ ( 'Short warning', 'NCrystalWarning' ),
                                         ( 'Python one', 'UserWarning' ),
                                         ( 'Short warning', 'NCrystalWarning' ) ] )
    long = ' '.join( f'word{i}' for i in range( 60 ) )
    show( 'wrapped', [ ( long, 'NCrystalWarning' ) ] )
    show( 'several lines of text', [ ( 'First line\n\n  Second line  \n'
                                       + long, 'UserWarning' ) ] )
    show( 'truncated', [ ( ' '.join( [ long ] * 3 ), 'NCrystalWarning' ) ] )
    show( 'truncated, nwmax=3', [ ( long, 'NCrystalWarning' ) ], nwmax = 3 )
    show( 'long category (at least 50 columns of text)',
          [ ( long, 'AVeryLongWarningCategoryNameFromSomeOtherPackage' ) ] )
    show( 'unbreakable word', [ ( 'x' * 100 + ' end', 'NCrystalWarning' ) ] )

if __name__ == '__main__':
    main()
