
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

# For checking some file that must be kept completely synchronised.

def get_content( relpath ):
    from .dirs import reporoot
    print(f'    {relpath}')
    return reporoot.joinpath(relpath).read_text()

def check_same( reffile, *otherfiles ):
    assert len(otherfiles) > 0
    print('  Checking for same contents:')
    ref = get_content( reffile )
    for o in otherfiles:
        if get_content(o) != ref:
            print()
            raise SystemExit(f'ERROR: Content of {o} and {reffile} differs')

def check_same_declarations( reffile, otherfile, typenames ):
    #Check that the C typedef'ed structs with the given names are the same in
    #both files, except for comments and whitespace (which may differ):
    import re
    print('  Checking for same C declarations:')
    def extract( relpath ):
        content = get_content( relpath )
        content = re.sub( r'/\*.*?\*/', ' ', content, flags = re.DOTALL )
        content = re.sub( r'//[^\n]*', ' ', content )
        res = {}
        for tn in typenames:
            m = re.findall( r'typedef\s+struct\s*\{[^}]*\}\s*%s\s*;'%tn,
                            content )
            if len(m) != 1:
                raise SystemExit(f'ERROR: Did not find exactly one'
                                 f' declaration of {tn} in {relpath}')
            res[tn] = ' '.join( re.sub( r'([(){};,*])', r' \1 ',
                                        m[0] ).split() )
        return res
    ref = extract( reffile )
    other = extract( otherfile )
    for tn in typenames:
        if ref[tn] != other[tn]:
            print()
            raise SystemExit(f'ERROR: Declaration of {tn} differs in'
                             f' {otherfile} and {reffile}')

def main():
    check_same( 'README.md',
                'ncrystal_metapkg/README.md' )

    check_same( 'ncrystal_core/README.md',
                'ncrystal_core/empty_pypkg/README.md' )

    check_same( 'LICENSE',
                'ncrystal_core/LICENSE',
                'ncrystal_core/empty_pypkg/LICENSE',
                'ncrystal_python/LICENSE',
                'ncrystal_metapkg/LICENSE',
                'ncrystal_pypluginmgr/LICENSE',
                'ncrystal_verify/LICENSE',
                'examples/plugin/LICENSE',
                'examples/plugin_dataonly/LICENSE')

    check_same( 'examples/ncrystal_example_cpp.cc',
                'examples/downstream_cmake/main.cc' )

    #NB: For ncrystal-verify we combine the more robust
    #ncrystal_core/app_test/main.cc with the CMake file from the downstream
    #example.
    check_same( 'ncrystal_core/app_test/main.cc',
                'ncrystal_verify/extra/data/downstream_cmake_main.cc',
                'tests/src/app_cpp/main.cc',
               )
    check_same( 'examples/downstream_cmake/CMakeLists.txt',
                'ncrystal_verify/extra/data/downstream_cmake_CMakeLists.txt' )

    #The NCrystalVAPIType2V1 class is tested in tests/src/app_virtcapicxx:
    for fn in ( 'NCrystalVAPIType2V1.hh', 'NCrystalVAPIType2V1.cc' ):
        check_same( f'examples/virtualcapi_project/src/{fn}',
                    f'tests/src/app_virtcapicxx/{fn}' )

    #The NCrystalVAPIType2V1 class, which applications copy into their
    #projects, contains a copy of the declarations from ncvirtapi.h:
    check_same_declarations(
        'ncrystal_core/include/NCrystal/virtualapi/ncvirtapi.h',
        'examples/virtualcapi_project/src/NCrystalVAPIType2V1.cc',
        [ 'ncrystal_vapi_error_t', 'ncrystal_vapi_type2_v1_t' ] )

if __name__=='__main__':
    main()
