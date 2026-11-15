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

#Tests NCrystal's handling of optional dependencies (_imports.py), using fake
#entries and fake modules.

def main():
    import pathlib
    import sys

    import NCrystalDev._depdb as depdb
    import NCrystalDev._imports as imp
    from NCTestUtils.common import ensure_error
    d = pathlib.Path('fakepypath').absolute()
    for name, content in [
            ( 'ncfake_old', '__version__ = "1.2.3"\n' ),
            ( 'ncfake_ok', '__version__ = "1.10.0"\n' ),
            ( 'ncfake_novers', '' ),
            ( 'ncfake_broken', 'raise ImportError("broken on purpose")\n' ),
            ( 'ncfake_sub', '__version__ = "2.0"\n' ) ]:
        d.joinpath(name).mkdir(parents=True)
        d.joinpath(name,'__init__.py').write_text(content)
    d.joinpath('ncfake_sub','extra.py').write_text('')
    sys.path.insert(0,str(d))
    for name, minv, modname in [ ( 'old', '1.10', 'ncfake_old' ),
                                 ( 'ok', '1.9.2', 'ncfake_ok' ),
                                 ( 'novers', '1.0', 'ncfake_novers' ),
                                 ( 'broken', '1.0', 'ncfake_broken' ),
                                 ( 'missing', '1.0', 'ncfake_missing' ),
                                 ( 'sub', '2', 'ncfake_sub.extra' ),
                                 ( 'nomin', None, 'ncfake_missing' ) ]:
        depdb.deps[f'fake-{name}'] = ( minv, modname,
                                       f'fake-{name}', f'fake-{name}-conda' )
    for name in sorted(depdb.deps):
        if not name.startswith('fake-'):
            continue
        mod, problem = imp.find_optional( name )
        print(f'{name}: module={mod.__name__ if mod else None}'
              f' problem={problem!r}')
        assert imp.find_optional( name ) == ( mod, problem )#cached
        if problem:
            print('  ',imp.optional_error_msg( name, problem ))
            with ensure_error(ImportError,
                              imp.optional_error_msg( name, problem,
                                                      'for testing' )):
                imp.import_optional( name, purpose = 'for testing' )
            with ensure_error(SystemExit,
                              'ERROR: '+imp.optional_error_msg( name,
                                                                problem )):
                imp.import_optional( name, sysexit = True )
        else:
            assert imp.import_optional( name ) is mod

if __name__ == '__main__':
    main()
