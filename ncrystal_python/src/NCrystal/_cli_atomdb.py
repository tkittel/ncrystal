
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

"""Command-line tool for browsing NCrystal's database of isotopes and natural
elements."""

from ._cliimpl import cli_entry_point, create_ArgumentParser, print


def climod_metadata():
    return {
        'displaygroup': 'misc',
        'displayorder': 5,
        'descr': "Browse the database of isotopes and natural elements."
    }

#Examples shown in the help (and run by the tests):
_examples = [
    ( 'He B10', 'natural He and its isotopes, and B10' ),
    ( '"absxs>1000"', 'strong absorbers' ),
    ( 'Gd --ncmat', 'output in NCMAT @ATOMDB format' ),
]
_cond_examples = [
    ( '"absxs>1000"', 'strong absorbers' ),
    ( 'A=0 "cohsl<0"', 'natural elements with negative cohsl' ),
    ( 'He A!=0', 'only the isotopes of He' ),
]

_patterns_str = """Each PATTERN can be:

  - an element (e.g. He, also selecting its isotopes),
  - an isotope (e.g. B10, with D and T meaning H2 and H3),
  - a condition, such as "absxs>1000" or A==0 (see --cond-help)."""

def _create_parser( progname ):
    import argparse

    from ._cli_browse import _examples_str, _helpw, _wrap
    def h( text ):
        return _wrap( text, _helpw )
    parser = create_ArgumentParser(
        prog = progname,
        usage = '%(prog)s [options] [PATTERN ...]',
        description = ( _wrap( 'Browse NCrystal\'s database of isotopes and'
                               ' natural elements, with their masses,'
                               ' scattering lengths, and cross sections.' )
                        + '\n\n' + _patterns_str ),
        epilog = ( 'examples:\n' + _examples_str( progname, _examples )
                   + '\n\n' + _wrap( 'Use --cond-help for the list of'
                                     ' properties.' ) ),
        formatter_class = argparse.RawTextHelpFormatter )
    #PATTERNs are explained in the description:
    parser.add_argument('pattern', type=str, nargs='*', metavar='PATTERN',
                        help=argparse.SUPPRESS)
    parser._optionals.title = 'other options'
    parser.add_argument('--cond-help', action='store_true',
                        help=h('Show the properties usable in conditions and'
                               ' exit.'))
    parser.add_argument('--json', action='store_true',
                        help=h('Output all properties of the selected'
                               ' entries as JSON.'))
    parser.add_argument('--ncmat', action='store_true',
                        help=h('Output in NCMAT @ATOMDB format.'))
    return parser

def parseArgs( progname, arglist, return_parser = False ):
    parser = _create_parser( progname )
    if return_parser:
        return parser
    args = parser.parse_intermixed_args( arglist )
    from . import atomdata as ad
    if args.cond_help:
        from ._cli_browse import _cond_help_str
        print( _cond_help_str( _cond_examples,
                               ad.atomdb_fields_doc( conditions = True ),
                               yesno = False ) )
        parser.exit()
    if 'help' in args.pattern:
        parser.print_help()
        parser.exit()
    if args.json and args.ncmat:
        parser.error('Do not specify both --json and --ncmat.')
    from ._cli_browse import _split_patterns
    args.elements, args.conditions = _split_patterns( args.pattern )
    from ._conditions import _Condition
    from .exceptions import NCBadInput
    try:
        names = [ n for n,d in ad.atomdb_fields_doc() ]
        for c in args.conditions:
            _Condition( c, names, ad._atomdb_kinds )
    except NCBadInput as e:
        parser.error( str(e) )
    for p in args.elements:
        if not ad.AtomDBBrowser( p ):
            parser.error(f'Unknown element or isotope "{p}" (specify elements'
                         ' or isotopes, e.g. He B10)')
    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs( progname, [], return_parser = True )

@cli_entry_point
def main( progname, arglist ):
    args = parseArgs( progname, arglist )
    from .atomdata import AtomDBBrowser
    a = AtomDBBrowser( *args.elements, where = args.conditions )
    if args.json:
        print( a.to_json(), end = '' )
        return
    if not a:
        print('No matching entries found.')
        return
    if args.ncmat:
        print( '@ATOMDB' )
        for ll in a.atomdb_lines():
            print( f'  {ll}' )
        return
    print( a.table(), end = '' )
