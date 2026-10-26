
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

"""Command-line tool for browsing and searching the available data."""

from ._cliimpl import cli_entry_point, create_ArgumentParser, print


def climod_metadata():
    return {
        'displaygroup': 'main',
        'displayorder': 15,
        'descr': "Browse and search available data."
    }

#Examples shown in the help (and run by the tests):
_examples = [
    ( 'Al_sg225', 'data with "al_sg225" in the name' ),
    ( 'Al', 'name has "al", or material has Al' ),
    ( 'B10', 'name has "b10", or material has B10' ),
    ( 'D', 'name has "d", or material has D' ),
    ( '"absxs>100"', 'strong absorbers' ),
    ( 'sg==225 "density<3"', 'all conditions must be fulfilled' ),
    ( 'scatknl=yes', 'materials with scattering kernels' ),
    ( 'state!=solid', 'liquids, gases, or unknown state' ),
    ( 'state=liquid H', 'liquids with "h" in name, or with H' ),
    ( '-s metal', 'header comments with "metal"' ),
    ( '-x Al_sg225.ncmat', 'show the full content' ),
]
_cond_examples = [
    ( 'sg=225', 'space group 225' ),
    ( '"absxs>100"', 'strong absorbers' ),
    ( 'state=liquid nelements=2', 'liquids with two elements' ),
    ( 'scatknl=yes', 'data with full scattering kernels' ),
    ( 'state!=solid', 'liquids, gases, or unknown state' ),
    ( 'crystal=no', 'non-crystalline materials' ),
]

_patterns_str = """Each PATTERN can be:

  - a substring of the data name (e.g. parts of the file name),
  - the name of an element or isotope in the materials (e.g. Al or B10),
  - a condition, such as "absxs>100" or sg==225 (see --cond-help).

Substrings are case-insensitive, but element and isotope names must be written
exactly (e.g. B for boron)."""

def _create_parser( progname, full ):
    #Parser for -h (full=False, one-line option help) or --full-help:
    import argparse
    def h( short, details = None, advanced = False ):
        #Advanced options are only documented in --full-help:
        if not full and advanced:
            return argparse.SUPPRESS
        return _wrap( ( details or short ) if full else short, _helpw )
    descr = ( _wrap( 'Browse and search the data available to NCrystal from'
                     ' built-in sources or plugins.' )
              + '\n\n' + _patterns_str )
    if full:
        txt = """
        Data from all sources is listed, e.g. files in the standard data
        library, in-memory data, or data created on-demand like
        "solid::B4C/2.52gcm3". Files found relative to the current directory
        are only listed with "-f relpath". Element names and conditions
        require all candidate materials to be loaded.

        Listed names include their source (like "stdlib::Al_sg225.ncmat").
        When several sources provide the same name, the first one listed is
        used when the name is given without its source.
        """
        descr += '\n\n' + '\n\n'.join( _wrap(e) for e in
                                         txt.strip().split('\n\n') )
    epilog = ( 'examples:\n' + _examples_str( progname, _examples ) + '\n\n'
               + _wrap( 'Use --cond-help for the list of properties,'
                        + ( ' and "-f help" for the data sources.' if full else
                            ' and --full-help for all options.' ) ) )
    parser = create_ArgumentParser( prog = progname,
                                    usage = ( None if full else
                                              '%(prog)s [options]'
                                              ' [PATTERN ...]' ),
                                    description = descr,
                                    epilog = epilog,
                                    formatter_class
                                    = argparse.RawTextHelpFormatter )
    parser.add_argument('--full-help', action='store_true',
                        help=h('Show detailed help and exit.'))
    parser.add_argument('--cond-help', action='store_true',
                        help=h('Show the properties usable in conditions and'
                               ' exit.'))
    #PATTERNs are explained in the description:
    parser.add_argument('pattern', type=str, nargs='*', metavar='PATTERN',
                        help=( h('Substring, element, or condition (see'
                                 ' above).') if full else argparse.SUPPRESS ))
    if not full:
        parser._optionals.title = 'other options'
    gsel = parser.add_argument_group('selecting data') if full else parser
    gsel.add_argument('-s','--search', type=str, action='append',
                      default=[], metavar='WORD',
                      help=h('Only show data with WORD in names or comments.',
                             'Only show data which contain WORD in their'
                             ' name or NCMAT header comments (case'
                             '-insensitive). Can be specified multiple'
                             ' times, in which case all words must be'
                             ' present. Matching comment lines are shown.'))
    gsel.add_argument('-f','--factory', type=str, default=None,
                      metavar='SRC',
                      help=h('Only show data from the source SRC (e.g.'
                             ' "stdlib" or "virtual", see "-f help").',
                             advanced=True))
    gout = parser.add_argument_group('output') if full else parser
    gout.add_argument('--json', action='store_true',
                      help=h('Output all information about the selected'
                             ' data (including physics properties) as'
                             ' JSON.',
                             advanced=True))
    gout.add_argument('--color','--colour', type=str,
                      default='auto', metavar='WHEN',
                      choices=sorted(_color_choices),
                      help=h('Highlight --search hits: "always", "never",'
                             ' or "auto" (default, only on a terminal).',
                             advanced=True))
    goth = parser.add_argument_group('other') if full else parser
    goth.add_argument('-x','--extract', type=str, default=None,
                      metavar='ENTRY',
                      help=h('Print the full content of ENTRY.',
                             'Print the full content of ENTRY (e.g. a'
                             ' file name).'))
    return parser

_helpw = 54#width of option help texts (in 80 columns)

def _wrap( text, width = 79 ):
    import textwrap
    return textwrap.fill( ' '.join(text.split()), width = width,
                          break_on_hyphens = False )

def _examples_str( progname, examples ):
    #Commands with a "# description" after them, or below them if too long:
    out = []
    for e,d in examples:
        cmd = f'  {progname} {e}'
        out.append( f'{cmd.ljust(40)}  # {d}' if len(cmd) <= 40
                    else f'{cmd}\n{"":42}# {d}' )
    return '\n'.join( out )

def _docs_str( pd ):
    #Aligned list of (name,description), where consecutive names with the
    #same description are shown together:
    import textwrap
    groups = []
    for n,d in pd:
        if groups and groups[-1][1] == d:
            groups[-1][0].append( n )
        else:
            groups.append( ( [n], d ) )
    w = max( len(n[0]) for n,d in groups if len(n) == 1 )
    out = []
    for names,d in groups:
        ll = textwrap.wrap( d, width = 77 - w - 4 )
        if len(names) > 1:
            out.append( '  ' + ', '.join( names ) )
            out += [ ' '*(w+4)+e for e in ll ]
        else:
            out.append( f'  {names[0].ljust(w)}  {ll[0]}' )
            out += [ ' '*(w+4)+e for e in ll[1:] ]
    return '\n'.join(out)

def _cond_help_str( examples, docs ):
    txt = ( 'Conditions compare a property with a value, like "absxs>100"'
            ' (conditions with < or > need quotes). Text comparisons ignore'
            ' case. Yes/no values can also be given as true/false, on/off,'
            ' or 1/0.' )
    import textwrap
    return ( 'conditions:\n'
             + ''.join( f'  {ll}\n' for ll in textwrap.wrap( txt, 77 ) )
             + '\nexamples:\n'
             + '\n'.join( f'  {e.ljust(30)}{d}' for e,d in examples )
             + '\n\nproperties:\n' + _docs_str( docs ) )

def _factory_help_str():
    return """--factory SRC (or -f SRC) only shows data from the source SRC.

Common sources:

  stdlib    the standard data library of NCrystal
  stdpath   files in directories listed in NCRYSTAL_DATA_PATH
  relpath   files found relative to the current directory (only shown on
            explicit request)
  virtual   in-memory data (e.g. registered with registerInMemoryFileData)
  solid, gasmix, freegas
            data created on demand (e.g. "solid::B4C/2.52gcm3")

Plugins can provide additional sources.

Listed names include their source (e.g. "stdlib::Al_sg225.ncmat"), which can be
used to request data from a specific source."""

def _split_patterns( patterns ):
    #Returns ( name patterns, conditions ):
    from ._conditions import _is_condition
    names = [ p for p in patterns if not _is_condition( p ) ]
    return names, [ p for p in patterns if _is_condition( p ) ]

def parseArgs( progname, arglist, return_parser = False ):
    parser = _create_parser( progname, full = return_parser )
    if return_parser:
        return parser
    #Like grep, a plain --color (without "=WHEN") means --color=auto:
    arglist = [ ( '--color=auto' if a in ('--color','--colour') else a )
                for a in arglist ]
    args = parser.parse_intermixed_args( arglist )
    from . import browse as nb
    if args.full_help:
        _create_parser( progname, full = True ).print_help()
        parser.exit()
    if args.cond_help:
        print( _cond_help_str( _cond_examples,
                               nb.physics_props_doc( conditions = True ) ) )
        parser.exit()
    if 'help' in args.pattern:
        parser.print_help()
        parser.exit()
    if args.factory == 'help':
        print( _factory_help_str() )
        parser.exit()
    if args.extract and ( args.pattern or args.search or args.factory
                          or args.json ):
        parser.error('--extract can not be combined with other options.')
    args.names, args.conditions = _split_patterns( args.pattern )
    #Validate conditions up front, so problems are reported as usage errors:
    from ._conditions import _Condition
    from .exceptions import NCBadInput
    try:
        for c in args.conditions:
            _Condition( c, nb._propnames(), nb._propkinds )
    except NCBadInput as e:
        parser.error( str(e) )
    return args

class _Progress:
    #Progress indicator on stderr, only shown on a terminal and only if the
    #work takes more than a second.
    def __init__( self, what ):
        import sys
        import time
        self.__t0 = time.time()
        self.__what = what
        self.__shown = False
        self.__enabled = ( hasattr(sys.stderr,'isatty')
                           and sys.stderr.isatty() )
    def update( self, n, ntotal ):
        if not self.__enabled:
            return
        import sys
        import time
        t = time.time()
        if not self.__shown and t - self.__t0 < 1.0:
            return
        if self.__shown and t - self.__tlast < 0.1:
            return#limit update rate
        self.__shown, self.__tlast = True, t
        sys.stderr.write(f'\r{self.__what}: {n}/{ntotal}')
        sys.stderr.flush()
    def done( self ):
        if self.__shown:
            import sys
            sys.stderr.write('\r\x1b[K')
            sys.stderr.flush()

#Same options and synonyms as GNU grep:
_color_choices = { 'always' : True, 'yes' : True, 'force' : True,
                   'never' : False, 'no' : False, 'none' : False,
                   'auto' : None, 'tty' : None, 'if-tty' : None }

def _use_color( when ):
    #Decide whether to use ANSI color codes. NB: These env vars are general
    #conventions (not NCrystal specific), so are never namespaced.
    v = _color_choices[when]
    if v is not None:
        return v
    import os
    import sys
    env = os.environ
    if env.get('NO_COLOR'):
        return False
    if any( env.get(k,'0') not in ('','0')
            for k in ('FORCE_COLOR','CLICOLOR_FORCE') ):
        return True
    if env.get('TERM') == 'dumb':
        return False
    from ._common import _builtin_print, get_ncrystal_print_fct
    if get_ncrystal_print_fct() is not _builtin_print():
        return False#output redirected (e.g. captured via cli.run)
    if not ( hasattr(sys.stdout,'isatty') and sys.stdout.isatty() ):
        return False
    if sys.platform == 'win32' and not ( 'WT_SESSION' in env
                                         or 'TERM' in env ):
        return False#classic Windows consoles might not handle ANSI codes
    return True

def _match_sgr():
    #SGR code for matches, default and GREP_COLORS handling like GNU grep.
    import os
    sgr = '01;31'
    for part in os.environ.get('GREP_COLORS','').split(':'):
        k,_,v = part.partition('=')
        if k in ('mt','ms') and v and all( c.isdigit() or c==';' for c in v ):
            sgr = v
    return sgr

class _Highlighter:
    #Wraps all occurrences of the search words (case-insensitive) in color
    #codes.
    def __init__( self, words, enabled ):
        import re
        self.__res = ( [ re.compile( re.escape(w), re.IGNORECASE )
                         for w in words if w ] if enabled else [] )
        self.__start = f'\x1b[{_match_sgr()}m\x1b[K'
        self.__end = '\x1b[m\x1b[K'
    def __call__( self, s ):
        if not self.__res or not s:
            return s
        spans = sorted( m.span() for r in self.__res for m in r.finditer(s)
                        if m.end() > m.start() )
        merged = []
        for a,b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max( merged[-1][1], b )
            else:
                merged.append( [a,b] )
        out, pos = [], 0
        for a,b in merged:
            out += [ s[pos:a], self.__start, s[a:b], self.__end ]
            pos = b
        out.append( s[pos:] )
        return ''.join(out)

def create_argparser_for_sphinx( progname ):
    return parseArgs( progname, [], return_parser = True )

def _browser( args ):
    #Create DataBrowser with selection according to the arguments:
    from . import browse as nb
    progress = _Progress( 'Loading materials' )
    b = nb.DataBrowser( factory = args.factory, _progress = progress.update )
    args.all_browser = b#for suggestions
    try:
        sel = b
        if args.names:
            #Names match as substrings, and element and isotope names also
            #select data containing them:
            labels = nb._atom_labels()
            keys = { e.fullkey for e in b.match( *args.names ) }
            elems = [ n for n in args.names if n in labels ]
            if elems:
                keys.update( e.fullkey for e in b.with_elements( *elems ) )
            sel = b._filter( lambda e : e.fullkey in keys )
        sel = sel.search( *args.search ).where( *args.conditions )
    finally:
        progress.done()
    return sel

def _print_no_matches( args ):
    print('No matching data found.')
    sugg = _suggestions( args )
    if sugg:
        print(f'Did you mean: {", ".join(sugg)}?')

def _suggestions( args ):
    #Suggest similar names if name patterns alone matched nothing:
    patterns = args.names
    if not patterns:
        return []
    b = args.all_browser
    if b.match( *patterns ):
        return []
    res = []
    for p in patterns:
        res += [ n for n in _similar_names( b, p ) if n not in res ]
    return res

def _similar_names( entries, pattern ):
    #Names of entries similar to the pattern:
    import difflib
    def stem( n ):
        return n.rsplit('.',1)[0].lower() if '.' in n else n.lower()
    #Compare with full names (without extension), and with their parts
    #(e.g. "diamond" in "C_sg227_Diamond.ncmat"):
    cands = {}
    for e in entries:
        st = stem(e.name)
        cands.setdefault( st, e.name )
        for part in st.split('_'):
            if len(part) >= 4:
                cands.setdefault( part, e.name )
    res = []
    for m in difflib.get_close_matches( stem(pattern), list(cands), n = 3,
                                        cutoff = 0.7 ):
        if cands[m] not in res:
            res.append( cands[m] )
    return res

def _listing( sel, words, hl ):
    #One full key per line, followed by the header comment lines containing
    #any of the search words:
    lw = [ w.lower() for w in words ]
    out = []
    for e in sel:
        out.append( hl( e.fullkey ) )
        out += [ '    | ' + hl( ll.strip() ) for ll in ( e.comments or [] )
                 if lw and any( w in ll.lower() for w in lw ) ]
    return ''.join( ll + '\n' for ll in out )

@cli_entry_point
def main( progname, arglist ):
    args = parseArgs( progname, arglist )
    if args.extract:
        from .core import createTextData
        print( createTextData( args.extract ).rawData, end='' )
        return
    sel = _browser( args )
    if args.json:
        print( sel.to_json(), end = '' )
        return
    if not sel:
        _print_no_matches( args )
        return
    hl = _Highlighter( args.search, _use_color( args.color ) )
    print( _listing( sel, args.search, hl ), end = '' )

#FIXME: Not happy with this, way too verbose in the default listing.
