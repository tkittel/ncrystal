
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

# Monitoring of cross section curves: a test evaluates the cross sections of a
# list of cfg-strings on an energy grid and compares them with the values in a
# single reference file per test. The file holds one block per cfg-string
# ("@cfg <cfgstr>" followed by "<energy> <xs>" lines, or a single "@error <msg>"
# line for cfg-strings expected to fail), so changes show up as readable diffs.
#
# Run with --update to (re)write the reference file, and with --plot to see
# the curves. Extra arguments select a subset of cfg-strings by exact match.

import numpy as np
from NCrystalDev._numpy import _np_geomspace, _np_linspace
from NCrystalDev.exceptions import NCException

import NCTestUtils.dirs as testdirs
import NCTestUtils.enable_fpe  # noqa F401


class XSMonitor:

    def __init__( self, *, testname, cfglist, matloadfct, egridgenfct,
                  rdtolfct, datadir = 'xsref' ):
        """Monitor cross sections of the cfg-strings in cfglist. The
        rdtolfct(cfgstr) returns the acceptable relative difference, which
        also decides the number of digits stored in the reference file."""
        self.__do_plot, self.__do_update, self.__select = _parse_sysargv()
        self.__testname = testname
        self.__cfgs = list( cfglist )
        assert len(set(self.__cfgs)) == len(self.__cfgs), 'duplicate cfgs'
        assert all( '\n' not in c for c in self.__cfgs )
        self.__load = matloadfct
        self.__egrid = egridgenfct
        self.__rdtol = rdtolfct
        if not self.__do_plot:
            from NCTestUtils.env import ncsetenv
            ncsetenv('FAKEPYPLOT','1')
        d = testdirs.get_named_test_data_dir( datadir,
                                              for_updates = self.__do_update )
        self.__reffile = d.joinpath(f'{testname}.txt')

    def __ndigits( self, cfg ):
        #Store three digits more than needed for the tolerance:
        import math
        return math.ceil( -math.log10( self.__rdtol(cfg) ) ) + 3

    def __evaluate( self, cfg ):
        try:
            mat = self.__load( cfg )
            evals = sorted({ float(f'{e:.7g}') for e in self.__egrid(mat) })
            evals = np.asarray( evals, dtype = float )
            return mat, evals, mat.scatter.xsect( evals ), None
        except NCException as e:
            msg = f'{e.__class__.__name__}: {e}'.replace('\n',' ')
            return None, None, None, msg

    def __format_block( self, cfg, evals, xsvals, errmsg ):
        lines = [ f'@cfg {cfg}' ]
        if errmsg is not None:
            lines.append( f'@error {errmsg}' )
        else:
            n = self.__ndigits( cfg )
            lines += [ f'{e:.7g} {x:.{n}g}' for e, x in zip(evals,xsvals) ]
        return lines

    def __read_ref( self ):
        if not self.__reffile.is_file():
            return {}
        res, cur = {}, None
        for line in self.__reffile.read_text().splitlines():
            if not line or line.startswith('#'):
                continue
            if line.startswith('@cfg '):
                cur = line[5:]
                assert cur not in res
                res[cur] = { 'error' : None, 'data' : [] }
            elif line.startswith('@error '):
                res[cur]['error'] = line[7:]
            else:
                e, x = line.split()
                res[cur]['data'].append( ( float(e), float(x) ) )
        return res

    def __write_ref( self, blocks ):
        import sys
        script = sys.argv[0].split('/')[-1] if sys.argv else '<test script>'
        out = [ ( f'# Reference cross sections [barn/atom] vs. energy [eV] for'
                  f' the test "{self.__testname}".' ),
                f'# Regenerate with: {script} --update' ]
        for cfg in self.__cfgs:
            out += blocks[cfg]
        self.__reffile.write_text( '\n'.join(out) + '\n' )
        print(f'Wrote: {self.__reffile}')

    def run( self ):
        from .common import calc_reldiff
        ref = self.__read_ref()
        sel = self.__select
        if sel and not sel.issubset(self.__cfgs):
            unknown = ', '.join(sorted(sel-set(self.__cfgs)))
            raise SystemExit(f'ERROR: Selected unknown cfg-strings: {unknown}')
        if not self.__do_update:
            missing = [ c for c in self.__cfgs if c not in ref ]
            excess = [ c for c in ref if c not in self.__cfgs ]
            if missing or excess:
                for c in missing:
                    print(f'ERROR: no reference data for: {c}')
                for c in excess:
                    print(f'ERROR: obsolete reference data for: {c}')
                raise SystemExit('ERROR: reference file out of sync with the'
                                 ' test list (rerun with --update).')
        blocks, issues = {}, []
        for cfg in self.__cfgs:
            if sel and cfg not in sel:
                if cfg in ref:
                    r = ref[cfg]
                    blocks[cfg] = ( [f'@cfg {cfg}']
                                    + ([f'@error {r["error"]}'] if r['error']
                                       else [ f'{e:.7g} {x:.{self.__ndigits(cfg)}g}'
                                              for e,x in r['data'] ]) )
                continue
            print(f'=====> {cfg}')
            mat, evals, xsvals, errmsg = self.__evaluate( cfg )
            blocks[cfg] = self.__format_block( cfg, evals, xsvals, errmsg )
            if self.__do_update:
                continue
            r = ref[cfg]
            if errmsg is not None or r['error'] is not None:
                if errmsg != r['error']:
                    issues.append( ( None, cfg ) )
                    print(f'ERROR: expected error {r["error"]!r}, got {errmsg!r}')
                continue
            #Compare against the values as they would be stored:
            stored = [ tuple(map(float,ln.split())) for ln in blocks[cfg][1:] ]
            ref_e = np.asarray([ e for e,x in r['data'] ])
            ref_x = np.asarray([ x for e,x in r['data'] ])
            cur_e = np.asarray([ e for e,x in stored ])
            cur_x = np.asarray([ x for e,x in stored ])
            rdtol = self.__rdtol( cfg )
            if cur_e.shape != ref_e.shape or not np.array_equal(cur_e,ref_e):
                issues.append( ( None, cfg ) )
                print('ERROR: energy grid differs from reference')
                continue
            rda = calc_reldiff( cur_x, ref_x )
            if rda.max() > rdtol:
                issues.append( ( rda.max(), cfg ) )
                print(f'ERROR: cross sections differ by up to {rda.max():g}'
                      f' (tolerance {rdtol:g}):')
                for e, xr, x, rd in zip( cur_e, ref_x, cur_x, rda ):
                    flag = '  <-- problem' if rd > rdtol else ''
                    print(f'   E={e:<12.7g} ref={xr:<16.10g} new={x:<16.10g}'
                          f' reldiff={rd:.3g}{flag}')
            if self.__do_plot:
                _plot( cfg, mat, cur_e, cur_x, ref_e, ref_x, rdtol )
        if self.__do_update:
            for c in self.__cfgs:
                if c not in blocks:
                    raise SystemExit(f'ERROR: no data for {c}')
            self.__write_ref( blocks )
            return
        if not issues:
            print('ALL OK')
            return
        print()
        print('ERROR - ISSUES DETECTED:')
        for rd, cfg in sorted( issues, key = lambda x : (x[0] is None, x[0] or 0) ):
            print(f'  {cfg}: ' + ( f'reldiff {rd:g}' if rd is not None
                                   else 'error or grid mismatch' ))
        print('Rerun with --plot or --update (optionally with selected cfgs).')
        raise SystemExit(1)

def _plot( cfg, mat, evals, xsvals, ref_evals, ref_xsvals, rdtol ):
    import matplotlib.pyplot as plt

    from .common import calc_reldiff
    elux = set(evals) | set(_np_linspace(evals[0],evals[-1],5000))
    elux |= set(_np_geomspace(evals[0],evals[-1],5000))
    elux = np.asarray(sorted(elux))
    fig, (ax, axdiff) = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax.set_title(cfg)
    ax.plot(elux,mat.scatter.xsect(elux),'-')
    ax.plot(evals,xsvals,'o',alpha=0.5,label='current')
    ax.plot(ref_evals,ref_xsvals,'d',alpha=0.5,label='ref')
    ax.legend().set_draggable(True)
    ax.loglog()
    ax.grid()
    ax.set_ylabel('XS [barn/atom]')
    axdiff.plot(evals,calc_reldiff(xsvals,ref_xsvals),label='observed difference')
    axdiff.loglog()
    axdiff.grid()
    axdiff.axhline(rdtol,color='red',ls=':',label='test tolerance')
    axdiff.set_ylabel('XS relative difference')
    axdiff.set_xlabel('Neutron energy [eV]')
    axdiff.legend().set_draggable(True)
    plt.show()

def _parse_sysargv(args=None):
    import sys
    if args is None:
        args = sys.argv[1:]
    args = list(args)
    do_plot, do_update = False, False
    while '--plot' in args:
        args.remove('--plot')
        do_plot = True
    while '--update' in args:
        args.remove('--update')
        do_update = True
    return do_plot, do_update, set(args)
