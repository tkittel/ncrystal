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


"""TEMPORARY tool: textual migration of .ncmat files to NCMAT v8, keeping all
comments (used to migrate the stdlib files for release 5.0.0).

  ncmat_migrate_v8.py IN.ncmat OUT.ncmat [--som STATE]

Changes: header -> v8; @SPACEGROUP -> explicit setting and @ATOMPOSITIONS ->
@ATOMSITES (via the NCrystal "util atomsites relaxed" query); vdos densities
get an explicit final 0 (one bin beyond the last point, as NCrystal appends
for older files anyway); non-crystalline files without @STATEOFMATTER get one
(as NCrystal classifies the original file, or --som). Prints notes to stderr.
"""

import importlib
import importlib.util
import re
import sys
from fractions import Fraction

#Works with an installed NCrystal or (in the simplebuild environment)
#NCrystalDev:
_modname = ( 'NCrystalDev' if importlib.util.find_spec('NCrystalDev')
             else 'NCrystal' )
NC = importlib.import_module( _modname )
evaluate_query = importlib.import_module( f'{_modname}.misc' ).evaluate_query
_impl = importlib.import_module( f'{_modname}._ncmatimpl' )
_standardise_vdos = _impl._standardise_vdos
_fmt_value_rtol = importlib.import_module( f'{_modname}._fmt' ).fmt_value_rtol
_uses_rhombohedral_axes = _impl._uses_rhombohedral_axes
formatVectorForNCMAT = _impl.formatVectorForNCMAT

def note(*a):
    print('NOTE:', *a, file=sys.stderr)

def split_sections(lines):
    """[(name or None for preamble, [lines])]"""
    out = [[None, []]]
    for line in lines:
        m = re.match(r'^@([A-Z_]+)\s*(#.*)?$', line.strip())
        if m and not line.lstrip().startswith('#'):
            out.append([m.group(1), [line]])
        else:
            out[-1][1].append(line)
    return out

def data_tokens(line):
    return line.split('#', 1)[0].split()

def parse_num(tok):
    if '/' in tok:
        return float(Fraction(tok))
    return float(tok)

def fmt_coord(v, orig_tokens):
    v = v % 1.0
    if abs(v - 1.0) < 1e-12:
        v = 0.0
    for den in (1, 2, 3, 4, 6, 8, 12, 16, 24):
        n = round(v * den)
        if abs(v - n / den) < 1e-12:
            if n == 0:
                return '0'
            f = Fraction(n, den)
            return str(f.numerator) if f.denominator == 1 else f'{f.numerator}/{f.denominator}'
    for t in orig_tokens:
        try:
            if parse_num(t) == v:
                return t
        except ValueError:
            pass
    return _fmt_value_rtol(v, 1e-14)

def migrate(text, som_override=None, origpath=None):
    lines = text.splitlines()
    m = re.match(r'^NCMAT\s+v(\d+)(.*)$', lines[0])
    assert m, 'bad header'
    oldver = int(m.group(1))
    lines[0] = 'NCMAT v8' + m.group(2)
    secs = split_sections(lines[1:])
    names = [s[0] for s in secs]
    assert 'DEBYETEMPERATURE' not in names, 'DEBYETEMPERATURE not handled'
    is_crystal = 'CELL' in names
    # --- space group and sites ---
    if 'SPACEGROUP' in names and 'ATOMPOSITIONS' in names:
        sgsec = secs[names.index('SPACEGROUP')]
        sgtoks = [t for e in sgsec[1][1:] for t in data_tokens(e)]
        assert len(sgtoks) == 1, sgtoks
        sgnum = int(sgtoks[0].split(':')[0])
        # cell alpha for rhombohedral detection:
        info = NC.createInfo(origpath) if origpath else None
        alpha = info.structure_info['alpha'] if info else 90.0
        spec = f'{sgnum}:r' if _uses_rhombohedral_axes(sgnum, alpha) else str(sgnum)
        apsec = secs[names.index('ATOMPOSITIONS')]
        pos, orig_tokens = [], []
        for line in apsec[1][1:]:
            t = data_tokens(line)
            if not t:
                continue
            assert len(t) == 4, line
            pos.append((t[0], *[parse_num(x) for x in t[1:]]))
            orig_tokens += t[1:]
        q = ['util', 'atomsites', 'relaxed', spec]
        q += [f'{lbl},{x!r},{y!r},{z!r}' for lbl, x, y, z in pos]
        r = evaluate_query(q)
        setting = r['spacegroup_setting']
        nsites = sum(n for _, _, n in r['sites'])
        assert nsites == len(pos), (nsites, len(pos))
        # keep comment lines of the sections (not trailing comments):
        sgcomments = [e for e in sgsec[1][1:] if e.strip().startswith('#')]
        sgsec[1] = [sgsec[1][0]] + sgcomments + [f'  {setting}']
        apcomments = [e for e in apsec[1][1:] if e.strip().startswith('#')]
        sitelines = [f'  {lbl} ' + ' '.join(fmt_coord(c, orig_tokens) for c in xyz)
                     for lbl, xyz, _n in r['sites']]
        apsec[0] = 'ATOMSITES'
        apsec[1] = ['@ATOMSITES'] + apcomments + sitelines
        if setting != str(sgnum):
            note(f'spacegroup {sgtoks[0]} -> {setting}')
        note(f'{len(pos)} positions -> {len(r["sites"])} sites')
    elif 'SPACEGROUP' in names:
        note('SPACEGROUP without ATOMPOSITIONS?!')
    # --- vdos final zero ---
    for s in secs:
        if s[0] != 'DYNINFO':
            continue
        body = s[1]
        toks = [t for e in body[1:] for t in data_tokens(e)]
        if 'vdos_density' not in toks:
            continue
        # locate keyword blocks by line ranges:
        kwline = {}
        cur = None
        for i, line in enumerate(body[1:], start=1):
            t = data_tokens(line)
            if t and re.match(r'^[a-z_]+$', t[0]) and not re.match(r'^[0-9.eE+-]+$', t[0]):
                cur = t[0]
                kwline.setdefault(cur, []).append(i)
            elif cur:
                kwline.setdefault(cur, []).append(i)
        def vals(kw, body, kwline):
            v = []
            for i in kwline[kw]:
                t = data_tokens(body[i])
                if t and t[0] == kw:
                    t = t[1:]
                for x in t:
                    if 'r' in x:
                        a, b = x.split('r')
                        v += [parse_num(a)] * int(b)
                    else:
                        v.append(parse_num(x))
            return v
        eg = vals('vdos_egrid', body, kwline)
        dens = vals('vdos_density', body, kwline)
        orig_eg_tokens = [t for i in kwline['vdos_egrid'] for t in data_tokens(body[i]) if t != 'vdos_egrid']
        if dens[-1] == 0.0:
            continue
        neg, ndens = _standardise_vdos(eg, dens)
        indent = re.match(r'^(\s*)', body[kwline['vdos_density'][0]]).group(1) or '  '
        def rewrite(kw, newvals, body, kwline, indent, orig_eg_tokens):
            idx = kwline[kw]
            if kw == 'vdos_egrid' and len(newvals) == 2:
                #full precision for the computed upper value (the loader
                #computes the same double for older files):
                newl = [f'{indent}vdos_egrid {orig_eg_tokens[0]} {newvals[1]!r}']
            else:
                newl = formatVectorForNCMAT(kw, newvals, indent).rstrip('\n').split('\n')
            # drop comment-only lines inside the block? keep them after
            comm = [body[i] for i in idx if body[i].strip().startswith('#')]
            for i in reversed(idx):
                del body[i]
            for j, nl in enumerate(newl + comm):
                body.insert(idx[0] + j, nl)
        # rewrite the later keyword first (indices of the earlier stay valid)
        order = sorted(['vdos_egrid', 'vdos_density'], key=lambda k: -kwline[k][0])
        for kw in order:
            rewrite(kw, neg if kw == 'vdos_egrid' else ndens,
                    body, kwline, indent, orig_eg_tokens)
            # recompute indices for remaining keyword
            kwline2 = {}
            cur = None
            for i, line in enumerate(body[1:], start=1):
                t = data_tokens(line)
                if t and re.match(r'^[a-z_]+$', t[0]):
                    cur = t[0]
                    kwline2.setdefault(cur, []).append(i)
                elif cur:
                    kwline2.setdefault(cur, []).append(i)
            kwline = kwline2
        note('vdos: appended final 0 density')
    # --- state of matter ---
    if not is_crystal and 'STATEOFMATTER' not in names:
        som = som_override
        if som is None and origpath:
            st = NC.createInfo(origpath).stateOfMatter()
            som = {NC.Info.StateOfMatter.Solid: 'solid',
                   NC.Info.StateOfMatter.Liquid: 'liquid',
                   NC.Info.StateOfMatter.Gas: 'gas'}.get(st)
        if som is None:
            raise RuntimeError('unknown state of matter (use --som)')
        # insert before first DYNINFO:
        i = names.index('DYNINFO')
        secs.insert(i, ['STATEOFMATTER', ['@STATEOFMATTER', f'  {som}']])
        note(f'added @STATEOFMATTER {som}')
    out = [lines[0]]
    for s in secs:
        out += s[1]
    return '\n'.join(out) + '\n', oldver

if __name__ == '__main__':
    a = sys.argv[1:]
    som = None
    if '--som' in a:
        i = a.index('--som')
        som = a[i+1]
        del a[i:i+2]
    src, dst = a
    with open(src) as fh:
        txt = fh.read()
    new, oldver = migrate(txt, som, src)
    with open(dst, 'w') as fh:
        fh.write(new)
