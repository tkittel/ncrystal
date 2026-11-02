
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

"""Implementation of the built-in tests."""

__all__ = ['test', 'test_all', 'test_cmdline']
from ._common import print as _nc_print


def test( verbose = False ):
    """Quick test that NCrystal works as expected in the current installation."""
    _actualtest( verbose = verbose )
    if verbose!='quiet':
        _nc_print("Tests completed succesfully")

def test_cmdline( verbose = False ):
    """Quick test that NCrystal command-line scripts are present and working in
    the current installation."""
    _actual_test_cmdline( verbose = verbose )
    if verbose!='quiet':
        _nc_print("Tests completed succesfully")

def test_all( verbose = False ):
    """Run both test() and test_cmdline()."""
    _actualtest( verbose = verbose )
    _actual_test_cmdline( verbose = verbose )
    if verbose!='quiet':
        _nc_print("Tests completed succesfully")

def _get_prfct( verbose ):
    if verbose and verbose != 'quiet':
        def prt(  *a, **kw ):
            kw['flush'] = True
            return _nc_print('::NCrystalTest::',*a,**kw)
        return prt
    else:
        return lambda *a, **kw : None

def _actualtest( verbose ):
    prfct = _get_prfct( verbose )
    prfct('starting standard Python-API testing')
    try:
        import numpy as _np
    except ImportError:
        _np = None
    from . import core as NC
    from .constants import wl2ekin
    from .plugins import hasFactory

    def require(b):
        if not b:
            raise RuntimeError('check failed')
    def flteq(a,b,rtol=1.0e-6,atol=1.0e-6):
        return abs(a-b) <= 0.5 * rtol * (abs(a) + abs(b)) + atol
    def require_flteq(a,b):
        if not flteq(a,b):
            raise RuntimeError(f'check failed ({a:.16g} != {b:.16g}, diff {a-b:g})')
        return True
    require(hasFactory('stdncmat'))
    from . import _fmt as nc_fmt
    require( nc_fmt.prettyFmtValue(0.25) == '1/4' )
    require( nc_fmt.prettyFmtValue(0.25 + 1e-12) != '1/4' )
    require( nc_fmt.prettyFmtValue(0.25 + 1e-15) == '1/4' )

    _cfgstr='stdlib::Al_sg225.ncmat;dcutoff=1.4'
    prfct(f'Trying to createInfo("{_cfgstr}")')
    al = NC.createInfo(_cfgstr)
    prfct('Verifying loaded Info object')
    require(al.hasTemperature() and require_flteq(al.getTemperature(),293.15))
    require_flteq(al.getXSectFree(),1.39667)
    require_flteq(al.getXSectAbsorption(),0.231)
    require_flteq(al.getDensity(),2.69864547673)
    require_flteq(al.getNumberDensity(),0.06023238256131625)
    require(al.hasDebyeTemperature())

    require(al.hasStructureInfo())
    si=al.getStructureInfo()
    require( si['spacegroup'] == 225 )
    require_flteq(si['a'],4.04958)
    require_flteq(si['b'],4.04958)
    require_flteq(si['c'],4.04958)
    require( si['alpha'] == 90.0 )
    require( si['beta'] == 90.0 )
    require( si['gamma'] == 90.0 )
    require( si['n_atoms'] == 4 )
    require_flteq(si['volume'],66.4094599932)
    require( al.hasHKLInfo() )
    require( al.nHKL() == 3 )
    require_flteq(al.hklDLower(),1.4)
    require( al.hklDUpper() > 1e36 )
    expected_hkl = { 0  : (1, 1, 1, 8, 2.3380261031049243, 1.7731590373262052),
                     1  : (2, 0, 0, 6, 2.02479, 1.7317882793764163),
                     2  : (2, 2, 0, 12, 1.4317427394787092, 1.5757351418107877) }
    for idx,hkl in enumerate(al.hklList()):
        h,k,l_,mult,dsp,fsq = hkl
        require(idx<len(expected_hkl))
        e = expected_hkl[idx]
        require( list(e)[0:4] == [h,k,l_,mult] )
        require_flteq(dsp, e[4])
        require_flteq(fsq, e[5])

    _cfgstr2='stdlib::Al_sg225.ncmat;dcutoff=1.4;incoh_elas=0;inelas=0'
    prfct(f'Trying to createScatter("{_cfgstr2}")')
    #We do all createScatter... here with independent RNG, for reproducibility
    #and to avoid consuming random numbers from other streams.
    alpc = NC.createScatterIndependentRNG(_cfgstr2)
    prfct('Verifying loaded Scatter object')
    require( alpc.name == 'PowderBragg' )
    require( isinstance(alpc.name,str) )
    require( alpc.refCount() in (1,2) )
    require( type(alpc.refCount()) is int )
    require( alpc.isNonOriented() )
    #_nc_print(alpc.xsect(wl=4.0))
    require_flteq(1.632435821586171,alpc.crossSectionIsotropic(wl2ekin(4.0)) )
    require_flteq(1.632435821586171,alpc.crossSection(wl2ekin(4.0),(1,0,0)))
    require( alpc.crossSectionIsotropic(wl2ekin(5.0)) == 0.0 )

    require( alpc.rngSupportsStateManipulation() )
    require(alpc.getRNGState()=='a79fd777407ba03b3d9d242b2b2a2e58b067bd44')

    alpc.setRNGState('deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    require(alpc.getRNGState()=='deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    alpc_clone = alpc.clone()
    require(alpc.getRNGState()=='deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    require(alpc_clone.getRNGState()=='e0fd16d42a2aced7706cffa08536d869b067bd44')
    alpc_clone2 = alpc_clone.clone(for_current_thread=True)
    require(alpc_clone2.getRNGState()=='cc762bb1160a0be514300da860f6d160b067bd44')
    alpc_clone3 = alpc_clone.clone(rng_stream_index = 12345 )
    require(alpc_clone3.getRNGState()=='3a20660a10fd581bd7cddef8fc3f32a2b067bd44')

    #Pick Nickel at 1.2 angstrom, to also both vdos + incoherent-elastic + coherent-elastic:
    _cfgstr3='stdlib::Ni_sg225.ncmat;dcutoff=0.6;vdoslux=2'
    _seed=2543577
    prfct(f'Trying to createScatterIndependentRNG("{_cfgstr3}",{_seed})')
    nipc = NC.createScatterIndependentRNG(_cfgstr3,_seed)
    prfct('Verifying loaded Scatter object')
    nipc_testwl = 1.2
    #_nc_print(nipc.xsect(wl=nipc_testwl),nipc.xsect(wl=5.0))
    require_flteq(16.76474410391571,nipc.xsect(wl=nipc_testwl))
    require_flteq(16.76474410391571,nipc.xsect(wl=nipc_testwl,direction=(1,0,0)))
    require_flteq(5.958467463288343,nipc.xsect(wl=5.0))

    require( nipc.name == 'ProcComposition' )

    expected = [ ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, -0.25360938329540783 ),
                 ( 0.056808478892590906, 0.8391056916029317 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, 0.8260541809964751 ),
                 ( 0.056808478892590906, 0.5783836868120302 ),
                 ( 0.056808478892590906, 0.07228896531453344 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.07976617780753162, -0.9254365194024979 ),
                 ( 0.056808478892590906, -0.3915665520281999 ),
                 ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, -0.8299176791481262 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, -0.16636790113624955 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.0440778421354693, -0.9838232729247124 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, 0.07228896531453344 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, 0.8826297782791799 ),
                 ( 0.056808478892590906, -0.5655123710317247 ),
                 ( 0.056808478892590906, 0.8632864868136316 ),
                 ( 0.056808478892590906, 0.07228896531453344 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, -0.15963879335683306 ),
                 ( 0.056808478892590906, 0.8260541809964751 ),
                 ( 0.056808478892590906, 0.07228896531453344 ),
                 ( 0.056808478892590906, -0.5655123710317247 ) ]

    if not _np:
        ekin,mu=[],[]
        for i in range(30):
            _ekin,_mu=nipc.sampleScatterIsotropic(wl2ekin(nipc_testwl))
            mu += [_mu]
            ekin += [_ekin]
    else:
        ekin,mu = nipc.sampleScatterIsotropic(wl2ekin(nipc_testwl),repeat=30)

    for i in range(len(ekin)):
        #print ( f'    ( {ekin[i]}, {mu[i]} ),');continue
        require_flteq(ekin[i],expected[i][0])
        require_flteq(mu[i],expected[i][1])

    expected = [ ( 0.056808478892590906, (-0.15963879335683306, -0.8670635323578942, -0.4719282641574095) ),
                 ( 0.056808478892590906, (-0.15963879335683306, -0.9853505670542051, -0.05999763046596856) ),
                 ( 0.056808478892590906, (0.733956301770001, -0.5862953549271661, -0.34287884723752143) ),
                 ( 0.056808478892590906, (0.36219866365374176, -0.21569056961469935, -0.9068019112390149) ),
                 ( 0.025631780579978305, (-0.9832870620764568, -0.09474117383381278, 0.1554691722935501) ),
                 ( 0.02355674950366303, (0.21067658085461424, -0.05628246202432177, 0.9759342512423109) ),
                 ( 0.056808478892590906, (0.3042167239859003, -0.1587592770142396, -0.9392804037182945) ),
                 ( 0.056808478892590906, (0.7672386484388488, 0.43335192552206203, 0.4728117648579081) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.7905105193171594, -0.6081672667471253) ),
                 ( 0.056808478892590906, (0.2855453142094114, -0.8544734300516215, -0.43398044986942114) ),
                 ( 0.056808478892590906, (0.07228896531453344, -0.39741541395284924, -0.914787021249449) ),
                 ( 0.056808478892590906, (0.3503361232953429, -0.9197415233860485, -0.17703144035425164) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.5148338899674894, -0.22930450453313347) ),
                 ( 0.041661065632835934, (-0.387441912842209, -0.6769735855838138, 0.625775941208165) ),
                 ( 0.056808478892590906, (-0.3915665520281999, 0.9125393050471081, -0.11810017814079642) ),
                 ( 0.056808478892590906, (0.5361444826572666, 0.7644444598215647, 0.358013633210089) ),
                 ( 0.056808478892590906, (-0.2350943803753518, 0.9708260859803082, -0.047194735893790966) ),
                 ( 0.056808478892590906, (0.7680722413286334, -0.4990243565808846, -0.4012975500042876) ),
                 ( 0.056808478892590906, (-0.3181901153046516, 0.11568324371218267, 0.9409423136657422) ),
                 ( 0.05346787094824556, (0.033236780550364886, 0.6578034043965023, -0.7524559771727637) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.7288673223223796, -0.6808279752950511) ),
                 ( 0.0729796424467488, (-0.45159038644090005, -0.8154315022452432, -0.36212924215012593) ),
                 ( 0.03140005709326457, (-0.2865595397907464, 0.9521950198823131, -0.1058691374586427) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.49528807392249224, 0.8627604580230847) ),
                 ( 0.056808478892590906, (-0.3915665520281999, 0.24454393258214657, 0.8870591301429717) ),
                 ( 0.056808478892590906, (0.5361444826572666, -0.23795484435765296, 0.8098929471003542) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.6839258999305058, -0.7224342513352359) ),
                 ( 0.056808478892590906, (-0.37568216004953503, -0.45939628236594093, 0.8048713998949573) ),
                 ( 0.056808478892590906, (0.3042167239859003, 0.6691911433212484, -0.6779641572735884) ),
                 ( 0.049137704161612276, (-0.1946367760174378, -0.9283828962854157, 0.31657814723420286) ) ]

    for i in range(30):
        out_ekin,outdir = nipc.sampleScatter(wl2ekin(nipc_testwl),(1.0,0.0,0.0))
        #print ( f'    ( {out_ekin}, {outdir} ),');continue
        require_flteq(out_ekin,expected[i][0])
        require_flteq(outdir[0],expected[i][1][0])
        require_flteq(outdir[1],expected[i][1][1])
        require_flteq(outdir[2],expected[i][1][2])

    _cfgstr4="""stdlib::Ge_sg227.ncmat;dcutoff=0.5;mos=40.0arcsec
                            ;dir1=@crys_hkl:5,1,1@lab:0,0,1
                            ;dir2=@crys_hkl:0,-1,1@lab:0,1,0""".replace(' ','').replace('\n','')
    _seed4 = 3453455
    prfct(f'Trying to createScatterIndependentRNG("{_cfgstr4}",{_seed4})')
    gesc = NC.createScatterIndependentRNG(_cfgstr4,_seed4)
    prfct('Verifying loaded Scatter object')
    require_flteq(591.026403658,gesc.crossSection(wl2ekin(1.540),( 0., 1., 1. )))
    require_flteq(1.66765659374,gesc.crossSection(wl2ekin(1.540),( 1., 1., 0. )))
    prfct('standard Python-API testing done')


class CallInspector:

    @property
    def name( self ):
        return self.__name

    @property
    def the_real_inspected_object( self ):
        return self.__realobj

    def __init__( self, name, *, realobj = None, subfcts = None, fmtcall = None ):
        self.__name = name
        self.__realobj = realobj
        self.__subfcts = {}
        self.__fmtcall = fmtcall or _fmtcall
        for sf in (subfcts or []):
            if isinstance(sf,str):
                self.__subfcts[sf] = None
            else:
                assert len(sf)==2 and isinstance(sf,tuple)
                self.__subfcts[sf[0]] = sf[1]

    def __call__( self, *a, **kwargs ):
        return self.__getattrimpl('__call__')(*a,**kwargs)

    def __getattr__( self, attrname ):
        return self.__getattrimpl(attrname)

    def __getattrimpl( self, attrname ):
        if attrname not in self.__subfcts:
            raise RuntimeError(f'Not allowed to access attribute {attrname} of {self.__name} objects.')
        sf_kwargs = self.__subfcts[attrname] or {}
        thefmtcall = sf_kwargs.get('fmtcall') or self.__fmtcall
        def wrapper( *args, **kwargs ):
            _n = f'{self.__name}.{attrname}' if attrname!='__call__' else self.__name
            _nc_print(f"CALLING {thefmtcall(_n,args,kwargs)}")
            res = getattr(self.__realobj,attrname)(*args,**kwargs) if self.__realobj else None
            return CallInspector( name = f'ResultOf[{_n}(..)]',
                                  realobj = res, **sf_kwargs ) if sf_kwargs else None
        return wrapper

_fmtvalue_default_ndigits = 3
def _fmtvalue( x, *, ndigits = _fmtvalue_default_ndigits ):
    import numbers
    if isinstance(x,numbers.Real) and not isinstance(x,numbers.Integral):
        fmtstring = f'%.{ndigits}g'
        s=fmtstring%x#only 3 to avoid too many spurious test issues due to FPE irrep.
        return s+'.0' if ( s.isdigit() or s[0]=='-' and s[1:].isdigit() ) else s
    return repr(x)

def _fmtcall(fctname,args=(),kwargs=None):
    import numbers
    def _fmt(a):
        if isinstance(a,numbers.Real):
            return _fmtvalue(a)
        def pruneaddr(s):
            while ' object at 0x' in s:
                _ = s.split(' object at 0x',1)
                while _[1] and _[1][0].isalnum():
                    _[1] = _[1][1:]
                s = ' object at SNIPADDR'.join(_)
            return s
        return f'Object[{a.name}]' if isinstance(a,CallInspector) else pruneaddr(repr(a))
    ll = [ _fmt(a) for a in args ]
    if kwargs is None:
        kwargs = {}
    ll += [ f'{k}={_fmt(v)}' for k,v in sorted(kwargs.items()) ]
    a=','.join(ll)
    return f'{fctname}({a})'

def _create_pyplot_inspector( pass_calls_to_real_plt ):
    import copy
    if pass_calls_to_real_plt:
        import matplotlib.pyplot as realplt
    else:
        realplt = None
    def shorten( x ):
        if hasattr(x,'shape') and len(x.shape)==2:
            return f'Array(shape={x.shape},content={shorten(x.flatten())})'
        if isinstance(x,str) or not hasattr(x,'__len__'):
            return x
        if len(x)==0:
            return x
        #Fewer digits to guard against annoying test errors:
        maxxval = max(x)
        def _fmtthislistval(val):
            if val==0.0:
                return '0.0'
            if abs(val)<abs(maxxval)*1e-13:
                return 'TINYVAL'
                #return _fmtvalue( val, ndigits = max(1,(_fmtvalue_default_ndigits-2)) )
            elif abs(val)<abs(maxxval)*1e-8:
                return 'SMALLVAL'
                #return _fmtvalue( val, ndigits = max(1,(_fmtvalue_default_ndigits-1)) )
            else:
                return _fmtvalue( val )
        if len(x) <= 10:
            return [_fmtthislistval(e) for e in x]
        else:
            return [_fmtthislistval(e) for e in x[0:3]]+['...']+[_fmtthislistval(e) for e in x[-3:]]
    def _create_shortening_fmtcall( nargs_to_shorten, kwargs_to_shorten = None ):
        def fmtcall_pltplot( name, args, kwargs ):
            plot_args, plot_kwargs = args, kwargs
            if nargs_to_shorten>0:
                plot_args = [ shorten(e) for e in args[0:nargs_to_shorten] ] + list(args[nargs_to_shorten:])
            if any( e in kwargs for e in (kwargs_to_shorten or [])):
                plot_kwargs = copy.deepcopy(kwargs)
                for xy in ('x','y','yerr'):
                    if xy in plot_kwargs:
                        plot_kwargs[xy] = shorten(plot_kwargs[xy])
            return _fmtcall(name, plot_args,plot_kwargs)
        return fmtcall_pltplot

    csf = _create_shortening_fmtcall

    fill_between = ( 'fill_between', { 'fmtcall': csf(3)} )
    plot = ( 'plot', { 'fmtcall': csf(2,('x','y'))} )
    bar = ( 'bar', { 'fmtcall': csf(1,('x',))} )
    errorbar = ( 'errorbar', { 'fmtcall': csf(3,('x','y','yerr',))} )
    get_figure = ( 'get_figure',{ 'subfcts': ['suptitle','colorbar'] } )
    pcolormesh = ( 'pcolormesh',{ 'subfcts': ['set_clim'], 'fmtcall': csf(3)} )
    axisfcts = [ pcolormesh,fill_between, plot, bar, errorbar,
                 'semilogx','semilogy',get_figure,
                 'set_title','set_xticks','set_xlabel','set_ylabel',
                 'set_ylim','set_xlim','legend','semilogx',
                 'semilogy', 'grid' ]

    return CallInspector( name = 'plt', realobj = realplt,
                          subfcts = [ 'show','subplots','savefig',
                                      'suptitle','tight_layout','close',
                                      ( 'gca',{ 'subfcts': axisfcts } ) ] )

def _create_pdfpages_inspector( real_pdfpages ):
    return CallInspector( name = 'PdfPages',
                          realobj = real_pdfpages,
                          subfcts = [ ('__call__',{'subfcts': ['savefig','close']}) ] )

def _run_cmd( cmd ):
    import subprocess
    import sys
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        p = subprocess.Popen( cmd, shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE )
        output = p.communicate()[0]
        ok = p.returncode==0
    except OSError:
        ok = False
    sys.stdout.flush()
    sys.stderr.flush()
    return ok, output

def _test_cmdline_script_availablity( prfct ):
    import shutil
    import subprocess

    from .cli import cli_tool_list, cli_tool_lookup
    for t in cli_tool_list():
        c = cli_tool_lookup(t)['shellcmd']
        prfct(f"Testing availability of command: {c}")
        if not shutil.which(c):
            raise RuntimeError(f'Command {c} not found!')
        ev = subprocess.run( [c,'--help'], capture_output = True,
                             check = False )
        if ev.returncode != 0:
            raise RuntimeError(f'Command "{c} --help" did not run succesfully!')

def _actual_test_cmdline( verbose ):
    import shlex

    from .cli import cli_tool_lookup
    prfct = _get_prfct( verbose )
    prfct('starting testing of cmd-line utilities')
    _test_cmdline_script_availablity( prfct )
    cmds = ['ncrystal-config --help',
            'ncrystal-config -s',
            'nctool --version',
            'nctool --help',
            'nctool --test',
            'ncrystal_endf2ncmat --help',
            'ncrystal_hfg2ncmat --help',
            'ncrystal_ncmat2cpp --help',
            'ncrystal_ncmat2hkl --help',
            'ncrystal_cif2ncmat --help',
            'ncrystal_vdos2ncmat --help',
            'ncrystal-config --show cmakedir']
    for cmd in cmds:
        cmd = shlex.split(cmd)
        cmd[0] = cli_tool_lookup( cmd[0] )['shellcmd']
        cmd = shlex.join(cmd)
        prfct('Trying to run:',cmd)
        ok, _output = _run_cmd(cmd)
        if not ok:
            raise RuntimeError(f'Command failed: {cmd}')
    prfct('testing of cmd-line utilities done')
