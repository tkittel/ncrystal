
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
    require( si['n_positions'] == 4 )
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
    _cfgstr3='stdlib::Ni_sg225.ncmat;dcutoff=0.6;vdoslux=2002;knllux=3'
    _seed=2543577
    prfct(f'Trying to createScatterIndependentRNG("{_cfgstr3}",{_seed})')
    nipc = NC.createScatterIndependentRNG(_cfgstr3,_seed)
    prfct('Verifying loaded Scatter object')
    nipc_testwl = 1.2
    #_nc_print(nipc.xsect(wl=nipc_testwl),nipc.xsect(wl=5.0))
    require_flteq(16.75748060461658,nipc.xsect(wl=nipc_testwl))
    require_flteq(16.75748060461658,nipc.xsect(wl=nipc_testwl,direction=(1,0,0)))
    require_flteq(5.96098519807626,nipc.xsect(wl=5.0))

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
                 ( 0.08378876515740953, -0.666237855402724 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, -0.16636790113624955 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.03556468458983686, -0.1546394403801546 ),
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
                 ( 0.056808478892590906, -0.5655123710317247 ),
                 ( 0.056808478892590906, -0.15963879335683306 ),
                 ( 0.056808478892590906, 0.8560205329926969 ),
                 ( 0.056808478892590906, -0.15963879335683306 ),
                 ( 0.08829603560816682, 0.7786896724702138 ) ]

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

    expected = [ ( 0.056808478892590906, (0.22033427485529333, -0.5075901639685345, 0.8329495979748045) ),
                 ( 0.05577002875736579, (0.712545541391594, 0.701489077419671, -0.013853725270081366) ),
                 ( 0.06210445905099917, (-0.6334954702075752, 0.3584528519655151, 0.6857076943889953) ),
                 ( 0.056808478892590906, (-0.3915665520281999, 0.12100356351333735, -0.9121588529142393) ),
                 ( 0.056808478892590906, (0.36219866365374176, -0.6430608594398735, 0.6747479967394456) ),
                 ( 0.056808478892590906, (-0.3915665520281999, -0.9196076837100372, 0.0315807440413816) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.539797243436807, 0.16202909009269678) ),
                 ( 0.056808478892590906, (0.36219866365374176, 0.9243915901826631, -0.11963409232744321) ),
                 ( 0.056808478892590906, (0.5361444826572666, 0.7644444598215647, 0.358013633210089) ),
                 ( 0.056808478892590906, (-0.2350943803753518, 0.9708260859803082, -0.047194735893790966) ),
                 ( 0.056808478892590906, (0.7680722413286334, -0.4990243565808846, -0.4012975500042876) ),
                 ( 0.056808478892590906, (-0.3181901153046516, 0.11568324371218267, 0.9409423136657422) ),
                 ( 0.056808478892590906, (0.11073398686077585, 0.009122144428738927, 0.9938082162243068) ),
                 ( 0.056808478892590906, (-0.5655123710317247, 0.6187696974560152, 0.5452704097237153) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.5615367611788097, -0.04807240271772545) ),
                 ( 0.031168296780048017, (-0.6176396811075454, -0.6503742443577937, 0.44219290654351134) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.7022089766740898, -0.7082915067767521) ),
                 ( 0.056808478892590906, (0.7680722413286334, 0.15586048220243717, -0.6211059025544946) ),
                 ( 0.056808478892590906, (-0.4534627002970746, 0.378312081625071, 0.8070015788929976) ),
                 ( 0.056808478892590906, (-0.6647130752435024, 0.46171929462796124, -0.5873430178086598) ),
                 ( 0.056808478892590906, (-0.5655123710317247, 0.28806649539030227, 0.7727958672532003) ),
                 ( 0.056808478892590906, (0.8260541809964751, -0.2793761820760892, 0.4894726130712837) ),
                 ( 0.056808478892590906, (0.3042167239859003, 0.6691911433212484, -0.6779641572735884) ),
                 ( 0.07427614208411026, (-0.680523862628922, -0.1250613759344416, 0.7219743240877728) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.7112385472736513, -0.6955613653518439) ),
                 ( 0.056808478892590906, (0.7680722413286334, -0.6392725088477443, 0.03735895517704262) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.2489220191768986, -0.9631737706493833) ),
                 ( 0.056808478892590906, (-0.15963879335683306, 0.933317892436944, -0.321610272399119) ),
                 ( 0.056808478892590906, (-0.34161995982477755, -0.4695432196737952, -0.8141406315297645) ),
                 ( 0.056808478892590906, (0.36219866365374176, 0.17637999530837417, -0.9152607417028558) ) ]

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
