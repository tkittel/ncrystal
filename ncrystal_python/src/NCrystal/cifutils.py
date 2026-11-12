
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

"""Utilities for converting local or online CIF data to NCMAT data.

Two utility classes (CIFSource and CIFLoader) are included, and along with the
NCMATComposer they are used to implement the related cmd-line scripts (in
particular ncrystal_cif2ncmat), but they can of course also be used directly via
the Python API below by expert users.

"""

__all__ = [
    'CIFLoader',
    'CIFSource',
    'produce_validation_plot',
    'produce_validation_plots',
]

from . import _common as _nc_common
from . import _sgsymimpl as _nc_sgsymimpl
from . import core as _nc_core
from ._cifcomposerimpl import _impl_create_ncmat_composer
from ._cifloadimpl import _actual_init_gemmicif
from ._cifsrcimpl import _cod_get_cifdata, _mp_get_cifdata


class CIFSource:

    """
    Generalised CIF source, either in the form of text data, a file path, a
    database ID (in the form of a string like "codid::9008460" or "mpid::87)",
    for either the Crystallography Open Database
    (https://www.crystallography.net/cod/) or the Materials Project
    (https://www.materialsproject.org/), an existing CIFSource object, or any
    object like CIFLoader which has a cifsrc property which is a CIFSource
    object.
    """

    def __init__( self, data, *, allow_fail = False, name = None ):
        """Initialise from data in various formats (see class
        description). Unless allow_fail=True, an unrecognised input format will
        result in an NCBadInput exception being thrown. If allow_fail=True the
        .invalid property can be used to check if the loading failed.

        The name argument can be used to assign a name to otherwise anonymous
        text data. It can also be used to override the name.

        """
        self.__codid = None
        self.__mpid = None
        self.__fp = None
        self.__textdata = None
        self.__name_override = name

        def check_status():
            if ( not allow_fail ) and self.invalid:
                raise _nc_core.NCBadInput('Could not detect CIF source type')

        #Check if data is CIFSource object or has one as a .cifsrc property (e.g. LoadedCIF, CIFAnalyser)
        o = getattr( data, 'cifsrc', None ) or data
        if isinstance(o,CIFSource):
            self.__codid = o.__codid
            self.__mpid = o.__mpid
            self.__fp = o.__fp
            self.__textdata = o.__textdata
            self.__name_override = o.__name_override
            return
        def _setfp( pth ):
            pth = pth.decode() if hasattr(pth,'decode') else pth
            self.__fp = _nc_common._lookup_existing_file( pth )
            if self.__fp is None:
                #Try to look up via NCrystal's TextData infrastructure:
                from .core import NCFileNotFound, createTextData
                try:
                    td = createTextData(pth)
                    tdname = td.dataSourceName
                except NCFileNotFound:
                    td = None
                    tdname = None
                if td is not None:
                    self.__textdata = td.rawData
                    if tdname and not self.__name_override:
                        self.__name_override = tdname
                    self.__fp = None
            if not allow_fail and self.__fp is None:
                raise _nc_core.NCBadInput('Could not detect CIF source type'
                                          ' (tried to load from path, but'
                                          f' could not locate "{pth}")')

        if hasattr(data,'__fspath__'):
            _setfp(data)
            check_status()
            return
        if hasattr( data, 'startswith' ):
            str_to_sb = ( lambda s : s ) if not hasattr(data,'decode') else ( lambda s : s.encode() )

            if data.startswith(str_to_sb('codid::')) and data[7:].isdigit():
                self.__codid = int(data[7:])
                check_status()
                return
            if data.startswith(str_to_sb('mpid::')) and data[6:].isdigit():
                self.__mpid = int(data[6:])
                check_status()
                return
        if isinstance( data, str ):
            if '\n' in data.strip():
                self.__textdata = data
            else:
                _setfp(data)
            check_status()
            return
        if isinstance( data, bytes ):
            if b'\n' in data.strip():
                self.__textdata = data.decode()
            else:
                _setfp(data)
            check_status()
            return
        check_status()

    @property
    def invalid( self ):
        """Check if instance is invalid (only possible if constructed with
        allow_fail=True). """
        return all(e is None for e in (self.__codid,self.__mpid,
                                       self.__fp,self.__textdata))

    @property
    def codid( self ):
        """
        None or a database ID (integer) for the Crystallography Open Database
        (https://www.crystallography.net/cod/).
        """
        return self.__codid

    @property
    def mpid( self ):
        """
        None or a database ID (integer) for the Materials Project
        (https://www.materialsproject.org/).
        """
        return self.__mpid

    @property
    def textdata( self ):
        """
        None or a string containing raw CIF data.
        """
        return self.__textdata

    @property
    def name( self ):
        """A short string describing the data. This might be a file name
        (without directory part), or a special string like 'mpid::xyz'
        'codid::xyz'. Returns None if no name is available.
        """
        if self.__name_override:
            return self.__name_override
        if self.__mpid:
            return f'mpid::{self.__mpid}'
        if self.__codid:
            return f'codid::{self.__codid}'
        if self.__fp:
            return self.__fp.name

    @property
    def filepath( self ):
        """
        None or a pathlib.Path to an on-disk CIF file.
        """
        return self.__fp

    @property
    def is_remote( self):
        """True if either .codid or .mpid is not None"""
        return any( ( e is not None ) for e in (self.__codid,self.__mpid) )

    def load_data( self, quiet = False, mp_apikey = None ):
        """Try to load the data and return it as a string. Depending on the
        source, this might result in reading a file or querying an online
        database. In case the data is defined via a materials project ID (mpid),
        a materials project API key must be provided either via the mp_apikey
        parameter, or in the MATERIALSPROJECT_USER_API_KEY environment variable.
        """
        if self.invalid:
            raise _nc_core.NCDataLoadError( 'Can not load_data'
                                            ' from invalid CIF source' )
        if self.__textdata is not None:
            return self.__textdata
        if self.__fp is not None:
            if not quiet:
                _nc_common.print(f"Loading data from file {self.__fp.name}")
            return self.__fp.read_text()
        if self.__mpid is not None:
            return _mp_get_cifdata( self.__mpid, quiet = quiet, apikey = mp_apikey )
        assert self.__codid
        return _cod_get_cifdata( self.__codid, quiet = quiet )

class CIFLoader:

    """Class which is used to load and analyse CIF data, primarily so it can be
    used to transfer the crystal structure into an NCMATComposer object, with
    the intention of producing NCMAT data and new materials for usage with
    NCrystal.

    Behind the scenes the class uses third-party projects gemmi and spglib to
    parse the CIF data and process the symmetries within.

    The perhaps biggest issue to understand when basing NCrystal materials on
    CIF data, is that usually CIF data does not contain any or only incomplete
    information about the dynamics of the atoms. Thus, for a complete
    description of a high-quality material, one might want to pair up the usage
    of a CIFLoader, with a PhononDOSAnalyser from the NCrystal.vdos
    module.

    However, some CIF data contain information about atomic displacements
    (either "Uiso" or the anisotropic equivalent), which NCrystal is able to
    convert into a Debye temperature in order to provide a temperature-dependent
    model of dynamics and atomic displacement. However, such "Uiso" information
    is only useful if one knows at which temperature it was obtained, and that
    information is almost always absent in most CIF data out there. Thus, in the
    last step - when using the .create_ncmat(..) or .create_ncmat_composer(..)
    methods of the CIFLoader object, one should ideally use the uiso_temperature
    parameter to provide this information to NCrystal (usually after a call to
    print(cifloader.raw_cif_textdata) to manually inspect the meta-data, read
    any linked journal papers, etc.).

    TODO: Add link here to online tutorial once it is ready.

    """

    def __init__( self, cifsrc, quiet = False, mp_apikey = None,
                  refine_with_spglib = True, merge_equiv = True,
                  override_spacegroup = None, *, tol = None,
                  tol_symprec = None, tol_uij = None, tol_occu = None,
                  permissive = False ):
        """Initialise from cifsrc which can either be a CIFSource object, or any
        sort of data which can be used to initialise a CIFSource object (cf. the
        CIFSource class documentation)

        If refine_with_spglib=False, the structure will not be verified and
        refined with spglib. It is strongly recommended to not disable this
        refinement, as it is needed to ensure correctness and consistency of the
        resulting material.

        The override_spacegroup parameter can in some cases be used to select a
        particular setting of the detected spacegroup. For instance, spacegroup
        number 70 is available in two settings, and loading a CIF file with
        those might result in a warning that both "F d d d:1" and "F d d d:2"
        are available. Trying again with override_spacegroup="F d d d:1" can be
        used to explicitly select one of them (it might require some
        investigation of course to determine which one is right). A list of
        space group settings can be found at:

        https://cci.lbl.gov/sginfo/hall_symbols.html

        Setting merge_equiv=False will prevent usage of the same NCMAT label for
        multiple lists of equivalent atomic positions. However, note that this
        might affect the symmetry and space group of the structure, and is in
        general not recommended.

        If quiet=True, no informative messages will be emitted.

        Tolerances: CIF data might give symmetry equivalent sites separately
        (e.g. with space group P1), in which case the symmetry is detected by
        spglib. Sites are then only considered equivalent if their positions
        (within tol_symprec, in Aa) as well as their compositions, occupancies
        (within the relative tolerance tol_occu), Uiso values and the
        eigenvalues of their anisotropic displacement tensors (both within the
        relative tolerance tol_uij) agree. Additionally, sites of the same
        composition, with occupancies and Uiso values agreeing within the
        tolerances, are given the same NCMAT label (unless merge_equiv=False).
        Unless given, the tolerances follow from the overall tolerance tol:
        tol_uij and tol_occu equal tol, and tol_symprec scales with tol, being
        0.01Aa for the default tol.

        CIF data with inconsistent information (e.g. symmetry operations not
        matching the space group symbol, site multiplicities not matching the
        symmetry, or anisotropic displacements of equivalent sites which are
        not each other's symmetry images) is rejected. With permissive=True,
        some such data is loaded anyway, with warnings, if the information
        needed can be recovered (e.g. by discarding anisotropic
        displacements). Use this with care.

        Finally, the mp_apikey can be used to specify an API key in case the
        cifsrc is an entry for the materialsproject.org. To avoid having to
        specify the key here, one can instead provide the key in the environment
        variable MATERIALSPROJECT_USER_API_KEY.
        """
        tols = _nc_sgsymimpl._resolve_tolerances( tol,
                                                  tol_symprec = tol_symprec,
                                                  tol_uij = tol_uij,
                                                  tol_occu = tol_occu )
        with _nc_common.WarningSpy( block = quiet ) as warnlist:
            load, gemmi, struct = _actual_init_gemmicif( cifsrc,
                                                         quiet = quiet,
                                                         mp_apikey = mp_apikey,
                                                         refine_with_spglib = refine_with_spglib,
                                                         merge_equiv = merge_equiv,
                                                         override_spacegroup = override_spacegroup,
                                                         tols = tols,
                                                         permissive = permissive )
        assert len(load)==12
        self.__cifsrc = load['cifsrc']
        self.__actual_codid = load['actual_codid']
        self.__actual_mpid = load['actual_mpid']
        self.__cif_chemformula = load['cif_chemformula']
        self.__cifdata = load['cifdata']
        self.__cif_raw = load['cif_raw']
        self.__cellsg = load['cellsg']
        self.__atoms = load['atoms']
        self.__sites = load['sites']
        self._site2atom = load['site2atom']
        self._sgsetting = load['sgsetting']
        self.__cifdescr = load['cifdescr']
        self.__warnlist = tuple( warnlist )

        #Hidden support for additional processing in derived class:
        _ = getattr( self, '_process_raw_gemmi_objects', None )
        if _:
            _( gemmi, struct )

    @property
    def cifsrc( self ):
        """The CIFSource object on which everything is based."""
        return self.__cifsrc

    @property
    def actual_codid( self ):
        """None or the database ID of the Crystallography Open Database. This
        might be present even if .cifsrc.codid is None, since the CIF data might
        have been downloaded manually, but provide the ID in its metadata."""
        return self.__actual_codid

    @property
    def actual_mpid( self ):
        """None or the database ID of the Material Project. This might be
        present even if .cifsrc.mpid is None, since the CIF data might have
        been downloaded manually, but provide the ID in its metadata.
        """
        return self.__actual_mpid

    @property
    def raw_cif_textdata( self ):
        """The raw CIF input data, as a string."""
        return self.__cifdata

    @property
    def raw_cif_chemformula( self ):
        """The chemical formula encoded in the CIF data (usually in the
        _chemical_formula_sum field) as a string."""
        return self.__cif_chemformula

    @property
    def raw_cif_dict( self ):
        """The CIF data parsed into a dictionary."""
        return self.__cif_raw

    @property
    def atoms( self ):
        """A tuple of all atoms in the crystal structure, including positions in
        the unit cell. The format of each atom is essentially a dictionary
        similar to:
               {'cif_labels': ('Al',),
                'composition': ((1.0, 'Al'),),
                'equivalent_positions': ((0.0, 0.0, 0.0),
                                         (0.0, 0.5, 0.5),
                                         (0.5, 0.0, 0.5),
                                         (0.5, 0.5, 0.0)),
                'occupancy': 1.0,
                'uiso': None}

        Here 'uiso' (Uiso) provides information about isotropic displacements
        (or the Ueq of anisotropic ones), and occupancy is the site
        occupancy. The 'cif_labels' are the corresponding labels in the CIF
        data, and the 'composition' + 'equivalent_positions' should be
        self-explanatory. Each atom corresponds to one or more of the .sites.
        """
        return self.__atoms

    @property
    def sites( self ):
        """A tuple of the crystallographic sites of the structure, one per
        orbit of the space group (in the setting of the final unit cell). The
        format of each site is a dictionary similar to:
               {'cif_labels': ('Al',),
                'composition': ((1.0, 'Al'),),
                'occupancy': 1.0,
                'position': (0.0, 0.0, 0.0),
                'multiplicity': 4,
                'uiso': 0.0123,
                'aniso': {'u11': 0.0123, 'u22': 0.0123, 'u33': 0.0123,
                          'u12': 0.0, 'u13': 0.0, 'u23': 0.0}}

        Here 'aniso' is either None or the anisotropic displacement tensor in
        the CIF convention (relative to the normalised reciprocal axes) of the
        final unit cell, 'uiso' its Ueq or the Uiso value of the CIF data, and
        the 'cif_labels' those of the CIF sites corresponding to the site.
        """
        return self.__sites

    @property
    def cellsg( self ):
        """
        Get information about unit cell parameters and space group, as a
        dictionary with a format like:

           {'a': 4.0496, 'b': 4.0496,'c': 4.0496,
            'alpha': 90,'beta': 90,'gamma': 90,
            'spacegroup': {'hm': 'Fm-3m', 'number': 225} }

        """
        return self.__cellsg

    @property
    def extracted_description( self ):
        """
        Description of the material based on the meta-data found in the CIF
        data. Returned as a list of strings.
        """
        return self.__cifdescr

    @property
    def warnings( self ):
        """
        List of warnings emitted during processing of the CIF data, as a
        sequence of tuples of two strings: (warning_message, warning_category).
        """
        return self.__warnlist

    def create_ncmat( self, *args, **kwargs ):
        """
        Return NCMAT data corresponding to the currently held crystal
        parameters. All parameters are forwarded (same as calling
        create_ncmat_composer(*args,**kwargs).create_ncmat())
        """
        c = self.create_ncmat_composer(*args,**kwargs)
        return c.create_ncmat()

    #FIXME: Discuss whether assuming uiso_temperature=300K when it is not
    #given (so the CIF Uiso values still set the msds) would be better than
    #the fallback to a Debye temperature of 300K.
    def create_ncmat_composer( self, *,
                               uiso_temperature = None,
                               skip_dyninfo = False,
                               top_comments = None,
                               fallback_debye_temp = 300.0,
                               quiet = False,
                               remap = None,
                               no_formula_check = False ):
        """Setup and return an NCMATComposer object, based on the currently held
        crystal parameters.

        The most important thing to be aware of is that, if you wish to take
        advantage of any Uiso/biso values found in the CIF data, you must
        provide the uiso_temperature parameter value (kelvin), to indicate the
        temperature at which these are valid (most likely after inspecting the
        .raw_cif_textdata manually for hints about the correct value).

        If quiet=True, no informative messages will be emitted.

        The fallback_debye_temp parameter can be used to change the default
        Debye temperature assigned to atoms, for which Uiso/biso information was
        not available.

        The other parameters are more esoteric:

        If skip_dyninfo=True, no ncmat.set_dyninfo_... calls will be performed,
        presumably because the caller will perform such calls subsequently
        (i.e. with a PhononDOSAnalyser).

        The remap parameter can be used to remap elements and isotopes found in
        the CIF data to any desired composition (cf. the
        NCMATComposer.set_composition for the allowed composition syntax. For
        instance, if one wishes to take a CIF file containing hydrogen atoms and
        deuterate it partly, one could used remap = { 'H' : '0.95 D 0.05 H', }.

        Finally, unless no_formula_check=True, an exception is raised if the
        chemical formula of the resulting material is incompatible with any
        chemical formula indicated in the input CIF meta data. This check is
        normally nice to keep, since it can catch many cases where the wrong
        setting of a spacegroup was used.
        """
        return _impl_create_ncmat_composer( self,
                                            uiso_temperature = uiso_temperature,
                                            skip_dyninfo = skip_dyninfo,
                                            top_comments = top_comments,
                                            fallback_debye_temp = fallback_debye_temp,
                                            #TODO? #merge_equiv = merge_equiv,
                                            quiet = quiet,
                                            remap = remap,
                                            no_formula_check = no_formula_check )


def produce_validation_plots( files, verbose_lbls = True, pdf_target = None,
                              **plot_kwargs ):
    """Function which can produce several of the validation plots resulting
    from usage of the produce_validation_plot(..) function, and possibly even
    embed them in a PDF file. Any plot_kwargs will be passed along to the
    produce_validation_plot(..) function.
    """
    from .plot import _import_matplotlib_pdfpages, _import_matplotlib_plt

    if pdf_target:
        pdfpages = _import_matplotlib_pdfpages()
        the_real_pdfpages = ( pdfpages.the_real_inspected_object
                              if hasattr(pdfpages,'the_real_inspected_object')
                              else None )
        already_pdfpages = isinstance(pdf_target,the_real_pdfpages) if the_real_pdfpages else False
        pdf = pdf_target if already_pdfpages else pdfpages(pdf_target)
        assert pdf
    else:
        pdf, already_pdfpages = None, False

    if pdf:
        plt = _import_matplotlib_plt()
    for f in files:
        produce_validation_plot( f, verbose_lbls = verbose_lbls, do_show = not pdf,
                                 line_width_scale = 0.5, **plot_kwargs )
        if pdf:
            pdf.savefig()
            plt.close()

    if pdf and not already_pdfpages:
        pdf.close()
        _nc_common.print(f"created {pdf_target}")

def produce_validation_plot( data_or_file, verbose_lbls = True, line_width_scale = 1,
                             quiet = False, xlabel = None, legend_args = None, do_newfig = True,
                             do_show = True, do_legend = True, do_grid=True, do_tight_layout=True ):

    """Function which can produce validation plots for NCMAT data defining
    crystalline materials like the ones on
    https://github.com/mctools/ncrystal/wiki/Data-library. The plot will compare
    the powder Bragg diffraction cross section curve of the NCMAT file, with
    those obtained by replacing the crystal structure with a crystal structure
    loaded directly from online database entries that might be mentioned in the
    comments of the NCMAT data. This makes it easy to verify that the crystal
    structure is still compatible (or not) with those crystal structures
    mentioned in the comments.

    Currently two online databases are supported, with recognised URLs in a
    format like either of the following:

       https://www.crystallography.net/cod/9008460.html
       https://www.materialsproject.org/materials/mp-87

    """

    def lookupAndProduce( cifsrc, *,dynamics, remap):
        cifloader = CIFLoader( cifsrc, quiet = quiet )
        ncmat = cifloader.create_ncmat_composer( uiso_temperature = None,#dynamics provided separately here
                                                 remap = remap,
                                                 skip_dyninfo = True )
        ncmat.transfer_dyninfo_objects( dynamics )
        return ncmat.create_ncmat()

    from .plot import _import_matplotlib_plt
    plt = _import_matplotlib_plt()
    if do_newfig:
        plt.subplots()

    def multcreate( data ):
        return _nc_core.directMultiCreate(data,cfg_params='comp=bragg')
    _file = None
    fn = None
    if hasattr( data_or_file, '__fspath__' ):
        pass
    elif isinstance( data_or_file, _nc_core.TextData ):
        mc = multcreate( data_or_file )
        contentiterable = data_or_file
        fn = data_or_file.dataSourceName
    elif isinstance(data_or_file, (bytes, str)):
        _ = data_or_file.decode() if isinstance( data_or_file, bytes ) else data_or_file
        if '\n' in _ or _.startswith('NCMAT'):
            contentiterable = _.splitlines()
            mc = multcreate(_)
            fn = 'Anonymous NCMAT data'
    if not fn:
        fn = data_or_file
        contentiterable = _nc_core.createTextData(data_or_file)
        mc = multcreate(contentiterable)

    if hasattr(fn,'__fspath__'):
        import pathlib
        fn = pathlib.Path(fn).name
    elif '/' in fn:
        fn = fn.split('/')[-1]
    if not quiet:
        _nc_common.print(f'Attempting to validate {fn}')

    if not mc.info.isCrystalline():
        if not quiet:
            _nc_common.print(f"Ignoring non-crystalline material {fn}")
        return

    def _extractID(s,pattern):
        if pattern not in s:
            return None
        ll=s.split(pattern)[1:]
        while ll:
            e,ll = ll[0],ll[1:]
            d=''
            while e and e[0].isdigit():
                d+=e[0]
                e=e[1:]
            if d and int(d)>0:
                yield int(d)

    import re as _re
    _re_atomdbspecs = _re.compile(r"\[ *with ([a-zA-Z ]+)->([ a-zA-Z0-9+-\.]+) *\]")
    def _extractAtomDBSpec(s):
        #Look for remapping specs like "[with H->D]" and return in @ATOMDB format
        #(i.e. "H is D").
        if '->' not in s:
            return
        m = _re_atomdbspecs.search(s)
        return ' '.join('{} is {}'.format(*m.groups()).split()) if m else None

    ids = []
    for ll in contentiterable:
        atomdb = _extractAtomDBSpec(ll)
        for e in _extractID(ll,'materialsproject.org/materials/mp-'):
            ids.append( {'dbtype': 'mp','entryid': e,'atomdb': atomdb})
        for e in _extractID(ll,'crystallography.net/cod/'):
            ids.append( {'dbtype': 'cod','entryid': e,'atomdb': atomdb})

    #order-preserving remove duplicates:
    _seen = set()
    newids=[]
    for d in ids:
        key = tuple(sorted((k,v) for k,v in d.items()))
        if key in _seen:
            continue
        newids.append(d)
        _seen.add(key)
    ids = newids

    def _atomdb_to_remap( atomdb):
        _atomdb = [ ' '.join(e.strip().split())
                    for e in (atomdb or '').replace(':',' ').split('@') ]
        _atomdb = [ e for e in _atomdb if e ]
        ll = []
        for c in _atomdb:
            p=c.replace(':',' ').split()
            if not len(p)>=3 or p[1]!='is':
                raise _nc_core.NCBadInput(f'invalid atomdb remap syntax in "{c}"')
            ll.append( (p[0],' '.join(p[2:]) ) )
        return ll

    dynamics = mc.info
    cmps = [ (fn, mc ) ]
    for _ in ids:
        dbname = _['dbtype']
        entryid = _['entryid']
        atomdb = _['atomdb']
        lpargs={'remap': _atomdb_to_remap(atomdb),'dynamics': dynamics}
        if dbname=='mp':
            lpargs['cifsrc']=f'mpid::{int(entryid)}'
        else:
            assert dbname=='cod'
            lpargs['cifsrc']=f'codid::{int(entryid)}'
        out = lookupAndProduce( **lpargs )
        lbl=f'{dbname}-{entryid}'
        if atomdb:
            lbl = f'{lbl} with replacement "{atomdb}"'
        cmps.append( ( lbl, multcreate(out) ) )
    from .constants import ekin2wl
    wlmax = max(1.05*ekin2wl(mc.scatter.domain()[0]) for _,mc in cmps)
    from ._numpy import _np_linspace
    wls = _np_linspace(0.001,wlmax,10000)
    if len(cmps)<=1:
        _nc_common.warn(f"No online DB IDs found in {fn}")
        return

    sgno_all = set()
    natoms_all = set()

    col_ordered = [_nc_common._palette_Few.get(k,k) for k in
                   ('black',
                    'blue',
                    'orange',
                    'green',
                    'red',
                    'brown',
                    'purple',
                    'yellow',
                    'pink',
                    'gray')]

    def fix_lbl_for_plt(lbl):
        return lbl.replace('_',"$"+'\\'+"mathrm{\\_}$")

    for i,(lbl,mc) in enumerate(cmps):
        lbl = str(lbl)
        if '/' in lbl:
            lbl=lbl.split('/')[-1]
        elif lbl.startswith('mp-'):
            lbl = 'Materials Project entry '+lbl[3:]
        elif lbl.startswith('cod-'):
            lbl = 'Crystallography Open Database entry '+lbl[4:]
        lw=4 if i>0 else 2
        lw *= line_width_scale

        si = mc.info.structure_info if mc.info.hasStructureInfo() else None
        if not quiet:
            _nc_common.print(f'Structure[{lbl}]:',si)
        if si and verbose_lbls:
            sgno=si.get('spacegroup',None)
            lbl += ( f' (SG-{sgno or "unspecified"},'
                     f' {si["n_atoms_average"]:g} atoms/cell,'
                     f' Vcell={si["volume"]:g})' )
            sgno_all.add(sgno)
            natoms_all.add(si['n_atoms_average'])

        plt.plot(wls,mc.scatter.xsect(wl=wls),
                 label=fix_lbl_for_plt(lbl),
                 linewidth=lw,
                 alpha=0.5,#was: 0.5 if i==0 else 0.5,
                 color = col_ordered[i] if i<len(col_ordered) else None,
                 dashes = [] if i==0 else [4 if len(cmps)>2 else 2,2]+[2,2]*i)
    if len(natoms_all)>1:
        _nc_common.warn('WARNING validation detected different number of atoms/cell in tested curves')
    if len(sgno_all)>1:
        _nc_common.warn('WARNING validation detected different space group numbers in tested curves')

    plt.xlim(0.0)
    plt.ylim(0.0)
    plt.title('NB: This compares the crystal structure (space group, lattice, atom positions). Phonons/dynamics always taken from .ncmat file',fontsize=6)
    plt.xlabel(xlabel or 'Neutron wavelength ({})'.format(b'\xc3\x85'.decode()))
    plt.ylabel('Coherent elastic cross section (barn)')
    if do_legend:
        plt.legend(loc='best',handlelength=5,**(legend_args or {}))
    if do_grid:
        plt.grid()
    if do_tight_layout:
        plt.tight_layout()
    if do_show:
        plt.show()


####################### NEW GEMMI STUFF ############################################
