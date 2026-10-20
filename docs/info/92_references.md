# Bibliography

## NCrystal publications (citation requests)

- X.-X. Cai, T. Kittelmann, "NCrystal: A library for thermal neutron
  transport", Comput. Phys. Commun. 246 (2020) 106851,
  doi:10.1016/j.cpc.2019.07.015 (arXiv:1901.08890). Primary reference.
- T. Kittelmann, X.-X. Cai, "Elastic neutron scattering models for
  NCrystal", Comput. Phys. Commun. 267 (2021) 108082,
  doi:10.1016/j.cpc.2021.108082 (arXiv:2012.04294). Elastic physics.
- X.-X. Cai, T. Kittelmann, E. Klinkby, J. I. Marquez Damian,
  "Rejection-based sampling of inelastic neutron scattering", J. Comput.
  Phys. 380 (2019) 400, doi:10.1016/j.jcp.2018.11.043. Inelastic physics.

## Related NCrystal-ecosystem publications

- K. Ramic et al., "NJOY+NCrystal: An open-source tool for creating
  thermal neutron scattering libraries with mixed elastic support", NIM A
  1027 (2022) 166227.
- K. Ramic et al., "Advances in nuclear data and software development for
  the HighNESS project", Nucl. Sci. Eng. 198 (2024) 74.
- N. Rizzi et al., nanodiamond SANS models and benchmark, Nucl. Sci. Eng.
  (2023), doi:10.1080/00295639.2023.2196926.
- J. R. Granada, D. D. DiJulio, J. I. Marquez Damian, G. Muhrer, liquid
  helium scattering kernel, NIM A 1053 (2023) 168284.
- S. Xu et al., "Impact of extinction effects on neutron transmission in
  solid beryllium metal", J. Appl. Cryst. 58 (2025) 1957,
  doi:10.1107/S1600576725007939.
- T. Kittelmann, D. D. DiJulio, S. Xu, J. I. Marquez Damian, "Revisiting
  Becker-Coppens (1974): updated recipes for estimating extinction factors
  in spherical crystallites", Acta Cryst. A82 (2026) 163,
  doi:10.1107/S2053273326001245.
- J. I. Marquez Damian et al., "Open source pipelines for thermal
  scattering data generation: from NCrystal to ENDF-6", ND2025.
- Grammer et al., "An NCrystal extension for thermal neutron scattering in
  MCNP", NIM A (2026), doi:10.1016/j.nima.2025.171079.
- Jalil et al., "Validation of ENDF/B-VIII.0 and NCrystal thermal
  scattering laws for graphite in criticality safety analysis using
  OpenMC", Nuclear Analysis (2026), doi:10.1016/j.nucana.2026.100219.
- D. D. DiJulio, J. I. Marquez Damian, T. Kittelmann, G. Muhrer, "A liquids
  plugin for NCrystal", ICANS XXV (2026).
- G. Romanelli et al., J. Phys.: Condens. Matter 33 (2021) 285901,
  doi:10.1088/1361-648X/abfc13 (hydrogen functional groups).

## Theory

- A. Sjolander, "Multi-phonon processes in slow neutron scattering by
  crystals", Arkiv for Fysik 14 (1958) 315.
- G. Placzek, L. Van Hove, Nuovo Cimento 1 (1955) 233 (incoherent
  approximation).
- R. J. Glauber, Phys. Rev. 98 (1955) 1692 (Debye-model displacements).
- W. E. Lamb, Phys. Rev. 55 (1939) 190 (impulse approximation).
- G. M. Borgonovi, Gulf General Atomic report GA-9950 (1970) (SCT).
- T. M. Sutton, T. H. Trumbull, C. R. Lubitz, Proc. M&C 2009 (SCT formula
  catalogue).
- V. F. Sears, Neutron News 3 (1992) 26; Acta Cryst. A53 (1997) 35.
- G. L. Squires, "Introduction to the theory of thermal neutron
  scattering" (free gas, chapter 4.5).
- B. Hammouda, "The SANS toolbox", NIST (2012).
- W. A. Dollase, J. Appl. Cryst. 19 (1986) 267 (March-Dollase texture);
  T. M. Sabine, Aust. J. Phys. 38 (1985) 507 and Acta Cryst. A44 (1988) 368
  (extinction); P. J. Becker, P. Coppens, Acta Cryst. A30 (1974) 129.
- O. Zimmer, Phys. Rev. C 93 (2016) 035503 (paramagnetic scattering).

## Data and formats

- ENDF-6 Formats Manual (ENDF-102), sections 7.4 (thermal scattering);
  ENDF/B-VIII.0 thermal sublibrary; R. E. MacFarlane et al., NJOY2016
  (LA-UR-17-20093); LA-12639-MS (1994).
- U. Shmueli, S. R. Hall, R. W. Grosse-Kunstleve, International Tables for
  Crystallography Vol. B (2010), Appendix A1.4.2, Table A1.4.2.7; S. R.
  Hall, Acta Cryst. A37 (1981) 517.
- W. J. Huang et al., "The AME2016 atomic mass evaluation (I)", Chinese
  Physics C 41 (2017) 030002; the ATI Vienna scattering-length table
  (Landolt-Bornstein I/16A, ADNDT 49 (1991) 65).
- A. Togo, I. Tanaka, Scripta Materialia 108 (2015) 1 (phonon database);
  J. I. Marquez Damian et al., Ann. Nucl. Energy 65 (2014) 280 (CAB water
  kernels).
- M. Boin, J. Appl. Cryst. 45 (2012) 603 (nxslib); G. Schnabel et al.,
  arXiv:2312.08249 (endf-parserpy).
- spglib, gemmi, ASE, SgInfo (symmetry cross-checks); the Crystallography
  Open Database and the Materials Project.

## Numerics

- F. N. Fritsch, R. E. Carlson, SIAM J. Numer. Anal. 17 (1980) 238
  (PCHIP); W. J. Cody, W. Waite (1980) (range reduction); Blackman and
  Vigna, xoroshiro128+; NIST DLMF 5.11 (Stirling series).

## NCrystal documentation and resources

- Repository: https://github.com/mctools/ncrystal (`docs/ncmat_doc.md`,
  `docs/data_browsing.md`, `docs/devel_*.md`, `CHANGELOG`).
- Wiki: https://github.com/mctools/ncrystal/wiki (Home, Using-NCrystal,
  Included-physics, CfgRefDoc, NCMAT-format, Data-library, Plugins,
  CuratedPlugins, ncmat2endf, minimc, release announcements 2.0.0, 2.1.0,
  2.5.0, 2.7.0, 3.0.0, 3.1.0, 3.2.0).
- Notebooks: https://github.com/mctools/ncrystal-notebooks.
- Extra data and tools: https://github.com/mctools/ncrystal-extra; Geant4
  interface: https://github.com/mctools/ncrystal-geant4.
