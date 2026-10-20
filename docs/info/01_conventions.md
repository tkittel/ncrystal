# Conventions: variables, units, normalisations, constants

## Units

- Energy in eV, wavelength in angstrom (Aa), temperature in kelvin (cfg
  strings also accept C and F), lengths in angstrom in data and in metres
  in MiniMC geometry, cross sections in barn per atom, scattering lengths
  in fm in NCMAT data (internally $\sqrt{\mathrm{barn}} = 10$ fm), densities
  in g/cm3, kg/m3 or atoms/Aa3, angles in radians internally with `deg`
  and `arcmin` accepted in cfg strings, masses in Daltons (amu).
- Conversions (`NCDefs.hh`): $E = h^2/(2m_n\lambda^2)$ (`wl2ekin`,
  `ekin2wl`, `k2ekin`, `ekin2k`, `ksq2ekin`), $k = 2\pi/\lambda$; the
  constants are CODATA values (`const_neutron_mass_amu` = 1.00866491588,
  `constant_hbar`, `constant_boltzmann`, ...). Thermal reference: 2200 m/s = 0.0253 eV
  = 1.798 Aa; $kT$ at 293.15 K = 25.26 meV.

## Kinematics

$E$ and $E'$ are the incident and outgoing kinetic energies, $\mu =
\cos\theta$ the cosine of the scattering angle, $\vec Q = \vec k_f - \vec
k_i$ the wave-vector transfer with $Q^2 = k^2 + k'^2 - 2kk'\mu$, and
$\hbar\omega = E - E'$ the energy transferred to the sample (positive for
energy loss of the neutron). In the kernel chapters the dimensionless
variables are

$$\alpha = \frac{E + E' - 2\mu\sqrt{EE'}}{kT} = \frac{\hbar^2 Q^2}{2m_n kT}, \qquad \beta = \frac{E'-E}{kT},$$

so $\beta < 0$ is downscatter, and $\alpha$ uses the neutron mass (ENDF's
$\alpha_E = \alpha/A$ with $A = M/m_n$). The accessible region at energy
$E$ is $\beta \ge -E/kT$, $\alpha_-(\beta) \le \alpha \le \alpha_+(\beta)$,
$\alpha_\pm = (2E + kT\beta \pm 2\sqrt{E(E+kT\beta)})/kT$, equivalently
$(\alpha-\beta)^2 \le 4\alpha E/kT$. Elastic scattering has $\beta = 0$,
$\alpha = 2(E/kT)(1-\mu)$.

## Scattering function and cross sections

$$\frac{d^2\sigma}{dE'd\Omega} = \frac{\sigma_b}{4\pi kT}\sqrt{\frac{E'}{E}}S(\alpha,\beta), \qquad \sigma(E) = \frac{\sigma_b kT}{4E}\iint_{D(E)}S\,d\alpha\,d\beta,$$

with the *asymmetric* $S$ obeying $S(\alpha,\beta) = e^{-\beta}S(\alpha,-\beta)$
and $S_{sym} = e^{\beta/2}S$. Cross sections are per atom of the material
(composition-averaged); a phase's process is scaled by the atom-count
weights in multiphase materials. Bound cross sections $\sigma_b =
\sigma_{coh} + \sigma_{inc} = 4\pi(\bar b^2 + \text{variance term}) +
\sum f_i\sigma_{inc,i}$; free cross sections $\sigma_{free} = \sigma_b
(A/(1+A))^2$. Macroscopic cross section $\Sigma = n_d\sigma$ with the
number density $n_d$ in atoms/Aa3 (so $\Sigma$ in 1/cm needs the factor
$10^{-24}$ cm2/barn times $10^{24}$ Aa3/cm3).

## Crystallography

Conventional unit cell of the space-group setting, fractional coordinates
in $[0,1)$, Miller indices $hkl$ with $d_{hkl}$ from the reciprocal metric,
$|\vec\tau_{hkl}| = 2\pi/d_{hkl}$, Bragg condition $\lambda \le 2d$ and
$\lambda = 2d\sin\theta_B$ (scattering angle $2\theta_B$), structure factor
$F = \sum_i b_i e^{-W_i}e^{i\vec\tau\cdot\vec r_i}$ with
$W = \frac12 Q^2\langle u^2\rangle$ and $\langle u^2\rangle$ the
one-dimensional isotropic mean-squared displacement ($U_{iso}$;
$B_{iso} = 8\pi^2 U_{iso}$). Multiplicities count both $hkl$ and $-h-k-l$.
Mosaicity `mos` is a FWHM; `dir1`/`dir2` map crystal to laboratory
directions.

## Dynamics

Phonon energies $E = \hbar\omega$ in eV; the VDOS $\rho(E)$ unit-normalised
over $[0, E_{max}]$ with $\rho \propto E^2$ below its first point; Debye
temperature $T_D$ with cut-off $kT_D$; effective temperature $T_{eff}$ from
$kT_{eff} = \int\rho(E)(E/2)\coth(E/2kT)dE$; Sjolander's $G_n$ normalised
to unit area; $x = 2W = \alpha\cdot\texttt{alpha2x}$ with
$\texttt{alpha2x} = 2m_n kT\langle u^2\rangle/\hbar^2$.
