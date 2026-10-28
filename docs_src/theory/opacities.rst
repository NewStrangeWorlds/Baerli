Opacities
=========

Bärli evaluates the gas opacity by **opacity sampling**: the model's spectral grid is a subset of
the native high-resolution wavenumber grid on which the molecular and atomic cross sections were
pre-computed (HELIOS-k output), and at every selected point the cross section is taken as it is,
without averaging, binning or :math:`k`-distributions. The radiative transfer is then
monochromatic at every grid point. Multi-species mixing is exact, since the absorption
coefficients of all species are simply summed at each wavenumber, in contrast to correlated-*k*
methods that need random-overlap or pre-mixed tables. The accuracy depends on how well the point
density resolves the line structure (:doc:`../user_guide/spectral_grid`). The user-facing
species list is described in :doc:`../user_guide/opacities`.

Absorption coefficient
----------------------

At every level the absorption and scattering coefficients (cm\ :sup:`-1`) are assembled from the
configured opacity species:

.. math::

   \kappa^{\mathrm{abs}}_{\tilde\nu} = \sum_{s} \sigma_{s,\tilde\nu}(p, T)\, n_s
   + \sum_{\mathrm{pairs}} \sigma_{ab,\tilde\nu}(T)\, n_a n_b
   + \sum_{\mathrm{continua}} \alpha_{\tilde\nu}(T, \vec n), \qquad
   \kappa^{\mathrm{sca}}_{\tilde\nu} = \sum_s \sigma^{\mathrm{Ray}}_{s,\tilde\nu}\, n_s ,

with the line cross sections :math:`\sigma_s` in cm\ :sup:`2` per molecule interpolated in
pressure and temperature, the tabulated collision-induced absorption (CIA) coefficients
:math:`\sigma_{ab}` in cm\ :sup:`5` multiplied by the number densities of both partners, the
built-in continua (H\ :sup:`-`, CO\ :sub:`2` CIA, the MT_CKD water continuum) returning absorption
coefficients directly, and the Rayleigh cross sections.

Tabulated cross sections
~~~~~~~~~~~~~~~~~~~~~~~~

Each species' data folder holds one binary file per :math:`(p, T)` point and a ``filelist.dat``
that maps pressure and temperature to the files (:doc:`../user_guide/opacities`). A file is loaded
on first use, sampled onto the grid and unloaded again. HELIOS-k opacities in cm\ :sup:`2`
g\ :sup:`-1` are converted to cm\ :sup:`2` per molecule with the species mass from the file list;
CIA tables (mass 0, stored in log) are used as they are. The sampled values are floored at
:math:`10^{-200}` and stored as :math:`\log_{10}\sigma`.

For a given :math:`(p, T)` the four neighbouring tabulated points are located and
:math:`\log_{10}\sigma` is interpolated **bilinearly in** :math:`(\log_{10} p, T)`; outside the
tabulated range the nearest boundary point is used (clamping, no extrapolation). The same scheme
is applied to the band integrals of the band-closure correction.

Band closure
~~~~~~~~~~~~

The sampled points miss the Doppler-narrow line cores. For every line species and every 30
cm\ :sup:`-1` band the difference between the native band integral of the cross section and its
trapezoidal representation on the sampled points,

.. math::

   \Delta_b = n_s \left[\int_b \sigma\, \dd\tilde\nu - \sum_{k \in b} w_k\, \sigma_k\right] \ge 0,

together with the peak coefficient of the band, is handed to the temperature solver, which adds
the missing absorption to the local radiative-equilibrium residual with a frozen band-mean
radiation field and an escape-probability weight (:doc:`temperature_iteration`). Rayleigh
scattering, CIA and the continua are spectrally smooth and already exact on the sampled grid.

Collision-induced absorption
----------------------------

H\ :sub:`2`-H\ :sub:`2`, H\ :sub:`2`-He and H-He CIA are read from tabulated files of the same
layout as the line cross sections (one placeholder pressure, many temperatures) and multiplied
by the product of the two partner number densities.

CO\ :sub:`2`-CO\ :sub:`2`
~~~~~~~~~~~~~~~~~~~~~~~~~

The built-in ``CO2-CIA`` species combines the analytic roto-translational CIA of Gruszka &
Borysow (1997, 1998) for 0 to 500 cm\ :sup:`-1` with the dimer band of Baranov et al. (2004) for
1000 to 2000 cm\ :sup:`-1` (the port added by R. Wordsworth to the LMD generic GCM), as

.. math::

   \alpha(\tilde\nu) = \left[k_1(T, \tilde\nu) + k_2(T, \tilde\nu)\right] \rho_{\mathrm{amg}}^2,
   \qquad
   \rho_{\mathrm{amg}} = \frac{273.15}{T}\, \frac{p_{\mathrm{CO_2}}}{1.01325\ \mathrm{bar}},

with the temperature clamped to 100 to 800 K. There is no CO\ :sub:`2` CIA between 500 and
1000 cm\ :sup:`-1` or above 2000 cm\ :sup:`-1`.

Water-vapour continuum
~~~~~~~~~~~~~~~~~~~~~~

The built-in ``H2O-CIA`` species is the MT_CKD self and foreign continuum (version 4.2), with the
tabulated reference coefficients :math:`C_s`, :math:`C_f` and the self-continuum temperature
exponent :math:`n_s` at :math:`T_{\mathrm{ref}} = 296` K and :math:`p_{\mathrm{ref}} = 1.013`
bar:

.. math::

   \alpha(\tilde\nu) = \left[C_s(\tilde\nu) \left(\frac{T_{\mathrm{ref}}}{T}\right)^{n_s(\tilde\nu)}
   x_{\mathrm{H_2O}}\, \rho_r + C_f(\tilde\nu)\, (1 - x_{\mathrm{H_2O}})\, \rho_r\right]
   R(\tilde\nu, T)\, n_{\mathrm{H_2O}}, \qquad
   \rho_r = \frac{p}{p_{\mathrm{ref}}}\, \frac{T_{\mathrm{ref}}}{T},

with the radiation term :math:`R = \tilde\nu \tanh(h c \tilde\nu / 2 k T)` and the data limited
to 0 to 20000 cm\ :sup:`-1`.

H\ :sup:`-` bound-free and free-free
------------------------------------

The ``H-`` species uses the polynomial fits of John (1988). The bound-free photo-detachment cross
section for :math:`0.125 \le \lambda \le 1.6419` micron,

.. math::

   \sigma_{\mathrm{bf}}(\lambda) = 10^{-18}\, \lambda^3 \left(\frac{1}{\lambda} - \frac{1}{\lambda_0}\right)^{3/2}
   \sum_{n=1}^{6} C_n \left(\frac{1}{\lambda} - \frac{1}{\lambda_0}\right)^{(n-1)/2}
   \ \mathrm{cm^2}, \qquad \lambda_0 = 1.6419\ \micron,

and the free-free coefficient

.. math::

   \kappa_{\mathrm{ff}}(\lambda, T) = 10^{-29} \sum_{n} \left(\frac{5040}{T}\right)^{(n+1)/2}
   \left(\lambda^2 A_n + B_n + \frac{C_n}{\lambda} + \frac{D_n}{\lambda^2} + \frac{E_n}{\lambda^3} + \frac{F_n}{\lambda^4}\right)
   \ \mathrm{cm^4\,dyn^{-1}}

with two coefficient sets for :math:`\lambda > 0.3645` micron and :math:`0.1823 \le \lambda \le
0.3645` micron. The absorption coefficient is
:math:`\alpha = \sigma_{\mathrm{bf}}\, n_{\mathrm{H^-}} + \kappa_{\mathrm{ff}}\, P_e\, n_{\mathrm{H}}`
with the electron pressure :math:`P_e = n_e k T`; the number densities of H\ :sup:`-`, H and
electrons come from the equilibrium chemistry.

Rayleigh scattering
-------------------

Rayleigh cross sections are built in for H\ :sub:`2`, He, H, CO, CO\ :sub:`2`, CH\ :sub:`4`,
H\ :sub:`2`\ O, N\ :sub:`2` and O\ :sub:`2`, from the refractive index :math:`n` at the reference
number density :math:`n_{\mathrm{ref}}` and the King correction factor :math:`F_K`,

.. math::

   \sigma_{\mathrm{Ray}}(\tilde\nu) = \frac{24\pi^3 \tilde\nu^4}{n_{\mathrm{ref}}^2}
   \left(\frac{n^2 - 1}{n^2 + 2}\right)^2 F_K ,

with the standard Sellmeier and Cauchy dispersion formulae for each gas (the forms tabulated by
Sneep & Ubachs 2005 and Thalman et al. 2014), the IAPWS refractive index of water, and the
low-energy expansion of Lee & Kim (2004) for atomic hydrogen,
:math:`\sigma = 8.41 \times 10^{-25} x^4 + 3.37 \times 10^{-24} x^6 + 4.71 \times 10^{-22} x^{14}`
cm\ :sup:`2` with :math:`x = 0.0912\,\micron / \lambda`. The cross sections of the molecules are
set to zero beyond 2 micron, where they are negligible; H\ :sub:`2`, He and H have no cut-off. The
tables are temperature and pressure independent and computed once.

Clouds
------

The cloud modules provide per-layer optical depths, single-scattering albedos and asymmetry
parameters that the radiative transfer adds to the gas layer by layer; the cloud absorption also
enters the total absorption coefficient that weights the local-equilibrium residual of the
solvers (:doc:`clouds`).
