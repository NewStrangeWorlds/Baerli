Atmospheric structure
=====================

This chapter describes the discretisation of the atmosphere, the hydrostatic structure and the
thermodynamic quantities the other components rely on.

State variables
---------------

The atmosphere is a column of ``nb_grid_points`` pressure **levels**. The following quantities
are defined on every level and exposed through :attr:`model.atmosphere`:

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Quantity
     - Unit
     - Set by
   * - :math:`p_i`, ``pressure``
     - bar
     - the grid construction
   * - :math:`T_i`, ``temperature``
     - K
     - the initial profile, then the temperature solver
   * - :math:`n_{i,s}`, ``number_densities``
     - cm\ :sup:`-3`
     - the chemistry modules; index 0 of the species axis is the total number density
       (``Total``)
   * - :math:`\bar\mu_i`, ``mean_molecular_weight``
     - g mol\ :sup:`-1`
     - the chemistry modules
   * - :math:`\rho_i`, ``mass_density``
     - g cm\ :sup:`-3`
     - the structure calculation
   * - :math:`z_i`, ``altitude``
     - cm
     - the structure calculation (:math:`z_0 = 0`)
   * - :math:`H_{p,i}`, ``scale_height``
     - cm
     - the structure calculation
   * - :math:`K_{zz,i}`, ``kzz``
     - cm\ :sup:`2` s\ :sup:`-1`
     - the eddy diffusion module
   * - ``convective``
     - 0/1
     - the convection scheme or the corrector (diagnostic flag)
   * - condensed phase
     - formula units cm\ :sup:`-3`
     - the equilibrium chemistry with condensation

Mixing ratios are not stored; they are the ratio :math:`n_{i,s}/n_{i,\mathrm{tot}}`
(:func:`baerli.mixing_ratios`). Heat capacities, adiabatic gradients and latent heats are
evaluated on demand from the composition and temperature (see below). The convective flux is
not a stored field either; it lives inside the mixing-length corrector
(:doc:`convection`).

The pressure grid
-----------------

The levels are equidistant in :math:`\log_{10} p` between the bottom pressure
``boundary_pressures[0]`` (level 0) and the top pressure ``boundary_pressures[1]``
(level :math:`N-1`):

.. math::

   \log_{10} p_i = \log_{10} p_{\mathrm{bot}} - i\,\frac{\log_{10} p_{\mathrm{bot}} - \log_{10} p_{\mathrm{top}}}{N-1},
   \qquad i = 0, \ldots, N-1 .

**Index 0 is the bottom** of the atmosphere, the highest pressure. Every array in Bärli, in the
Python interface and in the output files follows this convention. A **layer** :math:`j` is the
interval between levels :math:`j` and :math:`j+1`; layer quantities (optical depths, cloud
properties, the "links" of the convection schemes) have :math:`N-1` entries.

Hydrostatic structure
---------------------

The structure is recomputed from the current composition and temperature at every outer
iteration and in every trial evaluation of the Newton solvers that recomputes the opacity.

**Mass density** from the ideal gas law, with :math:`p` converted from bar to dyn cm\ :sup:`-2`
and :math:`R = 8.314472 \times 10^7` erg K\ :sup:`-1` mol\ :sup:`-1`:

.. math::

   \rho_i = \frac{\bar\mu_i \, (10^6\, p_i)}{R\, T_i} .

**Altitude** by trapezoidal integration of hydrostatic equilibrium
:math:`\dd z = -\dd p/(\rho g)` upward from :math:`z_0 = 0`:

.. math::

   z_i = z_{i-1} + \frac{1}{2}\left(\frac{1}{\rho_i\, g_i} + \frac{1}{\rho_{i-1}\, g_i}\right)
   10^6 \,(p_{i-1} - p_i) .

With constant gravity :math:`g_i = g` is the ``surface_gravity``. With ``variable_gravity``
(brown dwarfs and gas planets with a ``radius`` :math:`R_0`) the gravity of each step is
evaluated at the radius of the lower level,

.. math::

   g_i = g_0 \left(\frac{R_0}{R_0 + z_{i-1}}\right)^2 .

**Scale height**

.. math::

   H_{p,i} = \frac{R\, T_i}{\bar\mu_i\, g} .

The scale height, the mixing length of the convection scheme and the eddy diffusion coefficient
always use the bottom gravity, also with ``variable_gravity``; only the altitude integration
sees the altitude dependence. Terrestrial planets are always plane-parallel with constant
gravity.

Thermodynamic data
------------------

Heat capacities
~~~~~~~~~~~~~~~

Molar heat capacities come from the NASA 7-coefficient polynomials of McBride, Zehe & Gordon
(2002),

.. math::

   \frac{c_p}{R} = a_1 + a_2 T + a_3 T^2 + a_4 T^3 + a_5 T^4,

with two temperature ranges (200 to 1000 K and 1000 to 6000 K) and the temperature clamped to
200 to 6000 K. Coefficients are tabulated for H\ :sub:`2`, He, H, H\ :sub:`2`\ O, CO,
CO\ :sub:`2`, CH\ :sub:`4`, NH\ :sub:`3`, H\ :sub:`2`\ S, N\ :sub:`2`, HCN, C\ :sub:`2`\ H\
:sub:`2` and OH. Monatomic species (H, He, C, O, Fe, Ca, Ti, Na, K, V, Mn, Si, Cr and their
ions, electrons) use :math:`c_p/R = 5/2`; every other species falls back to the diatomic value
:math:`c_p/R = 7/2`. The fallback is a poor approximation for Ar and for the triatomics
SO\ :sub:`2`, SiO\ :sub:`2`, CS\ :sub:`2`, N\ :sub:`2`\ O and O\ :sub:`3`, which does not matter
for H\ :sub:`2`/He or N\ :sub:`2`/O\ :sub:`2`/CO\ :sub:`2` dominated atmospheres at the per-cent
level.

The mean specific heat capacity of the mixture is mole-fraction weighted,

.. math::

   \bar c_p = \frac{R \sum_s x_s\, (c_p/R)_s}{\sum_s x_s\, M_s}
   \qquad [\mathrm{erg\,g^{-1}\,K^{-1}}], \qquad x_s = \frac{n_s}{n_{\mathrm{tot}}} .

Dry adiabatic gradient
~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \nabla_{\mathrm{ad}} \equiv \left(\frac{\dd \ln T}{\dd \ln p}\right)_{\mathrm{ad}}
   = \frac{R}{\bar c_{p,\mathrm{molar}}} = \frac{1}{\sum_s x_s\, (c_p/R)_s},

the ideal-gas result :math:`(\gamma - 1)/\gamma` with the mixture's molar heat capacity.

Water: saturation vapour pressure and latent heat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The saturation vapour pressure of water follows Murphy & Koop (2005): their Eq. 7 over ice for
:math:`T \le 273.15` K,

.. math::

   \ln e_s[\mathrm{Pa}] = 9.550426 - \frac{5723.265}{T} + 3.53068 \ln T - 0.00728332\, T,

and their Eq. 10 over liquid water above,

.. math::

   \ln e_s[\mathrm{Pa}] = 54.842763 - \frac{6763.22}{T} - 4.210 \ln T + 0.000367\, T
   + \tanh\!\left[0.0415 (T - 218.8)\right]
   \left(53.878 - \frac{1331.22}{T} - 9.44523 \ln T + 0.014025\, T\right).

The latent heat follows from the Clausius-Clapeyron relation applied to the same fit,

.. math::

   L(T) = \frac{\dd \ln e_s}{\dd T}\, \frac{R}{M_{\mathrm{H_2O}}}\, T^2
   \qquad [\mathrm{erg\,g^{-1}}],

so that it is the sublimation enthalpy below the triple point and the vaporisation enthalpy
above.

Moist adiabatic gradient
~~~~~~~~~~~~~~~~~~~~~~~~

The moist (saturated pseudo-adiabatic) gradient is the standard closed form

.. math::

   \nabla_m = \nabla_{\mathrm{ad}}\,
   \frac{1 + \dfrac{L\, w_s}{R_d\, T}}{1 + \dfrac{L^2\, w_s}{c_{p,d}\, R_v\, T^2}},
   \qquad
   w_s = \frac{M_{\mathrm{H_2O}}}{\bar\mu}\, \frac{e_s}{p - e_s},

with :math:`R_d = R/\bar\mu` and :math:`c_{p,d}` the gas constant and heat capacity of the
mixture, :math:`R_v = R/M_{\mathrm{H_2O}}` and :math:`w_s` the saturation mass mixing ratio.
Where :math:`e_s \ge p` the dry gradient is returned.

Two properties of this formulation matter for its use in the convection schemes:

* It assumes the parcel is **saturated with water at the local temperature and pressure**,
  independent of the actual water abundance of the composition. In a subsaturated, cold
  stratosphere the saturated pseudo-adiabat becomes very flat, and a moist scheme without a
  pressure ceiling would read such a stratosphere as moist-unstable; this is why the
  ``min_pressure`` of ``mlt_moist`` is load-bearing (:doc:`convection`).
* Only water condensation is treated; there are no NH\ :sub:`3` or CH\ :sub:`4` moist adiabats.

The energy balance of the three objects
---------------------------------------

The structure described above is the same for all objects; they differ in the boundary
conditions and in the flux every level has to carry in equilibrium:

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * -
     - BrownDwarf
     - GasPlanet
     - TerrestrialPlanet
   * - Target flux :math:`F_\star`
     - :math:`\sigma \Teff^4`
     - :math:`\sigma \Tint^4`
     - 0
   * - Lower boundary
     - semi-infinite (diffusion)
     - semi-infinite (diffusion)
     - surface with albedo and temperature
   * - Upper boundary
     - no incident radiation
     - stellar beam :math:`\mu_\star S`
     - stellar beam :math:`\mu_\star S`
   * - Flux scale of the residual
     - :math:`F_\star`
     - :math:`\mu_\star S + F_\star`
     - :math:`\max(1, \mu_\star S)`
   * - Mixing length near the bottom
     - :math:`\alpha H_p`
     - :math:`\alpha H_p`
     - Blackadar wall law
   * - Default ``min_pressure``
     - :math:`10^{-4}` bar
     - :math:`10^{-4}` bar
     - :math:`10^{-3}` bar

In all cases the converged state satisfies, at every level,

.. math::

   F_{\mathrm{rad}}(p) + F_{\mathrm{conv}}(p) = F_\star,

where :math:`F_{\mathrm{rad}}` is the net radiative flux (thermal plus stellar, upward
positive) and :math:`F_{\mathrm{conv}}` the convective flux of the mixing-length scheme (zero
with the adjustment schemes, which instead enforce the adiabatic gradient directly). For a
terrestrial planet the bottom level carries the surface energy balance
:math:`F_{\mathrm{net}}(p_0) = 0`. See :doc:`temperature_iteration` for how the solvers enforce
this.
