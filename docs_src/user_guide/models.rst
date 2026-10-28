Model classes
=============

The three object classes share the pressure grid, the opacity sources and the pluggable
components; they differ in their boundary conditions and in the object-specific physical
parameters. All three are built the same way:

.. code-block:: python

   model = baerli.<Class>(grid, <object parameters>, <shared configuration>)
   model.initialize(<profile spec>)          # or model.initialize_from_file(...)
   converged = model.compute()

Shared configuration
--------------------

Every class accepts these keyword arguments after the grid:

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - Argument
     - Default
     - Meaning
   * - ``nb_grid_points``
     - ``100``
     - Number of pressure **levels**. Layers are the intervals between levels, so there are
       ``nb_grid_points - 1`` layers.
   * - ``boundary_pressures``
     - ``[1e2, 1e-6]``
     - ``[bottom, top]`` in bar. The levels are spaced uniformly in :math:`\ln p`.
   * - ``opacity_path``
     - the grid's
     - Root directory of the cross-section data; defaults to the path the spectral grid was
       built from.
   * - ``opacity_species``
     - ``[]``
     - List of ``(species symbol, data folder)`` pairs (:doc:`opacities`).
   * - ``chemistry``
     - **required**
     - List of chemistry module specs, applied in order (:doc:`chemistry`).
   * - ``clouds``
     - none
     - List of cloud module specs; their optical properties are summed (:doc:`clouds`).
   * - ``cloud_ramp``
     - ``"default"``
     - How the cloud opacity is switched on during the iteration (:doc:`clouds`).
   * - ``radiative_transfer``
     - ``("disort", {nb_streams: 4})``
     - Radiative transfer solver (:doc:`radiative_transfer`).
   * - ``convection``
     - ``"mlt_dry"``
     - Convection scheme (:doc:`convection`).
   * - ``kzz``
     - ``"mlt"``
     - Eddy diffusion profile (:doc:`convection`).
   * - ``solver``
     - ``"ratio_ul"``
     - Temperature correction scheme and its convergence settings (:doc:`solvers`).

The same list is available at run time as the string ``baerli.model_config_doc`` and in the
docstrings of the classes.

BrownDwarf
----------

A self-luminous, semi-infinite atmosphere: no irradiation and no surface. The internal flux
:math:`\sigma \Teff^4` escapes through the diffusion lower boundary and anchors the deep
temperature.

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Argument
     - Default
     - Meaning
   * - ``effective_temperature``
     - required
     - K. Sets the internal flux :math:`\Fint = \sigma \Teff^4` that every level must carry.
   * - ``surface_gravity``
     - required
     - cm s\ :sup:`-2` at the bottom of the domain.
   * - ``radius``
     - ``0.0``
     - cm, the radius at the bottom of the domain. Only used with ``variable_gravity``.
   * - ``variable_gravity``
     - ``False``
     - Let the gravity decrease with altitude as :math:`g(r) = g_0 (R_0/r)^2`; otherwise
       plane-parallel with constant gravity.

The convection scheme's ``min_pressure`` defaults to :math:`10^{-4}` bar for this class.

GasPlanet
---------

An irradiated, semi-infinite atmosphere without a surface. The top is illuminated by a stellar
beam, the deep atmosphere is anchored to the internal heat flux :math:`\sigma \Tint^4`. In
equilibrium the net **total** flux (thermal plus stellar) equals :math:`\Fint` at every level.

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Argument
     - Default
     - Meaning
   * - ``internal_temperature``
     - required
     - K. Sets :math:`\Fint = \sigma \Tint^4`.
   * - ``surface_gravity``
     - required
     - cm s\ :sup:`-2`.
   * - ``instellation``
     - required
     - Incident stellar flux at the planet in erg cm\ :sup:`-2` s\ :sup:`-1`, integrated over
       wavelength, measured perpendicular to the beam. The stellar spectrum is scaled to this
       value.
   * - ``zenith_angle``
     - required
     - **Cosine** of the stellar zenith angle :math:`\mu_\star`. The flux entering the
       atmosphere is :math:`\mu_\star S`. ``0.5`` approximates a global average.
   * - ``stellar_spectrum``
     - required
     - ``("tabulated", {file})`` or ``("blackbody", {temperature})``
       (:doc:`radiative_transfer`).
   * - ``radius``, ``variable_gravity``
     - ``0.0``, ``False``
     - As for the brown dwarf.

The residuals and the convergence metric of an irradiated planet are normalised by the flux the
atmosphere actually transports, :math:`\mu_\star S + \Fint`, not by :math:`\Fint` alone, which
can be tiny (:doc:`solvers`).

TerrestrialPlanet
-----------------

An irradiated atmosphere on top of a solid surface, plane-parallel with constant gravity. There
is no internal heat flux; the surface energy balance :math:`\Fnet(\mathrm{surface}) = 0`
closes the problem and the net flux vanishes at every level in equilibrium.

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Argument
     - Default
     - Meaning
   * - ``surface_gravity``
     - required
     - cm s\ :sup:`-2`.
   * - ``instellation``
     - required
     - Incident stellar flux, erg cm\ :sup:`-2` s\ :sup:`-1` (see the gas planet).
   * - ``zenith_angle``
     - required
     - Cosine of the stellar zenith angle.
   * - ``stellar_spectrum``
     - required
     - ``("tabulated", {file})`` or ``("blackbody", {temperature})``.
   * - ``surface``
     - required
     - ``"blackbody"``, ``("simple", {albedo, wavelength_switch})`` or
       ``("variable_albedo", {file})`` (:doc:`radiative_transfer`).

The surface temperature is a result of the model and available as
``model.surface_temperature``. With the Newton-type solvers it is identical to the temperature
of the bottom level (the surface is tied to the lowest level, which keeps the surface energy
balance as the flux equation of that level); with the relaxation solvers it is updated by the
surface's own energy balance after each temperature step. The convection scheme's
``min_pressure`` defaults to :math:`10^{-3}` bar for this class.

Setting the initial state
-------------------------

.. py:method:: initialize(profile, chemistry=None)
   :noindex:

   Initialise temperature and composition from an analytic profile spec. ``profile`` is one of
   ``adiabat``, ``milne``, ``const`` or ``guillot`` with their parameters (:ref:`initial-profiles`).
   ``chemistry`` optionally lists the chemistry modules used for the **initial** composition
   instead of the model's own (e.g. prescribed mixing ratios for a start far from equilibrium).
   Returns the model.

.. py:method:: initialize_from_file(filename)
   :noindex:

   Restart from an output file written by :func:`baerli.save_model`: temperature, number
   densities and mean molecular weight are taken from the file, which must have the same
   number of levels. Returns the model.

.. py:method:: initialize_from_arrays(temperature, number_densities, mean_molecular_weight)
   :noindex:

   The compiled core's restart entry point: plain lists on this model's grid.

After either call the hydrostatic structure, the eddy diffusion profile and the cloud modules
are evaluated once, so the model is ready for :meth:`compute` or for diagnostic evaluations.

Running the iteration
---------------------

.. py:method:: compute(tolerance=0.0, max_iterations=0, warm_start=False)
   :noindex:

   Iterate to radiative-convective equilibrium; returns ``True`` on convergence.

   ``tolerance`` and ``max_iterations`` override the solver spec's ``convergence_threshold``
   and ``max_iterations`` for this call only (values ``<= 0`` keep the spec's values).

   ``warm_start=True`` declares the current profile a converged state with a slightly changed
   problem (e.g. the previous pass of a coupling loop with a new composition): the corrector
   then skips its easy-start homotopy (the mixing-length ramp), and a cloud ramp applies the
   full cloud immediately.

   The call can end early without convergence when the iteration **stalls**: when the corrector
   proposes no temperature change at all for three iterations while the residual is above the
   threshold, or when the residual has not improved by 1 % over the last 25 iterations. Both
   cases are reported and return ``False``; an outer loop can simply carry on.

Several :meth:`compute` calls can be chained. The solver's internal state is reset at every
call, but the atmosphere is not: a second call continues from the converged profile, which is
what the kinetics coupling and parameter continuation rely on.

Handing in external data
------------------------

.. py:method:: set_composition(species, mixing_ratios)
   :noindex:

   Hand a composition to the ``external`` chemistry module: a list of species symbols and a
   table ``mixing_ratios[level][k]`` on this model's grid (index 0 = bottom). Unknown symbols
   are ignored. The supplied species are rescaled per level to the budget the preceding
   chemistry modules gave them, everything else is left untouched (:doc:`chemistry`).

.. py:method:: set_chemistry_enabled(type, enabled)
   :noindex:

   Enable or disable every chemistry module of the given type (``"quench"``, ``"external"``,
   ...); returns the number of modules affected. A disabled module passes the composition
   through unchanged.

.. py:method:: set_cloud_optical_properties(pressure, wavelength, optical_depth, single_scattering, asymmetry)
   :noindex:

   Hand a table of cloud optical properties to the ``fixed`` cloud module(s) configured without
   a file (:doc:`clouds`). ``pressure`` are levels in bar in any order, ``wavelength`` in
   micron, the three tables have shape ``[len(pressure)][len(wavelength)]``.

Results and diagnostics
-----------------------

.. py:attribute:: atmosphere
   :noindex:

   Read-only view of the atmospheric state (:class:`baerli.Atmosphere`): ``pressure`` (bar),
   ``temperature`` (K), ``altitude`` and ``scale_height`` (cm), ``mass_density``
   (g cm\ :sup:`-3`), ``mean_molecular_weight`` (amu), ``number_densities[level][species]``
   (cm\ :sup:`-3`, species order of :func:`baerli.species_symbols`), ``convective`` (1 in
   convective levels), ``kzz`` (cm\ :sup:`2` s\ :sup:`-1`), and the condensed phase
   ``condensate_symbols``, ``condensate_weights``, ``condensate_densities[level][k]``.

.. py:attribute:: radiation_field
   :noindex:

   The radiation field of the last radiative transfer solve
   (:class:`baerli.RadiativeTransferOutput`): the wavelength-integrated ``flux_total`` (net),
   ``flux_up_total``, ``flux_down_total`` and ``mean_intensity_total`` per level, the spectral
   ``flux``, ``flux_up``, ``flux_down``, ``mean_intensity`` (``[wavenumber][level]``), the
   ``spectrum`` (upward flux at the top of the atmosphere per micron), the
   ``flux_divergence`` :math:`\dd\Fnet/\dd p` per level, the ``net_heating`` per level and,
   when a Newton solver requested them, the Jacobians ``net_flux_jacobian`` and
   ``net_heating_jacobian`` (:doc:`../theory/radiative_transfer`).

.. py:attribute:: surface_temperature
   :noindex:

   Terrestrial planets only: the surface temperature in K.

.. py:method:: cloud_optical_depth(wavelength=1.0)
   :noindex:

   Summed cloud optical depth per layer (``nb_grid_points - 1`` values, 0 = bottom) at the
   spectral point closest to ``wavelength`` (micron), without the ramp factor.

.. py:method:: cloud_optical_depth_per_module(wavelength=1.0)
   :noindex:

   The same per cloud module.

.. py:method:: cloud_optics(wavelength=1.0)
   :noindex:

   The combined ``(optical_depth, single_scattering, asymmetry)`` per layer the radiative
   transfer actually uses, including the ramp factor.

.. py:method:: cloud_diagnostics()
   :noindex:

   Per cloud module a dict of per-layer diagnostics; the Ackerman & Marley module reports
   ``condensate_mmr``, ``column_density``, ``r_g``, ``r_eff``, ``r_w``, ``alpha`` and ``fsed``,
   prescribed modules report nothing.

.. py:method:: cloud_ramp_scale()
   :noindex:

   The current factor on the cloud optical depths (1 = full cloud).

.. py:method:: cloud_module_types()
   :noindex:

   The type keywords of the configured cloud modules, in order.

.. py:method:: recompute_chemistry(lagged)
   :noindex:

   Diagnostic: the composition the chemistry modules give at the current temperature on the
   committed path (``lagged=False``, e.g. including FastChem's rainout) or on the Newton's trial
   path (``lagged=True``, the depleted gas phase). Returns ``number_densities[level][species]``
   and leaves the model unchanged.

.. py:method:: eval_radiation(compute_jacobian)
   :noindex:

   Brown dwarfs: one radiative transfer solve at the current state (chemistry, structure and
   opacity recomputed) with or without the temperature Jacobians; returns the net flux per
   level.

.. py:method:: eval_forward(temperature, compute_jac, recompute_opacity=True)
   :noindex:

   Terrestrial planets: one forward solve at a given temperature profile; returns the net flux
   per level and fills the Jacobians in ``radiation_field`` if requested. Useful for
   finite-difference checks of the analytic Jacobian.

Helper functions
----------------

.. py:function:: baerli.mixing_ratios(model, species=None)
   :noindex:

   ``(symbols, table)`` with ``table[level][k]`` the volume mixing ratio of ``symbols[k]`` on
   the model's grid (default: every species Bärli knows, without the ``Total`` pseudo-species).
   The table has the layout :meth:`set_composition` expects.

.. py:function:: baerli.species_symbols()
   :noindex:

   The list of all chemical species symbols, in the order of the species index used by
   ``number_densities``.

.. _pressure-grid:

The pressure grid
-----------------

The ``nb_grid_points`` levels are distributed uniformly in :math:`\ln p` between
``boundary_pressures[0]`` (bottom, level 0) and ``boundary_pressures[1]`` (top, the last
level). Layer quantities such as cloud optical depths refer to the interval between level
:math:`j` and level :math:`j+1`. The altitude of each level follows from hydrostatic
equilibrium with the current composition and temperature (:doc:`../theory/atmosphere`); the
bottom level sits at altitude zero.

For self-luminous objects the bottom pressure should be deep enough that the lowest levels are
convective and optically thick at all wavelengths (100 to 1000 bar); for terrestrial planets it
is the surface pressure. Resolution matters for the Newton solvers: the radiative-convective
boundary of a terrestrial planet is only properly resolved with about 200 levels or more, and
an under-resolved boundary shows up as a kink in the flux profile.
