Configuration parameters
========================

Every parameter of a Bärli model on one page. Module specs are written as ``type {parameter =
default, ...}``; a parameter without default is required. Units: pressures in bar, temperatures
in K, gravity in cm s\ :sup:`-2`, fluxes in erg cm\ :sup:`-2` s\ :sup:`-1`, wavelengths and
particle radii in micron, :math:`K_{zz}` in cm\ :sup:`2` s\ :sup:`-1`.

Spectral grid
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Factory
     - Arguments
   * - ``constant_resolution``
     - ``opacity_path, resolution, wavelength_min, wavelength_max, wavenumber_file=""``
   * - ``constant_wavenumber_step``
     - ``opacity_path, step [cm^-1], wavelength_min, wavelength_max, wavenumber_file=""``
   * - ``constant_wavelength_step``
     - ``opacity_path, step [micron], wavelength_min, wavelength_max, wavenumber_file=""``
   * - ``covering``
     - ``opacity_path, wavelength_min, wavelength_max, temperature_min, temperature_max,
       nb_points, nb_temperatures=0, stellar_temperature=0, nb_points_stellar=0, wavenumber_file=""``

Object parameters
-----------------

.. list-table::
   :header-rows: 1
   :widths: 24 20 56

   * - Argument
     - Class
     - Meaning
   * - ``effective_temperature``
     - BrownDwarf
     - K; internal flux :math:`\sigma \Teff^4`.
   * - ``internal_temperature``
     - GasPlanet
     - K; internal flux :math:`\sigma \Tint^4`.
   * - ``surface_gravity``
     - all
     - cm s\ :sup:`-2`.
   * - ``radius = 0.0``
     - BrownDwarf, GasPlanet
     - cm; bottom-of-domain radius for ``variable_gravity``.
   * - ``variable_gravity = False``
     - BrownDwarf, GasPlanet
     - Gravity decreasing with altitude.
   * - ``instellation``
     - GasPlanet, TerrestrialPlanet
     - erg cm\ :sup:`-2` s\ :sup:`-1`, normal to the beam.
   * - ``zenith_angle``
     - GasPlanet, TerrestrialPlanet
     - Cosine of the stellar zenith angle.
   * - ``stellar_spectrum``
     - GasPlanet, TerrestrialPlanet
     - Module spec (below).
   * - ``surface``
     - TerrestrialPlanet
     - Module spec (below).
   * - ``nb_grid_points = 100``
     - all
     - Number of pressure levels.
   * - ``boundary_pressures = [1e2, 1e-6]``
     - all
     - ``[bottom, top]`` in bar.
   * - ``opacity_path = grid's``
     - all
     - Cross-section data root.
   * - ``opacity_species = []``
     - all
     - List of ``(symbol, folder)`` pairs.
   * - ``chemistry``
     - all
     - List of module specs (required).
   * - ``clouds = []``
     - all
     - List of module specs.
   * - ``cloud_ramp = "default"``
     - all
     - Module spec.
   * - ``radiative_transfer = disort``
     - all
     - Module spec.
   * - ``convection = "mlt_dry"``
     - all
     - Module spec.
   * - ``kzz = "mlt"``
     - all
     - Module spec.
   * - ``solver = "ratio_ul"``
     - all
     - Module spec.

Chemistry modules
-----------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``equilibrium`` (``eq``)
     - ``parameter_file``, ``metallicity = 1.0``, ``c_to_o = 0.5``,
       ``condensation = "none" | "equilibrium" | "rainout"``, ``condensation_acceleration = true``
   * - ``isoprofile`` (``iso``)
     - ``<species symbol> = <mixing ratio>``, ... (at least one)
   * - ``fixed`` (``fix``)
     - ``file``
   * - ``manabe_wetherald`` (``mw``)
     - ``surface_rh = 0.77``
   * - ``quench`` (``q``)
     - ``metallicity = 1.0``, ``relax = true``
   * - ``external`` (``ext``)
     - none

Cloud modules
-------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``fixed`` (``fix``)
     - ``file`` (optional; without it the table comes from :meth:`set_cloud_optical_properties`)
   * - ``mie_slab`` (``mie``)
     - ``condensate``, ``particle_radius`` [micron], ``sigma = 1.0``, ``top_pressure``,
       ``bottom_pressure``, ``optical_depth`` with ``reference_wavelength = 1.0`` **or**
       ``column_density`` [cm\ :sup:`-2`], ``mie_resolution = 200``, ``mie_radii = 120``
   * - ``ackerman_marley`` (``am``)
     - ``condensate``, ``fsed = 1.0``, ``fsed_mode = "const" | "exp"``, ``b = 1.0``,
       ``eps = 1e-2``, ``sigma = 2.0``, ``kzz_min = 1e5``,
       ``mixing_length = "scale_height" | "lapse_ratio"``, ``relax = 0.5``, ``max_change = 0.3``
       [dex], ``tolerance = 0.05`` [dex], ``gate = 1e-2``, ``seed = 1e-2``, ``bulk_density``
       [g cm\ :sup:`-3`, table default], ``molecule_diameter = 2.827e-8`` [cm],
       ``lj_temperature = 59.7`` [K], ``mie_resolution = 200``, ``mie_radii = 120``

Cloud ramp
----------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``default``
     - ``tau_start = 1.0``, ``gate = 1e-2``, ``factor = 2.0`` (> 1), ``backoff = 2.0`` (> 0)
   * - ``none``
     - none

Radiative transfer
------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``disort``
     - ``nb_streams = 4`` (total; 4, 6, 8, 10, 12, 14, 16, 32, 64)
   * - ``adding_doubling`` (``ad``)
     - ``nb_streams = 2`` (per hemisphere, at least 2)

Stellar spectrum
----------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``blackbody``
     - ``temperature`` [K]
   * - ``tabulated``
     - ``file``

Surface
-------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``blackbody``
     - none
   * - ``simple``
     - ``albedo``, ``wavelength_switch`` [micron]
   * - ``variable_albedo``
     - ``file``

Convection
----------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``mlt_dry``, ``mlt_moist``
     - ``alpha = 1.0``, ``min_pressure = 1e-4`` (brown dwarf, gas planet) or ``1e-3``
       (terrestrial); require ``solver = "ratio_ul"``
   * - ``dry``, ``moist``
     - ``min_pressure`` (same defaults), ``max_sweeps = 10``
   * - ``none``
     - none

Eddy diffusion
--------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``mlt``
     - ``scaling = "velocity" | "flux"``, ``radiative = "constant" | "power_law" | "fixed"``,
       ``slope`` (required for ``power_law``), ``value`` (required for ``fixed``),
       ``min = 1e4``, ``relax = 0.5``, ``tolerance = 0.05`` [dex]
   * - ``constant``
     - ``value``
   * - ``power_law``
     - ``value``, ``pressure = 1.0`` [bar], ``slope``

Solver
------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - all
     - ``max_iterations = 100``, ``convergence_threshold = 1e-4``
   * - ``ratio_ul`` (``rce``), ``flux_divergence`` (``flux``)
     - no further parameters
   * - ``ptc`` (``lin``)
     - ``max_change = 0.1``
   * - ``time_stepping`` (``ts``)
     - ``gamma = 0.5``, ``ng_interval = 10``, ``max_change = 0.1``
   * - ``time_stepping_lre`` (``ts_lre``)
     - ``gamma = 0.5``, ``ng_interval = 10``, ``max_change = 0.1``, ``lre_fraction = 0.5``
   * - ``helios`` (``hel``)
     - ``ng_interval = 0``, ``step_init = 10``, ``step_grow = 1.1``, ``step_shrink = 1.5``,
       ``adapt_interval = 20``, ``step_exponent = 0.1``, ``max_step = 500``,
       ``stencil = "backward" | "centered" | "forward"``, ``residual = "flux" | "heating"``

Initial profile (``initialize``)
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Type
     - Parameters
   * - ``adiabat``
     - ``surface_temperature``, ``stratosphere_temperature``
   * - ``milne``
     - ``kappa_ross`` [cm\ :sup:`2` g\ :sup:`-1`], ``effective_temperature`` (default: the
       object's)
   * - ``const``
     - ``temperature``
   * - ``guillot``
     - ``kappa_ir`` [cm\ :sup:`2` g\ :sup:`-1`], ``t_irr``, ``gamma``, ``t_int`` (default: the
       object's), ``mode = "isotropic"`` with ``f = 0.25`` or ``mode = "beam"`` with ``mu``
       (default: the object's zenith angle)
   * - file restart
     - :meth:`initialize_from_file(filename)`; in a config file ``{type = "file", file = "..."}``

Method arguments
----------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Call
     - Arguments
   * - ``compute``
     - ``tolerance = 0.0``, ``max_iterations = 0`` (``<= 0`` keeps the solver spec),
       ``warm_start = False``
   * - ``initialize``
     - ``profile``, ``chemistry = None``
   * - ``set_composition``
     - ``species``, ``mixing_ratios[level][k]``
   * - ``set_chemistry_enabled``
     - ``type``, ``enabled``
   * - ``set_cloud_optical_properties``
     - ``pressure``, ``wavelength``, ``optical_depth``, ``single_scattering``, ``asymmetry``
   * - ``KineticsCoupling``
     - ``model``, ``chem``, ``species = None``, ``rce_tolerance = (1e-2, 1e-5)``,
       ``dead_band = 0.05``, ``composition_tolerance = 0.05``, ``relax = 0.5``,
       ``tracked_species = None``, ``tracked_floor = 1e-10``, ``max_passes = 12``,
       ``seed_continuation = False``, ``workdir = None``, ``verbose = True``

Config file sections
--------------------

``[grid]`` (a factory ``type`` and its arguments), ``[model]`` (``type`` and the object
parameters, module specs as inline tables), ``[initial_profile]``, optional
``initial_chemistry`` list, ``[output]`` with ``file``.
