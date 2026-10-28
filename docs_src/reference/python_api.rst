Python API
==========

The ``baerli`` package is the only user interface. Its public names:

.. code-block:: python

   import baerli

   # classes
   baerli.SpectralGrid, baerli.BrownDwarf, baerli.GasPlanet, baerli.TerrestrialPlanet
   baerli.Atmosphere, baerli.RadiativeTransferOutput, baerli.KineticsCoupling

   # functions
   baerli.species_symbols, baerli.mixing_ratios, baerli.spec_to_dict
   baerli.save_model, baerli.load_model_data, baerli.load_temperature
   baerli.load_init_arrays, baerli.load_output_config
   baerli.load_config, baerli.build_model, baerli.run
   baerli.load_vulcan, baerli.neovulcan_dir, baerli.interpolate_log

   # submodules
   baerli.clouds          # helpers for prescribed cloud tables
   baerli.model_config_doc  # the shared configuration as a string

Spectral grid
-------------

.. py:class:: baerli.SpectralGrid

   A subset of the native wavenumber list of the opacity data (:doc:`../user_guide/spectral_grid`).
   Built with one of the factories; the raw constructors of the compiled class are not meant
   for direct use.

   .. py:classmethod:: constant_resolution(opacity_path, resolution, wavelength_min, wavelength_max, wavenumber_file="")
   .. py:classmethod:: constant_wavenumber_step(opacity_path, step, wavelength_min, wavelength_max, wavenumber_file="")
   .. py:classmethod:: constant_wavelength_step(opacity_path, step, wavelength_min, wavelength_max, wavenumber_file="")
   .. py:classmethod:: covering(opacity_path, wavelength_min, wavelength_max, temperature_min, temperature_max, nb_points, nb_temperatures=0, stellar_temperature=0.0, nb_points_stellar=0, wavenumber_file="")
   .. py:classmethod:: from_config(config)

      Build a grid from ``{"type": <factory name>, <its keyword arguments>}``.

   .. py:attribute:: opacity_path
   .. py:attribute:: wavenumber_list

      Wavenumbers in cm\ :sup:`-1`, increasing.

   .. py:attribute:: wavelength_list

      The corresponding wavelengths in micron (decreasing).

   .. py:attribute:: nb_spectral_points
   .. py:attribute:: nb_spectral_points_full

      Number of points of the grid and of the native list of the opacity data.

   .. py:attribute:: config

      The factory name and arguments the grid was built from.

   .. py:method:: wavelength_to_wavenumber(wavelength)
   .. py:method:: wavenumber_to_wavelength(wavenumber)

Model classes
-------------

.. py:class:: baerli.BrownDwarf(grid, *, effective_temperature, surface_gravity, radius=0.0, variable_gravity=False, **config)
.. py:class:: baerli.GasPlanet(grid, *, internal_temperature, surface_gravity, instellation, zenith_angle, stellar_spectrum, radius=0.0, variable_gravity=False, **config)
.. py:class:: baerli.TerrestrialPlanet(grid, *, surface_gravity, instellation, zenith_angle, stellar_spectrum, surface, **config)

   The three object classes (:doc:`../user_guide/models`). ``**config`` is the shared
   configuration: ``nb_grid_points``, ``boundary_pressures``, ``opacity_path``,
   ``opacity_species``, ``clouds``, ``cloud_ramp``, ``chemistry``, ``radiative_transfer``,
   ``convection``, ``solver``, ``kzz``. The model keeps a reference to the grid (``model.grid``)
   and the configuration (``model.config``).

   .. py:method:: initialize(profile, chemistry=None)
   .. py:method:: initialize_from_file(filename)
   .. py:method:: initialize_from_arrays(temperature, number_densities, mean_molecular_weight)
   .. py:method:: compute(tolerance=0.0, max_iterations=0, warm_start=False)
   .. py:method:: set_composition(species, mixing_ratios)
   .. py:method:: set_chemistry_enabled(type, enabled)
   .. py:method:: set_cloud_optical_properties(pressure, wavelength, optical_depth, single_scattering, asymmetry)
   .. py:method:: cloud_optical_depth(wavelength=1.0)
   .. py:method:: cloud_optical_depth_per_module(wavelength=1.0)
   .. py:method:: cloud_optics(wavelength=1.0)
   .. py:method:: cloud_diagnostics()
   .. py:method:: cloud_ramp_scale()
   .. py:method:: cloud_module_types()
   .. py:method:: recompute_chemistry(lagged)
   .. py:method:: eval_radiation(compute_jacobian)

      :class:`BrownDwarf` only.

   .. py:method:: eval_forward(temperature, compute_jac, recompute_opacity=True)

      :class:`TerrestrialPlanet` only.

   .. py:attribute:: atmosphere
   .. py:attribute:: radiation_field
   .. py:attribute:: surface_temperature

      :class:`TerrestrialPlanet` only.

   .. py:attribute:: config
   .. py:attribute:: grid

.. py:class:: baerli.Atmosphere

   Read-only view of the atmospheric state; all arrays run from the bottom level (index 0) to
   the top.

   .. py:attribute:: nb_grid_points
   .. py:attribute:: pressure

      bar

   .. py:attribute:: temperature

      K

   .. py:attribute:: altitude

      cm above the bottom level

   .. py:attribute:: scale_height

      cm

   .. py:attribute:: mass_density

      g cm\ :sup:`-3`

   .. py:attribute:: mean_molecular_weight

      amu

   .. py:attribute:: number_densities

      ``[level][species]`` in cm\ :sup:`-3`, species order of :func:`baerli.species_symbols`

   .. py:attribute:: convective

      1 in convective levels

   .. py:attribute:: kzz

      cm\ :sup:`2` s\ :sup:`-1`

   .. py:attribute:: condensate_symbols
   .. py:attribute:: condensate_weights

      g mol\ :sup:`-1`

   .. py:attribute:: condensate_densities

      ``[level][k]`` formula units per cm\ :sup:`3`; in rainout mode the amount condensing at
      each level

.. py:class:: baerli.RadiativeTransferOutput

   The radiation field of the last solve (:doc:`../user_guide/models`).

   .. py:attribute:: spectrum
   .. py:attribute:: flux_total
   .. py:attribute:: flux_up_total
   .. py:attribute:: flux_down_total
   .. py:attribute:: mean_intensity_total
   .. py:attribute:: flux
   .. py:attribute:: flux_up
   .. py:attribute:: flux_down
   .. py:attribute:: mean_intensity
   .. py:attribute:: flux_divergence
   .. py:attribute:: net_heating
   .. py:attribute:: net_flux_jacobian
   .. py:attribute:: net_heating_jacobian

Module specs
------------

.. py:function:: baerli.spec_to_dict(spec)

   Normalise a module spec (``"type"``, ``(type, {params})`` or ``{"type": ...}``) to a dict.

.. py:data:: baerli.model_config_doc

   The shared configuration (keyword arguments, module types and parameters) as a string, also
   appended to the class docstrings.

Species and composition
-----------------------

.. py:function:: baerli.species_symbols()

   All chemical species symbols in the order of the species index.

.. py:function:: baerli.mixing_ratios(model, species=None)

   ``(symbols, table)`` with the volume mixing ratios ``table[level][k]`` of ``symbols[k]``.

Output files
------------

.. py:function:: baerli.save_model(filename, model, grid=None, attrs=None)
.. py:function:: baerli.load_model_data(filename)

   The file as an :class:`xarray.Dataset`.

.. py:function:: baerli.load_temperature(filename)
.. py:function:: baerli.load_init_arrays(filename)

   ``{"temperature", "number_densities", "mean_molecular_weight"}`` as plain lists.

.. py:function:: baerli.load_output_config(filename)

   ``{"grid": ..., "model": ...}``.

Config files
------------

.. py:function:: baerli.load_config(path)

   Read a TOML, YAML or JSON file into a dict.

.. py:function:: baerli.build_model(config)

   ``(model, grid)`` from a config dict (file or output-file ``config``); the model is
   initialised if the config has an ``initial_profile``.

.. py:function:: baerli.run(config)

   Build, initialise, compute and save (if ``[output].file`` is given); ``config`` may be a dict
   or a file name. Returns ``(model, grid)``. Also the entry point of ``python -m baerli``.

Kinetics coupling
-----------------

.. py:class:: baerli.KineticsCoupling(model, chem, species=None, rce_tolerance=(1e-2, 1e-5), dead_band=0.05, composition_tolerance=0.05, relax=0.5, tracked_species=None, tracked_floor=1e-10, max_passes=12, seed_continuation=False, workdir=None, verbose=True)

   See :doc:`../user_guide/kinetics_coupling`.

   .. py:method:: run()

      Alternate until both codes agree; returns ``True`` on convergence.

   .. py:attribute:: history

      One dict per pass.

.. py:function:: baerli.load_vulcan(config_path, base_dir=None, cfg_overrides=None, regenerate_chemistry=True)
.. py:function:: baerli.neovulcan_dir()
.. py:function:: baerli.interpolate_log(p_from, values, p_to, log_values=False)

   Interpolate a profile (1D or ``(n, k)``) in :math:`\log_{10} p` with constant extension;
   with ``log_values`` on :math:`\log_{10}` of the values (floored at :math:`10^{-300}`).

Cloud helpers
-------------

.. py:function:: baerli.clouds.slab(optical_depth, top_pressure, bottom_pressure, single_scattering=0.0, asymmetry=0.0, wavelength=None)

   The keyword arguments of :meth:`set_cloud_optical_properties` for a cloud deck between two
   pressures with a uniform optical depth per unit :math:`\ln p`. Scalars give a grey deck;
   sequences give one value per ``wavelength`` (micron).

.. py:function:: baerli.clouds.write_table(filename, table, comment=None)

   Write such a table in the file format of the ``fixed`` cloud module.
