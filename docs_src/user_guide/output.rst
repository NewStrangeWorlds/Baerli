Output files
============

:func:`baerli.save_model` writes a model to a NetCDF4 file (through xarray with the
``h5netcdf`` engine). The file carries the atmospheric structure, the radiation field, the
spectrum, the cloud and condensate diagnostics, and the complete configuration the model was
built from, so a run can be inspected, replotted and reproduced from the file alone.

.. code-block:: python

   baerli.save_model("output.nc", model)            # grid defaults to the model's
   baerli.save_model("output.nc", model, grid, attrs={"comment": "test run"})

Variables
---------

.. list-table::
   :header-rows: 1
   :widths: 26 22 12 40

   * - Variable
     - Dimensions
     - Units
     - Meaning
   * - ``pressure``
     - ``grid_point``
     - bar
     - Level pressures, bottom to top.
   * - ``temperature``
     - ``grid_point``
     - K
     -
   * - ``altitude``
     - ``grid_point``
     - cm
     - Height above the bottom level.
   * - ``scale_height``
     - ``grid_point``
     - cm
     - Pressure scale height :math:`k_B T / (\bar m g)`.
   * - ``mass_density``
     - ``grid_point``
     - g cm\ :sup:`-3`
     -
   * - ``mean_molecular_weight``
     - ``grid_point``
     - amu
     -
   * - ``convective``
     - ``grid_point``
     -
     - 1 in convective levels, 0 otherwise.
   * - ``kzz``
     - ``grid_point``
     - cm\ :sup:`2` s\ :sup:`-1`
     - Eddy diffusion coefficient.
   * - ``number_densities``
     - ``grid_point``, ``species``
     - cm\ :sup:`-3`
     - All species in the order of the coordinate ``species_name``.
   * - ``condensate_densities``
     - ``grid_point``, ``condensate``
     - cm\ :sup:`-3`
     - Only with FastChem condensation: formula units per cm\ :sup:`3`; in rainout mode the
       amount condensing *at* each level. Coordinate ``condensate_symbol``.
   * - ``cloud_optical_depth``
     - ``layer``
     -
     - Summed cloud optical depth per layer at 1 micron (without the ramp factor).
   * - ``cloud<k>_<name>``
     - ``layer``
     - various
     - Per-module diagnostics of the microphysical cloud modules (``condensate_mmr``,
       ``column_density`` in cm\ :sup:`-2`, ``r_g``, ``r_eff``, ``r_w`` in micron, ``alpha``,
       ``fsed``).
   * - ``flux_total``
     - ``grid_point``
     - erg cm\ :sup:`-2` s\ :sup:`-1`
     - Wavelength-integrated **net** flux (upward positive).
   * - ``flux_up_total``, ``flux_down_total``
     - ``grid_point``
     - erg cm\ :sup:`-2` s\ :sup:`-1`
     - Integrated upward and downward fluxes; the stellar beam is folded into the downward
       flux.
   * - ``mean_intensity_total``
     - ``grid_point``
     - erg cm\ :sup:`-2` s\ :sup:`-1`
     - Integrated mean intensity :math:`4\pi J`.
   * - ``flux_divergence``
     - ``grid_point``
     - erg cm\ :sup:`-2` s\ :sup:`-1` bar\ :sup:`-1`
     - :math:`\dd \Fnet / \dd p`.
   * - ``spectrum``
     - ``wavenumber_index``
     - erg cm\ :sup:`-2` s\ :sup:`-1` micron\ :sup:`-1`
     - Upward flux at the top of the atmosphere **per micron** (the ``units`` attribute written
       by the current ``save_model`` still says per cm\ :sup:`-1`).
   * - ``wavenumber``, ``wavelength``
     - ``wavenumber_index``
     - cm\ :sup:`-1`, micron
     - The spectral grid.

Attributes
----------

* ``config``: the full configuration as a JSON string, ``{"grid": grid.config, "model":
  model.config}``. The model part contains the class under ``type``, every keyword argument,
  and the ``initial_profile`` spec.
* Every scalar entry of the model configuration is also stored as a flat attribute with the
  nested keys joined by dots (``solver.convergence_threshold``, ``chemistry`` lists are not
  flattened), for quick inspection with ``ncdump -h``.
* ``surface_temperature`` for terrestrial planets.
* Any extra ``attrs`` passed to :func:`baerli.save_model`.

Reading files back
------------------

.. code-block:: python

   ds = baerli.load_model_data("output.nc")        # an xarray.Dataset
   ds["temperature"].plot(y="pressure", yscale="log", yincrease=False)

   T = baerli.load_temperature("output.nc")         # numpy array, bottom -> top
   config = baerli.load_output_config("output.nc")  # {"grid": ..., "model": ...}

   # rebuild the model from the file's configuration (initialised from its initial profile)
   model, grid = baerli.build_model(config)

   # or restart a new model from the saved profile and composition
   model.initialize_from_file("output.nc")

``build_model`` reproduces the model as it was *configured*; ``initialize_from_file`` loads
the *converged* profile. Combine both to continue a run with a changed parameter:

.. code-block:: python

   config = baerli.load_output_config("output.nc")
   config["model"]["solver"]["convergence_threshold"] = 1e-6
   config["initial_profile"] = {"type": "file", "file": "output.nc"}
   model, grid = baerli.build_model(config)
   model.compute(warm_start=True)
