Configuring a model
===================

A Bärli model is configured once, in the constructor of its object class, plus one call that
sets the initial state. This page describes the configuration mechanism; the individual
parameters are documented in the chapters that follow and collected in
:doc:`../reference/parameters`.

Keyword arguments
-----------------

All constructor arguments are keyword-only, apart from the spectral grid that comes first.
They fall into two groups:

* **Object parameters** are plain numbers or strings: ``effective_temperature``,
  ``surface_gravity``, ``instellation``, ``zenith_angle``, ``nb_grid_points``,
  ``boundary_pressures``, ``opacity_path``, ``opacity_species``, ``radius``,
  ``variable_gravity``. They are listed per object in :doc:`models`.
* **Module specs** select and parametrise a pluggable component: ``chemistry``, ``clouds``,
  ``cloud_ramp``, ``radiative_transfer``, ``convection``, ``kzz``, ``solver``,
  ``stellar_spectrum``, ``surface``, and the ``profile`` argument of
  :meth:`~baerli.BrownDwarf.initialize`.

.. _module-specs:

Module specs
------------

A module spec names a *type* and gives that type's *named parameters*. Three equivalent
spellings are accepted everywhere:

.. code-block:: python

   convection = "mlt_dry"                                  # a bare type name
   convection = ("mlt_dry", {"alpha": 1.0})                # a (type, {parameters}) tuple
   convection = {"type": "mlt_dry", "alpha": 1.0}          # a dict with a "type" key

The families ``chemistry`` and ``clouds`` take a **list** of specs that are applied in order
(a single spec is also accepted and treated as a one-element list). Every other family takes
exactly one spec.

Parameter values are handed to the compiled core as strings and converted by the module that
consumes them. Numbers, booleans (``True``/``False``, ``"yes"``/``"no"``, ``"on"``/``"off"``,
``1``/``0``) and strings are accepted; lists are joined with commas (used e.g. by
``mie_radii``). Integers must be whole numbers.

Validation
~~~~~~~~~~

Each module reads exactly the parameters it understands and then **rejects every key it was not
asked for**. The error lists the accepted parameters:

.. code-block:: text

   Unknown parameter for convection 'mlt_dry': 'max_sweeps'. Accepted parameters: alpha, min_pressure

This is deliberate: a knob that belongs to one scheme can never be passed silently to another,
so a misspelt or misplaced parameter is caught at construction time. Missing required
parameters and values of the wrong kind (``"abc"`` for a number, ``1.5`` for an integer) raise
the same kind of error.

Unknown type names raise an error that lists the available types of the family. Some families
accept short aliases (e.g. ``eq`` for ``equilibrium``); the documentation uses the long names.

Defaults
~~~~~~~~

Omitted components get these defaults (the Python layer records them in ``model.config`` for
provenance):

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Family
     - Default
     - Note
   * - ``radiative_transfer``
     - ``("disort", {"nb_streams": 4})``
     -
   * - ``convection``
     - ``"mlt_dry"``
     - requires the ``ratio_ul`` solver
   * - ``kzz``
     - ``"mlt"``
     - eddy diffusion from the mixing-length convection
   * - ``solver``
     - ``"ratio_ul"``
     -
   * - ``cloud_ramp``
     - ``"default"``
     - only active with cloud modules
   * - ``clouds``
     - none
     -
   * - ``chemistry``
     - **required**
     - at least one module

Pairing rules
~~~~~~~~~~~~~

A few combinations are checked at construction time:

* The mixing-length convection schemes ``mlt_dry`` and ``mlt_moist`` require the ``ratio_ul``
  solver. The other solvers work with ``dry``, ``moist`` and ``none``.
* The ``ackerman_marley`` cloud module requires an ``equilibrium`` chemistry module with
  ``condensation = "rainout"`` and a FastChem parameter file that lists the condensate data.
* ``quench`` must follow ``equilibrium`` in the chemistry list; ``external`` must be last.
* ``set_composition`` requires an ``external`` chemistry module, ``set_cloud_optical_properties``
  a ``fixed`` cloud module without a file.

Configuration record and provenance
-----------------------------------

The Python classes remember how they were built: ``grid.config`` holds the factory name and its
arguments, ``model.config`` the model type and every keyword argument (module specs normalised
to dicts, tuples to lists, numpy scalars to Python numbers), plus the ``initial_profile`` spec
once :meth:`initialize` or :meth:`initialize_from_file` has been called.
:func:`baerli.save_model` writes both into the output file as the JSON attribute ``config``, and
:func:`baerli.build_model` rebuilds a model from it (:doc:`output`).

.. _config-files:

Config files
------------

The same configuration can be written as a TOML, YAML or JSON file. The sections mirror the
Python interface:

.. code-block:: toml

   [grid]
   type = "constant_resolution"          # a SpectralGrid factory name + its arguments
   opacity_path = "/data/opacities/"
   resolution = 1000.0
   wavelength_min = 0.3
   wavelength_max = 100.0

   [model]
   type = "BrownDwarf"                   # BrownDwarf | GasPlanet | TerrestrialPlanet
   effective_temperature = 1000.0        # ... followed by the class's keyword arguments
   surface_gravity = 31622.8
   opacity_species = [["H2O", "Molecules/1H2-16O__POKAZATEL_e2b"], ["CH4", "..."]]
   chemistry = [{type = "equilibrium", parameter_file = "fastchem_parameters.dat"}]
   radiative_transfer = {type = "disort", nb_streams = 4}
   convection = "mlt_dry"
   solver = {type = "ratio_ul", max_iterations = 200, convergence_threshold = 1e-5}

   [initial_profile]
   type = "milne"
   kappa_ross = 1e-2

   [output]
   file = "output_brown_dwarf.nc"

* ``[grid]`` is a :class:`~baerli.SpectralGrid` factory name (``constant_resolution``,
  ``constant_wavenumber_step``, ``constant_wavelength_step``, ``covering``) and its arguments.
* ``[model]`` holds ``type`` (``BrownDwarf``, ``GasPlanet``, ``TerrestrialPlanet`` or their
  snake_case spellings) and the class's keyword arguments. Module specs are inline tables with
  a ``type`` entry or bare type names; ``opacity_species`` is a list of two-element lists.
* ``[initial_profile]`` is the spec of :meth:`initialize`; ``{type = "file", file = "run.nc"}``
  restarts from a saved model instead. An optional ``initial_chemistry`` list overrides the
  chemistry used for the initial composition.
* ``[output]`` names the file written after :meth:`compute`.

The ``model`` section of an output file's ``config`` attribute has the same structure, so an
output file can be rerun unchanged.

Run a config file from the shell or from Python:

.. code-block:: sh

   python -m baerli configs/earth.toml            # writes [output].file
   python -m baerli configs/earth.toml -o run.nc  # overrides the output file

.. code-block:: python

   config = baerli.load_config("configs/earth.toml")
   model, grid = baerli.build_model(config)       # built and initialised
   model.compute()

   # or everything at once
   model, grid = baerli.run("configs/earth.toml")

Relative paths in a configuration (``parameter_file``, data files, opacity folders) are
resolved against the current working directory, not against the config file.
