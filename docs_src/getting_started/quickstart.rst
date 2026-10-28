Quickstart
==========

This page runs through a complete model: a brown dwarf, then the two irradiated objects. The
full scripts are in the repository root (``run_brown_dwarf.py``, ``run_gas.py``,
``run_earth.py``, ``run_venus.py``) and are discussed in :doc:`../user_guide/examples`.

A brown dwarf
-------------

**1. The spectral grid.** Every model needs a :class:`baerli.SpectralGrid`, a subset of the
native wavenumber list of the opacity data. The common choice is a grid of constant resolving
power:

.. code-block:: python

   import baerli

   opacity_path = "/media/data/opacity_data/helios-k/"

   grid = baerli.SpectralGrid.constant_resolution(
       opacity_path, resolution=1000.0, wavelength_min=0.3, wavelength_max=100.0)

**2. The opacity sources.** A list of ``(species symbol, data folder)`` pairs. The symbol is one
of Bärli's chemical species (:func:`baerli.species_symbols`), the folder is relative to
``opacity_path``; ``"Rayleigh"`` and ``"none"`` select built-in continuum sources
(:doc:`../user_guide/opacities`):

.. code-block:: python

   opacity_species = [
       ("CIA-H2-H2", "CIA/H2-H2"),
       ("CIA-H2-He", "CIA/H2-He"),
       ("K",         "Alkali_Allard/K"),
       ("Na",        "Alkali_Allard/Na"),
       ("H2O",       "Molecules/1H2-16O__POKAZATEL_e2b"),
       ("CH4",       "Molecules/12C-1H4__YT34to10_e2b"),
       ("NH3",       "Molecules/14N-1H3__CoYuTe_e2b"),
       ("CO",        "Molecules/12C-16O__Li2015_e2b"),
       ("CO2",       "Molecules/12C-16O2__CDSD_4000_e2b"),
   ]

**3. The model.** Physical parameters are plain keyword arguments; the pluggable components
are module specs (:ref:`module-specs`):

.. code-block:: python

   model = baerli.BrownDwarf(
       grid,
       effective_temperature=1000.0,      # K
       surface_gravity=10**4.5,           # cm/s^2
       nb_grid_points=100,
       boundary_pressures=[1e2, 1e-6],    # bar, bottom -> top
       opacity_species=opacity_species,
       chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters.dat",
                                       metallicity=1.0, c_to_o=0.5))],
       radiative_transfer=("disort", dict(nb_streams=4)),
       convection="mlt_dry",
       solver=("ratio_ul", dict(max_iterations=200, convergence_threshold=1e-5)))

**4. The initial profile.** The temperature iteration needs a starting profile and an initial
composition. The grey Milne-Eddington solution is a good start for self-luminous objects
(the effective temperature defaults to the model's):

.. code-block:: python

   model.initialize(("milne", dict(kappa_ross=1e-2)))

Alternatively restart from a previous output file:

.. code-block:: python

   model.initialize_from_file("output_brown_dwarf.nc")

**5. Iterate and save.**

.. code-block:: python

   converged = model.compute()          # True on convergence
   baerli.save_model("output_brown_dwarf.nc", model)

The iteration prints one line per step: the largest relative and absolute temperature change
and the level where it occurred, the flux residual, the bottom and top temperatures and the
number of convective levels (:doc:`../user_guide/solvers`). After convergence the results are
available on the model object:

.. code-block:: python

   import numpy as np

   atm = model.atmosphere
   rf = model.radiation_field

   pressure = np.array(atm.pressure)            # bar, bottom -> top
   temperature = np.array(atm.temperature)      # K
   wavelength = np.array(grid.wavelength_list)  # micron
   spectrum = np.array(rf.spectrum)             # erg/cm^2/s/cm^-1 at the top of the atmosphere

An irradiated gas planet
------------------------

A :class:`baerli.GasPlanet` adds the stellar irradiation: the incident flux at the planet,
the cosine of the zenith angle and the stellar spectrum. The internal temperature replaces the
effective temperature. A Guillot (2010) profile is the natural start:

.. code-block:: python

   model = baerli.GasPlanet(
       grid,
       internal_temperature=100.0,        # K
       surface_gravity=2500.0,            # cm/s^2
       instellation=2.48e8,               # erg/cm^2/s
       zenith_angle=0.5,                  # cos(60 deg)
       stellar_spectrum=("tabulated", dict(file="data/stellar_spectra/spectrum_sun.dat")),
       boundary_pressures=[1e3, 1e-6],
       opacity_species=opacity_species,
       chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters.dat"))],
       radiative_transfer=("adding_doubling", dict(nb_streams=2)),
       convection="mlt_dry",
       solver=("ratio_ul", dict(convergence_threshold=1e-5)))

   model.initialize(("guillot", dict(kappa_ir=1e-2, t_irr=1446.0, gamma=0.4, f=0.25)))
   model.compute()

A terrestrial planet
--------------------

A :class:`baerli.TerrestrialPlanet` has a surface instead of an internal heat flux. The
composition of a modern Earth is prescribed from a file, water follows a relative-humidity
profile if desired, and the convection scheme is the moist mixing-length theory:

.. code-block:: python

   model = baerli.TerrestrialPlanet(
       grid,
       surface_gravity=980.0,
       instellation=1361e3 * 0.5,         # solar constant, 0.5 for a fast rotator
       zenith_angle=0.5,
       stellar_spectrum=("tabulated", dict(file="data/stellar_spectra/spectrum_sun.dat")),
       surface=("simple", dict(albedo=0.3, wavelength_switch=2.0)),
       boundary_pressures=[1.0, 1e-5],
       opacity_species=[("N2", "Rayleigh"), ("O2", "Rayleigh"),
                        ("H2O", "Molecules/H2O_HITRAN_cut25"),
                        ("CO2", "Molecules/12C-16O2__CDSD_4000_e2b"),
                        ("CO2-CIA", "none"), ("H2O-CIA", "none")],
       chemistry=[("fixed", dict(file="data/Earth/earth_standard_composition.dat")),
                  ("manabe_wetherald", dict(surface_rh=0.77))],
       radiative_transfer=("adding_doubling", dict(nb_streams=2)),
       convection=("mlt_moist", dict(min_pressure=1e-2)),
       solver=("ratio_ul", dict(convergence_threshold=1e-5)))

   model.initialize(("adiabat", dict(surface_temperature=288.0, stratosphere_temperature=160.0)))
   model.compute()
   print(model.surface_temperature)

Running from a config file
--------------------------

The same model can be written as a TOML, YAML or JSON file whose sections mirror the keyword
arguments (:ref:`config-files`) and run from the command line:

.. code-block:: sh

   python -m baerli configs/earth.toml -o output_earth.nc

Where to go next
----------------

* :doc:`../user_guide/configuration` for the configuration mechanism and config files.
* :doc:`../user_guide/models` for all constructor arguments and model methods.
* :doc:`../user_guide/solvers` if a model does not converge.
* :doc:`../theory/overview` for the physics.
