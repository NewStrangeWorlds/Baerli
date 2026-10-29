Bärli
=====

**Bärli** (ASCII: *Baerli*), developed at the  University of Bern, computes one-dimensional 
radiative-convective equilibrium atmospheres of brown dwarfs, gas giant planets and terrestrial planets.
In Swiss German,  **Bärli** means "little bear", a nod to the fact that it inherits some of its code from 
the BeAR retrieval codebase.

A compiled C++ core does the heavy lifting: opacity sampling of line-by-line cross sections, 
multi-stream radiative transfer with analytic temperature Jacobians, equilibrium and disequilibrium 
chemistry, mixing-length convection, condensate clouds, and a family of Newton-type and relaxation 
solvers for the temperature structure. Everything is driven from the Python package ``baerli``.

.. code-block:: python

   import baerli

   grid = baerli.SpectralGrid.constant_resolution(
       "/data/opacities/", resolution=1000.0, wavelength_min=0.3, wavelength_max=100.0)

   model = baerli.BrownDwarf(
       grid,
       effective_temperature=1000.0,      # K
       surface_gravity=10**4.5,           # cm/s^2
       opacity_species=[("CIA-H2-H2", "CIA/H2-H2"),
                        ("H2O", "Molecules/1H2-16O__POKAZATEL_e2b")],
       chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters.dat"))],
       convection="mlt_dry",
       solver=("ratio_ul", dict(convergence_threshold=1e-5)))

   model.initialize(("milne", dict(kappa_ross=1e-2)))
   model.compute()
   baerli.save_model("output_brown_dwarf.nc", model)

The documentation is organised in four parts:

* :doc:`Getting started <getting_started/overview>` explains what Bärli models, how to build
  it and how to run a first model.
* The :doc:`user guide <user_guide/configuration>` is the practical reference of every model
  option: the objects, the spectral grid, opacities, chemistry, clouds, convection, radiative
  transfer, solvers, output files and the kinetics coupling.
* The :doc:`theory <theory/overview>` chapters document the physics and the numerics that are
  actually implemented, with the equations and the literature they come from.
* The :doc:`reference <reference/parameters>` collects every configuration parameter in one
  place, the Python API and the bibliography.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   getting_started/overview
   getting_started/installation
   getting_started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: User guide

   user_guide/configuration
   user_guide/models
   user_guide/spectral_grid
   user_guide/opacities
   user_guide/chemistry
   user_guide/clouds
   user_guide/convection
   user_guide/radiative_transfer
   user_guide/solvers
   user_guide/output
   user_guide/kinetics_coupling
   user_guide/examples

.. toctree::
   :maxdepth: 2
   :caption: Theory

   theory/overview
   theory/atmosphere
   theory/radiative_transfer
   theory/opacities
   theory/chemistry
   theory/clouds
   theory/convection
   theory/temperature_iteration

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/parameters
   reference/python_api
   reference/bibliography
   reference/glossary

Indices
-------

* :ref:`genindex`
* :ref:`search`
