Overview
========

What Bärli computes
-------------------

Bärli solves the one-dimensional radiative-convective equilibrium (RCE) of a planetary or
substellar atmosphere. Given the physical parameters of an object (gravity, internal or
effective temperature, stellar irradiation, surface properties) and a set of opacity sources,
it iterates the temperature profile :math:`T(p)` on a fixed pressure grid until the net energy
flux through every layer is constant: radiation carries the flux in radiative layers, convection
transports whatever radiation cannot in convectively unstable layers. The converged model
provides the temperature structure, the chemical composition, the cloud structure, the radiation
field inside the atmosphere and the emergent spectrum.

Three object classes cover the usual cases (see :doc:`../user_guide/models`):

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Object
     - Lower boundary
     - Upper boundary
     - Flux constraint
   * - :class:`BrownDwarf`
     - semi-infinite, internal flux :math:`\sigma \Teff^4`
     - no irradiation
     - :math:`\Fnet(p) = \sigma \Teff^4`
   * - :class:`GasPlanet`
     - semi-infinite, internal flux :math:`\sigma \Tint^4`
     - stellar beam (instellation, zenith angle)
     - :math:`\Fnet(p) = \sigma \Tint^4` (thermal + stellar)
   * - :class:`TerrestrialPlanet`
     - solid surface, surface energy balance
     - stellar beam (instellation, zenith angle)
     - :math:`\Fnet(p) = 0`, :math:`\Fnet(\mathrm{surface}) = 0`

Architecture
------------

Bärli is a C++17 library wrapped with pybind11. The compiled core ``baerli._baerli`` is placed
inside the pure-Python package ``baerli`` and is never used directly: the package adds
configuration bookkeeping, config files and the output format.

.. figure:: /_static/architecture.svg
   :width: 100%
   :alt: Components of a Bärli model and the order in which an iteration evaluates them.

   The components of a model and the order in which one iteration of the temperature solver
   evaluates them. Each box with a type keyword is a pluggable *module* selected by a
   :ref:`module spec <module-specs>`.

Every model is assembled from pluggable **modules**. Each module family has a fixed set of
types, and every type accepts only its own named parameters:

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Family
     - Types
     - Chapter
   * - ``chemistry`` (list)
     - ``equilibrium``, ``isoprofile``, ``fixed``, ``manabe_wetherald``, ``quench``, ``external``
     - :doc:`../user_guide/chemistry`
   * - ``clouds`` (list)
     - ``fixed``, ``mie_slab``, ``ackerman_marley``
     - :doc:`../user_guide/clouds`
   * - ``cloud_ramp``
     - ``default``, ``none``
     - :doc:`../user_guide/clouds`
   * - ``radiative_transfer``
     - ``disort``, ``adding_doubling``
     - :doc:`../user_guide/radiative_transfer`
   * - ``convection``
     - ``mlt_dry``, ``mlt_moist``, ``dry``, ``moist``, ``none``
     - :doc:`../user_guide/convection`
   * - ``kzz``
     - ``mlt``, ``constant``, ``power_law``
     - :doc:`../user_guide/convection`
   * - ``solver``
     - ``ratio_ul``, ``flux_divergence``, ``ptc``, ``time_stepping``, ``time_stepping_lre``, ``helios``
     - :doc:`../user_guide/solvers`
   * - ``stellar_spectrum``
     - ``tabulated``, ``blackbody``
     - :doc:`../user_guide/radiative_transfer`
   * - ``surface``
     - ``blackbody``, ``simple``, ``variable_albedo``
     - :doc:`../user_guide/radiative_transfer`
   * - initial profile
     - ``adiabat``, ``milne``, ``const``, ``guillot``, file restart
     - :doc:`../user_guide/solvers`

The spectral grid is a separate object shared by models (:doc:`../user_guide/spectral_grid`),
and the opacity sources are a list of species and data folders (:doc:`../user_guide/opacities`).

Dependencies
------------

The build fetches the external components at configure time and pins their versions
(:doc:`installation`):

* `DisORT++ <https://github.com/NewStrangeWorlds/DisORT>`_, a C++ discrete-ordinate radiative
  transfer solver with analytic flux and mean-intensity temperature Jacobians;
* `AddingDoublingRT <https://github.com/NewStrangeWorlds/AddingDoublingRT>`_, an adding-doubling
  radiative transfer solver with the same Jacobian interface;
* `FastChem <https://github.com/NewStrangeWorlds/FastChem>`_, the equilibrium chemistry
  including its condensed phase;
* `LX-MIE <https://github.com/NewStrangeWorlds/LX-MIE>`_, Mie efficiencies for the cloud optics;
* `neoVULCAN <https://github.com/NewStrangeWorlds/neoVULCAN>`_ (optional), the chemical kinetics
  code of the :doc:`kinetics coupling <../user_guide/kinetics_coupling>`;
* Eigen (through DisORT++) and pybind11.

Units
-----

Bärli works in cgs units throughout. Pressures are given in **bar** in every user-facing
interface (constructor arguments, module parameters, output), temperatures in K, gravities in
cm s\ :sup:`-2`, fluxes in erg cm\ :sup:`-2` s\ :sup:`-1`, number densities in cm\ :sup:`-3`,
eddy diffusion coefficients in cm\ :sup:`2` s\ :sup:`-1`, lengths in cm except for wavelengths
and particle radii, which are in micron. Spectra are tabulated per wavenumber (cm\ :sup:`-1`).

Pressure grids run from the **bottom** (index 0, highest pressure) to the **top** of the
atmosphere; every per-level array in the output and in the Python interface follows this order.
