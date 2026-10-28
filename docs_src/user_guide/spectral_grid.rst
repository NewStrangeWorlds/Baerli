The spectral grid
=================

Every model is built on a :class:`baerli.SpectralGrid`, a subset of the native wavenumber list
of the opacity data. The radiative transfer is solved monochromatically at each point and all
spectral integrals are trapezoidal sums over the grid, so the grid controls both the accuracy and
the cost of a model (:doc:`../theory/opacities`).

Factories
---------

All factories take the opacity path first and the wavelength window in micron:

.. code-block:: python

   grid = baerli.SpectralGrid.constant_resolution(
       opacity_path, resolution=1000.0, wavelength_min=0.3, wavelength_max=100.0)

   grid = baerli.SpectralGrid.constant_wavenumber_step(
       opacity_path, step=1.0, wavelength_min=0.3, wavelength_max=100.0)      # step in cm^-1

   grid = baerli.SpectralGrid.constant_wavelength_step(
       opacity_path, step=0.01, wavelength_min=0.3, wavelength_max=100.0)     # step in micron

   grid = baerli.SpectralGrid.covering(
       opacity_path, wavelength_min=0.15, wavelength_max=100.0,
       temperature_min=100.0, temperature_max=350.0, nb_temperatures=30,
       nb_points=10000, stellar_temperature=5772.0, nb_points_stellar=10000)

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Factory
     - Point selection
   * - ``constant_resolution``
     - Constant resolving power :math:`R = \lambda/\Delta\lambda`: from each selected point the
       next native point with :math:`\lambda \le \lambda_{\mathrm{last}}(1 - 1/R)`.
   * - ``constant_wavenumber_step``
     - The next native point with :math:`\tilde\nu \ge \tilde\nu_{\mathrm{last}} + \Delta\tilde\nu`.
   * - ``constant_wavelength_step``
     - The next native point with :math:`\lambda \le \lambda_{\mathrm{last}} - \Delta\lambda`.
   * - ``covering``
     - The composite-Planck covering distribution (below).

The two window edges are always included. The optional ``wavenumber_file`` argument names a file
with the native wavenumber list (first entry the number of points, then the wavenumbers); without
it the HELIOS-k grid of 0 to 90000 cm\ :sup:`-1` in 0.01 cm\ :sup:`-1` steps is assumed.

The covering grid
-----------------

For an atmosphere with a wide temperature range, or an irradiated one, a grid of constant
resolving power spends many points where neither the thermal nor the stellar radiation carries
energy. The ``covering`` factory places the points where they matter, after Helling & Jørgensen
(1998): the point density is proportional to the *covering curve*

.. math::

   c(\tilde\nu) = \max_k \frac{B_{\tilde\nu}(T_k)}{T_k^4}

over a set of temperatures :math:`T_k` spanning ``temperature_min`` to ``temperature_max``
(``nb_temperatures`` of them; 0 gives steps of about 500 K). Walking along the native grid, a
point is selected whenever the accumulated measure :math:`\int c\, \dd\tilde\nu` since the last
selected point exceeds a threshold, and the threshold is found by bisection so that about
``nb_points`` points result. With ``stellar_temperature`` and ``nb_points_stellar`` a second set
is computed with the star's Planck function alone and the two sets are merged, so the final count
is at most ``nb_points + nb_points_stellar``. Set ``temperature_max`` to about the temperature of
the deepest layer.

.. list-table::
   :header-rows: 1
   :widths: 26 14 60

   * - Parameter
     - Default
     - Meaning
   * - ``wavelength_min``, ``wavelength_max``
     - required
     - Window in micron.
   * - ``temperature_min``, ``temperature_max``
     - required
     - Bounds of the thermal covering temperatures in K.
   * - ``nb_points``
     - required
     - Target number of thermal points.
   * - ``nb_temperatures``
     - ``0``
     - Number of covering temperatures; 0 selects steps of about 500 K.
   * - ``stellar_temperature``
     - ``0``
     - K; 0 omits the stellar set.
   * - ``nb_points_stellar``
     - ``0``
     - Target number of stellar points; 0 omits the stellar set.

Choosing the resolution
-----------------------

* Self-luminous objects: :math:`R = 1000` over 0.3 to 100 micron (about 6000 points) is the
  working resolution of the brown dwarf examples; :math:`R = 300` over 0.5 to 30 micron is
  enough for cloud experiments.
* Irradiated gas planets: :math:`R = 10000` in the examples, or a covering grid with a stellar
  set.
* Terrestrial planets: a covering grid with 10000 thermal and 10000 stellar points over 0.15 to
  100 micron, or :math:`R \ge 5000` for a constant-resolution grid. Below about :math:`R = 1000`
  the water-vapour window regions are not resolved and the Newton can stall on a near-null
  mode of the flux Jacobian.
* The band-closure correction of the ``ratio_ul`` solver recovers most of the line-core opacity
  that any sampled grid misses, so the converged temperature profile is far less sensitive to the
  resolution than the sampled opacity itself. Very dense grids (more than about 30000 points)
  are not recommended for the Newton solvers.

Attributes
----------

``grid.wavenumber_list`` (cm\ :sup:`-1`, ascending), ``grid.wavelength_list`` (micron,
descending), ``grid.nb_spectral_points`` and ``grid.opacity_path``. The spectrum and the spectral
radiation field of a model are defined on this list, in this order. ``grid.config`` records the
factory and its arguments for provenance, and :meth:`SpectralGrid.from_config` rebuilds a grid
from such a dict.
