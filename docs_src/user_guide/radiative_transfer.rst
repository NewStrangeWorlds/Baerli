Radiative transfer, stellar spectra and surfaces
================================================

Three module families describe how radiation enters, traverses and leaves the atmosphere: the
``radiative_transfer`` solver, the ``stellar_spectrum`` of an irradiated object and the
``surface`` of a terrestrial planet. The physics is in :doc:`../theory/radiative_transfer`.

Radiative transfer solvers
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - Type
     - Parameter
     - Meaning
   * - ``disort`` (default)
     - ``nb_streams = 4``
     - DISORT++ discrete-ordinate solver. ``nb_streams`` is the **total** number of streams over
       both hemispheres; allowed values 4, 6, 8, 10, 12, 14, 16, 32, 64.
   * - ``adding_doubling``
     - ``nb_streams = 2``
     - AddingDoublingRT matrix-operator solver. ``nb_streams`` is the number of Gauss-Legendre
       quadrature points **per hemisphere** (at least 2; 2, 4, 8, 16 and 32 have specialised
       kernels).

Short alias: ``ad`` for ``adding_doubling``. Both solvers provide the analytic temperature
Jacobians the Newton solvers need, both handle thermal emission, a Lambertian surface and the
diffusion lower boundary of a semi-infinite atmosphere, and both are parallel over the spectral
points.

.. code-block:: python

   radiative_transfer=("disort", dict(nb_streams=4))            # brown dwarf reference
   radiative_transfer=("adding_doubling", dict(nb_streams=2))   # fast choice for planets

Choosing between them: with two streams per hemisphere the adding-doubling solver is the faster
choice for terrestrial and gas-planet models at high spectral resolution; DISORT with four streams
is the reference for self-luminous objects. The two agree to the accuracy of their angular
quadrature, and the TOA flux error reported at convergence (about :math:`10^{-4}` of the
transported flux for the adding-doubling solver) is the thin-layer accuracy floor of the backend.

Stellar spectra
---------------

Gas planets and terrestrial planets take a ``stellar_spectrum`` spec together with the scalar
``instellation`` and ``zenith_angle`` constructor arguments. The module provides the **shape** of
the spectrum; its normalisation and the incidence geometry are properties of the planet.

.. list-table::
   :header-rows: 1
   :widths: 18 22 60

   * - Type
     - Parameter
     - Meaning
   * - ``blackbody``
     - ``temperature`` (K, required)
     - A Planck spectrum :math:`\pi B_{\tilde\nu}(T_\star)`.
   * - ``tabulated``
     - ``file`` (required)
     - A tabulated spectrum, two whitespace-separated columns ``wavelength [micron]`` and
       ``F_lambda`` in any unit; lines starting with ``#`` are comments; at least two points.
       Converted to per wavenumber, interpolated linearly onto the model grid, zero outside the
       file's range.

In both cases the spectrum is rescaled so that its trapezoidal integral over the **model grid**
equals ``instellation`` (erg cm\ :sup:`-2` s\ :sup:`-1`); the absolute units of a file are
irrelevant. A grid that truncates part of the stellar spectrum therefore redistributes the
missing energy into the covered range; choose ``wavelength_min`` small enough (0.15 to 0.3
micron for a solar-type star).

The repository ships ``data/stellar_spectra/spectrum_sun.dat`` (Gueymard 2004) and
``data/stellar_spectra/spectrum_m4_5v.dat`` (AD Leo, M3.5e).

Irradiation geometry
~~~~~~~~~~~~~~~~~~~~

``instellation`` is the bolometric flux at the planet's orbit measured normal to the beam,
:math:`S = \sigma T_\star^4 (R_\star / d)^2`, and ``zenith_angle`` is :math:`\mu_0 =
\cos\theta_0`. The direct flux entering the top of the atmosphere is :math:`\mu_0 S`. Bärli does
not compute the instellation from a stellar radius and an orbital distance; the user provides it,
together with any redistribution factor. The common global-average approximation uses half the
full instellation with :math:`\mu_0 = 0.5`, as in the Earth and Venus examples. The
redistribution factor ``f`` of the Guillot initial profile is a separate quantity.

Surfaces
--------

Terrestrial planets take a ``surface`` spec. The surface is Lambertian with a spectral albedo
:math:`A_{\tilde\nu}` and emissivity :math:`1 - A_{\tilde\nu}`.

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Type
     - Parameters
     - Albedo
   * - ``blackbody``
     - none
     - 0 everywhere.
   * - ``simple``
     - ``albedo``, ``wavelength_switch`` (micron), both required
     - ``albedo`` for :math:`\lambda <` ``wavelength_switch``, 0 at longer wavelengths (a grey
       reflector in the stellar range, a black body in the thermal infrared).
   * - ``variable_albedo``
     - ``file`` (required)
     - Tabulated: two columns ``wavelength [micron]`` and ``reflectance``, ``#`` comments,
       at least two points, interpolated linearly in wavelength with the end values held
       outside the range.

The repository ships ``data/Earth/earth_spectral_surface_reflection.dat``, a global-average
surface reflectance of the Earth from 0.3 to 25 micron.

The surface temperature is a result of the model (``model.surface_temperature``). With the Newton
solvers it is the temperature of the bottom level, which carries the surface energy balance as
its flux equation; with the relaxation solvers it is stepped explicitly with a 50 cm slab-ocean
heat capacity after each temperature update (:doc:`../theory/temperature_iteration`).

.. code-block:: python

   model = baerli.TerrestrialPlanet(
       grid,
       surface_gravity=980.0,
       instellation=1361e3 * 0.5,
       zenith_angle=0.5,
       stellar_spectrum=("tabulated", dict(file="data/stellar_spectra/spectrum_sun.dat")),
       surface=("variable_albedo", dict(file="data/Earth/earth_spectral_surface_reflection.dat")),
       ...)
