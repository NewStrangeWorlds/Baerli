Clouds
======

The ``clouds`` option takes a list of cloud module specs whose optical properties are summed;
``cloud_ramp`` controls how the cloud opacity is switched on during the iteration. The physics is
in :doc:`../theory/clouds`.

.. code-block:: python

   # prescribed grey slab handed over at run time
   model = baerli.BrownDwarf(grid, ..., clouds=[("fixed", {})], cloud_ramp="none")
   model.set_cloud_optical_properties(**baerli.clouds.slab(
       2.0, top_pressure=0.01, bottom_pressure=1.0, single_scattering=0.3, asymmetry=0.6))

   # a Mie deck of 1 micron forsterite particles with tau = 1 at 1 micron
   clouds=[("mie_slab", dict(condensate="Mg2SiO4", particle_radius=1.0, sigma=1.5,
                             top_pressure=0.1, bottom_pressure=1.0, optical_depth=1.0))]

   # Ackerman & Marley forsterite and iron clouds on the rainout chemistry
   model = baerli.BrownDwarf(
       grid, ...,
       chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters_condensates.dat",
                                       condensation="rainout"))],
       clouds=[("ackerman_marley", dict(condensate="Mg2SiO4", fsed=2.0, sigma=2.0)),
               ("ackerman_marley", dict(condensate="Fe", fsed=2.0, sigma=2.0))],
       kzz=("mlt", dict(scaling="flux", radiative="constant")))

Module types
------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Type
     - Model
   * - ``fixed``
     - Prescribed layer optical depth, single-scattering albedo and asymmetry parameter over
       pressure and wavelength, from a file or handed over with
       :meth:`set_cloud_optical_properties`. Mapped conservatively onto the model layers.
   * - ``mie_slab``
     - A deck of Mie particles of one condensate between two pressures, uniform in
       :math:`\ln p`, lognormal size distribution; optics from LX-MIE.
   * - ``ackerman_marley``
     - Ackerman & Marley (2001) eddy-sedimentation cloud of one condensate on the rainout
       condensation of the ``equilibrium`` chemistry; particle sizes from the fall-velocity
       balance with the model's :math:`K_{zz}`.

Short aliases: ``fix``, ``mie``, ``am``.

``fixed``
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 14 66

   * - Parameter
     - Default
     - Meaning
   * - ``file``
     - none
     - Table file. Without it the module does nothing until a table arrives through
       :meth:`set_cloud_optical_properties`.

**File format.** Whitespace-separated columns ``pressure[bar] wavelength[micron] optical_depth
single_scattering asymmetry``, one row per (pressure, wavelength) pair, ``#`` starts a comment.
The rows must cover the full Cartesian grid of the distinct pressures and wavelengths; at least
two pressure levels are needed. Pressures are **levels** in any order; the row at :math:`p_k`
holds the optical depth of the input layer between :math:`p_k` and the next level towards lower
pressure (the optical depth of the lowest-pressure row is ignored), plus that layer's
:math:`\omega` and :math:`g`. Constraints: :math:`\tau \ge 0`, :math:`0 \le \omega \le 1`,
:math:`-1 \le g \le 1`. Everything is linear in wavelength and constant outside the table's
range; a single wavelength gives a grey cloud.

.. code-block:: text

   # pressure[bar]  wavelength[micron]  optical_depth  single_scattering  asymmetry
   1.00000000e+00  1.00000000e+00  2.00000000e+00  3.00000000e-01  6.00000000e-01
   1.00000000e-02  1.00000000e+00  0.00000000e+00  3.00000000e-01  6.00000000e-01

:func:`baerli.clouds.slab` builds the table of a deck with a uniform optical depth per unit
:math:`\ln p` between two pressures (grey, or per wavelength), and :func:`baerli.clouds.write_table`
writes any such table in the file format. :meth:`set_cloud_optical_properties` hands the table
to every file-less ``fixed`` module and takes effect at the next outer iteration.

``mie_slab``
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Parameter
     - Default
     - Meaning
   * - ``condensate``
     - required
     - A name from the condensate table below, or the path to a refractive-index file in the
       same format.
   * - ``particle_radius``
     - required
     - Geometric mean radius :math:`r_g` in micron.
   * - ``sigma``
     - ``1.0``
     - Geometric standard deviation of the lognormal distribution; 1 is monodisperse.
   * - ``top_pressure``, ``bottom_pressure``
     - required
     - Deck boundaries in bar, ``0 < top < bottom``.
   * - ``optical_depth``
     - none
     - Extinction optical depth of the deck at ``reference_wavelength``. Alternative to
       ``column_density``; exactly one of the two must be given.
   * - ``reference_wavelength``
     - ``1.0``
     - Micron.
   * - ``column_density``
     - none
     - Particle column in cm\ :sup:`-2`.
   * - ``mie_resolution``
     - ``200``
     - Resolving power of the coarse wavelength grid of the Mie table.
   * - ``mie_radii``
     - ``120``
     - Number of radii of the Mie table (:math:`10^{-4}` to 542 micron, log-spaced).

``ackerman_marley``
~~~~~~~~~~~~~~~~~~~

Requires ``chemistry = [("equilibrium", {parameter_file: "fastchem_parameters_condensates.dat",
condensation: "rainout"})]``; the model refuses to build otherwise.

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Parameter
     - Default
     - Meaning
   * - ``condensate``
     - required
     - A name from the condensate table that has a FastChem counterpart.
   * - ``fsed``
     - ``1.0``
     - Sedimentation efficiency.
   * - ``fsed_mode``
     - ``"const"``
     - ``"const"`` or ``"exp"`` (Rooney et al. 2022: decays from ``fsed`` at the cloud base to
       ``eps`` over ``b`` mixing lengths).
   * - ``b``
     - ``1.0``
     - Scale of the exponential ``fsed`` in mixing lengths.
   * - ``eps``
     - ``1e-2``
     - Asymptotic ``fsed`` of the exponential mode.
   * - ``sigma``
     - ``2.0``
     - Geometric standard deviation of the lognormal size distribution.
   * - ``kzz_min``
     - ``1e5``
     - Floor on the layer :math:`K_{zz}` in cm\ :sup:`2` s\ :sup:`-1`.
   * - ``mixing_length``
     - ``"scale_height"``
     - ``"scale_height"`` or ``"lapse_ratio"`` (scale height times
       :math:`\max(0.1, \nabla/\nabla_{\mathrm{ad}})`).
   * - ``relax``
     - ``0.5``
     - Log-space relaxation weight of the structure between outer iterations (1 = none).
   * - ``max_change``
     - ``0.3``
     - Cap in dex on the change of mass column and radius per update (0 = none).
   * - ``tolerance``
     - ``0.05``
     - Dead band in dex: smaller changes are not committed.
   * - ``gate``
     - ``1e-2``
     - The structure is recomputed only while the corrector's residual is below this.
   * - ``seed``
     - ``1e-2``
     - Fraction of its target a newly cloudy layer starts at.
   * - ``bulk_density``
     - table value
     - Condensate density in g cm\ :sup:`-3`; required for a refractive-index file path.
   * - ``molecule_diameter``
     - ``2.827e-8``
     - Gas molecule diameter in cm for the viscosity (H\ :sub:`2`).
   * - ``lj_temperature``
     - ``59.7``
     - Lennard-Jones well depth over :math:`k_B` in K (H\ :sub:`2`).
   * - ``mie_resolution``, ``mie_radii``
     - ``200``, ``120``
     - As for ``mie_slab``.

The module reports its state once per iteration (``A&M <condensate>: base ..., top ..., max q_c
..., r_g ..., tau(1um) ..., committed/held``) and exposes per-layer diagnostics through
:meth:`cloud_diagnostics` and the output file: ``condensate_mmr``, ``column_density``,
``mass_column``, ``r_g``, ``r_eff``, ``r_w``, ``alpha``, ``fsed``.

Condensates
-----------

Names accepted by ``condensate``, with the bulk density used by default and the FastChem symbol
through which the Ackerman & Marley module reads the condensation. Names without a FastChem
symbol can only be used by ``mie_slab``. The refractive indices are the compilation of Kitzmann
& Heng (2018) in ``data/refractive_indices/``; the source of each data set is given in the
file header.

.. list-table::
   :header-rows: 1
   :widths: 22 36 16 26

   * - Name
     - Refractive-index file
     - :math:`\rho_p` (g cm\ :sup:`-3`)
     - FastChem symbol
   * - ``Al2O3``
     - Al2O3.dat (amorphous)
     - 3.987
     - Al2O3(s,l)
   * - ``C``
     - C.dat (graphite)
     - 2.26
     - C(s)
   * - ``CaTiO3``
     - CaTiO3.dat
     - 3.987
     - CaTiO3(s)
   * - ``Cr``
     - Cr.dat
     - 7.15
     - Cr(s,l)
   * - ``Fe``
     - Fe.dat
     - 7.875
     - Fe(s,l)
   * - ``FeO``
     - FeO.dat
     - 5.74
     - FeO(s,l)
   * - ``FeS``
     - FeS.dat
     - 4.84
     - FeS(s,l)
   * - ``Fe2O3``
     - Fe2O3.dat
     - 5.24
     - Fe2O3(s)
   * - ``Fe2SiO4``
     - Fe2SiO4.dat
     - 4.39
     - Fe2SiO4(s)
   * - ``H2O_ice``
     - H2O_s.dat
     - 0.93
     - H2O(s,l)
   * - ``KCl``
     - KCl.dat
     - 1.99
     - KCl(s,l)
   * - ``Mg2SiO4``
     - Mg2SiO4_amorph_sol-gel.dat
     - 3.214
     - Mg2SiO4(s,l)
   * - ``MgAl2O4``
     - MgAl2O4.dat
     - 3.64
     - MgAl2O4(s,l)
   * - ``MgFeSiO4``
     - MgFeSiO4_amorph_glass.dat
     - 3.8
     - none
   * - ``MgO``
     - MgO.dat
     - 3.58
     - MgO(s,l)
   * - ``MgSiO3``
     - MgSiO3_amorph_glass.dat
     - 3.192
     - MgSiO3(s,l)
   * - ``MgSiO3_sol-gel``
     - MgSiO3_amorph_sol-gel.dat
     - 3.192
     - MgSiO3(s,l)
   * - ``Mg08Fe02SiO3``
     - Mg08Fe02SiO3_amorph_glass.dat
     - 3.3
     - none
   * - ``Mg05Fe05SiO3``
     - Mg05Fe05SiO3_amorph_glass.dat
     - 3.5
     - none
   * - ``Mg04Fe06SiO3``
     - Mg04Fe06SiO3_amorph_glass.dat
     - 3.6
     - none
   * - ``Mg08Fe12SiO4``
     - Mg08Fe12SiO4_amorph_glass.dat
     - 3.9
     - none
   * - ``MnS``
     - MnS.dat
     - 4.0
     - MnS(s)
   * - ``Na2S``
     - Na2S.dat
     - 1.856
     - Na2S(s,l)
   * - ``NaCl``
     - NaCl.dat
     - 2.165
     - NaCl(s,l)
   * - ``SiC``
     - SiC.dat
     - 3.21
     - SiC(s)
   * - ``SiO``
     - SiO.dat
     - 2.13
     - SiO(s)
   * - ``SiO2``
     - SiO2_amorph.dat
     - 2.648
     - SiO2(s,l)
   * - ``SiO2_alpha``
     - SiO2_alpha.dat
     - 2.648
     - SiO2(s,l)
   * - ``TiC``
     - TiC.dat
     - 4.93
     - TiC(s,l)
   * - ``TiO2``
     - TiO2_anatase.dat
     - 4.25
     - TiO2(s,l)
   * - ``tholin``
     - Titan_tholin.dat
     - 1.0
     - none
   * - ``ZnS``
     - ZnS.dat
     - 4.04
     - none

The cloud ramp
--------------

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Type
     - Behaviour
   * - ``default``
     - Nothing until the corrector's easy start is over, then the cloud enters with at most
       ``tau_start`` of optical depth at 1 micron and grows by ``factor`` on every iteration
       whose residual is below ``gate``; a step that lifts the residual above ``backoff`` is
       undone and the growth factor reduced to its square root. A warm start applies the full
       cloud at once.
   * - ``none``
     - The full cloud from the first iteration.

.. list-table::
   :header-rows: 1
   :widths: 20 14 66

   * - Parameter
     - Default
     - Meaning
   * - ``tau_start``
     - ``1.0``
     - Total cloud optical depth at 1 micron at onset.
   * - ``gate``
     - ``1e-2``
     - Residual below which the scale grows.
   * - ``factor``
     - ``2.0``
     - Growth factor per step (must be > 1).
   * - ``backoff``
     - ``2.0``
     - Residual above which the last step is undone.

The ramp scale is available as :meth:`cloud_ramp_scale`; :meth:`cloud_optics` includes it, the
module diagnostics do not.

Practical notes
---------------

* Cloudy brown dwarfs: use ``kzz = ("mlt", {scaling: "flux", radiative: "constant"})``. The
  default velocity scaling gives :math:`K_{zz}` of :math:`10^9` to :math:`10^{10}`
  cm\ :sup:`2` s\ :sup:`-1` above the convective zone and hence 100 to 500 micron particles.
  Do not use ``radiative = "fixed"`` with a low value under a deck; the particle size then flips
  by a large factor whenever a detached convective patch appears or vanishes and the iteration
  cycles.
* A converged cloudy run typically takes 100 to 200 iterations; the example in
  ``run_brown_dwarf_clouds.py`` converges in about 110 with the constant :math:`K_{zz}`
  extension (optical depths of 4 to 5 per deck at 1 micron).
* The cloud decks change the convective structure: cloud heating creates detached convective
  patches above the main zone, which the mixing-length scheme handles without any special
  treatment.
* Restarting from an output file restores temperature and composition but not the condensate
  increments; the Ackerman & Marley module keeps its committed structure until the first
  chemistry call of the new run.
* The tests in ``tests/`` document the expected behaviour of every module (conservative
  remapping, Mie cross sections against LX-MIE, the ``fsed`` scaling of the suspended
  condensate).
