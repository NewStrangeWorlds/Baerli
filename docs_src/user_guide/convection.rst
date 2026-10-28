Convection and eddy diffusion
=============================

Two module families control vertical energy transport and vertical mixing: ``convection``
selects the convection scheme, ``kzz`` the eddy diffusion profile that the quench chemistry,
the microphysical clouds and the kinetics coupling consume. The physics is described in
:doc:`../theory/convection`.

Convection schemes
------------------

.. list-table::
   :header-rows: 1
   :widths: 16 44 40

   * - Type
     - Model
     - Works with
   * - ``mlt_dry`` (default)
     - Mixing-length convection with the dry adiabatic neutrality gradient. The convective flux
       is a smooth function of the super-adiabaticity and is solved for inside the Newton
       corrector; no hard adjustment.
     - ``ratio_ul`` solver only
   * - ``mlt_moist``
     - As ``mlt_dry`` with the saturated (water) pseudo-adiabat as neutrality gradient.
     - ``ratio_ul`` solver only
   * - ``dry``
     - Hard dry-adiabatic convective adjustment after each temperature step (relaxation
       solvers), or adiabatic slaving of the convective zone (Newton solvers).
     - all solvers
   * - ``moist``
     - Hard moist-adiabatic adjustment with the water pseudo-adiabat.
     - all solvers
   * - ``none``
     - No convection: a purely radiative atmosphere.
     - all solvers

Parameters of ``mlt_dry`` and ``mlt_moist``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 16 64

   * - Parameter
     - Default
     - Meaning
   * - ``alpha``
     - ``1.0``
     - Mixing length in pressure scale heights, :math:`\lambda = \alpha H_p` (away from the
       surface of a terrestrial planet, where the Blackadar wall law takes over). The
       converged temperature profile is insensitive to it; the ``mlt`` eddy diffusion profile
       is not.
   * - ``min_pressure``
     - object default
     - Pressure ceiling in bar: no convective flux above it. Defaults to :math:`10^{-4}` bar
       for brown dwarfs and gas planets and :math:`10^{-3}` bar for terrestrial planets. For
       ``mlt_moist`` the ceiling must exclude the subsaturated stratosphere (the Earth example
       uses :math:`10^{-2}` bar), because the pseudo-adiabat is evaluated for a saturated
       parcel.

Parameters of ``dry`` and ``moist``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 16 64

   * - Parameter
     - Default
     - Meaning
   * - ``min_pressure``
     - object default
     - Pressure ceiling in bar; the adjustment sweeps stop at the first level above it.
   * - ``max_sweeps``
     - ``10``
     - Maximum number of bottom-to-top adjustment sweeps per iteration; a sweep without any
       adjustment ends the loop early.

Choosing a scheme
~~~~~~~~~~~~~~~~~

* Use the default ``mlt_dry`` with ``ratio_ul`` for brown dwarfs and gas planets. It converges
  at least as fast as the adjustment schemes everywhere and gives kink-free
  radiative-convective boundaries.
* Use ``mlt_moist`` for terrestrial planets with a condensable water inventory, with a
  ``min_pressure`` in the lower stratosphere.
* Use ``dry``/``moist`` with the ``flux_divergence``, ``ptc``, ``helios`` and time-stepping
  solvers. The mixing-length schemes are rejected at construction time with these solvers,
  because the flux law is consumed inside the ``ratio_ul`` corrector and the other schemes
  rely on the hard adjustment.
* ``none`` is useful as a test bed for the radiative solvers.

Eddy diffusion (``kzz``)
------------------------

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Type
     - Profile
   * - ``mlt`` (default)
     - Self-consistent profile from the mixing-length flux law in convective layers, extended
       into the radiative layers by a prescription and relaxed between iterations.
   * - ``constant``
     - The same value at every level.
   * - ``power_law``
     - :math:`K_{zz} = \mathrm{value}\,(p/\mathrm{pressure})^{\mathrm{slope}}`.

Parameters of ``mlt``
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 18 16 66

   * - Parameter
     - Default
     - Meaning
   * - ``scaling``
     - ``"velocity"``
     - ``"velocity"``: :math:`K_{zz} = \lambda\, w_{\mathrm{conv}}` (AGNI's default).
       ``"flux"``: the Ackerman & Marley / Charnay et al. scaling
       :math:`K_{zz} = \frac{H_p}{3} (\lambda/H_p)^{4/3} (R_{\mathrm{spec}} F_{\mathrm{conv}} /
       \rho c_p)^{1/3}` with the convective flux of the flux law. The velocity scaling gives
       :math:`10^9` to :math:`10^{10}` cm\ :sup:`2` s\ :sup:`-1` in the convective zone of a
       brown dwarf and, continued upward, 100 to 500 micron cloud particles; the flux scaling
       is the usual choice for cloudy runs.
   * - ``radiative``
     - ``"constant"``
     - Extension above each convective zone, anchored at the zone's maximum within one
       pressure decade below its top: ``"constant"`` continues that value; ``"power_law"``
       scales it as :math:`(p/p_{\mathrm{rcb}})^{\mathrm{slope}}` and requires ``slope``
       (negative grows upward, for irradiated planets; positive decays to ``min``, for
       self-luminous objects); ``"fixed"`` uses ``value`` above every zone.
   * - ``slope``
     - required for ``power_law``
     - Exponent of the radiative power law.
   * - ``value``
     - required for ``fixed``
     - cm\ :sup:`2` s\ :sup:`-1`, must be positive.
   * - ``min``
     - ``1e4``
     - Floor in cm\ :sup:`2` s\ :sup:`-1`; also the value when nothing is convective.
   * - ``relax``
     - ``0.5``
     - Initial Aitken weight of the relaxation between iterations, in :math:`(0, 1]`.
   * - ``tolerance``
     - ``0.05``
     - Dead band in dex: a raw profile that differs from the committed one by less than this
       everywhere is not committed. ``0`` disables the band.

Parameters of ``constant`` and ``power_law``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 18 16 66

   * - Parameter
     - Default
     - Meaning
   * - ``value``
     - required
     - cm\ :sup:`2` s\ :sup:`-1`, positive. For ``power_law`` the value at the reference
       pressure.
   * - ``pressure``
     - ``1.0``
     - ``power_law`` only: reference pressure in bar.
   * - ``slope``
     - required
     - ``power_law`` only: the exponent.

Examples
~~~~~~~~

.. code-block:: python

   # cloudy brown dwarf: flux-scaled MLT Kzz continued above the convective zone
   kzz=("mlt", dict(scaling="flux", radiative="constant"))

   # self-luminous object with a decaying radiative-zone Kzz
   kzz=("mlt", dict(radiative="power_law", slope=2.0, min=1e5))

   # irradiated planet whose envelope is mixed by the circulation (Parmentier et al. 2013)
   kzz=("power_law", dict(value=5e8, pressure=1.0, slope=-0.5))

Notes
~~~~~

* With cloud modules, avoid ``radiative = "fixed"`` with a low value: the cloud decks of a
  brown dwarf sit where the cloud heating creates detached convective patches, and a step of
  several decades in :math:`K_{zz}` between a patch and the fixed value flips the particle size
  by a large factor whenever a patch appears or vanishes, so the coupled iteration cycles.
* :math:`K_{zz}` is updated once per outer iteration from the committed profile and is frozen
  inside the Newton's trial evaluations.
* The ``mlt`` relaxation state is kept across :meth:`compute` calls of the same model.
* The resulting profile is stored in ``model.atmosphere.kzz`` and in the output file.
