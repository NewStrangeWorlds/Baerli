Solvers and initial profiles
============================

The ``solver`` spec selects the temperature-correction scheme and its convergence settings; the
``profile`` argument of :meth:`initialize` selects the starting profile. The numerics are
described in :doc:`../theory/temperature_iteration`.

Solver types
------------

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Type
     - Scheme
     - Convection
   * - ``ratio_ul`` (default)
     - Newton on the collocated ratio residual with the full Unsöld-Lucy correction, NLEQ-ERR
       damping. Recommended for every object class: the only scheme whose converged root is
       free of a grid-scale checkerboard.
     - ``mlt_dry``, ``mlt_moist``, ``dry``, ``moist``, ``none``
   * - ``flux_divergence``
     - Newton on clima's heat-capacity-weighted flux-difference residual (MINPACK ``hybrj``).
       Conserves flux exactly per level at the cost of a sawtooth in the root.
     - ``dry``, ``moist``, ``none``
   * - ``ptc``
     - The legacy linearised flux Newton with pseudo-transient continuation.
     - ``dry``, ``moist``, ``none``
   * - ``time_stepping``
     - Explicit relaxation on the radiative heating rate with the local radiative timescale.
     - ``dry``, ``moist``, ``none``
   * - ``time_stepping_lre``
     - As above, blended with a local-radiative-equilibrium update.
     - ``dry``, ``moist``, ``none``
   * - ``helios``
     - HELIOS' per-level adaptive pseudo-time relaxation.
     - ``dry``, ``moist``, ``none``

Short aliases are accepted: ``rce``, ``flux``, ``lin``, ``ts``, ``ts_lre``, ``hel``.

Parameters
----------

Common to all solvers:

.. list-table::
   :header-rows: 1
   :widths: 26 14 60

   * - Parameter
     - Default
     - Meaning
   * - ``max_iterations``
     - ``100``
     - Outer iterations, i.e. radiative transfer solves with recomputed opacities.
   * - ``convergence_threshold``
     - ``1e-4``
     - Threshold on the scheme's residual (dimensionless, see below).

Relaxation schemes (``time_stepping``, ``time_stepping_lre``; ``max_change`` also for ``ptc``):

.. list-table::
   :header-rows: 1
   :widths: 26 14 60

   * - Parameter
     - Default
     - Meaning
   * - ``gamma``
     - ``0.5``
     - Relaxation factor on the local radiative timescale and on the bottom anchor step.
   * - ``ng_interval``
     - ``10``
     - Ng acceleration every *n* iterations; ``0`` disables it.
   * - ``max_change``
     - ``0.1``
     - Largest relative temperature change per level and iteration; ``0`` disables the clip.
   * - ``lre_fraction``
     - ``0.5``
     - ``time_stepping_lre`` only: weight of the local-radiative-equilibrium update.

``helios``:

.. list-table::
   :header-rows: 1
   :widths: 26 14 60

   * - Parameter
     - Default
     - Meaning
   * - ``ng_interval``
     - ``0``
     - Ng acceleration (off by default).
   * - ``step_init``
     - ``10``
     - Initial per-level temperature step in K.
   * - ``step_grow``
     - ``1.1``
     - Factor applied to a level's step when it moves monotonically.
   * - ``step_shrink``
     - ``1.5``
     - Divisor applied when a level oscillates.
   * - ``adapt_interval``
     - ``20``
     - Iterations between step adaptations (at least 2).
   * - ``step_exponent``
     - ``0.1``
     - Exponent :math:`q` in :math:`\Delta T = s\, \mathrm{sgn}(r)\, |r|^q`.
   * - ``max_step``
     - ``500``
     - Cap on :math:`|\Delta T|` per level and iteration in K.
   * - ``stencil``
     - ``"backward"``
     - Flux-difference stencil: ``"backward"`` (recommended), ``"centered"``, ``"forward"``.
   * - ``residual``
     - ``"flux"``
     - ``"flux"`` (flux difference) or ``"heating"`` (local radiative heating of the slab).

The Newton schemes ``ratio_ul`` and ``flux_divergence`` take no scheme parameters beyond the
common two; their step control is internal. Parameters of one scheme passed to another are
rejected at construction.

The convergence residual
------------------------

The iteration converges when the scheme's residual is below ``convergence_threshold``. What that
residual is depends on the scheme:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Solver
     - Residual tested
   * - ``ratio_ul``
     - :math:`\max(\max_i |g_i|,\ \max|\Delta T/T|\text{ of the inner solve})`: the ratio rows
       are relative local-radiative-equilibrium errors, the zone and surface rows flux errors
       relative to the flux scale. Set to 1 after a large move of the convective mask and floored
       at 0.5 while the mixing-length ramp runs.
   * - ``flux_divergence``
     - :math:`\max(\text{flux-conservation error}/F_{\mathrm{norm}},\ \max|\Delta T/T|)`.
   * - ``ptc``
     - Self-luminous: the largest :math:`|F_i - F_\star|/F_{\mathrm{scale}}` outside the
       optically thin skin; terrestrial with a convective zone: :math:`\max|\Delta T/T|`.
   * - ``helios``
     - :math:`\max |F_i - F_\star|/F_{\mathrm{scale}}` over the non-convective levels.
   * - ``time_stepping``, ``time_stepping_lre``
     - :math:`\max_i |\Delta T_i/T_i|` of the outer iteration.

The flux-conservation error of ``ratio_ul`` has a floor of about :math:`10^{-3}` for a
multi-stream, sampled-opacity radiative transfer, so a threshold of :math:`10^{-5}` on the
printed residual does not translate into :math:`10^{-5}` flux conservation. The TOA net-flux
error printed at convergence for self-luminous objects is the honest energy-balance error bar
and is about :math:`10^{-4}` of the transported flux for the adding-doubling backend, invariant
to the stream count and grid resolution.

Reading the iteration log
-------------------------

Each outer iteration prints one line. Brown dwarf:

.. code-block:: text

   iter     max|dT/T|    max|dT|   @ lv    dF/F(TOA)  flux_divergence     T_bot     T_top        N_conv
       1    1.2345e-01     45.678     12   8.7654e-02     1.2345e+03      2456       312      31

Gas planet (``dF/F(s+i)`` instead of ``dF/F(TOA)``) and terrestrial planet (``F_net(TOA)``,
``conv_resid`` and ``T_surf``) differ only in a few columns:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Column
     - Meaning
   * - ``max|dT/T|``
     - Largest relative temperature change of this iteration.
   * - ``max|dT|``, ``@ lv``
     - The signed change at the level where :math:`|\Delta T|` is largest, and that level
       (0 = bottom).
   * - ``dF/F(TOA)``, ``dF/F(s+i)``, ``conv_resid``
     - The scheme's convergence residual when it provides one (all Newton schemes and ``helios``);
       otherwise the TOA flux error relative to the flux scale (self-luminous) or ``-1``
       (terrestrial). This is the number tested against the threshold.
   * - ``flux_divergence``
     - Brown dwarfs: :math:`\dd F/\dd p` at the level of the largest change.
   * - ``F_net(TOA)``
     - Terrestrial: the TOA net flux relative to the incident flux.
   * - ``T_bot``, ``T_top``, ``T_surf``
     - Bottom-level, top-level and (terrestrial) surface temperatures in K.
   * - ``N_conv``
     - Number of levels flagged convective.
   * - ``[Ng]``
     - Appended when the Ng acceleration fired.

Additional lines report the cloud ramp (``cloud ramp: scale ...``), the Ackerman & Marley
modules (``A&M <condensate>: ...``), a stall (``Stalled: ...``) and the final state (``Converged
after N iterations.`` with the TOA energy-balance diagnostic, or ``Warning: did not converge``).

.. _initial-profiles:

Initial profiles
----------------

:meth:`initialize(profile, chemistry=None)` sets temperature and composition. The profile spec
types and their parameters:

.. list-table::
   :header-rows: 1
   :widths: 14 40 46

   * - Type
     - Parameters
     - Use
   * - ``adiabat``
     - ``surface_temperature``, ``stratosphere_temperature`` (K, both required)
     - Terrestrial planets. The adiabat of the active convection scheme from the surface,
       reduced by 2 % and floored at the stratosphere temperature; lets the mixing-length
       corrector skip its ramp.
   * - ``milne``
     - ``kappa_ross`` (cm\ :sup:`2` g\ :sup:`-1`, required), ``effective_temperature``
       (default: the object's :math:`\Teff` or :math:`\Tint`; required for terrestrial planets)
     - Brown dwarfs. The grey Milne-Eddington solution.
   * - ``const``
     - ``temperature`` (K, required)
     - Isothermal start. Slow for irradiated gas planets.
   * - ``guillot``
     - ``kappa_ir`` (cm\ :sup:`2` g\ :sup:`-1`), ``t_irr`` (K), ``gamma`` (:math:`\kappa_{\mathrm{vis}}/\kappa_{\mathrm{ir}}`), all required;
       ``t_int`` (default: the object's :math:`\Tint`); ``mode = "isotropic"`` (default) with
       ``f = 0.25``, or ``mode = "beam"`` with ``mu`` (default: the object's zenith angle)
     - Irradiated gas planets. Guillot (2010) analytic profile; :math:`T_{\mathrm{irr}} =
       (S/\sigma)^{1/4}`.

The optional ``chemistry`` argument lists chemistry module specs used for the **initial**
composition only, e.g. ``[("isoprofile", dict(H2=0.85, He=0.15, H2O=1e-4))]`` to start an
equilibrium model from prescribed abundances. Restarts use
:meth:`initialize_from_file(filename)` (:doc:`output`).

Recommendations
---------------

* **Brown dwarfs**: ``ratio_ul`` with ``mlt_dry`` from a ``milne`` profile. The mixing-length
  ramp costs about fifteen iterations on a cold start; a converged model typically takes 20 to
  40 iterations.
* **Gas planets**: ``ratio_ul`` with ``mlt_dry`` from a ``guillot`` profile with
  :math:`T_{\mathrm{irr}}` from the instellation. The ramp is skipped automatically when the
  internal flux is below 1 % of the transported flux and the start is on the stable side.
* **Terrestrial planets**: ``ratio_ul`` with ``mlt_moist`` (or ``mlt_dry`` for a dry
  atmosphere) from an ``adiabat`` profile. Use about 200 levels or more; with 100 levels the
  radiative-convective boundary is under-resolved and shows up as a kink in the flux profile.
* **Thresholds**: :math:`10^{-4}` (default) to :math:`10^{-5}`. For a coupling loop start loose
  and tighten per pass.
* **Restarts and continuation**: after a small change of the problem (composition, a parameter)
  call :meth:`compute(warm_start=True)`; it skips the ramp and applies the full cloud at once.
* **Relaxation schemes** are kept for comparison and for cases where a Newton scheme fails; they
  need many more iterations (hundreds) and the convective adjustment schemes.

Troubleshooting
---------------

The iteration stalls at a floor above the threshold
   The flux-conservation floor of the radiative transfer is above the requested threshold.
   Relax the threshold to just above the floor, increase the spectral resolution (the floor
   decreases with resolution up to about :math:`R = 5000`), or increase the stream count.

The residual sits at 0.5 for many iterations
   The mixing-length ramp is running (cold start of a self-luminous object). This is expected
   and lasts about twelve to fifteen iterations.

``convection mlt_dry/mlt_moist requires solver ratio_ul``
   The mixing-length schemes only work with the ``ratio_ul`` corrector; switch to ``dry`` or
   ``moist`` for the other solvers.

The deep atmosphere of an irradiated planet crawls by a few kelvin per iteration
   The start is far from the root. Use a ``guillot`` profile, or restart from a converged
   neighbouring model with ``warm_start=True``.

A cloud ramp step is undone
   The step lifted the residual above ``backoff``; the ramp continues with a smaller growth
   factor. Harmless, but a sign that the deck is thick; see :doc:`clouds`.

Very high spectral resolution diverges
   With more than about 30000 spectral points the line-core heating aloft and the minimum
   damping of the Newton can interact badly; use the ``covering`` grid or a lower resolution
   with the band-closure correction, which recovers most of the missing line-core opacity.
