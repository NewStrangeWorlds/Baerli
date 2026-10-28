Kinetics coupling with neoVULCAN
================================

Bärli's own chemistry modules give the equilibrium composition and the quench approximation of
the most important disequilibrium species. For full chemical kinetics and photochemistry the
model is coupled to `neoVULCAN <https://github.com/NewStrangeWorlds/neoVULCAN>`_, a JAX-based
reimplementation of the VULCAN kinetics code. The coupling is a Python driver,
:class:`baerli.KineticsCoupling`, that alternates the two codes until both agree. The theory is
in :ref:`theory-kinetics-coupling`; this page is the practical guide.

How the coupling works
----------------------

Both codes solve a steady state: Bärli the radiative-convective equilibrium for a frozen
composition, the kinetics code the chemical steady state for a frozen atmosphere. The driver
alternates them (a Gauss-Seidel iteration), with the kinetics code on its own pressure grid:

.. code-block:: text

   Bärli.compute(tol)  ->  T(p), Kzz(p)  ->  kinetics to steady state  ->  mixing ratios
        ^                                                                       |
        +------- external chemistry module  <--  Aitken relaxation + dead band  <-+

* **Pass 0.** Bärli converges with its seed composition, normally the quench approximation
  (``chemistry = [equilibrium, quench, external]``). The kinetics code is initialised with
  Bärli's temperature, Kzz and composition interpolated onto its grid.
* **Every further pass.** The kinetics code runs to its steady state (warm-started from the
  previous pass). The returned mixing ratios are interpolated back onto Bärli's grid and
  compared with the composition last committed. If the change exceeds the **dead band** it is
  committed with **Aitken relaxation** on the logarithmic mixing ratios through
  :meth:`set_composition`; Bärli then re-converges with ``warm_start=True`` at a tolerance that
  tightens by a decade per pass.
* **Convergence** is declared when the raw composition change is below
  ``composition_tolerance`` at the final RCE tolerance and the RCE itself converged.

Relaxation and dead band exist because the undamped alternation oscillates and because every
committed composition change costs the Newton solver its quadratic convergence tail.

Setting up a coupled run
------------------------

.. code-block:: python

   import baerli
   from baerli import KineticsCoupling, load_vulcan, mixing_ratios

   model = baerli.BrownDwarf(
       grid, effective_temperature=1000.0, surface_gravity=10**4.5,
       opacity_species=opacity_species,
       chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters.dat")),
                  ("quench", dict(metallicity=1.0)),   # first-pass guess
                  ("external", {})],                   # receives neoVULCAN's composition
       convection="mlt_dry",
       solver=("ratio_ul", dict(max_iterations=200, convergence_threshold=1e-5)))
   model.initialize(("milne", dict(kappa_ross=1e-2)))

   chem, neovulcan_path = load_vulcan("configs/neovulcan_brown_dwarf.cfg")

   coupling = KineticsCoupling(
       model, chem,
       rce_tolerance=(1e-3, 1e-5),
       dead_band=0.05,
       composition_tolerance=0.05,
       workdir=neovulcan_path)

   converged = coupling.run()
   symbols, table = mixing_ratios(model, ["CO", "CH4", "H2O", "NH3", "N2", "HCN", "CO2"])

The chemistry list must end in ``("external", {})``. The ``kzz`` option of the model provides
the eddy diffusion profile handed to the kinetics code: the self-consistent ``mlt`` profile
for self-luminous objects, typically a ``power_law`` for irradiated planets whose radiative
envelope is mixed by the circulation rather than by convection.

Loading neoVULCAN
~~~~~~~~~~~~~~~~~

.. py:function:: baerli.load_vulcan(config_path, base_dir=None, cfg_overrides=None, regenerate_chemistry=True)
   :noindex:

   Create and initialise a neoVULCAN ``VulcanChemistry`` instance. ``config_path`` is
   neoVULCAN's TOML configuration; relative paths inside it are resolved against ``base_dir``,
   which defaults to the checkout fetched by the build (:func:`baerli.neovulcan_dir`). The
   initialisation runs from ``base_dir`` because neoVULCAN reads its thermodynamic data
   relative to it. Returns ``(chem, base_dir)``.

.. py:function:: baerli.neovulcan_dir()
   :noindex:

   The neoVULCAN checkout compiled into the build; raises if the build was configured with
   ``BAERLI_FETCH_NEOVULCAN=OFF``.

neoVULCAN reads its data relative to its base directory on every rate update, while Bärli's
FastChem parameter file lists paths relative to the Bärli directory, so the two cannot share a
working directory: pass the neoVULCAN directory as ``workdir`` and the driver switches
directories around every call into the kinetics code.

KineticsCoupling parameters
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 16 58

   * - Parameter
     - Default
     - Meaning
   * - ``model``
     - required
     - An initialised Bärli model whose chemistry list ends in ``external``.
   * - ``chem``
     - required
     - The initialised kinetics code instance.
   * - ``species``
     - all common
     - Species exchanged between the codes; default every name both codes know.
   * - ``rce_tolerance``
     - ``(1e-2, 1e-5)``
     - RCE convergence threshold of the first pass and of the final one; tightened by a decade
       per pass in between.
   * - ``dead_band``
     - ``0.05``
     - Composition changes below this (dex, over the tracked species and levels) are not
       committed. The kinetics code's own convergence noise is about 0.03 dex; with
       photochemistry on, use 0.1.
   * - ``composition_tolerance``
     - ``0.05``
     - The coupling is converged when the raw composition change falls below this (dex) at the
       final RCE tolerance.
   * - ``relax``
     - ``0.5``
     - Initial Aitken weight of the composition update, then adapted and clamped to
       :math:`[0.05, 1]`.
   * - ``tracked_species``
     - opacity species
     - Species whose change decides relaxation, dead band and convergence. Default: the
       exchanged species that are also opacity sources of the model, falling back to all
       exchanged species. Keeps trace radicals swinging across the floor from driving the Aitken
       weight.
   * - ``tracked_floor``
     - ``1e-10``
     - Mixing ratios below this in both old and new composition are ignored when measuring
       changes.
   * - ``max_passes``
     - ``12``
     -
   * - ``seed_continuation``
     - ``False``
     - Run pass 0 in two stages: first with every ``quench`` module disabled (pure
       equilibrium), then with it enabled from the converged profile. Rescues strongly
       irradiated planets started from a cold analytic profile, where the quench composition is
       so far from the answer that the Newton never recovers; it *hurts* self-luminous objects.
       No effect without a quench module.
   * - ``workdir``
     - ``None``
     - Directory to switch to around every call into the kinetics code.
   * - ``verbose``
     - ``True``
     - Print one line per pass, prefixed ``[coupling]``.

:meth:`run` returns ``True`` on convergence. ``coupling.history`` holds one dict per pass with
``pass_``, ``change_dex``, ``committed``, ``weight``, ``rce_tolerance``, ``rce_converged`` and
the kinetics code's convergence info.

The kinetics interface
----------------------

The kinetics code is accessed through a small duck-typed interface that neoVULCAN's
``VulcanChemistry`` satisfies; any code offering the same methods can be coupled:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Member
     - Meaning
   * - ``species``
     - List of species names.
   * - ``data_atm.pco``
     - Pressure grid in dyn cm\ :sup:`-2`, bottom to top.
   * - ``set_atmosphere(T=, P=, Kzz=)``
     - Temperature and pressure on the levels, Kzz on the ``nz - 1`` interfaces (geometric
       means of adjacent level pressures).
   * - ``set_composition(species, table)``
     - Seed composition ``table[level][k]``.
   * - ``run_to_convergence(warm_start=)``
     - Run to steady state.
   * - ``get_mixing_ratios()``
     - ``(nz, n_species)`` array.
   * - ``get_convergence_info()``
     - Dict with at least ``count`` and ``end_case``.

Profiles are transferred by interpolation in :math:`\log_{10} p` with constant extension
beyond the range (:func:`baerli.interpolate_log`); mixing ratios and Kzz are interpolated in
their logarithm.

The neoVULCAN configuration
---------------------------

neoVULCAN keeps its own TOML configuration (``configs/neovulcan_brown_dwarf.cfg`` and
``configs/neovulcan_gas_planet.cfg`` in the repository, derived from neoVULCAN's Jupiter
example). Only the network, the grid and the switches matter for the coupling: temperature,
pressure, :math:`K_{zz}` and the initial composition are supplied by Bärli and overwrite the
placeholders in the file before the first integration. Relative paths are resolved against the
neoVULCAN base directory.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Key
     - Role in the coupling
   * - ``[network] network``, ``atom_list``
     - The reaction network and element set; they determine ``chem.species`` and hence the
       exchanged species (an NCHO low-temperature network for the brown dwarf, an SNCHO
       photochemical network for the gas planet).
   * - ``[atmosphere] nz``, ``P_b``, ``P_t``
     - neoVULCAN's own pressure grid in dyn cm\ :sup:`-2`; ``P_b`` should match Bärli's bottom
       pressure. The coupling interpolates between the grids.
   * - ``gs``, ``Rp``
     - Gravity (cm s\ :sup:`-2`) and radius (cm); ``gs`` must equal Bärli's ``surface_gravity``.
   * - ``use_Kzz = true``, ``use_moldiff = true``, ``use_vm_mol = true``
     - Eddy and molecular diffusion; ``use_vm_mol`` is required for the upper atmosphere and
       needs the ``Ros2`` solver.
   * - ``r_star``, ``orbit_radius``, ``sflux_file``
     - Stellar flux for photochemistry. neoVULCAN takes the file as the flux at the **stellar
       surface** and scales it by :math:`(r_\star R_\odot / a)^2`, so ``r_star`` and
       ``orbit_radius`` (AU) must reproduce Bärli's ``instellation``: with the Gueymard solar
       spectrum, ``r_star = 1`` and ``orbit_radius = 0.0742`` give the :math:`2.48 \times 10^8`
       erg cm\ :sup:`-2` s\ :sup:`-1` of the gas-planet example.
   * - ``sl_angle_deg``, ``f_diurnal``
     - Solar zenith angle and the diurnal factor on the photolysis rates; match Bärli's
       ``zenith_angle``.
   * - ``[photochemistry] use_photo``
     - Photochemistry on (gas planet) or off (brown dwarf).
   * - ``[boundary_conditions] diff_esc``
     - Diffusion-limited escape of H at the top for the irradiated planet.
   * - ``[solver] ode_solver = "Ros2"``, ``count_max``, ``yconv_cri``, ``slope_cri``,
       ``yconv_min``
     - Integrator and the steady-state criterion; the loose branch of the criterion is the
       origin of the 0.03 dex (thermal) to 0.1 dex (photochemical) noise floor that bounds
       ``dead_band`` from below.
   * - ``[paths] output_dir``, ``[output]``, ``[plotting]``
     - Not used through the coupling API; neoVULCAN writes no output files. The only file side
       effect is the regenerated chemistry module in the neoVULCAN checkout at
       :func:`baerli.load_vulcan` time.

Practical notes
---------------

* Start irradiated planets from a Guillot profile. From an isothermal start the first pass is
  a long damped-Newton crawl far from the root.
* With photochemistry neoVULCAN stops at a residual change of several per cent, so ``dead_band``
  and ``composition_tolerance`` must stay above about 0.1 dex.
* For an irradiated planet it can be better to leave the ``quench`` module out entirely and let
  neoVULCAN handle the disequilibrium: with the module in the chain the temperature correction
  of every pass crawled towards its iteration cap in the hot-Jupiter test, without it no pass
  needed more than about twenty iterations.
* Each pass's Bärli call uses ``warm_start=True`` and therefore skips the mixing-length ramp;
  pass 0 (and pass 0b of the seed continuation) start cold.
