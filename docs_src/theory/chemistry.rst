Chemistry
=========

The composition of the atmosphere is produced by a list of **chemistry modules** applied in order
at every outer iteration; each module may overwrite the species it is responsible for and leaves
the others as they are. This chapter describes the models behind the six module types; the
parameters are in :doc:`../user_guide/chemistry`.

Composition representation
--------------------------

The composition is stored as number densities :math:`n_{i,s}` (cm\ :sup:`-3`) on every level
:math:`i` for the species :math:`s` of Bärli's species table (:func:`baerli.species_symbols`);
column 0 is the total number density from the ideal gas law,

.. math::

   n_{\mathrm{tot}} = \frac{10^6\, p\,[\mathrm{bar}]}{k_B T} .

Volume mixing ratios are :math:`f_s = n_s / n_{\mathrm{tot}}`. The species table is the only
composition vocabulary of the model: it maps a Bärli symbol (``H2O``) to a FastChem symbol
(``H2O1``) and a molecular weight, and the same symbols name the opacity species, the keys of the
``isoprofile`` module, the columns of a ``fixed`` file and the species handed to
:meth:`set_composition`. The 49 species are H, He, C, O, Fe, Fe\ :sup:`+`, Ca, Ti,
Ti\ :sup:`+`, H\ :sub:`2`, H\ :sub:`2`\ O, CO\ :sub:`2`, CO, CH\ :sub:`4`, HCN, NH\ :sub:`3`,
C\ :sub:`2`\ H\ :sub:`2`, N\ :sub:`2`, Na, K, H\ :sub:`2`\ S, H\ :sup:`-`, TiO, VO, FeH, SH,
MgO, AlO, CaO, CrH, MgH, CaH, TiH, OH, e\ :sup:`-`, V, V\ :sup:`+`, Mn, Si, Cr, Cr\ :sup:`+`,
SiO, SiO\ :sub:`2`, SO\ :sub:`2`, CS\ :sub:`2`, O\ :sub:`2`, Ar, N\ :sub:`2`\ O and O\ :sub:`3`.

The **mean molecular weight** is a separate per-level array. The prescribed modules compute it
from the Bärli species, :math:`\bar\mu = \sum_s f_s m_s`; the equilibrium module takes FastChem's
own value, which includes every species of its network (about 500) and not only the 49 Bärli
species. The post-processing modules (``quench``, ``external``) therefore update it
**differentially**, :math:`\bar\mu \leftarrow \bar\mu + \sum_{s\ \mathrm{touched}} (f_s^{\mathrm{new}}
- f_s^{\mathrm{old}})\, m_s`, so that the heavy trace species outside the Bärli list keep their
contribution.

Equilibrium chemistry (``equilibrium``)
---------------------------------------

Thermochemical equilibrium is computed with `FastChem <https://github.com/NewStrangeWorlds/FastChem>`_
(Stock et al. 2018, 2022; Kitzmann et al. 2024), fetched by the build and driven through its
library interface. FastChem solves the mass-action laws of about 500 gas species and 185
condensates with the equilibrium constants

.. math::

   \log K(T) = \frac{a_1}{T} + a_2 \ln T + a_3 + a_4 T + a_5 T^2

of its data files (``fastchem_data/logK.dat``, ``logK_condensates.dat``, fitted to the JANAF
tables and others) for the element abundances of ``element_abundances_solar.dat`` (Asplund et al.
2009, astronomical scale with :math:`\varepsilon_{\mathrm{H}} = 12`). The parameter file names
these data files and the solver tolerances. The temperature is floored at 100 K in every FastChem
call.

Element abundances
~~~~~~~~~~~~~~~~~~

Two parameters scale the solar abundances before each call: every element except H and He is
multiplied by the **metallicity** factor :math:`m`, and the **C/O ratio** is then imposed by
setting the oxygen abundance from carbon,

.. math::

   \varepsilon_X \to m\, \varepsilon_X \ (X \ne \mathrm{H, He}), \qquad
   \varepsilon_{\mathrm{O}} = \frac{\varepsilon_{\mathrm{C}}}{\mathrm{C/O}} .

Carbon keeps its scaled solar value. This is applied unconditionally: the default ``c_to_o =
0.5`` already differs from the solar value 0.55 of the shipped abundance file and lowers the
oxygen abundance by 10 %; pass ``c_to_o=0.55`` for a strictly solar composition.

Condensation
~~~~~~~~~~~~

Three modes exist. ``none`` is the gas phase only. ``equilibrium`` lets FastChem condense at
every level independently: the condensate stays where it forms, the :math:`f_{\mathrm{sed}} \to
0` limit of a cloud. ``rainout`` is FastChem's rainout approximation (Kitzmann et al. 2024, Sect.
3.6): the levels are solved from the bottom upward, the elements locked in condensates at a level
are removed from the gas phase handed to the next level, and the condensate density stored at a
level is the amount that condenses *at* that level from the vapour arriving from below, the
increment the Ackerman & Marley cloud module integrates (:doc:`clouds`). The per-element
cumulative degree of condensation :math:`\delta_j(i)` is recorded as well. Both condensing modes
need the condensate data file in the FastChem parameter file (``fastchem_parameters_condensates.dat``).

Operator splitting
~~~~~~~~~~~~~~~~~~

The rainout march is sequential over the levels and cannot be evaluated level by level in parallel,
so it runs only on the **committed** composition call once per outer iteration. The trial
evaluations of the Newton corrector solve the gas phase level by level in parallel (one FastChem
instance per thread) with the element abundances depleted by the committed degrees of
condensation, :math:`\varepsilon_j (1 - \delta_j(i))`. The committed gas phase itself is taken
from the same depleted gas-only solve, so that committed and trial compositions agree by
construction: FastChem's own rainout gas phase differs from it by its element-conservation
tolerance, a per-cent effect on the remainder of a 99 %-condensed element, enough to floor the
flux residual at :math:`10^{-3}` in the alkali depletion zone of a 1300 K brown dwarf. The split
is exact at convergence, where :math:`\delta` no longer changes.

The outer map :math:`\delta \to T \to \delta` contracts slowly (a factor 0.96 per iteration was
measured for the alkali depletion at the top of a brown dwarf). With
``condensation_acceleration`` (default on), every third committed call applies a Steffensen
extrapolation to :math:`x = \ln(1 - \delta)`, the gas remainder the opacity sees, wherever the
last two steps show a slow geometric tail (ratio between 0.5 and 0.95).

Prescribed compositions
-----------------------

``isoprofile``
   Constant volume mixing ratios at every level for the listed species. Sodium, if not listed, is
   slaved to potassium with the solar ratio :math:`n_{\mathrm{Na}} = 16.2181\, n_{\mathrm{K}}`.

``fixed``
   A table of mixing ratios over pressure (bottom to top), interpolated in :math:`\log_{10} p`:
   log-linearly where both bracketing values are positive, linearly otherwise, and clamped to the
   table's range.

``manabe_wetherald``
   The relative-humidity profile of Manabe & Wetherald (1967, their Eq. 2),

   .. math::

      \mathrm{RH}(p) = \max\!\left[0,\; \mathrm{RH}_0\, \frac{p/p_s - 0.02}{1 - 0.02}\right],

   which reaches zero at :math:`p/p_s = 0.02` (a dry stratosphere, the parameterisation's own
   cold trap), with :math:`p_s` the bottom level. The water mixing ratio is :math:`f_{\mathrm{H_2O}}
   = \mathrm{RH}\, e_s(T)/p` with the Murphy & Koop (2005) saturation pressure
   (:doc:`atmosphere`), clipped to 0.99, and subject to a **cold trap**: walking upward, the
   mixing ratio is the running minimum of the saturation value, so a warm (ozone-heated)
   stratosphere cannot re-inflate the water. Without the cold trap the fixed-humidity water-vapour
   feedback fills the stratosphere and drives an unphysical runaway. Only water is written; the
   module is meant to follow a dry-composition module such as ``fixed``.

Quench chemistry (``quench``)
-----------------------------

The quenching approximation of Zahnle & Marley (2014) for disequilibrium chemistry in
H\ :sub:`2`-dominated atmospheres: a species family stays in equilibrium where the chemical
timescale is shorter than the vertical mixing timescale and is frozen at its equilibrium value at
the quench point above it. The module is a post-processing step on the equilibrium composition and
must be listed after ``equilibrium``.

Timescales
~~~~~~~~~~

The mixing timescale uses the scale height as mixing length,

.. math::

   t_{\mathrm{mix}} = \frac{H^2}{K_{zz}}, \qquad H = \frac{k_B T}{\bar\mu\, m_u\, g},

with the model's eddy diffusion profile (:ref:`theory-kzz`) and the surface gravity. The chemical
timescales are the fits of Zahnle & Marley (2014), with :math:`p` in bar, :math:`T` in K and
:math:`m` the metallicity relative to solar:

.. math::

   t_{q1} = 1.5 \times 10^{-6}\; p^{-1} m^{-0.7} e^{42000/T}\ \mathrm{s}, \qquad
   t_{q2} = 40\; p^{-2} e^{25000/T}\ \mathrm{s}, \qquad
   t_{\mathrm{CO}} = \left(t_{q1}^{-1} + t_{q2}^{-1}\right)^{-1}

(their Eqs. 12 to 14, evaluated in log space);

.. math::

   t_{\mathrm{NH_3}} = 1.0 \times 10^{-7}\; p^{-1} e^{52000/T}\ \mathrm{s}, \qquad
   t_{\mathrm{HCN}} = 1.5 \times 10^{-4}\; p^{-1} m^{-0.7} e^{36000/T}\ \mathrm{s}, \qquad
   t_{\mathrm{CO_2}} = 1.0 \times 10^{-10}\; p^{-0.5} e^{38000/T}\ \mathrm{s}

(Eqs. 32, 40 and 44). The families are {CO, CH\ :sub:`4`, H\ :sub:`2`\ O}, {NH\ :sub:`3`,
N\ :sub:`2`}, {HCN} and {CO\ :sub:`2`}; water is frozen together with CO and CH\ :sub:`4`
because oxygen conservation ties it to their quench.

Quench point and frozen profile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With :math:`D(i) = \log_{10} t_{\mathrm{chem}}(T_i, p_i) - \log_{10} t_{\mathrm{mix}, i}`, the
quench pressure is the root of :math:`D` interpolated linearly in :math:`\log p` between the last
equilibrium level and the first quenched level,

.. math::

   \log_{10} p_q = \log_{10} p_{i-1} + \left(\log_{10} p_i - \log_{10} p_{i-1}\right)
   \frac{D_{i-1}}{D_{i-1} - D_i} .

The continuous placement matters because the module is re-evaluated at every trial temperature of
the Newton: a quench point jumping between levels would make the composition and the flux residual
discontinuous in :math:`T`. The frozen value is the equilibrium mixing ratio interpolated
log-linearly at :math:`p_q`, applied to every level above. A family that is already quenched at
the bottom of the grid is frozen at the bottom level with a warning to extend the grid deeper.

CO\ :sub:`2` is treated last and follows the quenched CO, H\ :sub:`2`\ O and H\ :sub:`2` in
quasi-equilibrium below its own quench point (Zahnle & Marley Sect. 6, Eq. 43),

.. math::

   f_{\mathrm{CO_2}} = f^{\mathrm{eq}}_{\mathrm{CO_2}}\,
   \frac{f_{\mathrm{CO}}\, f_{\mathrm{H_2O}}}{f^{\mathrm{eq}}_{\mathrm{CO}}\, f^{\mathrm{eq}}_{\mathrm{H_2O}}},

with the equilibrium constant implied by the equilibrium composition itself rather than taken from
a fit, so that it reduces exactly to the equilibrium CO\ :sub:`2` wherever CO and water are
unquenched.

Whatever the quenched species gain or lose is taken from or given to H\ :sub:`2`, as in the
picaso and Exo-REM implementations, and the mean molecular weight is updated differentially.

Relaxation
~~~~~~~~~~

The outer loop "solve :math:`T` for a frozen composition, then move the quench points" is a
Picard iteration whose gain exceeds unity at the photospheric CO/CH\ :sub:`4` transition of a
warm Jupiter (measured: the CO quench pressure flipping between 2 and 11 bar and one level
swinging by :math:`\pm 184` K every outer iteration). The quench pressures are therefore relaxed
between committed evaluations with a scalar Aitken weight per family on :math:`\ln p_q`, starting
at 0.5 and clamped to :math:`[0.05, 1]`, with a dead band of 0.02 dex. In lagged (trial)
evaluations the committed quench levels and frozen values are re-applied on top of the trial
equilibrium composition instead of being re-derived.

The fits are valid for H\ :sub:`2`-dominated atmospheres of roughly solar to a few times solar
metallicity without photochemistry; the module refuses to run without H\ :sub:`2`.

External composition (``external``)
-----------------------------------

The ``external`` module holds mixing-ratio profiles for a subset of species supplied by the driver
through :meth:`set_composition`, typically from a kinetics code, and overwrites those species in
the composition of the preceding modules. The supplied species are **rescaled per level to the
budget** the preceding modules gave the same set of species,

.. math::

   f_k \leftarrow f_k^{\mathrm{ext}}\, \frac{\sum_k f_k^{\mathrm{prev}}}{\sum_k f_k^{\mathrm{ext}}},

so that the total stays conserved and species the external code does not know (TiO, VO, Na, K,
...) keep their share: the kinetics code normalises to one over its own network, FastChem's trace
species outside the network keep theirs. The module is a no-op until a composition is set, and it
must be listed last.

.. _theory-kinetics-coupling:

Coupling to chemical kinetics
-----------------------------

Bärli and a kinetics code such as neoVULCAN both solve a steady state: the radiative-convective
equilibrium for a frozen composition on one side, the chemical steady state for a frozen
atmosphere on the other. :class:`baerli.KineticsCoupling` alternates them in a Gauss-Seidel
fashion (:doc:`../user_guide/kinetics_coupling`). Temperature and :math:`K_{zz}` are
interpolated onto the kinetics grid in :math:`\log p` (mixing ratios and :math:`K_{zz}` in their
logarithm); :math:`K_{zz}` is given on the interfaces :math:`\sqrt{p_i p_{i+1}}` of that grid.

The composition returned by the kinetics code is committed to the ``external`` module with
**Aitken relaxation** on the logarithmic mixing ratios of the tracked species: with the residual
:math:`r = \ln y_{\mathrm{raw}} - \ln y_{\mathrm{old}}` and its previous value,

.. math::

   \omega \leftarrow \mathrm{clamp}\!\left(-\omega\,
   \frac{\langle r^{\mathrm{prev}}, r - r^{\mathrm{prev}}\rangle}{\langle r - r^{\mathrm{prev}}, r - r^{\mathrm{prev}}\rangle},
   \; 0.05, \; 1\right), \qquad
   \ln y_{\mathrm{new}} = \ln y_{\mathrm{old}} + \omega\, r ,

and only if the largest tracked change exceeds a **dead band** (0.05 dex by default, above the
kinetics code's own convergence noise). The undamped alternation oscillates, and every committed
composition change costs the Newton its quadratic convergence tail. The RCE tolerance is loose on
the first passes and tightens by a decade per pass; the coupling is converged when the raw
composition change falls below the tolerance at the final RCE tolerance.

For a strongly irradiated planet started from a cold analytic profile the quench composition can
be so far from the converged one that the Newton never recovers; the optional **seed
continuation** then converges pass 0 with pure equilibrium chemistry first and switches the quench
module on afterwards (HELIOS forces equilibrium chemistry on its first iteration for the same
reason). The same device hurts a self-luminous object, where switching the quench module on at a
converged profile is an order-of-magnitude opacity change.

Limitations
-----------

* The chemistry is lagged in the temperature Newton: trial evaluations recompute the gas phase at
  the trial temperature, but the Jacobian is Planck-only.
* Module order is not enforced beyond the quench module's check for H\ :sub:`2`.
* The species table spells three FastChem names differently from the shipped data file
  (Fe\ :sup:`+`, Ti\ :sup:`+`, FeH), so the equilibrium module leaves these three species at
  zero.
* The quench timescales use the surface gravity at all levels and are order-of-magnitude fits.
