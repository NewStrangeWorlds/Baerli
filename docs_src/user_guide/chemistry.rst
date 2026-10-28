Chemistry
=========

The ``chemistry`` argument is a **list** of chemistry module specs that are applied in order at
every outer iteration; at least one module is required. Each module overwrites the species it is
responsible for and leaves the rest untouched, which is what makes the list composable: an
equilibrium module followed by the quench approximation followed by an external composition, or a
tabulated dry composition followed by a humidity profile. The models are described in
:doc:`../theory/chemistry`.

.. code-block:: python

   chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters.dat",
                                   metallicity=1.0, c_to_o=0.55))]

   chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters_condensates.dat",
                                   condensation="rainout"))]

   chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters.dat")),
              ("quench", dict(metallicity=1.0)),
              ("external", {})]

   chemistry=[("fixed", dict(file="data/Earth/earth_standard_composition.dat")),
              ("manabe_wetherald", dict(surface_rh=0.77))]

   chemistry=[("isoprofile", dict(CO2=0.965, N2=0.035))]

Module types
------------

.. list-table::
   :header-rows: 1
   :widths: 22 56 22

   * - Type
     - Model
     - Writes
   * - ``equilibrium``
     - Thermochemical equilibrium with FastChem, optionally with condensation.
     - every species FastChem knows, total density, mean molecular weight
   * - ``isoprofile``
     - Constant mixing ratios.
     - the listed species (plus Na slaved to K), total density
   * - ``fixed``
     - Mixing ratios tabulated over pressure.
     - the file's columns, total density
   * - ``manabe_wetherald``
     - Manabe & Wetherald (1967) relative-humidity water profile with a cold trap.
     - H\ :sub:`2`\ O only
   * - ``quench``
     - Zahnle & Marley (2014) quenching of CO/CH\ :sub:`4`/H\ :sub:`2`\ O, NH\ :sub:`3`/N\
       :sub:`2`, HCN and CO\ :sub:`2` with the model's :math:`K_{zz}`.
     - those seven species and H\ :sub:`2`
   * - ``external``
     - Composition handed in through :meth:`set_composition`, e.g. from a kinetics code.
     - the supplied species

Short aliases: ``eq``, ``iso``, ``fix``, ``mw``, ``q``, ``ext``.

``equilibrium``
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Parameter
     - Default
     - Meaning
   * - ``parameter_file``
     - required
     - FastChem parameter file (``fastchem_parameters.dat`` for the gas phase,
       ``fastchem_parameters_condensates.dat`` with the condensate data). Paths inside it are
       relative to the working directory.
   * - ``metallicity``
     - ``1.0``
     - Factor applied to the solar abundance of every element except H and He.
   * - ``c_to_o``
     - ``0.5``
     - Carbon-to-oxygen ratio, imposed by setting the oxygen abundance from carbon. The solar
       value with the shipped abundances is 0.55.
   * - ``condensation``
     - ``"none"``
     - ``"none"``, ``"equilibrium"`` (condensates stay where they form) or ``"rainout"``
       (FastChem's rainout: elements condensed at a level are removed from the gas above it; the
       per-level condensation increments feed the Ackerman & Marley cloud modules). Both
       condensing modes need the condensate data file in the parameter file.
   * - ``condensation_acceleration``
     - ``true``
     - Steffensen acceleration of the lagged depletion of the gas phase.

The condensed phase is available as ``model.atmosphere.condensate_symbols``,
``condensate_weights`` and ``condensate_densities`` and in the output file.

``isoprofile``
~~~~~~~~~~~~~~

Every key is a species symbol and its value the volume mixing ratio (``dict(H2=0.85, He=0.15,
H2O=1e-4)``); at least one species is required and unknown symbols are rejected. Sodium, if not
listed, is set to the solar ratio to potassium, :math:`16.2` times the potassium mixing ratio
(zero if potassium is not listed either).

``fixed``
~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Parameter
     - Default
     - Meaning
   * - ``file``
     - required
     - Composition table.

The file starts with a header line ``#p_bar  N2  O2  ...`` naming the pressure column and the
species columns (unknown species are skipped with a notice), followed by one row per pressure
with the volume mixing ratios, from **high to low pressure**, at least two rows, tab or space
separated; further ``#`` lines are comments. Mixing ratios are interpolated in :math:`\log p`
(log-linearly where positive) and clamped to the table's range. ``data/Earth/earth_standard_composition.dat``
is a modern-Earth example.

``manabe_wetherald``
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Parameter
     - Default
     - Meaning
   * - ``surface_rh``
     - ``0.77``
     - Relative humidity at the surface. The humidity falls linearly to zero at 2 % of the
       surface pressure, and the water mixing ratio cannot increase with height (cold trap).

List it after the dry-composition module whose water it replaces.

``quench``
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Parameter
     - Default
     - Meaning
   * - ``metallicity``
     - ``1.0``
     - Metallicity relative to solar in the timescale fits; should match the ``equilibrium``
       module's value.
   * - ``relax``
     - ``true``
     - Aitken relaxation of the quench pressures between outer iterations (dead band 0.02 dex).

The module needs an H\ :sub:`2`-dominated equilibrium composition from a preceding module, the
model's ``kzz`` profile and the surface gravity; it prints the quench pressures of the four
families whenever one moves by more than 10 %. It is valid for roughly solar to a few times solar
metallicity without photochemistry. For a kinetics coupling it provides the first-pass guess.

``external``
~~~~~~~~~~~~

No parameters. The module does nothing until :meth:`set_composition(species, mixing_ratios)`
hands it a table on the model's grid; the supplied species are then rescaled per level to the
budget the preceding modules gave the same species and overwrite them (:doc:`kinetics_coupling`).
List it last.

Enabling and disabling modules
------------------------------

:meth:`set_chemistry_enabled(type, enabled)` switches every module of a type off or on without
rebuilding the model; a disabled module passes the composition through. The kinetics coupling
uses this to converge with pure equilibrium chemistry before switching the quench module on.

Initial composition
-------------------

:meth:`initialize(profile, chemistry=...)` accepts a separate list of module specs for the
**initial** composition, e.g. prescribed mixing ratios to start an equilibrium model from a
composition far from equilibrium; the model's own modules take over at the first iteration.

Inspecting the composition
--------------------------

.. code-block:: python

   symbols, table = baerli.mixing_ratios(model, ["CO", "CH4", "H2O", "NH3"])
   nd = np.array(model.atmosphere.number_densities)      # [level][species], cm^-3
   mmw = np.array(model.atmosphere.mean_molecular_weight)
   trial = model.recompute_chemistry(lagged=True)         # the Newton's trial-path composition
