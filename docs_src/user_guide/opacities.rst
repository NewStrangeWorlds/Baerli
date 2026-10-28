Opacity sources
===============

The ``opacity_species`` argument lists the opacity sources of a model as ``(species symbol, data
folder)`` pairs. The symbol is one of Bärli's chemical species (:func:`baerli.species_symbols`),
whose number density from the chemistry multiplies the cross section; the folder is the
cross-section data directory relative to ``opacity_path`` (by default the path the spectral grid
was built from). Two reserved folder names select built-in sources. The physics is in
:doc:`../theory/opacities`.

.. code-block:: python

   opacity_species = [
       ("CIA-H2-H2", "CIA/H2-H2"),                          # tabulated CIA
       ("CIA-H2-He", "CIA/H2-He"),
       ("H2",        "Rayleigh"),                           # built-in Rayleigh scattering
       ("He",        "Rayleigh"),
       ("H2O",       "Molecules/1H2-16O__POKAZATEL_e2b"),   # tabulated line cross sections
       ("CH4",       "Molecules/12C-1H4__YT34to10_e2b"),
       ("K",         "Alkali_Allard/K"),
       ("H-",        "none"),                               # built-in continuum
   ]

Kinds of opacity species
------------------------

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Symbol
     - Folder
     - Source
   * - any species symbol (``H2O``, ``CH4``, ``TiO``, ``K``, ...)
     - a data folder
     - Tabulated line cross sections of that species (HELIOS-k format, below). The same symbol
       may appear several times with different folders, e.g. an infrared and a visible line
       list of O\ :sub:`3`.
   * - ``CIA-H2-H2``, ``CIA-H2-He``, ``CIA-H-He``
     - a data folder
     - Tabulated collision-induced absorption; the coefficient is multiplied by the number
       densities of both partners.
   * - ``H2``, ``He``, ``H``, ``CO``, ``CO2``, ``CH4``, ``H2O``, ``N2``, ``O2``
     - ``"Rayleigh"``
     - Built-in Rayleigh scattering cross sections.
   * - ``H-``
     - ignored (``"none"``)
     - Built-in H\ :sup:`-` bound-free and free-free continuum (John 1988); needs H\ :sup:`-`,
       H and electrons from the equilibrium chemistry.
   * - ``CO2-CIA``
     - ignored (``"none"``)
     - Built-in CO\ :sub:`2`-CO\ :sub:`2` collision-induced absorption (Gruszka & Borysow;
       Baranov et al. 2004), 0 to 500 and 1000 to 2000 cm\ :sup:`-1`.
   * - ``H2O-CIA``
     - ignored (``"none"``)
     - Built-in MT_CKD water-vapour self and foreign continuum, 0 to 20000 cm\ :sup:`-1`.

A symbol that is not in the species table, or a folder without data, stops the model at
construction with a list of the species and whether their data were found.

Data format
-----------

The tabulated cross sections follow the layout of `HELIOS-k <https://github.com/exoclime/HELIOS-K>`_
output, one directory per species:

.. code-block:: text

   <opacity_path>/Molecules/1H2-16O__POKAZATEL_e2b/
       filelist.dat
       Out_00000_42000_00050_n800.bin
       Out_00000_42000_00050_n766.bin
       ...

``filelist.dat`` starts with the species mass in amu (``18.0105646837``), optionally followed by
the token ``log`` if the binary files store :math:`\log_{10}` of the cross section (``0 log`` for
the CIA tables), and then one line per :math:`(p, T)` point: ``pressure [bar]  temperature [K]
file name``. Each binary file is a headerless array of little-endian 32-bit floats, one value per
point of the **native wavenumber grid**, in cm\ :sup:`2` g\ :sup:`-1` (converted to cm\ :sup:`2`
per molecule with the mass) or, for CIA with mass 0, in cm\ :sup:`5` as stored.

The native grid is the HELIOS-k grid from 0 to 90000 cm\ :sup:`-1` in steps of 0.01 cm\ :sup:`-1`
unless a ``wavenumber_file`` is given to the spectral grid factories (:doc:`spectral_grid`). The
cross-section files must be laid out on exactly that grid; a shorter file is padded with
negligible values. Outside the tabulated pressure and temperature range the nearest tabulated
point is used (no extrapolation), so the tables should cover the model's :math:`(p, T)` domain.

The repository does not include opacity data. The examples use a HELIOS-k tree with ExoMol and
HITRAN line lists under ``Molecules/``, the HITRAN CIA tables under ``CIA/`` and the Allard alkali
profiles under ``Alkali_Allard/``.

Built-in data
-------------

The CO\ :sub:`2` CIA table (``src/transport_coeff/data/CO2_dimer_data.bin``) and the MT_CKD
continuum table (``src/transport_coeff/data/mt_ckd_h2o_continuum.dat``) are part of the
repository and are loaded by paths relative to the current working directory, so models that use
``CO2-CIA`` or ``H2O-CIA`` must be run from the repository root.

Typical species lists
---------------------

* **Brown dwarfs and gas planets** (H\ :sub:`2`/He atmospheres): ``CIA-H2-H2``, ``CIA-H2-He``,
  the alkalis ``K`` and ``Na``, and the molecules ``H2O``, ``CH4``, ``NH3``, ``CO``, ``CO2``,
  ``H2S``; for hot atmospheres ``TiO``, ``VO``, ``FeH``, ``H-`` and the Rayleigh scattering of
  ``H2`` and ``He`` for the optical.
* **Terrestrial planets**: the Rayleigh scattering of the bulk gases (``N2``, ``O2`` or ``CO2``),
  ``H2O``, ``CO2``, ``CH4``, ``O3``, ``N2O`` line lists, and the ``CO2-CIA`` and ``H2O-CIA``
  continua, which carry a large part of the greenhouse effect of a dense CO\ :sub:`2` or a humid
  atmosphere.

Cloud optical properties are configured separately (:doc:`clouds`).
