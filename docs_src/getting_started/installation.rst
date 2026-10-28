Installation
============

Requirements
------------

* A C++17 compiler (GCC or Clang), CMake 3.14 or newer, OpenMP.
* Python 3.10 or newer with ``numpy`` and ``xarray`` (with the ``h5netcdf`` engine for the
  output files). ``tomli`` is needed for TOML config files on Python 3.10; ``pyyaml`` for YAML
  config files.
* Network access during the first CMake configure run: the dependencies are fetched with
  ``FetchContent`` (see below).
* Opacity data in the sampled cross-section format of `HELIOS-K
  <https://github.com/exoclime/HELIOS-K>`_ (:doc:`../user_guide/opacities`). The data are not
  part of the repository.
* For the kinetics coupling: `JAX <https://jax.readthedocs.io>`_ and the other Python
  dependencies of neoVULCAN.

Building
--------

.. code-block:: sh

   git clone <repository>
   cd <repository>
   mkdir -p build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   make -j

The build produces the static library ``baerli_lib`` and the Python extension module
``_baerli*.so``, which CMake places inside the ``baerli/`` package directory of the source
tree. Afterwards ``import baerli`` works from the repository root, or from anywhere with the
repository root on ``PYTHONPATH``:

.. code-block:: sh

   export PYTHONPATH=/path/to/repository:$PYTHONPATH
   python -c "import baerli; print(baerli.species_symbols()[:5])"

The compiler flag ``-march=native`` is set for the whole build, including the fetched
dependencies, because Bärli, DisORT++ and FastChem exchange Eigen types and must agree on their
SIMD alignment. Do not copy the compiled module to a machine with a different CPU architecture.

CMake options
-------------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Option
     - Default
     - Meaning
   * - ``CMAKE_BUILD_TYPE``
     - ``Release``
     - Use ``Release``; the radiative transfer is far too slow in a debug build.
   * - ``BUILD_PYTHON_BINDINGS``
     - ``ON``
     - Build the Python module (the only user interface). Also switches the whole build to
       position-independent code.
   * - ``BAERLI_FETCH_NEOVULCAN``
     - ``ON``
     - Fetch the neoVULCAN checkout at configure time and compile its location into the module
       (``baerli.neovulcan_dir()``). ``OFF`` skips it; the kinetics coupling then needs an
       explicit ``base_dir``.
   * - ``FETCHCONTENT_SOURCE_DIR_<NAME>``
     - unset
     - Use a local checkout instead of the pinned Git revision for a dependency
       (``DISORTPP``, ``ADDING_DOUBLING``, ``FASTCHEM``, ``NEOVULCAN``, ``LXMIE``, ``PYBIND11``),
       e.g. ``-DFETCHCONTENT_SOURCE_DIR_NEOVULCAN=/path/to/neoVULCAN``.

Dependencies fetched by CMake
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 48 30

   * - Dependency
     - Role
     - Pinned to
   * - DisORT++
     - discrete-ordinate radiative transfer (``radiative_transfer = "disort"``), the shared
       Eigen 3.4 dependency
     - branch ``main``
   * - AddingDoublingRT
     - adding-doubling radiative transfer (``radiative_transfer = "adding_doubling"``)
     - commit ``a397261``
   * - FastChem
     - equilibrium chemistry and condensation (``chemistry = "equilibrium"``); built as a
       library only
     - commit ``cf931e7``
   * - LX-MIE
     - Mie efficiencies for ``mie_slab`` and ``ackerman_marley`` clouds; only the library target
       in ``src/mie`` is built
     - commit ``38cc6d9``
   * - neoVULCAN
     - chemical kinetics for ``baerli.KineticsCoupling`` (Python/JAX, fetched only)
     - commit ``3bc981b``
   * - pybind11
     - Python bindings
     - ``v2.13.6``

The pinned revisions are in ``CMakeLists.txt``. The data the compiled core needs at run time
(refractive indices of the cloud condensates, the CO\ :sub:`2` dimer and MT_CKD water continuum
tables) are shipped in the repository and located through compile-time paths, so the source tree
must stay where it was built.

Data files in the repository
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Path
     - Content
   * - ``fastchem_parameters.dat``
     - FastChem parameter file: element abundances and the gas-phase equilibrium constants.
   * - ``fastchem_parameters_condensates.dat``
     - The same plus the condensate data file; required for ``condensation = "rainout"`` or
       ``"equilibrium"`` and hence for the microphysical cloud modules.
   * - ``fastchem_data/``
     - FastChem's element abundances, log *K* tables and condensate data.
   * - ``data/stellar_spectra/``
     - Tabulated stellar spectra (Sun after Gueymard 2004, AD Leo).
   * - ``data/Earth/``
     - A modern-Earth composition profile and the spectral surface reflectance used by the
       terrestrial examples.
   * - ``data/refractive_indices/``
     - Complex refractive indices of the condensates known to the cloud modules.
   * - ``configs/``
     - A TOML model configuration (``earth.toml``) and the neoVULCAN configurations of the
       coupled examples.
   * - ``run_*.py``
     - Complete example scripts (:doc:`../user_guide/examples`).

Running the tests
-----------------

The Python tests in ``tests/`` exercise the cloud modules against analytic expectations and
against the ``lxmie`` Python package. They build small models and therefore need the opacity data
root used by the examples (they are skipped without it):

.. code-block:: sh

   python -m pytest tests/

Building this documentation
---------------------------

The documentation sources live in ``docs_src`` of the documentation branch and build with
Sphinx and the furo theme:

.. code-block:: sh

   pip install -r docs_src/requirements.txt
   sphinx-build -b html docs_src docs_test/html
