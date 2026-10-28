Example scripts
===============

The repository root contains complete, runnable examples. They all expect the opacity data
under ``/media/data/opacity_data/helios-k/``; change ``opacity_path`` to your installation.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Script
     - What it shows
   * - ``run_brown_dwarf.py``
     - A 1000 K, log *g* = 4.5 brown dwarf with equilibrium chemistry, DISORT with four
       streams, mixing-length convection and the ``ratio_ul`` solver from a Milne start.
   * - ``run_brown_dwarf_clouds.py``
     - The same object at 1300 K, log *g* = 5, run clear and with Ackerman & Marley forsterite
       and iron clouds (``fsed = 2``) on FastChem's rainout condensation. Prints the temperature
       and band-flux differences and the per-module cloud diagnostics. ``--cloudy-only`` reuses
       the clear output. The docstring of the script records the Kzz choices that make the
       coupled cloud iteration converge.
   * - ``run_gas.py``
     - A warm Jupiter (:math:`\Tint = 100` K, irradiation :math:`2.5 \times 10^8` erg
       cm\ :sup:`-2` s\ :sup:`-1`) with TiO/VO opacities, adding-doubling radiative transfer,
       dry convective adjustment, started isothermal with a prescribed initial composition.
   * - ``run_earth.py``
     - Modern Earth: fixed composition from a file, optional Manabe-Wetherald humidity, CO\
       :sub:`2` and H\ :sub:`2`\ O continua, Rayleigh scattering, a composite-Planck covering
       grid, moist mixing-length convection, an adiabatic start.
   * - ``run_venus.py``
     - A 92 bar CO\ :sub:`2`/N\ :sub:`2` atmosphere with the ``flux_divergence`` solver and dry
       convective adjustment.
   * - ``run_brown_dwarf_neovulcan.py``
     - The brown dwarf coupled to neoVULCAN's thermal kinetics
       (:doc:`kinetics_coupling`).
   * - ``run_gas_neovulcan.py``
     - The warm Jupiter coupled to neoVULCAN with photochemistry, a Guillot start, a power-law
       Kzz and the seed continuation.
   * - ``configs/earth.toml``
     - The Earth model as a config file: ``python -m baerli configs/earth.toml``.

Plotting the output
-------------------

The output files are plain NetCDF and open with xarray:

.. code-block:: python

   import matplotlib.pyplot as plt
   import baerli

   ds = baerli.load_model_data("output_brown_dwarf.nc")

   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

   ax1.semilogy(ds["temperature"], ds["pressure"])
   ax1.invert_yaxis()
   ax1.set_xlabel("T [K]"); ax1.set_ylabel("p [bar]")
   convective = ds["convective"].values.astype(bool)
   ax1.scatter(ds["temperature"][convective], ds["pressure"][convective], s=6, c="C3")

   ax2.loglog(ds["wavelength"], ds["spectrum"])                                 # per micron
   ax2.set_xlabel("wavelength [micron]"); ax2.set_ylabel("F [erg/cm^2/s/micron]")
   ax2.set_xlim(0.5, 30)

   plt.tight_layout()
   plt.show()

The ``spectrum`` variable is the top-of-atmosphere upward flux per micron; the spectral
radiation field (``flux``, ``flux_up``, ``mean_intensity`` on the model object) is per
wavenumber, and converts with :math:`F_\lambda = F_{\tilde\nu}\, \tilde\nu^2 / 10^4`.
