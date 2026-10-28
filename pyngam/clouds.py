"""Helpers for the `clouds` model option.

A cloud module of type ``fixed`` takes a table of layer optical depth, single-scattering albedo and
asymmetry parameter over (pressure, wavelength), either from a file (``("fixed", {"file": ...})``)
or handed over at run time with ``model.set_cloud_optical_properties(**table)``. The functions here
build such tables.

Table convention (see src/cloud/fixed_cloud.h): ``pressure`` are LEVELS in bar, in any order; the
row at pressure ``p_k`` holds the optical depth of the input layer between ``p_k`` and the next
level towards lower pressure (the optical depth of the lowest-pressure row is ignored), plus that
layer's ``w`` and ``g``. The model maps the table conservatively onto its own layers (the
cumulative optical depth is interpolated in ln p), so the column optical depth is independent of
the model grid. Everything is linear in wavelength and constant outside the table's wavelength
range; a single-wavelength table is grey.
"""

import numpy as np


def slab(optical_depth, top_pressure, bottom_pressure, single_scattering=0.0, asymmetry=0.0,
         wavelength=None):
    """A cloud deck between two pressures with a uniform optical depth per unit ln p.

    Parameters
    ----------
    optical_depth : float or sequence
        Column optical depth of the deck; a sequence gives one value per ``wavelength``.
    top_pressure, bottom_pressure : float
        Cloud top and base in bar (``top_pressure < bottom_pressure``).
    single_scattering, asymmetry : float or sequence
        Single-scattering albedo (0 = pure absorber) and asymmetry parameter, scalar or per
        wavelength.
    wavelength : sequence, optional
        Wavelengths in micron of the per-wavelength values; omit for a grey deck.

    Returns the keyword arguments of ``model.set_cloud_optical_properties``.
    """
    if not 0 < top_pressure < bottom_pressure:
        raise ValueError("need 0 < top_pressure < bottom_pressure")

    if wavelength is None:
        wavelength = [1.0]
    wavelength = np.atleast_1d(np.asarray(wavelength, dtype=float))
    nw = len(wavelength)

    def per_wavelength(value, name):
        v = np.atleast_1d(np.asarray(value, dtype=float))
        if len(v) == 1:
            v = np.full(nw, v[0])
        if len(v) != nw:
            raise ValueError(f"{name} must be a scalar or have one value per wavelength")
        return v

    tau = per_wavelength(optical_depth, "optical_depth")
    ssa = per_wavelength(single_scattering, "single_scattering")
    g = per_wavelength(asymmetry, "asymmetry")

    return {
        "pressure": [float(bottom_pressure), float(top_pressure)],
        "wavelength": wavelength.tolist(),
        "optical_depth": [tau.tolist(), [0.0] * nw],
        "single_scattering": [ssa.tolist(), ssa.tolist()],
        "asymmetry": [g.tolist(), g.tolist()],
    }


def write_table(filename, table, comment=None):
    """Write a table (the dict returned by :func:`slab`, or one with the same keys) in the file
    format of the ``fixed`` cloud module: columns pressure [bar], wavelength [micron], optical
    depth, single-scattering albedo, asymmetry parameter."""
    pressure = np.asarray(table["pressure"], dtype=float)
    wavelength = np.asarray(table["wavelength"], dtype=float)
    tau = np.asarray(table["optical_depth"], dtype=float)
    ssa = np.asarray(table["single_scattering"], dtype=float)
    g = np.asarray(table["asymmetry"], dtype=float)

    with open(filename, "w") as f:
        if comment:
            for line in str(comment).splitlines():
                f.write(f"# {line}\n")
        f.write("# pressure[bar]  wavelength[micron]  optical_depth  single_scattering  asymmetry\n")
        for k, p in enumerate(pressure):
            for m, wl in enumerate(wavelength):
                f.write(f"{p:.8e}  {wl:.8e}  {tau[k, m]:.8e}  {ssa[k, m]:.8e}  {g[k, m]:.8e}\n")
