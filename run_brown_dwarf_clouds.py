"""Cloudy brown dwarf with the Ackerman & Marley cloud module.

Runs the same 1300 K, log g = 5 brown dwarf twice -- clear (equilibrium chemistry with FastChem's
rainout condensation, no cloud opacity) and with Ackerman & Marley forsterite and iron clouds
(f_sed = 2) -- and prints the temperature and spectral differences.

    python run_brown_dwarf_clouds.py                # both runs
    python run_brown_dwarf_clouds.py --cloudy-only  # reuse output_brown_dwarf_clear.nc

Notes (see doc/cloud_module_design.md, "Coupled-run debugging log"):
  * the cloud modules need the equilibrium chemistry with condensation="rainout" and a FastChem
    parameter file that lists the condensate data (fastchem_parameters_condensates.dat);
  * cloudy runs use the flux-scaled MLT Kzz, continued above the convective zone at the zone's
    value (radiative="constant"); the default velocity scaling gives 1e9-1e10 cm2/s everywhere
    above the RCB and hence 100-500 micron particles;
  * do NOT use radiative="fixed" with a low value here: the decks sit where the cloud heating
    creates detached convective patches, and a 3-dex Kzz step between a patch and the fixed
    value flips the particle size 30x whenever a patch appears or vanishes -- the coupled loop
    then cycles with a period of two cloud commits and never converges (measured 2026-10-27:
    T_bot swinging +-120 K, iron base jumping between 1.3 and 10 bar). A steeply decaying
    power law (slope 2, picaso-like) gives Kzz ~1e7 at the deck base and thick decks during the
    ramp; it converges too (214 iterations, with one ramp back-off) but the constant extension
    is faster (107 iterations, tau ~ 4-5 per deck);
  * the cloud opacity enters through the model-level ramp (cloud_ramp): nothing during the
    corrector's easy start, then an onset capped at tau = 1 at 1 micron that grows by `factor`
    (2) while the corrector is settled; a step that lifts the residual above `backoff` (2) is
    undone and the factor reduced (a doubling into a thick deck can throw the Newton out of
    its contraction region -- measured as a 90-iteration freeze at lambda_min);
  * the corrector freeze seen in early coupled runs was a hard temperature cut in the
    band-closure correction of the ratio_ul solver (fixed 2026-10-14 with a smooth taper).
"""
import sys
import time

import numpy as np

import pyngam


# --- Model configuration ---

opacity_path = "/media/data/opacity_data/helios-k/"

opacity_species = [
    ("CIA-H2-H2", "CIA/H2-H2"),
    ("CIA-H2-He",  "CIA/H2-He"),
    ("K",          "Alkali_Allard/K"),
    ("Na",         "Alkali_Allard/Na"),
    ("H2O",        "Molecules/1H2-16O__POKAZATEL_e2b"),
    ("CH4",        "Molecules/12C-1H4__YT34to10_e2b"),
    ("NH3",        "Molecules/14N-1H3__CoYuTe_e2b"),
    ("CO",         "Molecules/12C-16O__Li2015_e2b"),
    ("CO2",        "Molecules/12C-16O2__CDSD_4000_e2b"),
    ("H2S",        "Molecules/1H2-32S__AYT2_e2b"),
]

clouds = [
    ("ackerman_marley", dict(condensate="Mg2SiO4", fsed=2.0, sigma=2.0)),
    ("ackerman_marley", dict(condensate="Fe",      fsed=2.0, sigma=2.0)),
]

grid = pyngam.SpectralGrid.constant_resolution(
    opacity_path, resolution=300.0, wavelength_min=0.5, wavelength_max=30.0)


def build(cloud_modules):
    return pyngam.BrownDwarf(
        grid,
        effective_temperature=1300.0,      # K
        surface_gravity=10**5.0,           # cm/s^2
        nb_grid_points=100,
        boundary_pressures=[1e2, 1e-6],    # bar: bottom -> top
        opacity_species=opacity_species,
        chemistry=[("equilibrium", dict(parameter_file="fastchem_parameters_condensates.dat",
                                        condensation="rainout"))],
        clouds=cloud_modules,
        convection="mlt_dry",
        kzz=("mlt", dict(scaling="flux", radiative="constant")),
        solver=("ratio_ul", dict(max_iterations=200, convergence_threshold=1e-4)))


def run(cloud_modules, filename):
    model = build(cloud_modules)
    model.initialize(("milne", dict(kappa_ross=1e-2)))
    start = time.time()
    converged = model.compute()
    print(f"\n  {filename}: converged = {converged} after {time.time() - start:.0f} s")
    pyngam.save_model(filename, model, grid)
    return model


# --- Run ---

cloudy_only = "--cloudy-only" in sys.argv

if cloudy_only:
    clear = pyngam.load_model_data("output_brown_dwarf_clear.nc")
    pressure = clear["pressure"].values
    temperature_clear = clear["temperature"].values
    spectrum_clear = clear["spectrum"].values
    clear.close()
else:
    model_clear = run(None, "output_brown_dwarf_clear.nc")
    pressure = np.asarray(model_clear.atmosphere.pressure)
    temperature_clear = np.asarray(model_clear.atmosphere.temperature)
    spectrum_clear = np.asarray(model_clear.radiation_field.spectrum)

model_cloudy = run(clouds, "output_brown_dwarf_cloudy.nc")


# --- Comparison ---

wavelength = np.asarray(grid.wavelength_list)
log_p = np.log(pressure)[::-1]                   # increasing, for np.interp


def temperature_at(temperature, p):
    return np.interp(np.log(p), log_p, np.asarray(temperature)[::-1])


def band_flux(spectrum, lo, hi):
    mask = (wavelength > lo) & (wavelength < hi)
    return np.asarray(spectrum)[mask].mean()


temperature_cloudy = np.asarray(model_cloudy.atmosphere.temperature)
spectrum_cloudy = np.asarray(model_cloudy.radiation_field.spectrum)

print("\n   p [bar]   T_clear  T_cloudy")
for p in (30, 10, 3, 1, 0.3, 0.1, 0.01):
    print(f"  {p:8.2f}  {temperature_at(temperature_clear, p):8.1f}  {temperature_at(temperature_cloudy, p):8.1f}")

for name, lo, hi in (("J", 1.17, 1.33), ("H", 1.5, 1.8), ("K", 2.0, 2.4), ("L", 3.4, 4.1)):
    print(f"  F_{name} cloudy / clear = {band_flux(spectrum_cloudy, lo, hi) / band_flux(spectrum_clear, lo, hi):.3f}")

layer_pressure = np.sqrt(pressure[:-1] * pressure[1:])
for k, (kind, diag) in enumerate(zip(model_cloudy.cloud_module_types(), model_cloudy.cloud_diagnostics())):
    tau = np.asarray(model_cloudy.cloud_optical_depth_per_module(1.0)[k])
    q_c = np.asarray(diag["condensate_mmr"])
    base = np.argmax(q_c)
    print(f"  module {k} ({kind}, {clouds[k][1]['condensate']}): tau(1 um) = {tau.sum():.2f}, "
          f"base {layer_pressure[base]:.2f} bar, r_g(base) = {diag['r_g'][base]:.3f} um, "
          f"max q_c = {q_c.max():.2e}, ramp scale {model_cloudy.cloud_ramp_scale():.3g}")
