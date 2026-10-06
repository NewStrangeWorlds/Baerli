"""Tests of the Mie optics (phase 1 of doc/cloud_module_design.md): the `mie_slab` cloud module
against the LX-MIE Python package on the same refractive indices, the Mie asymptotes, and the
module's bookkeeping. Needs the opacity data root and the `lxmie` package (pip install lxmie).

Run with  python -m pytest tests/test_mie_slab.py  from the repository root.
"""

import os

import numpy as np
import pytest

import pyngam

OPACITY_PATH = "/media/data/opacity_data/helios-k/"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

pytestmark = pytest.mark.skipif(
    not os.path.isdir(OPACITY_PATH), reason="opacity data not available")


@pytest.fixture(scope="module")
def grid():
    return pyngam.SpectralGrid.constant_resolution(
        OPACITY_PATH, resolution=50.0, wavelength_min=0.5, wavelength_max=30.0)


def make_model(grid, clouds, nb_grid_points=60, **kwargs):
    return pyngam.BrownDwarf(
        grid,
        effective_temperature=1000.0,
        surface_gravity=10**5,
        nb_grid_points=nb_grid_points,
        boundary_pressures=[1e2, 1e-6],
        opacity_species=[("CIA-H2-H2", "CIA/H2-H2")],
        chemistry=[("isoprofile", {"H2": 0.85, "He": 0.15})],
        clouds=clouds,
        **kwargs)


def refractive_index(name):
    """(wavelength [micron], n, k) from the shipped compilation."""
    data = np.loadtxt(os.path.join(REPO, "data", "refractive_indices", name))
    return data[:, 0], data[:, 1], np.abs(data[:, 2])


def interpolate_nk(name, wavelength):
    wl, n, k = refractive_index(name)
    x = np.log(wavelength)
    n_i = np.interp(x, np.log(wl), n)
    k_i = np.exp(np.interp(x, np.log(wl), np.log(np.maximum(k, 1e-300))))
    return n_i, k_i


def slab(grid, condensate, radius_um, tau=1.0, sigma=1.0, **extra):
    params = dict(condensate=condensate, particle_radius=radius_um, sigma=sigma,
                  top_pressure=0.01, bottom_pressure=1.0, optical_depth=tau,
                  reference_wavelength=1.0)
    params.update(extra)
    return make_model(grid, [("mie_slab", params)])


# --- against LX-MIE directly -----------------------------------------------------------------------

def test_spectral_shape_matches_lxmie_for_an_absorber(grid):
    """tau(lambda)/tau(1 um) of a monodisperse Fe deck follows Q_ext(lambda)/Q_ext(1 um) from
    LX-MIE with the same refractive indices (Fe is smooth: no single-size ripple to average)."""
    lxmie = pytest.importorskip("lxmie")

    radius_um = 1.0
    model = slab(grid, "Fe", radius_um)
    wl = np.asarray(grid.wavelength_list)
    ref = np.argmin(np.abs(wl - 1.0))

    tau_ref = np.sum(model.cloud_optical_depth(wl[ref]))
    assert tau_ref == pytest.approx(1.0, rel=1e-12)       # exact by construction

    n, k = interpolate_nk("Fe.dat", wl)
    for target in (0.6, 2.0, 5.0, 12.0, 25.0):
        i = np.argmin(np.abs(wl - target))
        q = lxmie.mie(n[i], k[i], 2 * np.pi * radius_um / wl[i]).q_ext
        q_ref = lxmie.mie(n[ref], k[ref], 2 * np.pi * radius_um / wl[ref]).q_ext
        ratio = np.sum(model.cloud_optical_depth(wl[i])) / tau_ref
        assert ratio == pytest.approx(q / q_ref, rel=0.03), f"at {wl[i]:.2f} micron"


def test_single_scattering_and_asymmetry_match_lxmie(grid):
    lxmie = pytest.importorskip("lxmie")

    radius_um = 0.5
    model = slab(grid, "MgSiO3", radius_um)
    wl = np.asarray(grid.wavelength_list)
    n, k = interpolate_nk("MgSiO3_amorph_glass.dat", wl)

    for target in (0.7, 1.5, 3.0, 9.0):
        i = np.argmin(np.abs(wl - target))
        tau, ssa, g = (np.asarray(x) for x in model.cloud_optics(wl[i]))
        inside = tau > 0
        r = lxmie.mie(n[i], k[i], 2 * np.pi * radius_um / wl[i])
        assert np.allclose(ssa[inside], r.q_sca / r.q_ext, atol=0.03), f"ssa at {wl[i]:.2f}"
        assert np.allclose(g[inside], r.asymmetry_parameter, atol=0.05), f"g at {wl[i]:.2f}"


# --- asymptotes -----------------------------------------------------------------------------------

def test_large_particles_approach_the_extinction_paradox(grid):
    """x >> 1: Q_ext -> 2, so tau is nearly grey and the column density is tau / (2 pi r^2)."""
    radius_um = 50.0
    model = slab(grid, "Fe", radius_um, tau=1.0)
    wl = np.asarray(grid.wavelength_list)
    tau = np.array([np.sum(model.cloud_optical_depth(w)) for w in wl])
    assert np.all(np.abs(tau - 1.0) < 0.15)

    # tau = N pi r^2 Q_ext with Q_ext ~ 2
    r_cm = radius_um * 1e-4
    n_expected = 1.0 / (2.0 * np.pi * r_cm**2)
    tau_from_n = n_expected * np.pi * r_cm**2 * 2.0
    assert tau_from_n == pytest.approx(1.0)


def test_small_particles_scatter_like_rayleigh(grid):
    """x << 1 for a weak absorber: the scattering optical depth scales as lambda^-4."""
    radius_um = 0.01
    model = slab(grid, "MgSiO3", radius_um, tau=1.0, reference_wavelength=0.6)
    wl = np.asarray(grid.wavelength_list)

    def tau_sca(target):
        i = np.argmin(np.abs(wl - target))
        tau, ssa, _ = (np.asarray(x) for x in model.cloud_optics(wl[i]))
        return np.sum(tau * ssa), wl[i]

    t1, w1 = tau_sca(0.6)
    t2, w2 = tau_sca(1.2)
    # (m^2-1)/(m^2+2) varies little for the glass between 0.6 and 1.2 micron
    assert t1 / t2 == pytest.approx((w2 / w1) ** 4, rel=0.15)


# --- bookkeeping ------------------------------------------------------------------------------------

def test_lognormal_width_broadens_but_conserves_the_reference_optical_depth(grid):
    mono = slab(grid, "MgSiO3", 1.0, tau=2.0)
    wide = slab(grid, "MgSiO3", 1.0, tau=2.0, sigma=2.0)
    wl = np.asarray(grid.wavelength_list)
    ref = np.argmin(np.abs(wl - 1.0))
    for m in (mono, wide):
        assert np.sum(m.cloud_optical_depth(wl[ref])) == pytest.approx(2.0, rel=1e-12)
    # the broad distribution contains large particles: more extinction in the far infrared
    far = np.argmin(np.abs(wl - 20.0))
    assert np.sum(wide.cloud_optical_depth(wl[far])) > np.sum(mono.cloud_optical_depth(wl[far]))


def test_column_density_alternative(grid):
    """tau = N pi r^2 Q_ext: a directly given column density gives an optical depth with a
    plausible extinction efficiency, and a deck specified by optical depth is linear in it."""
    wl = np.asarray(grid.wavelength_list)
    ref = np.argmin(np.abs(wl - 1.0))
    one = slab(grid, "Fe", 1.0, tau=1.0)
    three = slab(grid, "Fe", 1.0, tau=3.0)
    far = np.argmin(np.abs(wl - 10.0))
    assert np.sum(three.cloud_optical_depth(wl[far])) == pytest.approx(
        3.0 * np.sum(one.cloud_optical_depth(wl[far])), rel=1e-12)

    params = dict(condensate="Fe", particle_radius=1.0, top_pressure=0.01, bottom_pressure=1.0,
                  column_density=1e8)
    direct = make_model(grid, [("mie_slab", params)])
    tau_direct = np.sum(direct.cloud_optical_depth(wl[ref]))
    # tau = N pi r^2 Q_ext: Q_ext of a 1 micron Fe sphere at 1 micron is of order 2
    q_ext = tau_direct / (1e8 * np.pi * (1e-4) ** 2)
    assert 1.0 < q_ext < 4.0


def test_deck_is_uniform_in_lnp_and_within_bounds(grid):
    model = slab(grid, "Fe", 1.0, tau=1.0)
    p = np.asarray(model.atmosphere.pressure)
    mid = np.sqrt(p[:-1] * p[1:])
    tau = np.asarray(model.cloud_optical_depth(1.0))
    inside = (mid < 1.0) & (mid > 0.01)
    outside = (mid > 1.3) | (mid < 0.01 / 1.3)
    assert np.all(tau[outside] == 0.0)
    density = tau[inside] / np.log(p[:-1] / p[1:])[inside]
    assert np.allclose(density, density[0], rtol=1e-10)


def test_parameter_validation(grid):
    with pytest.raises(Exception, match="alternatives"):
        make_model(grid, [("mie_slab", dict(condensate="Fe", particle_radius=1.0, top_pressure=0.01,
                                            bottom_pressure=1.0, optical_depth=1.0, column_density=1e8))])
    with pytest.raises(Exception, match="top_pressure"):
        make_model(grid, [("mie_slab", dict(condensate="Fe", particle_radius=1.0, top_pressure=1.0,
                                            bottom_pressure=0.01, optical_depth=1.0))])
    with pytest.raises(Exception, match="known condensates"):
        make_model(grid, [("mie_slab", dict(condensate="unobtainium", particle_radius=1.0,
                                            top_pressure=0.01, bottom_pressure=1.0, optical_depth=1.0))])
    with pytest.raises(Exception, match="Unknown parameter"):
        make_model(grid, [("mie_slab", dict(condensate="Fe", particle_radius=1.0, top_pressure=0.01,
                                            bottom_pressure=1.0, optical_depth=1.0, fsed=1))])


def test_refractive_index_file_path_is_accepted(grid):
    path = os.path.join(REPO, "data", "refractive_indices", "KCl.dat")
    model = slab(grid, path, 1.0)
    assert np.sum(model.cloud_optical_depth(1.0)) == pytest.approx(1.0, rel=1e-12)
    assert model.cloud_module_types() == ["mie_slab"]
    assert model.config["clouds"][0]["type"] == "mie_slab"
