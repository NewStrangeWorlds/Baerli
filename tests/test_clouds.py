"""Tests of the `clouds` model option (phase 0 of doc/cloud_module_design.md): the `fixed` module's
conservative mapping onto the model grid, the scattering-weighted combination of several modules,
the Python helpers, and the zero-cloud identity of the iteration loop.

Run with  python -m pytest tests/test_clouds.py  from the repository root (needs the opacity data
root below; the tests are skipped without it).
"""

import os

import numpy as np
import pytest

import baerli

OPACITY_PATH = "/media/data/opacity_data/helios-k/"

pytestmark = pytest.mark.skipif(
    not os.path.isdir(OPACITY_PATH), reason="opacity data not available")


@pytest.fixture(scope="module")
def grid():
    return baerli.SpectralGrid.constant_resolution(
        OPACITY_PATH, resolution=20.0, wavelength_min=0.5, wavelength_max=30.0)


def make_model(grid, clouds, nb_grid_points=60, **kwargs):
    return baerli.BrownDwarf(
        grid,
        effective_temperature=1000.0,
        surface_gravity=10**5,
        nb_grid_points=nb_grid_points,
        boundary_pressures=[1e2, 1e-6],
        opacity_species=[("CIA-H2-H2", "CIA/H2-H2"), ("CIA-H2-He", "CIA/H2-He")],
        chemistry=[("isoprofile", {"H2": 0.85, "He": 0.15})],
        clouds=clouds,
        cloud_ramp=kwargs.pop("cloud_ramp", "none"),
        **kwargs)


def layer_dlnp(model):
    p = np.asarray(model.atmosphere.pressure)
    return np.log(p[:-1] / p[1:])


def layer_midpoints(model):
    p = np.asarray(model.atmosphere.pressure)
    return np.sqrt(p[:-1] * p[1:])


# --- the fixed module ----------------------------------------------------------------------------

@pytest.mark.parametrize("nb_grid_points", [60, 121])
def test_slab_is_conservative_and_uniform_in_lnp(grid, nb_grid_points):
    model = make_model(grid, [("fixed", {})], nb_grid_points=nb_grid_points)
    table = baerli.clouds.slab(2.0, top_pressure=0.01, bottom_pressure=1.0,
                               single_scattering=0.3, asymmetry=0.6)
    model.set_cloud_optical_properties(**table)

    tau, ssa, g = (np.asarray(x) for x in model.cloud_optics(1.0))
    assert tau.shape == (nb_grid_points - 1,)

    # the column optical depth is reproduced on any model grid
    assert tau.sum() == pytest.approx(2.0, rel=1e-10)

    # nothing outside the deck; a uniform dtau/dlnp inside it
    mid = layer_midpoints(model)
    inside = (mid < 1.0) & (mid > 0.01)
    outside = (mid > 1.0 * 1.3) | (mid < 0.01 / 1.3)       # layers not straddling the edges
    assert np.all(tau[outside] == 0.0)
    density = tau[inside] / layer_dlnp(model)[inside]
    assert np.allclose(density, 2.0 / np.log(1.0 / 0.01), rtol=1e-10)

    # w and g of the deck where it exists, zero elsewhere
    assert np.allclose(ssa[tau > 0], 0.3)
    assert np.allclose(g[tau > 0], 0.6)
    assert np.all(ssa[tau == 0] == 0.0)

    # the diagnostics agree with the combined optics for a single module
    assert np.allclose(model.cloud_optical_depth(1.0), tau)
    assert model.cloud_module_types() == ["fixed"]


def test_multi_level_table_is_conservative(grid):
    """A table with several input layers of different optical depth."""
    model = make_model(grid, [("fixed", {})])
    pressure = [10.0, 1.0, 0.1, 0.01]            # three input layers
    tau_in = [[0.5], [2.0], [0.25], [0.0]]
    model.set_cloud_optical_properties(
        pressure=pressure, wavelength=[1.0], optical_depth=tau_in,
        single_scattering=[[0.0]] * 4, asymmetry=[[0.0]] * 4)

    tau = np.asarray(model.cloud_optical_depth(1.0))
    assert tau.sum() == pytest.approx(2.75, rel=1e-10)

    # the model's cumulative optical depth at every model level equals the table's cumulative
    # optical depth interpolated (linearly in ln p) to that level
    p = np.asarray(model.atmosphere.pressure)
    cumulative = np.concatenate([[0.0], np.cumsum(tau)])        # at the model levels
    expected = np.interp(-np.log(p), -np.log(pressure), np.cumsum([0.0, 0.5, 2.0, 0.25]))
    assert np.allclose(cumulative, expected, atol=1e-12)


def test_wavelength_interpolation(grid):
    model = make_model(grid, [("fixed", {})])
    table = baerli.clouds.slab([1.0, 3.0], top_pressure=0.01, bottom_pressure=1.0,
                               wavelength=[1.0, 10.0])
    model.set_cloud_optical_properties(**table)

    # linear between the table's wavelengths, constant outside (evaluated at the grid's own points)
    wl = np.asarray(grid.wavelength_list)
    for target in (0.6, 1.0, 5.5, 10.0, 25.0):
        k = np.argmin(np.abs(wl - target))
        expected = 1.0 + 2.0 * np.clip((wl[k] - 1.0) / 9.0, 0.0, 1.0)
        assert np.sum(model.cloud_optical_depth(wl[k])) == pytest.approx(expected, rel=1e-10)


def test_file_round_trip(grid, tmp_path):
    table = baerli.clouds.slab(1.5, top_pressure=0.05, bottom_pressure=0.5,
                               single_scattering=0.9, asymmetry=0.7)
    path = tmp_path / "deck.dat"
    baerli.clouds.write_table(path, table, comment="test deck")

    from_file = make_model(grid, [("fixed", {"file": str(path)})])
    from_file.initialize(("const", {"temperature": 1000.0}))
    from_file.compute(max_iterations=1)                         # the loop runs update()

    from_setter = make_model(grid, [("fixed", {})])
    from_setter.set_cloud_optical_properties(**table)

    for wl in (1.0, 5.0):
        a = np.asarray(from_file.cloud_optics(wl))
        b = np.asarray(from_setter.cloud_optics(wl))
        assert np.allclose(a, b, rtol=1e-12, atol=0.0)
    assert np.asarray(from_file.cloud_optical_depth(1.0)).sum() == pytest.approx(1.5, rel=1e-10)


# --- combining several modules -------------------------------------------------------------------

def test_two_modules_combine_with_scattering_weights(grid, tmp_path):
    """An absorbing and a scattering deck of equal optical depth: w = 1/2, g = the scatterer's."""
    absorber = baerli.clouds.slab(1.0, top_pressure=0.01, bottom_pressure=1.0,
                                  single_scattering=0.0, asymmetry=0.0)
    path = tmp_path / "absorber.dat"
    baerli.clouds.write_table(path, absorber)

    model = make_model(grid, [("fixed", {"file": str(path)}), ("fixed", {})])
    scatterer = baerli.clouds.slab(1.0, top_pressure=0.01, bottom_pressure=1.0,
                                   single_scattering=1.0, asymmetry=0.8)
    model.set_cloud_optical_properties(**scatterer)
    model.initialize(("const", {"temperature": 1000.0}))
    model.compute(max_iterations=1)

    tau, ssa, g = (np.asarray(x) for x in model.cloud_optics(1.0))
    assert tau.sum() == pytest.approx(2.0, rel=1e-10)
    assert np.allclose(ssa[tau > 0], 0.5)
    assert np.allclose(g[tau > 0], 0.8)

    per_module = np.asarray(model.cloud_optical_depth_per_module(1.0))
    assert per_module.shape == (2, model.atmosphere.nb_grid_points - 1)
    assert per_module.sum() == pytest.approx(2.0, rel=1e-10)
    assert model.cloud_module_types() == ["fixed", "fixed"]


# --- configuration -------------------------------------------------------------------------------

def test_unknown_parameter_is_rejected(grid):
    with pytest.raises(Exception, match="Unknown parameter"):
        make_model(grid, [("fixed", {"foo": 1})])


def test_unknown_type_is_rejected(grid):
    with pytest.raises(Exception, match="unknown"):
        make_model(grid, ["nimbus"])


def test_setter_needs_a_fixed_module(grid):
    model = make_model(grid, None)
    with pytest.raises(Exception, match="fixed"):
        model.set_cloud_optical_properties(**baerli.clouds.slab(1.0, 0.01, 1.0))


def test_bad_table_is_rejected(grid):
    model = make_model(grid, [("fixed", {})])
    table = baerli.clouds.slab(1.0, 0.01, 1.0)
    table["single_scattering"] = [[1.5], [1.5]]
    with pytest.raises(Exception, match="single-scattering"):
        model.set_cloud_optical_properties(**table)


def test_config_records_clouds(grid):
    model = make_model(grid, [("fixed", {})])
    assert model.config["clouds"] == [{"type": "fixed"}]
    model = make_model(grid, None)
    assert "clouds" not in model.config


# --- the iteration loop --------------------------------------------------------------------------

def test_empty_fixed_module_leaves_the_model_unchanged(grid):
    """Zero-cloud identity: a `fixed` module without a table is the same model as no clouds."""
    kwargs = dict(solver=("ratio_ul", dict(max_iterations=3, convergence_threshold=1e-5)))
    reference = make_model(grid, None, **kwargs)
    with_module = make_model(grid, [("fixed", {})], **kwargs)
    for m in (reference, with_module):
        m.initialize(("milne", dict(kappa_ross=1e-2)))
        m.compute()

    assert np.array_equal(reference.atmosphere.temperature, with_module.atmosphere.temperature)
    assert np.array_equal(reference.radiation_field.flux_total, with_module.radiation_field.flux_total)
    assert np.all(np.asarray(with_module.cloud_optical_depth(1.0)) == 0.0)


def test_cloud_changes_the_solution_and_is_saved(grid, tmp_path):
    kwargs = dict(solver=("ratio_ul", dict(max_iterations=3, convergence_threshold=1e-5)))
    clear = make_model(grid, None, **kwargs)
    cloudy = make_model(grid, [("fixed", {})], **kwargs)
    cloudy.set_cloud_optical_properties(**baerli.clouds.slab(
        5.0, top_pressure=0.1, bottom_pressure=1.0, single_scattering=0.5, asymmetry=0.7))
    for m in (clear, cloudy):
        m.initialize(("milne", dict(kappa_ross=1e-2)))
        m.compute()

    assert not np.array_equal(clear.atmosphere.temperature, cloudy.atmosphere.temperature)

    out = tmp_path / "cloudy.nc"
    baerli.save_model(str(out), cloudy, grid)
    ds = baerli.load_model_data(str(out))
    tau = ds["cloud_optical_depth"].values
    ds.close()
    assert tau.shape == (cloudy.atmosphere.nb_grid_points - 1,)
    assert tau.sum() == pytest.approx(5.0, rel=1e-10)
    config = baerli.load_output_config(str(out))
    assert config["model"]["clouds"] == [{"type": "fixed"}]
