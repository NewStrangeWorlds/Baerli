"""Tests of the FastChem condensation option and the Ackerman & Marley cloud module (phase 2 of
doc/cloud_module_design.md). Needs the opacity data root; FastChem's condensate data ships with
the repository (fastchem_data/logK_condensates.dat, fastchem_parameters_condensates.dat).

Run with  python -m pytest tests/test_ackerman_marley.py  from the repository root.
"""

import os

import numpy as np
import pytest

import baerli

OPACITY_PATH = "/media/data/opacity_data/helios-k/"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARAMS = os.path.join(REPO, "fastchem_parameters_condensates.dat")
PARAMS_GAS = os.path.join(REPO, "fastchem_parameters.dat")

pytestmark = pytest.mark.skipif(
    not os.path.isdir(OPACITY_PATH), reason="opacity data not available")


@pytest.fixture(scope="module")
def grid():
    return baerli.SpectralGrid.constant_resolution(
        OPACITY_PATH, resolution=20.0, wavelength_min=0.5, wavelength_max=30.0)


def make_model(grid, clouds, condensation="rainout", params=PARAMS, teff=1400.0, **kwargs):
    chem = {"parameter_file": params}
    if condensation:
        chem["condensation"] = condensation
    return baerli.BrownDwarf(
        grid,
        effective_temperature=teff,
        surface_gravity=10**5,
        nb_grid_points=60,
        boundary_pressures=[1e2, 1e-6],
        opacity_species=[("CIA-H2-H2", "CIA/H2-H2"), ("CIA-H2-He", "CIA/H2-He"),
                         ("H2O", "Molecules/1H2-16O__POKAZATEL_e2b")],
        chemistry=[("equilibrium", chem)],
        clouds=clouds,
        solver=("ratio_ul", dict(max_iterations=kwargs.pop("max_iterations", 2),
                                 convergence_threshold=1e-4)),
        **kwargs)


def species_index(symbol):
    return list(baerli.species_symbols()).index(symbol)


# --- the chemistry side ----------------------------------------------------------------------------

def test_rainout_provides_increments_and_depletes_the_gas(grid):
    model = make_model(grid, None)
    model.initialize(("milne", dict(kappa_ross=1e-2)))

    atm = model.atmosphere
    symbols = list(atm.condensate_symbols)
    assert "Mg2SiO4(s,l)" in symbols and "Fe(s,l)" in symbols
    weights = dict(zip(symbols, atm.condensate_weights))
    assert weights["Fe(s,l)"] == pytest.approx(55.845, rel=1e-3)
    assert weights["Mg2SiO4(s,l)"] == pytest.approx(140.69, rel=1e-3)
    dens = np.asarray(atm.condensate_densities)
    assert dens.shape == (60, len(symbols))
    assert np.all(dens >= 0)

    # something condenses somewhere in a 1400 K atmosphere spanning 100 bar to 1 microbar
    assert dens.sum() > 0

    # the gas phase is depleted above the condensation level: atomic Fe at the top is far below its
    # deep value relative to the total density
    nd = np.asarray(atm.number_densities)
    fe = nd[:, species_index("Fe")] / nd[:, species_index("Total")]
    assert fe[-1] < 1e-3 * fe[0]


def test_rainout_is_lagged_inside_the_iteration(grid):
    """The committed call runs the rainout; the trial evaluations run the depleted gas phase.
    Both paths must leave a consistent composition: after a short run the gas is still depleted
    and the condensate table is still present."""
    model = make_model(grid, None, max_iterations=2)
    model.initialize(("milne", dict(kappa_ross=1e-2)))
    model.compute()
    atm = model.atmosphere
    assert len(atm.condensate_symbols) > 0
    nd = np.asarray(atm.number_densities)
    fe = nd[:, species_index("Fe")] / nd[:, species_index("Total")]
    assert fe[-1] < 1e-3 * fe[0]


def test_condensation_without_condensate_data_is_rejected(grid):
    with pytest.raises(Exception, match="condensate data"):
        make_model(grid, None, params=PARAMS_GAS)


def test_unknown_condensation_mode_is_rejected(grid):
    with pytest.raises(Exception, match="condensation must be"):
        make_model(grid, None, condensation="snow")


# --- the Ackerman & Marley module ------------------------------------------------------------------

def test_condensate_without_increments_gives_no_cloud(grid):
    """In rainout equilibrium forsterite takes the Mg before enstatite can form, so an MgSiO3
    module is legitimately empty in this atmosphere."""
    model = make_model(grid, [("ackerman_marley", dict(condensate="MgSiO3"))])
    model.initialize(("milne", dict(kappa_ross=1e-2)))
    assert np.sum(model.cloud_optical_depth(1.0)) == 0.0


def test_cloud_forms_on_the_rainout_increments(grid):
    model = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", fsed=2.0, relax=1.0))])
    model.initialize(("milne", dict(kappa_ross=1e-2)))

    diag = model.cloud_diagnostics()[0]
    q_c = np.asarray(diag["condensate_mmr"])
    N = np.asarray(diag["column_density"])
    r_g = np.asarray(diag["r_g"])          # micron
    assert q_c.max() > 0 and N.max() > 0
    cloudy = N > 0
    assert np.all(r_g[cloudy] > 0)
    assert np.all((r_g[cloudy] >= 1e-4) & (r_g[cloudy] <= 542.0))   # inside the Mie table
    # at the cloud base (the layer with the largest column) the particles are sub-micron to micron
    base = np.argmax(N)
    assert 1e-3 < r_g[base] < 100.0
    assert np.allclose(np.asarray(diag["fsed"]), 2.0)

    tau = np.asarray(model.cloud_optical_depth(1.0))
    assert tau.sum() > 0 and np.all(tau[~cloudy] == 0)
    # the whole solar Mg budget condensing near 10 bar in a log g = 5 object gives a column of
    # ~0.2 g/cm^2 -- optically very thick at the base (tau ~ 1e3 for micron particles) but finite
    assert 1.0 < tau.sum() < 1e5
    # mass budget: the condensate mass mixing ratio is bounded by the solar Mg budget (~1e-3)
    assert q_c.max() < 1e-2

    # above the last condensing layer the suspended condensate decays like exp(-f_sed dz / H)
    inc = np.asarray(model.atmosphere.condensate_densities)[:, list(model.atmosphere.condensate_symbols).index("Mg2SiO4(s,l)")]
    top_level = np.max(np.nonzero(inc)[0])            # highest level with an increment
    j0 = top_level                                     # layer index just above it
    z = np.asarray(model.atmosphere.altitude); H = np.asarray(model.atmosphere.scale_height)
    for j in range(j0 + 1, min(j0 + 4, len(q_c))):
        if q_c[j] <= 0 or q_c[j-1] <= 0:
            continue
        dz = z[j+1] - z[j]; L = 0.5 * (H[j] + H[j+1])
        assert q_c[j] / q_c[j-1] == pytest.approx(np.exp(-2.0 * dz / L), rel=1e-8)


def test_larger_fsed_gives_larger_particles_and_less_column(grid):
    low = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", fsed=0.5, relax=1.0, max_change=0))])
    high = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", fsed=5.0, relax=1.0, max_change=0))])
    for m in (low, high):
        m.initialize(("milne", dict(kappa_ross=1e-2)))
    d_low, d_high = low.cloud_diagnostics()[0], high.cloud_diagnostics()[0]
    N_low, N_high = np.asarray(d_low["column_density"]), np.asarray(d_high["column_density"])
    rg_low, rg_high = np.asarray(d_low["r_g"]), np.asarray(d_high["r_g"])
    # where most of the condensate MASS sits the particles grow with f_sed (A&M eq. 13: r_g ~
    # f_sed^(1/alpha)); the column NUMBER density is dominated by the nanometre tail, so it is
    # not the right weight, and the thinner f_sed = 5 cloud has no tail where the other still has
    q_low, q_high = np.asarray(d_low["condensate_mmr"]), np.asarray(d_high["condensate_mmr"])
    main = (q_low > 0.1 * q_low.max()) & (q_high > 0)
    assert main.any()
    assert np.all(rg_high[main] > rg_low[main])
    assert q_high.sum() < q_low.sum()
    # (the optical depth at 1 micron need not decrease: small f_sed makes sub-micron particles that
    # extinguish less per unit mass at this wavelength)


def test_cloud_ramp_caps_the_first_onset(grid):
    """The model-level cloud ramp: no cloud while the corrector has not run, then the first onset
    carries at most tau_start of optical depth at 1 micron (convection off so there is no MLT
    easy start to wait for)."""
    model = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", relax=1.0, max_change=0))],
                       max_iterations=2, convection="none")
    model.initialize(("milne", dict(kappa_ross=1e-2)))
    assert model.cloud_ramp_scale() == 0.0
    full = np.sum(model.cloud_optical_depth(1.0))          # the module's optics, unscaled
    assert full > 1.0
    model.compute()
    scale = model.cloud_ramp_scale()
    assert 0.0 < scale <= 1.0
    tau_rt = np.sum(np.asarray(model.cloud_optics(1.0)[0]))   # what the RT saw last
    assert tau_rt == pytest.approx(scale * np.sum(model.cloud_optical_depth(1.0)), rel=1e-6)
    assert tau_rt <= 1.0 + 1e-9


def test_cloud_ramp_none_applies_the_full_cloud(grid):
    model = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", relax=1.0, max_change=0))],
                       cloud_ramp="none")
    model.initialize(("milne", dict(kappa_ross=1e-2)))
    assert model.cloud_ramp_scale() == 1.0
    assert np.sum(np.asarray(model.cloud_optics(1.0)[0])) == pytest.approx(
        np.sum(model.cloud_optical_depth(1.0)), rel=1e-12)


def test_runs_inside_the_iteration_loop(grid):
    model = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", fsed=2.0)),
                              ("ackerman_marley", dict(condensate="Fe", fsed=2.0))],
                       max_iterations=3)
    model.initialize(("milne", dict(kappa_ross=1e-2)))
    model.compute()
    assert model.cloud_module_types() == ["ackerman_marley", "ackerman_marley"]
    tau = np.asarray(model.cloud_optical_depth_per_module(1.0))
    assert tau.shape[0] == 2 and tau.sum() > 0
    diag = model.cloud_diagnostics()
    assert set(diag[0]) >= {"condensate_mmr", "column_density", "r_g", "r_eff", "fsed", "alpha", "r_w"}


def test_module_needs_condensation_chemistry(grid):
    with pytest.raises(Exception, match="rainout"):
        make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4"))], condensation=None,
                   params=PARAMS_GAS)


def test_restart_from_file_keeps_working(grid, tmp_path):
    """A restart restores the composition without running the chemistry; the cloud module must
    wait for the first iteration instead of failing."""
    model = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4"))])
    model.initialize(("milne", dict(kappa_ross=1e-2)))
    out = tmp_path / "state.nc"
    baerli.save_model(str(out), model, grid)
    restarted = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4"))], max_iterations=1)
    restarted.initialize_from_file(str(out))
    restarted.compute()
    assert np.sum(restarted.cloud_optical_depth(1.0)) > 0


def test_parameter_validation(grid):
    with pytest.raises(Exception, match="no counterpart in FastChem"):
        make_model(grid, [("ackerman_marley", dict(condensate="tholin"))])
    with pytest.raises(Exception, match="unknown condensate"):
        make_model(grid, [("ackerman_marley", dict(condensate="unobtainium"))])
    with pytest.raises(Exception, match="fsed_mode"):
        make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", fsed_mode="linear"))])
    with pytest.raises(Exception, match="Unknown parameter"):
        make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", supersaturation=0.1))])


def test_output_file_contains_condensates_and_diagnostics(grid, tmp_path):
    model = make_model(grid, [("ackerman_marley", dict(condensate="Mg2SiO4", relax=1.0, max_change=0))])
    model.initialize(("milne", dict(kappa_ross=1e-2)))
    out = tmp_path / "am.nc"
    baerli.save_model(str(out), model, grid)
    ds = baerli.load_model_data(str(out))
    assert "condensate_densities" in ds and "cloud0_column_density" in ds and "cloud0_r_g" in ds
    assert "Mg2SiO4(s,l)" in list(ds["condensate_symbol"].values)
    ds.close()
