/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#include "ackerman_marley_cloud.h"
#include "fall_velocity.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "../additional/exceptions.h"
#include "../additional/physical_const.h"
#include "../additional/thermodynamic_data.h"


namespace baerli {


AckermanMarleyCloud::AckermanMarleyCloud(
  SpectralGrid* spectral_grid_, const size_t nb_grid_points, const Parameters& parameters)
  : CloudModel(spectral_grid_, nb_grid_points), par(parameters)
{
  type = "ackerman_marley";

  const CondensateData* data = findCondensate(par.condensate);

  if (data == nullptr)
    throw InvalidInput("clouds",
      "ackerman_marley: unknown condensate '" + par.condensate + "'. Known: " + condensateNames() + "\n");
  if (data->fastchem_symbol.empty())
    throw InvalidInput("clouds",
      "ackerman_marley: condensate '" + par.condensate + "' has no counterpart in FastChem's "
      "condensate data, so its condensation increments cannot come from the equilibrium chemistry\n");

  fastchem_symbol = data->fastchem_symbol;
  bulk_density = par.bulk_density > 0 ? par.bulk_density : data->bulk_density;

  if (!(par.fsed > 0)) throw InvalidInput("clouds", "ackerman_marley: fsed must be > 0\n");
  if (!(par.sigma >= 1.0)) throw InvalidInput("clouds", "ackerman_marley: sigma must be >= 1\n");
  if (!(par.relax > 0 && par.relax <= 1)) throw InvalidInput("clouds", "ackerman_marley: relax must be in (0, 1]\n");
  if (!(par.seed > 0 && par.seed <= 1)) throw InvalidInput("clouds", "ackerman_marley: seed must be in (0, 1]\n");
  if (par.fsed_exponential && (!(par.b > 0) || !(par.eps > 0)))
    throw InvalidInput("clouds", "ackerman_marley: fsed_mode exp needs b > 0 and eps > 0\n");

  mie_table = std::make_unique<MieTable>(
    refractiveIndexFile(par.condensate), spectral_grid, par.mie_resolution, par.mie_radii);

  q_c.assign(nb_layers, 0.0);
  r_eff.assign(nb_layers, 0.0);
  fsed_layer.assign(nb_layers, 0.0);
  alpha_layer.assign(nb_layers, 0.0);
  r_w_layer.assign(nb_layers, 0.0);

  std::cout << "- Cloud: Ackerman & Marley " << par.condensate << " (" << fastchem_symbol
            << ", rho_p = " << bulk_density << " g/cm3), f_sed = " << par.fsed
            << (par.fsed_exponential ? " (exponential, b = " + std::to_string(par.b) + " H, eps = "
                                       + std::to_string(par.eps) + ")" : "")
            << ", sigma = " << par.sigma << ", Kzz floor " << std::scientific << std::setprecision(1)
            << par.kzz_min << std::defaultfloat << " cm2/s, mixing length "
            << (par.lapse_ratio_mixing_length ? "max(0.1, nabla/nabla_ad) H" : "H")
            << ", relax " << par.relax << ", cap " << par.max_change << " dex/update, dead band "
            << par.tolerance << " dex, gate " << par.gate << ", seed " << par.seed << "\n";
}



void AckermanMarleyCloud::update(const Atmosphere& atmosphere)
{
  const size_t nb_levels = atmosphere.nb_grid_points;

  // not initialised yet (constructor-time call): nothing to do, optics stay empty
  if (atmosphere.temperature.empty() || !(atmosphere.temperature[0] > 0)
      || atmosphere.mass_density.empty() || !(atmosphere.mass_density[0] > 0))
    return;

  // the chemistry has not run yet (e.g. a restart from a file restores the composition without
  // it): keep whatever is committed; the owning object checked at construction that a chemistry
  // module provides the condensed phase
  if (atmosphere.condensates.empty()) return;

  const int k = atmosphere.condensates.index(fastchem_symbol);
  if (k < 0)
  {
    std::string available;
    for (const auto& s : atmosphere.condensates.symbols) available += (available.empty() ? "" : ", ") + s;
    throw InvalidInput("clouds",
      "ackerman_marley: condensate " + fastchem_symbol + " is not in the chemistry's condensate list ("
      + available + ")\n");
  }

  const double molar_mass = atmosphere.condensates.weights[k];      // g/mol
  const std::vector<double>& z = atmosphere.altitude;
  const std::vector<double>& rho = atmosphere.mass_density;
  const std::vector<double>& T = atmosphere.temperature;
  const std::vector<double>& p = atmosphere.pressure;               // bar
  const std::vector<double>& mu = atmosphere.mean_molecular_weight;
  const std::vector<double>& H = atmosphere.scale_height;

  // --- condensation increments per layer (mass mixing ratio): the level above the layer, plus the
  //     bottom level for the bottom layer
  std::vector<double> dq(nb_layers, 0.0);
  constexpr double avogadro = constants::gas_constant / constants::boltzmann_k;   // 1/mol
  auto increment = [&](size_t level) {
    // formula units per cm^3 -> g/cm^3 -> mass mixing ratio
    return atmosphere.condensates.densities[level][k] * molar_mass / avogadro / rho[level];
  };
  for (size_t j = 0; j < nb_layers; ++j) dq[j] = increment(j + 1);
  dq[0] += increment(0);

  // --- layer quantities and the recursion for the suspended condensate
  std::vector<double> dz(nb_layers), L(nb_layers), T_l(nb_layers), p_l(nb_layers), rho_l(nb_layers),
                      mu_l(nb_layers), kzz_l(nb_layers), z_mid(nb_layers);
  for (size_t j = 0; j < nb_layers; ++j)
  {
    dz[j] = z[j+1] - z[j];
    z_mid[j] = 0.5 * (z[j] + z[j+1]);
    T_l[j] = 0.5 * (T[j] + T[j+1]);
    p_l[j] = std::sqrt(p[j] * p[j+1]) * 1e6;                        // dyn/cm^2
    rho_l[j] = 0.5 * (rho[j] + rho[j+1]);
    mu_l[j] = 0.5 * (mu[j] + mu[j+1]);
    const double kzz_level = kzz.size() == nb_levels ? 0.5 * (kzz[j] + kzz[j+1]) : 0.0;
    kzz_l[j] = std::max(kzz_level, par.kzz_min);

    L[j] = 0.5 * (H[j] + H[j+1]);
    if (par.lapse_ratio_mixing_length)
    {
      const double nabla = (std::log(T[j]) - std::log(T[j+1])) / (std::log(p[j]) - std::log(p[j+1]));
      const double nabla_ad = ThermodynamicData::adiabaticGradient(atmosphere.number_densities[j], T[j]);
      L[j] *= std::max(0.1, nabla_ad > 0 ? nabla / nabla_ad : 1.0);
    }
  }

  // cloud base: the first layer with an increment
  size_t base = nb_layers;
  for (size_t j = 0; j < nb_layers; ++j) if (dq[j] > 0) { base = j; break; }

  std::vector<double> q_c_new(nb_layers, 0.0);
  for (size_t j = 0; j < nb_layers; ++j)
  {
    fsed_layer[j] = par.fsed;
    if (par.fsed_exponential && base < nb_layers && j >= base)
      fsed_layer[j] = par.eps + (par.fsed - par.eps) * std::exp(-(z_mid[j] - z_mid[base]) / (par.b * L[j]));

    const double carried = j > 0 ? q_c_new[j-1] * std::exp(-fsed_layer[j] * dz[j] / L[j]) : 0.0;
    q_c_new[j] = carried + dq[j];
    if (q_c_new[j] < 1e-30) q_c_new[j] = 0.0;
  }

  // the suspended tail decays by orders of magnitude above the deck; below 1e-6 of the maximum it
  // carries no optical depth but would still produce (nanometre) particles with Kzz at its floor
  const double q_c_max = *std::max_element(q_c_new.begin(), q_c_new.end());
  for (auto& q : q_c_new) if (q < 1e-6 * q_c_max) q = 0.0;

  // --- particle sizes (A&M eqs. 13-18) and condensate mass columns per layer
  const std::vector<double>& radius_grid = mie_table->radii();
  std::vector<double> M_new(nb_layers, 0.0), rg_new(nb_layers, 0.0);
  const double l = 0.5 * std::log(par.sigma) * std::log(par.sigma);

  for (size_t j = 0; j < nb_layers; ++j)
  {
    q_c[j] = q_c_new[j];
    r_eff[j] = 0.0; alpha_layer[j] = 0.0; r_w_layer[j] = 0.0;
    if (q_c_new[j] <= 0) continue;

    const FallMedium medium = fallMedium(
      T_l[j], p_l[j], mu_l[j], gravity, par.molecule_diameter, par.lj_temperature);
    const double w_convect = kzz_l[j] / L[j];

    auto excess = [&](double r) { return fallVelocity(r, bulk_density, medium) - w_convect; };

    // r_w: the radius falling at the convective velocity scale
    double r_lo = 1e-8, r_hi = 1.0;
    double r_w;
    if (excess(r_lo) >= 0)
      r_w = r_lo;                                   // even the smallest particles fall faster
    else
    {
      while (excess(r_hi) < 0 && r_hi < 1e4) r_hi *= 10.0;
      r_w = brentRoot(excess, r_lo, r_hi, 1e-8);
    }

    // alpha: least-squares exponent of v(r) ~ w* (r/r_w)^alpha over the radius grid (virga)
    double num = 0.0, den = 0.0;
    for (const double r : radius_grid)
    {
      const double x = std::log(r / r_w);
      const double v = fallVelocity(r, bulk_density, medium);
      if (v <= 0) continue;
      num += std::log(v / w_convect) * x;
      den += x * x;
    }
    const double alpha = den > 0 ? std::max(num / den, 0.05) : 2.0;

    double rg = r_w * std::pow(fsed_layer[j], 1.0 / alpha) * std::exp(-(alpha + 6.0) * l);
    rg = std::min(std::max(rg, radius_grid.front()), radius_grid.back());

    rg_new[j] = rg;
    r_eff[j] = rg * std::exp(5.0 * l);
    M_new[j] = rho_l[j] * q_c_new[j] * dz[j];        // g/cm^2 of condensate in the layer
    alpha_layer[j] = alpha;
    r_w_layer[j] = r_w;
  }

  // --- while the corrector is far from converged (residual above the gate), keep the committed
  //     structure: a cloud that re-forms on every iteration from a Kzz that is itself still ramping
  //     (the MLT easy-start) keeps the Newton from ever settling (measured: 27 iterations at the
  //     MLT ramp's residual floor, then a stall). The first structure is always taken.
  const bool corrector_settled = residual < 0 || residual < par.gate;
  if (!mass_column.empty() && !warm_start && !corrector_settled)
  {
    ++nb_updates;
    structure_changed = false;
    std::cout << "  A&M " << par.condensate << ": held (corrector residual " << std::scientific
              << std::setprecision(2) << residual << " above the gate " << par.gate << ")"
              << std::defaultfloat << "\n";
    return;
  }

  // --- log-space relaxation, per-update cap and dead band on the committed MASS columns AND
  //     geometric radii (a decaying layer keeps its last radius; a layer the cloud enters takes
  //     the physical radius at once).
  // A layer the cloud enters starts at `seed` of its target and grows by at most the cap per
  // update; a layer the cloud leaves decays by the cap per update instead of vanishing at once.
  // Measured without this: every commit heated the layers below the deck, the base evaporated and
  // re-formed higher, and five commits later the cloud sat at 0.1 bar over a 2100 K inversion.
  // The radius must be relaxed as well: with the physical radius committed at once, the optical
  // depth (~ M/r) is free to jump even though M is capped. Measured on the forsterite + iron brown
  // dwarf with a fixed radiative Kzz of 1e6: a period-two cycle between the deck and the detached
  // convective patches it creates (patch -> Kzz 3e9 -> micron particles -> thin deck -> patch
  // gone -> Kzz 1e6 -> 0.2 um particles -> deck 10x thicker -> ...), tau(1 um) alternating 18 <->
  // 170 between commits until the corrector froze at a ratio residual of 160.
  const double cap = par.max_change > 0 ? std::pow(10.0, par.max_change) : 0.0;
  auto relaxed = [&](double old_value, double new_value, double floor) {
    if (!(new_value > 0))
      return (cap > 0 && old_value > floor) ? old_value / cap : 0.0;          // leaves: decay
    if (!(old_value > 0))
      return new_value * par.seed;                                             // enters: seed
    double v = std::exp(std::log(old_value) + par.relax * (std::log(new_value) - std::log(old_value)));
    if (cap > 0) v = std::min(std::max(v, old_value / cap), old_value * cap);
    return v;
  };

  std::string action;
  if (mass_column.empty())
  {
    mass_column = M_new;
    r_g = rg_new;
    action = "initial";
    structure_changed = true;
  }
  else
  {
    std::vector<double> M_cand(nb_layers), rg_cand(nb_layers);
    double max_change = 0.0;
    const double M_max = std::max(*std::max_element(M_new.begin(), M_new.end()),
                                  *std::max_element(mass_column.begin(), mass_column.end()));
    const double M_floor = 1e-6 * M_max;
    // the dead band judges only layers that matter optically: a tail layer creeping along the
    // 1e-6 floor would otherwise commit a "0.3 dex change" every update and kick the Newton
    // (measured: residual 3e-4 -> 1e-3 per commit for 15 iterations at the end of a run whose
    // tau had long been constant to 1 %)
    const double M_metric = 1e-3 * M_max;

    for (size_t j = 0; j < nb_layers; ++j)
    {
      M_cand[j] = relaxed(mass_column[j], M_new[j], M_floor);
      if (rg_new[j] > 0)
        rg_cand[j] = r_g[j] > 0 ? relaxed(r_g[j], rg_new[j], 0.0) : rg_new[j];
      else
        rg_cand[j] = M_cand[j] > 0 ? r_g[j] : 0.0;

      if (M_floor > 0 && (M_cand[j] > M_metric || mass_column[j] > M_metric))
        max_change = std::max(max_change,
          std::abs(std::log10((M_cand[j] + M_floor) / (mass_column[j] + M_floor))));
      if (r_g[j] > 0 && rg_cand[j] > 0 && (M_cand[j] > M_metric || mass_column[j] > M_metric))
        max_change = std::max(max_change, std::abs(std::log10(rg_cand[j] / r_g[j])));
    }

    if (max_change < par.tolerance)
    {
      action = "held (change " + std::to_string(max_change) + " dex)";
      structure_changed = false;
    }
    else
    {
      mass_column = M_cand;
      r_g = rg_cand;
      action = "committed (change " + std::to_string(max_change) + " dex)";
      structure_changed = true;
    }
  }
  ++nb_updates;

  computeOptics();

  // --- report
  size_t top = nb_layers;
  for (size_t j = nb_layers; j-- > 0;) if (q_c_new[j] > 0) { top = j; break; }
  double tau_ref = 0.0;
  const std::vector<double> tau = opticalDepth(1.0);
  for (const double t : tau) tau_ref += t;

  std::cout << "  A&M " << par.condensate << ": ";
  if (base == nb_layers)
    std::cout << "no condensation";
  else
  {
    double rg_min = 1e99, rg_max = 0.0, qc_max = 0.0;
    for (size_t j = 0; j < nb_layers; ++j)
      if (q_c_new[j] > 0) { rg_min = std::min(rg_min, rg_new[j]); rg_max = std::max(rg_max, rg_new[j]); qc_max = std::max(qc_max, q_c_new[j]); }
    std::cout << std::scientific << std::setprecision(2)
              << "base " << std::sqrt(p[base] * p[base+1]) << " bar, top " << std::sqrt(p[top] * p[top+1])
              << " bar, max q_c " << qc_max << ", r_g " << rg_min * 1e4 << "-" << rg_max * 1e4
              << " um, tau(1um) " << tau_ref << std::defaultfloat;
  }
  std::cout << ", " << action << "\n";
}



void AckermanMarleyCloud::computeOptics()
{
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();
  committed.assignZero(nb_spectral_points, nb_layers);

  const std::vector<double> column_density = columnDensity();

  for (size_t j = 0; j < nb_layers; ++j)
  {
    if (!(column_density[j] > 0) || !(r_g[j] > 0)) continue;

    std::vector<double> radius, number_density;
    lognormalQuadrature(r_g[j], par.sigma, column_density[j], radius, number_density);

    // keep the quadrature inside the table's radius range
    const double r_min = mie_table->radii().front(), r_max = mie_table->radii().back();
    for (auto& r : radius) r = std::min(std::max(r, r_min), r_max);

    std::vector<double> tau_ext, tau_sca, tau_sca_g;
    mie_table->integrate(radius, number_density, tau_ext, tau_sca, tau_sca_g);

    std::vector<double> ssa(tau_ext.size(), 0.0), g(tau_ext.size(), 0.0);
    for (size_t m = 0; m < tau_ext.size(); ++m)
    {
      ssa[m] = tau_ext[m] > 0 ? std::min(tau_sca[m] / tau_ext[m], 1.0) : 0.0;
      g[m] = tau_sca[m] > 0 ? tau_sca_g[m] / tau_sca[m] : 0.0;
    }

    const std::vector<double> tau_grid = mie_table->toSpectralGrid(tau_ext);
    const std::vector<double> ssa_grid = mie_table->toSpectralGrid(ssa);
    const std::vector<double> g_grid = mie_table->toSpectralGrid(g);

    for (size_t nu = 0; nu < nb_spectral_points; ++nu)
    {
      committed.optical_depth[nu][j] = tau_grid[nu];
      committed.single_scattering[nu][j] = ssa_grid[nu];
      committed.asymmetry[nu][j] = g_grid[nu];
    }
  }
}



// column number density per layer from the committed mass column and radius (A&M eq. 18)
std::vector<double> AckermanMarleyCloud::columnDensity() const
{
  std::vector<double> N(nb_layers, 0.0);
  if (mass_column.empty()) return N;
  const double l = 0.5 * std::log(par.sigma) * std::log(par.sigma);
  for (size_t j = 0; j < nb_layers; ++j)
    if (mass_column[j] > 0 && r_g[j] > 0)
      N[j] = 3.0 * mass_column[j] / (4.0 * constants::pi * bulk_density * r_g[j] * r_g[j] * r_g[j])
             * std::exp(-9.0 * l);
  return N;
}



std::map<std::string, std::vector<double>> AckermanMarleyCloud::diagnostics() const
{
  std::vector<double> rg_um(nb_layers, 0.0), reff_um(nb_layers, 0.0), rw_um(nb_layers, 0.0);
  for (size_t j = 0; j < nb_layers; ++j)
  {
    rg_um[j] = (r_g.empty() ? 0.0 : r_g[j]) * 1e4;
    reff_um[j] = r_eff[j] * 1e4;
    rw_um[j] = r_w_layer[j] * 1e4;
  }
  return {
    {"condensate_mmr", q_c},
    {"column_density", columnDensity()},
    {"mass_column", mass_column.empty() ? std::vector<double>(nb_layers, 0.0) : mass_column},
    {"r_g", rg_um},
    {"r_eff", reff_um},
    {"r_w", rw_um},
    {"alpha", alpha_layer},
    {"fsed", fsed_layer}};
}


}
