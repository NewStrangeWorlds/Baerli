/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#include "fall_velocity.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "../additional/physical_const.h"


namespace baerli {


FallMedium fallMedium(
  const double temperature, const double pressure_dyn, const double mean_molecular_weight,
  const double gravity, const double molecule_diameter, const double lj_temperature)
{
  FallMedium m;
  m.temperature = temperature;
  m.pressure = pressure_dyn;
  m.mean_molecular_weight = mean_molecular_weight;
  m.gravity = gravity;

  const double n = pressure_dyn / (constants::boltzmann_k * temperature);
  const double d2 = molecule_diameter * molecule_diameter;
  const double mass = mean_molecular_weight * constants::amu;

  m.mean_free_path = 1.0 / (std::sqrt(2.0) * n * constants::pi * d2);
  m.viscosity = 5.0 / 16.0 * std::sqrt(constants::pi * constants::boltzmann_k * temperature * mass)
                / (constants::pi * d2) / (1.22 * std::pow(temperature / lj_temperature, -0.16));

  return m;
}



double fallVelocity(const double radius, const double bulk_density, const FallMedium& medium)
{
  constexpr double cdrag = 0.45;

  const double rho_atmos = medium.pressure
    / ((constants::gas_constant / medium.mean_molecular_weight) * medium.temperature);
  const double drho = bulk_density - rho_atmos;

  const double knudsen = medium.mean_free_path / radius;
  const double beta_slip = 1.0 + 1.26 * knudsen;

  double vfall = beta_slip * (2.0 / 9.0) * drho * medium.gravity * radius * radius / medium.viscosity;

  double reynolds = 2.0 * radius * rho_atmos * vfall / medium.viscosity;

  if (reynolds > 1.0 && reynolds <= 1e3)
  {
    const double cd_nre2 = 32.0 * radius * radius * radius * drho * rho_atmos * medium.gravity
                           / (3.0 * medium.viscosity * medium.viscosity);
    const double x = std::log(cd_nre2);
    const double b0 = -0.318657e1, b1 = 0.992696, b2 = -0.153193e-2, b3 = -0.987059e-3,
                 b4 = -0.578878e-3, b5 = 0.855176e-4, b6 = -0.327815e-5;
    const double y = b0 + b1*x + b2*x*x + b3*x*x*x + b4*x*x*x*x + b5*x*x*x*x*x + b6*x*x*x*x*x*x;

    reynolds = std::exp(y);
    vfall = medium.viscosity * reynolds / (2.0 * radius * rho_atmos);
  }

  if (reynolds > 1e3)
    vfall = beta_slip * std::sqrt(8.0 * drho * radius * medium.gravity / (3.0 * cdrag * rho_atmos));

  return vfall;
}



double brentRoot(
  const std::function<double(double)>& f, double a, double b,
  const double rel_tol, const int max_iter)
{
  double fa = f(a), fb = f(b);

  if (fa * fb > 0)
    throw std::runtime_error("brentRoot: the root is not bracketed");

  if (std::abs(fa) < std::abs(fb)) { std::swap(a, b); std::swap(fa, fb); }

  double c = a, fc = fa, d = b - a;
  bool bisected = true;

  for (int it = 0; it < max_iter; ++it)
  {
    if (fb == 0.0 || std::abs(b - a) <= rel_tol * std::abs(b)) return b;

    double s;
    if (fa != fc && fb != fc)
      s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc))
          + c * fa * fb / ((fc - fa) * (fc - fb));
    else
      s = b - fb * (b - a) / (fb - fa);

    const double lo = std::min((3.0 * a + b) / 4.0, b), hi = std::max((3.0 * a + b) / 4.0, b);
    const bool cond1 = !(s > lo && s < hi);
    const bool cond2 = bisected && std::abs(s - b) >= std::abs(b - c) / 2.0;
    const bool cond3 = !bisected && std::abs(s - b) >= std::abs(c - d) / 2.0;
    const bool cond4 = bisected && std::abs(b - c) < rel_tol * std::abs(b);
    const bool cond5 = !bisected && std::abs(c - d) < rel_tol * std::abs(b);

    if (cond1 || cond2 || cond3 || cond4 || cond5)
    {
      s = 0.5 * (a + b);
      bisected = true;
    }
    else
      bisected = false;

    const double fs = f(s);
    d = c; c = b; fc = fb;

    if (fa * fs < 0) { b = s; fb = fs; }
    else             { a = s; fa = fs; }

    if (std::abs(fa) < std::abs(fb)) { std::swap(a, b); std::swap(fa, fb); }
  }

  return b;
}


}
