/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#include "fixed_cloud.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>

#include "../additional/exceptions.h"


namespace ngam {


FixedCloud::FixedCloud(
  SpectralGrid* spectral_grid_, const size_t nb_grid_points, const std::string& file)
  : CloudModel(spectral_grid_, nb_grid_points), file_name(file)
{
  type = "fixed";

  if (!file_name.empty()) readFile(file_name);
}



void FixedCloud::setTable(
  const std::vector<double>& pressure,
  const std::vector<double>& wavelength,
  const std::vector<std::vector<double>>& optical_depth,
  const std::vector<std::vector<double>>& single_scattering,
  const std::vector<std::vector<double>>& asymmetry)
{
  const size_t np = pressure.size();
  const size_t nw = wavelength.size();

  if (np < 2 || nw < 1)
    throw InvalidInput("fixed cloud",
      "the table needs at least two pressure levels and one wavelength\n");

  auto checkShape = [&](const std::vector<std::vector<double>>& t, const char* name) {
    if (t.size() != np)
      throw InvalidInput("fixed cloud",
        std::string(name) + " must have one row per pressure level\n");
    for (const auto& row : t)
      if (row.size() != nw)
        throw InvalidInput("fixed cloud",
          std::string(name) + " must have one column per wavelength\n");
  };
  checkShape(optical_depth, "optical_depth");
  checkShape(single_scattering, "single_scattering");
  checkShape(asymmetry, "asymmetry");

  for (const auto& row : optical_depth)
    for (const double t : row)
      if (!(t >= 0))
        throw InvalidInput("fixed cloud", "optical depths must be >= 0\n");
  for (const auto& row : single_scattering)
    for (const double w : row)
      if (!(w >= 0 && w <= 1))
        throw InvalidInput("fixed cloud", "single-scattering albedos must be in [0, 1]\n");
  for (const auto& row : asymmetry)
    for (const double g : row)
      if (!(g >= -1 && g <= 1))
        throw InvalidInput("fixed cloud", "asymmetry parameters must be in [-1, 1]\n");

  // sort: pressure decreasing (bottom -> top), wavelength increasing
  std::vector<size_t> ip(np), iw(nw);
  std::iota(ip.begin(), ip.end(), 0);
  std::iota(iw.begin(), iw.end(), 0);
  std::sort(ip.begin(), ip.end(), [&](size_t a, size_t b) { return pressure[a] > pressure[b]; });
  std::sort(iw.begin(), iw.end(), [&](size_t a, size_t b) { return wavelength[a] < wavelength[b]; });

  for (size_t k = 1; k < np; ++k)
    if (pressure[ip[k]] == pressure[ip[k-1]])
      throw InvalidInput("fixed cloud", "pressure levels must be distinct\n");
  for (size_t m = 1; m < nw; ++m)
    if (wavelength[iw[m]] == wavelength[iw[m-1]])
      throw InvalidInput("fixed cloud", "wavelengths must be distinct\n");

  table_pressure.resize(np);
  table_wavelength.resize(nw);
  table_tau.assign(np, std::vector<double>(nw, 0.0));
  table_ssa.assign(np, std::vector<double>(nw, 0.0));
  table_g.assign(np, std::vector<double>(nw, 0.0));

  for (size_t k = 0; k < np; ++k)
  {
    table_pressure[k] = pressure[ip[k]];
    for (size_t m = 0; m < nw; ++m)
    {
      table_tau[k][m] = optical_depth[ip[k]][iw[m]];
      table_ssa[k][m] = single_scattering[ip[k]][iw[m]];
      table_g[k][m] = asymmetry[ip[k]][iw[m]];
    }
  }
  for (size_t m = 0; m < nw; ++m) table_wavelength[m] = wavelength[iw[m]];

  // the top row carries no layer
  for (size_t m = 0; m < nw; ++m) table_tau[np-1][m] = 0.0;
}



void FixedCloud::readFile(const std::string& file)
{
  std::fstream in(file.c_str(), std::ios::in);

  if (in.fail())
    throw InvalidInput("fixed cloud", "unable to open cloud table " + file + "\n");

  struct Row { double p, wl, tau, w, g; };
  std::vector<Row> rows;

  std::string line;
  while (std::getline(in, line))
  {
    const size_t hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    std::istringstream s(line);
    Row r;
    if (s >> r.p >> r.wl >> r.tau >> r.w >> r.g) rows.push_back(r);
  }

  if (rows.empty())
    throw InvalidInput("fixed cloud", "no data rows in cloud table " + file + "\n");

  std::vector<double> pressure, wavelength;
  for (const auto& r : rows)
  {
    if (std::find(pressure.begin(), pressure.end(), r.p) == pressure.end()) pressure.push_back(r.p);
    if (std::find(wavelength.begin(), wavelength.end(), r.wl) == wavelength.end())
      wavelength.push_back(r.wl);
  }

  const size_t np = pressure.size(), nw = wavelength.size();
  if (rows.size() != np * nw)
    throw InvalidInput("fixed cloud",
      "cloud table " + file + " must cover the full (pressure, wavelength) grid: found "
      + std::to_string(rows.size()) + " rows for " + std::to_string(np) + " pressures x "
      + std::to_string(nw) + " wavelengths\n");

  std::vector<std::vector<double>> tau(np, std::vector<double>(nw, 0.0)), ssa = tau, g = tau;

  for (const auto& r : rows)
  {
    const size_t k = std::find(pressure.begin(), pressure.end(), r.p) - pressure.begin();
    const size_t m = std::find(wavelength.begin(), wavelength.end(), r.wl) - wavelength.begin();
    tau[k][m] = r.tau;
    ssa[k][m] = r.w;
    g[k][m] = r.g;
  }

  setTable(pressure, wavelength, tau, ssa, g);

  std::cout << "\n  Fixed cloud table " << file << ": " << np << " pressure levels, "
            << nw << " wavelengths\n";
}



namespace {

// linear interpolation in x on a sorted (increasing) abscissa, constant outside the range
double interpolateClamped(const std::vector<double>& x, const std::vector<double>& y, const double xq)
{
  if (xq <= x.front()) return y.front();
  if (xq >= x.back()) return y.back();

  const size_t hi = std::upper_bound(x.begin(), x.end(), xq) - x.begin();
  const size_t lo = hi - 1;
  const double f = (xq - x[lo]) / (x[hi] - x[lo]);

  return y[lo] + f * (y[hi] - y[lo]);
}

}



void FixedCloud::update(const Atmosphere& atmosphere)
{
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  committed.assignZero(nb_spectral_points, nb_layers);

  if (!hasTable()) return;

  const size_t np = table_pressure.size();
  const std::vector<double>& wl = spectral_grid->wavelength_list;

  // 1. interpolate the table in wavelength onto the spectral grid (constant outside its range)
  std::vector<std::vector<double>> tau(np, std::vector<double>(nb_spectral_points, 0.0));
  std::vector<std::vector<double>> ssa = tau, g = tau;

  for (size_t k = 0; k < np; ++k)
    for (size_t nu = 0; nu < nb_spectral_points; ++nu)
    {
      tau[k][nu] = interpolateClamped(table_wavelength, table_tau[k], wl[nu]);
      ssa[k][nu] = interpolateClamped(table_wavelength, table_ssa[k], wl[nu]);
      g[k][nu] = interpolateClamped(table_wavelength, table_g[k], wl[nu]);
    }

  // 2. abscissae in -ln p (increasing bottom -> top): the table's levels and layer midpoints, the
  //    model's levels and layer midpoints
  std::vector<double> x_table(np), x_table_mid(np - 1);
  for (size_t k = 0; k < np; ++k) x_table[k] = -std::log(table_pressure[k]);
  for (size_t k = 0; k + 1 < np; ++k) x_table_mid[k] = 0.5 * (x_table[k] + x_table[k+1]);

  const std::vector<double>& p = atmosphere.pressure;
  std::vector<double> x_model(p.size()), x_model_mid(nb_layers);
  for (size_t i = 0; i < p.size(); ++i) x_model[i] = -std::log(p[i]);
  for (size_t j = 0; j < nb_layers; ++j) x_model_mid[j] = 0.5 * (x_model[j] + x_model[j+1]);

  // 3. per spectral point: cumulative optical depth from the bottom, interpolated onto the model
  //    levels and differenced; w and g interpolated at the layer midpoints
  std::vector<double> cumulative(np), ssa_mid(np - 1), g_mid(np - 1), cum_model(p.size());

  for (size_t nu = 0; nu < nb_spectral_points; ++nu)
  {
    cumulative[0] = 0.0;
    for (size_t k = 1; k < np; ++k) cumulative[k] = cumulative[k-1] + tau[k-1][nu];
    for (size_t k = 0; k + 1 < np; ++k) { ssa_mid[k] = ssa[k][nu]; g_mid[k] = g[k][nu]; }

    for (size_t i = 0; i < p.size(); ++i)
      cum_model[i] = interpolateClamped(x_table, cumulative, x_model[i]);

    for (size_t j = 0; j < nb_layers; ++j)
    {
      const double t = std::max(0.0, cum_model[j+1] - cum_model[j]);
      committed.optical_depth[nu][j] = t;
      if (t > 0)
      {
        committed.single_scattering[nu][j] = interpolateClamped(x_table_mid, ssa_mid, x_model_mid[j]);
        committed.asymmetry[nu][j] = interpolateClamped(x_table_mid, g_mid, x_model_mid[j]);
      }
    }
  }
}


}
