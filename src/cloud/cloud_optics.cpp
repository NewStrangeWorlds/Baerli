/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#include "cloud_optics.h"
#include "condensate_data.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <mie.h>   // LX-MIE

#include "../additional/exceptions.h"
#include "../additional/physical_const.h"


namespace ngam {


#ifndef NGAM_DATA_DIR
#define NGAM_DATA_DIR "data"
#endif


std::string refractiveIndexFile(const std::string& condensate)
{
  if (const CondensateData* c = findCondensate(condensate))
    return std::string(NGAM_DATA_DIR) + "/refractive_indices/" + c->refractive_index_file;

  return condensate;
}



namespace {

// linear interpolation on an increasing abscissa, constant outside the range
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



MieTable::MieTable(
  const std::string& refractive_index_file,
  SpectralGrid* spectral_grid,
  const double mie_resolution,
  const size_t nb_radii,
  const double r_min,
  const double r_max,
  const unsigned int sub_bins)
  : file_name(refractive_index_file)
{
  if (mie_resolution <= 0 || nb_radii < 2 || r_min <= 0 || r_max <= r_min || sub_bins < 1)
    throw InvalidInput("MieTable", "invalid Mie table parameters\n");

  // the coarse wavelength grid: log-spaced at the given resolving power over the model's range
  const double wl_min = spectral_grid->wavelength_list.back();    // the list is decreasing
  const double wl_max = spectral_grid->wavelength_list.front();
  const size_t nb_wavelengths =
    static_cast<size_t>(std::ceil(mie_resolution * std::log(wl_max / wl_min))) + 1;

  coarse_wavelength.resize(nb_wavelengths);
  for (size_t m = 0; m < nb_wavelengths; ++m)
    coarse_wavelength[m] = wl_min * std::exp(
      std::log(wl_max / wl_min) * static_cast<double>(m) / static_cast<double>(nb_wavelengths - 1));

  // the radius grid (virga's get_r_grid)
  radius_grid.resize(nb_radii);
  log_radius.resize(nb_radii);
  for (size_t i = 0; i < nb_radii; ++i)
  {
    radius_grid[i] = r_min * std::pow(r_max / r_min, static_cast<double>(i) / static_cast<double>(nb_radii - 1));
    log_radius[i] = std::log(radius_grid[i]);
  }

  // interpolation weights spectral grid <- coarse grid
  const std::vector<double>& wl = spectral_grid->wavelength_list;
  grid_index.resize(wl.size());
  grid_weight.resize(wl.size());
  for (size_t k = 0; k < wl.size(); ++k)
  {
    const double x = std::log(std::min(std::max(wl[k], wl_min), wl_max) / wl_min)
                     / std::log(wl_max / wl_min) * static_cast<double>(nb_wavelengths - 1);
    size_t lo = static_cast<size_t>(std::floor(x));
    if (lo >= nb_wavelengths - 1) lo = nb_wavelengths - 2;
    grid_index[k] = lo;
    grid_weight[k] = std::min(std::max(x - static_cast<double>(lo), 0.0), 1.0);
  }

  std::vector<double> wavelength_data, n_data, k_data;
  readRefractiveIndex(refractive_index_file, wavelength_data, n_data, k_data);

  const auto start = std::chrono::steady_clock::now();
  buildTable(wavelength_data, n_data, k_data, sub_bins);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

  std::cout << "- Mie table " << refractive_index_file << ": " << nb_radii << " radii ("
            << std::scientific << std::setprecision(1) << r_min * 1e4 << " - " << r_max * 1e4
            << " micron), " << nb_wavelengths << " wavelengths (R = " << std::defaultfloat
            << mie_resolution << ", " << std::setprecision(3) << wl_min << " - " << wl_max
            << " micron; data " << wavelength_data.front() << " - " << wavelength_data.back()
            << " micron), " << std::setprecision(2) << seconds << " s\n";

  if (wavelength_data.front() > wl_min || wavelength_data.back() < wl_max)
    std::cout << "  NOTE: refractive index held constant outside the data range\n";
}



void MieTable::readRefractiveIndex(
  const std::string& file, std::vector<double>& wavelength,
  std::vector<double>& n, std::vector<double>& k) const
{
  std::fstream in(file.c_str(), std::ios::in);

  if (in.fail())
    throw InvalidInput("MieTable",
      "unable to open refractive index file " + file + "\n  (known condensates: "
      + condensateNames() + ")\n");

  std::string line;
  while (std::getline(in, line))
  {
    const size_t hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    std::istringstream s(line);
    double wl, nr, ni;
    if (s >> wl >> nr >> ni)
    {
      wavelength.push_back(wl);
      n.push_back(nr);
      k.push_back(std::abs(ni));
    }
  }

  if (wavelength.size() < 2)
    throw InvalidInput("MieTable", "refractive index file " + file + " has fewer than two rows\n");

  // sort by wavelength
  std::vector<size_t> idx(wavelength.size());
  for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return wavelength[a] < wavelength[b]; });

  std::vector<double> w2, n2, k2;
  for (const size_t i : idx) { w2.push_back(wavelength[i]); n2.push_back(n[i]); k2.push_back(k[i]); }
  wavelength = w2; n = n2; k = k2;
}



void MieTable::buildTable(
  const std::vector<double>& wavelength_data, const std::vector<double>& n_data,
  const std::vector<double>& k_data, const unsigned int sub_bins)
{
  const size_t nb_wavelengths = coarse_wavelength.size();
  const size_t nb_radii = radius_grid.size();

  // refractive index on the coarse grid: n linear in ln lambda, k linear in ln lambda of ln k
  // (k spans decades); a zero k is kept linear
  std::vector<double> log_wl_data(wavelength_data.size());
  for (size_t i = 0; i < wavelength_data.size(); ++i) log_wl_data[i] = std::log(wavelength_data[i]);

  std::vector<double> n_coarse(nb_wavelengths), k_coarse(nb_wavelengths);
  for (size_t m = 0; m < nb_wavelengths; ++m)
  {
    const double x = std::log(coarse_wavelength[m]);
    n_coarse[m] = interpolateClamped(log_wl_data, n_data, x);

    if (x <= log_wl_data.front()) { k_coarse[m] = k_data.front(); continue; }
    if (x >= log_wl_data.back())  { k_coarse[m] = k_data.back(); continue; }

    const size_t hi = std::upper_bound(log_wl_data.begin(), log_wl_data.end(), x) - log_wl_data.begin();
    const size_t lo = hi - 1;
    const double f = (x - log_wl_data[lo]) / (log_wl_data[hi] - log_wl_data[lo]);
    if (k_data[lo] > 0 && k_data[hi] > 0)
      k_coarse[m] = std::exp(std::log(k_data[lo]) + f * (std::log(k_data[hi]) - std::log(k_data[lo])));
    else
      k_coarse[m] = k_data[lo] + f * (k_data[hi] - k_data[lo]);
  }

  // bin edges (virga's get_r_grid: linearly centred bins with log-spaced widths)
  std::vector<double> dr(nb_radii);
  for (size_t i = 1; i + 1 < nb_radii; ++i)
    dr[i] = 2.0 * (radius_grid[i+1] - radius_grid[i]) * (radius_grid[i] - radius_grid[i-1])
            / (radius_grid[i+1] - radius_grid[i-1]);
  dr[0] = 2.0 * (radius_grid[1] - radius_grid[0]) - dr[1];
  dr[nb_radii-1] = 2.0 * (radius_grid[nb_radii-1] - radius_grid[nb_radii-2]) - dr[nb_radii-2];

  q_ext.assign(nb_radii, std::vector<double>(nb_wavelengths, 0.0));
  q_sca.assign(nb_radii, std::vector<double>(nb_wavelengths, 0.0));
  g.assign(nb_radii, std::vector<double>(nb_wavelengths, 0.0));

  #pragma omp parallel for schedule(dynamic, 1) collapse(2)
  for (size_t i = 0; i < nb_radii; ++i)
    for (size_t m = 0; m < nb_wavelengths; ++m)
    {
      // LX-MIE works with m = n - ik internally and (since 38cc6d9) normalises the sign of the
      // imaginary part itself, so the compilation's positive k can be passed as is
      const std::complex<double> refractive_index(n_coarse[m], k_coarse[m]);
      const double r_lo = std::max(radius_grid[i] - 0.5 * dr[i], 0.5 * radius_grid[i]);
      const double r_hi = radius_grid[i] + 0.5 * dr[i];

      double sum_ext = 0, sum_sca = 0, sum_sca_g = 0;

      for (unsigned int s = 0; s < sub_bins; ++s)
      {
        const double r = sub_bins == 1
          ? radius_grid[i]
          : r_lo + (r_hi - r_lo) * static_cast<double>(s) / static_cast<double>(sub_bins - 1);
        const double size_parameter = 2.0 * constants::pi * r / (coarse_wavelength[m] * 1e-4);

        const lxmie::MieResult mie = lxmie::Mie(refractive_index, size_parameter);

        sum_ext += mie.q_ext;
        sum_sca += mie.q_sca;
        sum_sca_g += mie.q_sca * mie.asymmetry_parameter;
      }

      q_ext[i][m] = sum_ext / sub_bins;
      q_sca[i][m] = std::min(sum_sca / sub_bins, q_ext[i][m]);
      g[i][m] = sum_sca > 0 ? sum_sca_g / sum_sca : 0.0;
    }
}



void MieTable::radiusWeights(const double radius, size_t& lower, double& weight) const
{
  const double x = std::log(radius);
  const size_t n = log_radius.size();

  if (x <= log_radius.front()) { lower = 0; weight = 0.0; return; }
  if (x >= log_radius.back())  { lower = n - 2; weight = 1.0; return; }

  const size_t hi = std::upper_bound(log_radius.begin(), log_radius.end(), x) - log_radius.begin();
  lower = hi - 1;
  weight = (x - log_radius[lower]) / (log_radius[hi] - log_radius[lower]);
}



void MieTable::integrate(
  const std::vector<double>& radius,
  const std::vector<double>& number_density,
  std::vector<double>& tau_ext,
  std::vector<double>& tau_sca,
  std::vector<double>& tau_sca_g) const
{
  const size_t nb_wavelengths = coarse_wavelength.size();

  tau_ext.assign(nb_wavelengths, 0.0);
  tau_sca.assign(nb_wavelengths, 0.0);
  tau_sca_g.assign(nb_wavelengths, 0.0);

  if (radius.size() != number_density.size())
    throw InvalidInput("MieTable::integrate", "radius and number density lists differ in size\n");

  for (size_t j = 0; j < radius.size(); ++j)
  {
    if (number_density[j] <= 0) continue;

    if (radius[j] < radius_grid.front() || radius[j] > radius_grid.back())
      throw InvalidInput("MieTable::integrate",
        "particle radius " + std::to_string(radius[j] * 1e4) + " micron outside the Mie table ("
        + std::to_string(radius_grid.front() * 1e4) + " - "
        + std::to_string(radius_grid.back() * 1e4) + " micron)\n");

    size_t lo; double f;
    radiusWeights(radius[j], lo, f);

    const double area = constants::pi * radius[j] * radius[j] * number_density[j];

    for (size_t m = 0; m < nb_wavelengths; ++m)
    {
      // efficiencies interpolated as ln Q in ln r (exact for the Rayleigh power laws, flat in
      // the geometric limit), the asymmetry parameter linearly
      const double qe = std::exp((1.0 - f) * std::log(q_ext[lo][m]) + f * std::log(q_ext[lo+1][m]));
      const double qs = (q_sca[lo][m] > 0 && q_sca[lo+1][m] > 0)
        ? std::exp((1.0 - f) * std::log(q_sca[lo][m]) + f * std::log(q_sca[lo+1][m]))
        : (1.0 - f) * q_sca[lo][m] + f * q_sca[lo+1][m];
      const double qsg = qs * ((1.0 - f) * g[lo][m] + f * g[lo+1][m]);

      tau_ext[m] += area * qe;
      tau_sca[m] += area * qs;
      tau_sca_g[m] += area * qsg;
    }
  }
}



std::vector<double> MieTable::toSpectralGrid(const std::vector<double>& coarse) const
{
  std::vector<double> out(grid_index.size());

  for (size_t k = 0; k < out.size(); ++k)
    out[k] = (1.0 - grid_weight[k]) * coarse[grid_index[k]] + grid_weight[k] * coarse[grid_index[k] + 1];

  return out;
}



double MieTable::atWavelength(const std::vector<double>& coarse, const double wavelength) const
{
  std::vector<double> log_wl(coarse_wavelength.size());
  for (size_t m = 0; m < log_wl.size(); ++m) log_wl[m] = std::log(coarse_wavelength[m]);

  return interpolateClamped(log_wl, coarse, std::log(wavelength));
}



void lognormalQuadrature(
  const double r_g, const double sigma, const double column_density,
  std::vector<double>& radius, std::vector<double>& number_density,
  const size_t nb_points)
{
  if (!(r_g > 0) || !(sigma >= 1.0))
    throw InvalidInput("lognormalQuadrature", "need r_g > 0 and sigma >= 1\n");

  if (sigma < 1.0 + 1e-6 || nb_points < 3)
  {
    radius = {r_g};
    number_density = {column_density};
    return;
  }

  const double ln_sigma = std::log(sigma);
  const double half_width = 4.0 * ln_sigma;
  const double dx = 2.0 * half_width / static_cast<double>(nb_points - 1);

  radius.resize(nb_points);
  number_density.resize(nb_points);

  double norm = 0.0;
  for (size_t j = 0; j < nb_points; ++j)
  {
    const double x = -half_width + dx * static_cast<double>(j);          // ln(r / r_g)
    radius[j] = r_g * std::exp(x);
    double w = std::exp(-0.5 * x * x / (ln_sigma * ln_sigma));
    if (j == 0 || j + 1 == nb_points) w *= 0.5;                          // trapezoid
    number_density[j] = w;
    norm += w;
  }

  for (auto& n : number_density) n *= column_density / norm;
}



std::vector<double> layerFractionsUniformInLnP(
  const std::vector<double>& pressure, const double p_bottom, const double p_top)
{
  const size_t nb_layers = pressure.size() - 1;
  std::vector<double> fraction(nb_layers, 0.0);

  if (!(p_top > 0) || !(p_bottom > p_top)) return fraction;

  const double x_bottom = -std::log(p_bottom), x_top = -std::log(p_top);
  const double width = x_top - x_bottom;

  for (size_t j = 0; j < nb_layers; ++j)
  {
    const double x_lo = -std::log(pressure[j]), x_hi = -std::log(pressure[j+1]);
    const double overlap = std::min(x_hi, x_top) - std::max(x_lo, x_bottom);
    if (overlap > 0) fraction[j] = overlap / width;
  }

  return fraction;
}


}
