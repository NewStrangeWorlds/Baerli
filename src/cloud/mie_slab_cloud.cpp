/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#include "mie_slab_cloud.h"
#include "condensate_data.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "../additional/exceptions.h"


namespace ngam {


MieSlabCloud::MieSlabCloud(
  SpectralGrid* spectral_grid_,
  const size_t nb_grid_points,
  const std::string& condensate_,
  const double particle_radius_,
  const double sigma_,
  const double top_pressure_,
  const double bottom_pressure_,
  const double optical_depth,
  const double reference_wavelength,
  const double column_density_,
  const double mie_resolution,
  const size_t mie_radii)
  : CloudModel(spectral_grid_, nb_grid_points)
  , condensate(condensate_)
  , particle_radius(particle_radius_)
  , sigma(sigma_)
  , top_pressure(top_pressure_)
  , bottom_pressure(bottom_pressure_)
{
  type = "mie_slab";

  if (!(particle_radius > 0))
    throw InvalidInput("clouds", "mie_slab: particle_radius must be > 0\n");
  if (!(sigma >= 1.0))
    throw InvalidInput("clouds", "mie_slab: sigma (geometric standard deviation) must be >= 1\n");
  if (!(top_pressure > 0) || !(bottom_pressure > top_pressure))
    throw InvalidInput("clouds", "mie_slab: need 0 < top_pressure < bottom_pressure\n");
  if (!(optical_depth > 0) && !(column_density_ > 0))
    throw InvalidInput("clouds",
      "mie_slab: give either optical_depth (with reference_wavelength) or column_density\n");
  if (optical_depth > 0 && column_density_ > 0)
    throw InvalidInput("clouds", "mie_slab: optical_depth and column_density are alternatives\n");

  mie_table = std::make_unique<MieTable>(
    refractiveIndexFile(condensate), spectral_grid, mie_resolution, mie_radii);

  // optical properties per unit column density: integrate the size distribution on the coarse
  // grid, then move to the spectral grid
  std::vector<double> radius, number_density;
  lognormalQuadrature(particle_radius, sigma, 1.0, radius, number_density);

  std::vector<double> tau_ext, tau_sca, tau_sca_g;
  mie_table->integrate(radius, number_density, tau_ext, tau_sca, tau_sca_g);

  std::vector<double> ssa_coarse(tau_ext.size(), 0.0), g_coarse(tau_ext.size(), 0.0);
  for (size_t m = 0; m < tau_ext.size(); ++m)
  {
    ssa_coarse[m] = tau_ext[m] > 0 ? std::min(tau_sca[m] / tau_ext[m], 1.0) : 0.0;
    g_coarse[m] = tau_sca[m] > 0 ? tau_sca_g[m] / tau_sca[m] : 0.0;
  }

  sigma_ext = mie_table->toSpectralGrid(tau_ext);
  single_scattering = mie_table->toSpectralGrid(ssa_coarse);
  asymmetry = mie_table->toSpectralGrid(g_coarse);

  if (optical_depth > 0)
  {
    // scale to the requested optical depth at the spectral point closest to the reference wavelength
    const std::vector<double>& wl = spectral_grid->wavelength_list;
    size_t best = 0;
    for (size_t k = 1; k < wl.size(); ++k)
      if (std::abs(wl[k] - reference_wavelength) < std::abs(wl[best] - reference_wavelength)) best = k;

    if (!(sigma_ext[best] > 0))
      throw InvalidInput("clouds", "mie_slab: zero extinction at the reference wavelength\n");

    column_density = optical_depth / sigma_ext[best];
  }
  else
    column_density = column_density_;

  std::cout << "- Cloud: Mie slab of " << condensate << ", r_g = " << std::defaultfloat
            << std::setprecision(4) << particle_radius * 1e4 << " micron, sigma = " << sigma
            << ", " << bottom_pressure << " - " << top_pressure << " bar, column density "
            << std::scientific << std::setprecision(3) << column_density << " cm^-2"
            << std::defaultfloat << "\n";
}



void MieSlabCloud::update(const Atmosphere& atmosphere)
{
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  committed.assignZero(nb_spectral_points, nb_layers);

  const std::vector<double> fraction =
    layerFractionsUniformInLnP(atmosphere.pressure, bottom_pressure, top_pressure);

  #pragma omp parallel for
  for (size_t nu = 0; nu < nb_spectral_points; ++nu)
    for (size_t j = 0; j < nb_layers; ++j)
      if (fraction[j] > 0)
      {
        committed.optical_depth[nu][j] = column_density * fraction[j] * sigma_ext[nu];
        committed.single_scattering[nu][j] = single_scattering[nu];
        committed.asymmetry[nu][j] = asymmetry[nu];
      }
}


}
