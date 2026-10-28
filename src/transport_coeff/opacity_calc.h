/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef OPACITY_CALC_H
#define OPACITY_CALC_H


#include <vector>
#include <string>
#include <cstdlib>

#include "transport_coeff.h"
#include "../spectral_grid/spectral_grid.h"
#include "../atmosphere/atmosphere.h"
#include "../additional/quadrature.h"



namespace baerli{


class OpacityCalculation {
  public:
    OpacityCalculation(
      const std::string& cross_section_file_path,
      SpectralGrid* spectral_grid_,
      Atmosphere* atmosphere_,
      const std::vector<std::string>& opacity_species_symbol,
      const std::vector<std::string>& opacity_species_folder)
      : spectral_grid(spectral_grid_)
      , atmosphere(atmosphere_)
      , transport_coeff(
          cross_section_file_path,
          spectral_grid_,
          opacity_species_symbol,
          opacity_species_folder)
      {
        const size_t nb_layers = atmosphere->nb_grid_points - 1;
        cloud_optical_depths.assign(spectral_grid->nbSpectralPoints(), std::vector<double>(nb_layers, 0.0));
        cloud_single_scattering.assign(spectral_grid->nbSpectralPoints(), std::vector<double>(nb_layers, 0.0));
        cloud_asym_param.assign(spectral_grid->nbSpectralPoints(), std::vector<double>(nb_layers, 0.0));
      }
    //gas opacities (absorption and scattering coefficients per level); the cloud arrays below are
    //NOT touched here -- the owning object writes them (GenericObject::applyClouds)
    void calculate();

    std::vector< std::vector<double> > absorption_coeff;
    std::vector< std::vector<double> > scattering_coeff;

    //band-closure correction for the local energy balance (clima_rce_correction):
    //band_correction[b][level] = sum over line species of n_s*(int_band sigma dnu - sampled
    //representation), i.e. the absorption the sampled grid misses (can be negative where the
    //sampling overcounts). Empty when disabled (CLIMA_BANDCORR=0). Units: cm^-1 * cm^-1(nu).
    std::vector< std::vector<double> > band_correction;
    //peak absorption coefficient per band per level [cm^-1] (max over line species): the core
    //optical depth to space built from this decides the escape-probability weight of the closure
    std::vector< std::vector<double> > band_peak_coeff;
    std::vector<double> band_wavenumber;   //band centres [cm^-1]

    //per-LAYER cloud optical properties [spectral point][layer], nb_layers = nb_grid_points-1,
    //layer j between levels j and j+1 (0 = bottom); the sum over the model's cloud modules
    //(cloud/cloud_model.h). Zero without clouds.
    std::vector< std::vector<double> > cloud_optical_depths;
    std::vector< std::vector<double> > cloud_single_scattering;
    std::vector< std::vector<double> > cloud_asym_param;

    //TOTAL absorption coefficient per level [spectral point][level] = gas + cloud absorption, the
    //latter from the layer optical depths above (tau_c (1-w_c) / dz, averaged onto the levels).
    //This is the kappa that weights the local radiative-equilibrium residual kappa (J - B) and its
    //Jacobian in the temperature correctors and the RT backends: with a cloud, local RE is
    //(kappa_gas + kappa_cloud)(J - B) = 0, and weighting with the gas alone drives the Newton to a
    //wrong root (measured: a residual of 1.6 with no admissible step once an Fe deck sat in the
    //deep atmosphere). Equal to absorption_coeff without clouds. Filled by
    //GenericObject::applyClouds after calculate().
    std::vector< std::vector<double> > absorption_coeff_total;

    //absorption_coeff_total from the current gas coefficients and cloud arrays
    void updateTotalAbsorption(const std::vector<double>& altitude)
    {
      const size_t nb_spectral_points = absorption_coeff.size();
      const size_t nb_levels = atmosphere->nb_grid_points;
      absorption_coeff_total = absorption_coeff;

      if (cloud_optical_depths.empty()) return;

      #pragma omp parallel for
      for (size_t nu = 0; nu < nb_spectral_points; ++nu)
      {
        const std::vector<double>& tau = cloud_optical_depths[nu];
        const std::vector<double>& ssa = cloud_single_scattering[nu];
        for (size_t i = 0; i < nb_levels; ++i)
        {
          // (the altitude grid is zero before the first structure calculation, e.g. for a
          // prescribed cloud applied at construction: skip, the loop recomputes it)
          double k = 0.0; int count = 0;
          if (i > 0 && tau[i-1] > 0 && altitude[i] > altitude[i-1])
          { k += tau[i-1] * (1.0 - ssa[i-1]) / (altitude[i] - altitude[i-1]); ++count; }
          if (i + 1 < nb_levels && tau[i] > 0 && altitude[i+1] > altitude[i])
          { k += tau[i] * (1.0 - ssa[i]) / (altitude[i+1] - altitude[i]); ++count; }
          if (count > 0) absorption_coeff_total[nu][i] += k / count;
        }
      }
    }

  private:
    SpectralGrid* spectral_grid;
    Atmosphere* atmosphere;
    TransportCoefficients transport_coeff;
};



inline void OpacityCalculation::calculate()
{
  const size_t nb_grid_points = atmosphere->nb_grid_points;
  const size_t nb_spectral_points = spectral_grid->nbSpectralPoints();

  //band-closure correction setup (DEFAULT ON since the v3 escape-probability weighting;
  //kill switch CLIMA_BANDCORR=0): trapezoid weights of the sampled list and the correction
  //band of each point, built once per grid
  static const bool bandcorr_enabled = [] {
    const char* e = std::getenv("CLIMA_BANDCORR");
    return !(e != nullptr && std::string(e) == "0"); }();

  const size_t nb_bands = bandcorr_enabled ? spectral_grid->nbCorrectionBands() : 0;
  BandCorrectionSpec band_spec;
  std::vector<double> point_weights;
  std::vector<int> band_of_point;

  if (nb_bands > 0)
  {
    point_weights = aux::trapezoidalWeights(spectral_grid->wavenumber_list);
    band_of_point.assign(nb_spectral_points, -1);
    for (size_t k=0; k<nb_spectral_points; ++k)
      band_of_point[k] = static_cast<int>(
        spectral_grid->wavenumber_list[k]/SpectralGrid::correction_band_width);

    band_spec.point_weights = &point_weights;
    band_spec.band_of_point = &band_of_point;
    band_spec.nb_bands = nb_bands;

    band_correction.assign(nb_bands, std::vector<double>(nb_grid_points, 0.0));
    band_peak_coeff.assign(nb_bands, std::vector<double>(nb_grid_points, 0.0));
    band_wavenumber.assign(nb_bands, 0.0);
    for (size_t b=0; b<nb_bands; ++b)
      band_wavenumber[b] = (b + 0.5)*SpectralGrid::correction_band_width;
  }
  else
  {
    band_correction.clear();
    band_peak_coeff.clear();
    band_wavenumber.clear();
  }

  //calculate gas transport coefficients
  absorption_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));
  scattering_coeff.assign(nb_spectral_points, std::vector<double>(nb_grid_points, 0.0));

  for (size_t i=0; i<nb_grid_points; ++i)
  {
    std::vector<double> absorption_coeff_level(nb_spectral_points, 0.0);
    std::vector<double> scattering_coeff_level(nb_spectral_points, 0.0);
    std::vector<double> band_correction_level;
    std::vector<double> band_peak_level;

    transport_coeff.calculate(
      atmosphere->temperature[i],
      atmosphere->pressure[i],
      atmosphere->number_densities[i],
      absorption_coeff_level,
      scattering_coeff_level,
      nb_bands > 0 ? &band_spec : nullptr,
      nb_bands > 0 ? &band_correction_level : nullptr,
      nb_bands > 0 ? &band_peak_level : nullptr);

    for (size_t j=0; j<nb_spectral_points; ++j)
    {
      absorption_coeff[j][i] = absorption_coeff_level[j];
      scattering_coeff[j][i] = scattering_coeff_level[j];
    }

    for (size_t b=0; b<nb_bands; ++b)
    {
      band_correction[b][i] = band_correction_level[b];
      band_peak_coeff[b][i] = band_peak_level[b];
    }
  }
}


}

#endif
