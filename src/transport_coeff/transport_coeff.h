/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef TRANSPORT_COEFF_H
#define TRANSPORT_COEFF_H

#include <vector>
#include <memory>

#include "opacity_species.h"
#include "../spectral_grid/spectral_grid.h"


namespace baerli{


class TransportCoefficients {
  public:
    TransportCoefficients(
      const std::string& cross_section_file_path,
      SpectralGrid* grid_ptr,
      const std::vector<std::string>& opacity_species_symbol,
      const std::vector<std::string>& opacity_species_folder);
    void calculate(
      const double temperature,
      const double pressure,
      const std::vector<double>& number_densities,
      std::vector<double>& absorption_coeff,
      std::vector<double>& scattering_coeff,
      const BandCorrectionSpec* band_spec = nullptr,
      std::vector<double>* band_correction = nullptr,
      std::vector<double>* band_peak_coeff = nullptr);

  private:
    std::string cross_section_file_path;
    SpectralGrid* spectral_grid = nullptr;

    std::vector<std::unique_ptr<OpacitySpecies>> gas_species;

    bool addOpacitySpecies(
      const std::string& species_symbol, const std::string& species_folder);
};


}

#endif
