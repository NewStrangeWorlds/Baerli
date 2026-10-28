/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef VAR_ALBEDO_SURFACE_H
#define VAR_ALBEDO_SURFACE_H

#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "generic_surface.h"
#include "../spectral_grid/spectral_grid.h"


namespace baerli{

//A surface model with wavelength-dependent albedo, read from a two-column
//data file (wavelength in micron, reflectance) and interpolated onto the
//spectral grid.
class VariableAlbedoSurface : public GenericSurface {
  public:
    VariableAlbedoSurface(
      SpectralGrid* spectral_grid_,
      const std::string& file_path)
      : GenericSurface(spectral_grid_)
      {
        readFile(file_path);
      }
    virtual ~VariableAlbedoSurface() {}
  protected:
    void readFile(const std::string& file_path) {
      std::ifstream file(file_path);

      if (!file)
        throw std::runtime_error(
          "VariableAlbedoSurface: could not open file: " + file_path);

      std::vector<double> wavelength_data;
      std::vector<double> albedo_data;

      std::string line;

      while (std::getline(file, line))
      {
        if (line.empty() || line[0] == '#')
          continue;

        std::istringstream iss(line);
        double wl, refl;

        if (!(iss >> wl >> refl))
          continue;

        wavelength_data.push_back(wl);
        albedo_data.push_back(refl);
      }

      if (wavelength_data.size() < 2)
        throw std::runtime_error(
          "VariableAlbedoSurface: need at least 2 data points in " + file_path);

      albedo = spectral_grid->interpolateToWavelengthGrid(
        wavelength_data, albedo_data, false, true);
    }
};


} // namespace baerli

#endif // VAR_ALBEDO_SURFACE_H
