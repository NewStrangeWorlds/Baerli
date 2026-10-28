/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#include "spectral_grid.h"


#include <iostream>
#include <string>
#include <fstream>
#include <vector>

namespace baerli{


//convert wavenumbers in cm-1 to wavelengths in microns
std::vector<double> SpectralGrid::wavenumberToWavelength(
  const std::vector<double>& wavenumbers)
{
  const size_t nb_wavenumbers = wavenumbers.size();

  if (nb_wavenumbers == 0) return std::vector<double>(0,0);

  std::vector<double> wavelengths(nb_wavenumbers, 0);

  for (size_t i=0; i<nb_wavenumbers; ++i)
    wavelengths[i] = wavenumberToWavelength(wavenumbers[i]);

  return wavelengths;
}



std::vector<double> SpectralGrid::wavelengthToWavenumber(
  const std::vector<double>& wavelengths)
{
  size_t nb_wavelengths = wavelengths.size();

   if (nb_wavelengths == 0) return std::vector<double>(0,0);

  std::vector<double> wavenumbers(nb_wavelengths, 0);

  for (size_t i=0; i<nb_wavelengths; ++i)
    wavenumbers[i] = wavelengthToWavenumber(wavelengths[i]);

  return wavenumbers;
}




}

