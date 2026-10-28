/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef SIMPLE_SURFACE_H
#define SIMPLE_SURFACE_H

#include <memory>
#include <vector>

#include "generic_surface.h"
#include "../spectral_grid/spectral_grid.h"


namespace baerli{

//A simple surface model with a constant albedo in the shortwave and zero albedo in the longwave.
//The switch between the two regimes is defined by wavelength_switch (in micron).
class SimpleSurface : public GenericSurface {
  public:
    SimpleSurface(
      SpectralGrid* spectral_grid_,
      const double albedo_shortwave,
      const double wavelength_switch_)
      : GenericSurface(spectral_grid_)
      {
        for (size_t i = 0; i < spectral_grid->nbSpectralPoints(); ++i)
        {
          if (spectral_grid->wavelength_list[i] < wavelength_switch_)
            albedo[i] = albedo_shortwave;
          else
            albedo[i] = 0.0;
        }
      }
    virtual ~SimpleSurface() {}
  protected:
};


} // namespace baerli

#endif // SIMPLE_SURFACE_H
