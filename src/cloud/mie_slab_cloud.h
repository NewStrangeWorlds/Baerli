/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#ifndef _mie_slab_cloud_h
#define _mie_slab_cloud_h

#include <memory>
#include <string>
#include <vector>

#include "cloud_model.h"
#include "cloud_optics.h"


namespace ngam {


// A prescribed deck of Mie-scattering particles of one condensate between two pressures, with the
// column number density distributed uniformly in ln p. The particle size distribution is lognormal
// (geometric mean radius, geometric standard deviation; sigma = 1 is monodisperse). The column is
// set either directly (column_density [cm^-2]) or through the extinction optical depth at a
// reference wavelength (optical_depth, reference_wavelength) -- the optical depth at the spectral
// point closest to the reference wavelength is then exactly the requested value.
//
// The optical properties come from the condensate's refractive index through a MieTable
// (cloud_optics.h). Since nothing here depends on the atmosphere's state but the pressure grid,
// the per-unit-column optical properties are computed once at construction.
class MieSlabCloud : public CloudModel {
  public:
    MieSlabCloud(
      SpectralGrid* spectral_grid_,
      const size_t nb_grid_points,
      const std::string& condensate,
      const double particle_radius,      // cm
      const double sigma,
      const double top_pressure,         // bar
      const double bottom_pressure,      // bar
      const double optical_depth,        // <= 0: use column_density
      const double reference_wavelength, // micron
      const double column_density,       // cm^-2
      const double mie_resolution,
      const size_t mie_radii);
    virtual ~MieSlabCloud() {}

    void update(const Atmosphere& atmosphere) override;

    double columnDensity() const { return column_density; }
    const MieTable& mieTable() const { return *mie_table; }

  private:
    std::string condensate;
    double particle_radius, sigma, top_pressure, bottom_pressure;
    double column_density = 0.0;

    std::unique_ptr<MieTable> mie_table;

    // per unit column density, on the spectral grid
    std::vector<double> sigma_ext, single_scattering, asymmetry;
};


}
#endif
