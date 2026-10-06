/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#ifndef _select_cloud_model_h
#define _select_cloud_model_h

#include <memory>
#include <string>
#include <vector>

#include "cloud_model.h"
#include "fixed_cloud.h"
#include "mie_slab_cloud.h"
#include "../config/module_params.h"
#include "../spectral_grid/spectral_grid.h"


namespace ngam {


// Cloud modules, selected by a module spec {type, named parameters}; the model's `clouds` option
// is a LIST of them, their optical properties are summed (cloud_model.h):
//
//   fixed (fix)      prescribed optical properties over (pressure, wavelength), see fixed_cloud.h
//                    file        path to the table (optional; without it the table is handed over
//                                from Python with model.set_cloud_optical_properties, no-op until
//                                then)
//   mie_slab (mie)   a deck of Mie-scattering particles of one condensate between two pressures,
//                    uniform in ln p, lognormal size distribution; see mie_slab_cloud.h
//                    condensate            name from condensate_data.h or a refractive-index file
//                    particle_radius       geometric mean radius [micron]
//                    sigma                 geometric standard deviation, default 1 (monodisperse)
//                    top_pressure, bottom_pressure   [bar]
//                    optical_depth         extinction optical depth of the deck at
//                    reference_wavelength  [micron], default 1 -- OR --
//                    column_density        [cm^-2]
//                    mie_resolution        resolving power of the Mie table's wavelength grid,
//                                          default 200
//                    mie_radii             radii in the Mie table (log-spaced 1e-4 - 542 micron),
//                                          default 120
namespace cloud_modules{
  enum id {fixed, mie_slab};
  const std::vector<std::string> description {"fixed", "mie_slab"};
  const std::vector<std::string> description_short {"fix", "mie"};
}


inline std::unique_ptr<CloudModel> selectCloudModel(
  const ModuleSpec& spec, SpectralGrid* spectral_grid, const size_t nb_grid_points)
{
  const auto module_id = static_cast<cloud_modules::id>(resolveModuleType(
    spec.type, cloud_modules::description, cloud_modules::description_short, "clouds"));

  ParamReader reader(spec, "clouds");
  std::unique_ptr<CloudModel> module;

  switch (module_id)
  {
    case cloud_modules::fixed :
      module = std::make_unique<FixedCloud>(
        spectral_grid, nb_grid_points, reader.getString("file", ""));
      break;

    case cloud_modules::mie_slab :
      module = std::make_unique<MieSlabCloud>(
        spectral_grid, nb_grid_points,
        reader.requireString("condensate"),
        reader.requireDouble("particle_radius") * 1e-4,      // micron -> cm
        reader.getDouble("sigma", 1.0),
        reader.requireDouble("top_pressure"),
        reader.requireDouble("bottom_pressure"),
        reader.getDouble("optical_depth", 0.0),
        reader.getDouble("reference_wavelength", 1.0),
        reader.getDouble("column_density", 0.0),
        reader.getDouble("mie_resolution", 200.0),
        static_cast<size_t>(reader.getInt("mie_radii", 120)));
      break;
  }

  reader.finish();
  module->type = cloud_modules::description[module_id];
  return module;
}


inline std::vector<std::unique_ptr<CloudModel>> selectCloudModels(
  const std::vector<ModuleSpec>& specs, SpectralGrid* spectral_grid, const size_t nb_grid_points)
{
  std::vector<std::unique_ptr<CloudModel>> modules;
  for (const auto& spec : specs)
    modules.push_back(selectCloudModel(spec, spectral_grid, nb_grid_points));
  return modules;
}


}
#endif
