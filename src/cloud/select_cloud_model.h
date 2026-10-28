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
#include "ackerman_marley_cloud.h"
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
//   ackerman_marley (am)   Ackerman & Marley (2001) eddy-sedimentation cloud of one condensate on
//                    the condensation increments of the equilibrium chemistry (needs
//                    condensation = rainout there); see ackerman_marley_cloud.h
//                    condensate            name from condensate_data.h with a FastChem counterpart
//                    fsed                  sedimentation efficiency, default 1
//                    fsed_mode             const (default) | exp (Rooney et al. 2022, with b, eps)
//                    b, eps                default 1 (scale heights), 1e-2
//                    sigma                 lognormal geometric standard deviation, default 2
//                    kzz_min               floor on Kzz [cm^2/s], default 1e5
//                    mixing_length         scale_height (default) | lapse_ratio
//                    relax                 log-space relaxation weight of the structure between
//                                          outer iterations, default 0.5 (1 = none)
//                    max_change            cap on the change of column density and r_g per outer
//                                          iteration [dex], default 0.3 (0 = none)
//                    tolerance             dead band [dex] on column density and r_g, default 0.05
//                    gate                  the structure is recomputed only while the corrector's
//                                          residual is below this, default 1e-2
//                    seed                  a layer the cloud enters starts at this fraction of its
//                                          target, default 1e-2
//                    (the onset of the cloud opacity itself is the model's `cloud_ramp` option)
//                    bulk_density          g/cm^3, default from the condensate table
//                    molecule_diameter, lj_temperature   gas properties for the viscosity, default
//                                          H2 (2.827e-8 cm, 59.7 K)
//                    mie_resolution, mie_radii   as for mie_slab
namespace cloud_modules{
  enum id {fixed, mie_slab, ackerman_marley};
  const std::vector<std::string> description {"fixed", "mie_slab", "ackerman_marley"};
  const std::vector<std::string> description_short {"fix", "mie", "am"};
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

    case cloud_modules::ackerman_marley :
    {
      AckermanMarleyCloud::Parameters par;
      par.condensate = reader.requireString("condensate");
      par.fsed = reader.getDouble("fsed", 1.0);
      const std::string fsed_mode = reader.getString("fsed_mode", "const");
      if (fsed_mode != "const" && fsed_mode != "exp")
        throw InvalidInput("clouds", "ackerman_marley: fsed_mode must be const or exp\n");
      par.fsed_exponential = (fsed_mode == "exp");
      par.b = reader.getDouble("b", 1.0);
      par.eps = reader.getDouble("eps", 1e-2);
      par.sigma = reader.getDouble("sigma", 2.0);
      par.kzz_min = reader.getDouble("kzz_min", 1e5);
      const std::string mixing_length = reader.getString("mixing_length", "scale_height");
      if (mixing_length != "scale_height" && mixing_length != "lapse_ratio")
        throw InvalidInput("clouds", "ackerman_marley: mixing_length must be scale_height or lapse_ratio\n");
      par.lapse_ratio_mixing_length = (mixing_length == "lapse_ratio");
      par.relax = reader.getDouble("relax", 0.5);
      par.max_change = reader.getDouble("max_change", 0.3);
      par.tolerance = reader.getDouble("tolerance", 0.05);
      par.gate = reader.getDouble("gate", 1e-2);
      par.seed = reader.getDouble("seed", 1e-2);
      par.bulk_density = reader.getDouble("bulk_density", 0.0);
      par.molecule_diameter = reader.getDouble("molecule_diameter", 2.827e-8);
      par.lj_temperature = reader.getDouble("lj_temperature", 59.7);
      par.mie_resolution = reader.getDouble("mie_resolution", 200.0);
      par.mie_radii = static_cast<size_t>(reader.getInt("mie_radii", 120));
      module = std::make_unique<AckermanMarleyCloud>(spectral_grid, nb_grid_points, par);
      break;
    }
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
