/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _select_radiative_transfer_h
#define _select_radiative_transfer_h

#include <memory>
#include <string>
#include <vector>

#include "radiative_transfer.h"
#include "discrete_ordinate.h"
#include "adding_doubling.h"
#include "../config/module_params.h"


namespace baerli {


// Radiative-transfer solvers, selected by a module spec {type, named parameters}:
//
//   disort            discrete-ordinate solver
//                     nb_streams    total number of streams, default 4
//   adding_doubling   adding-doubling solver with analytic temperature Jacobians
//                     nb_streams    quadrature points per hemisphere, default 2
namespace rt_modules{
  enum id {disort, adding_doubling};
  const std::vector<std::string> description {"disort", "adding_doubling"};
  const std::vector<std::string> description_short {"disort", "ad"};
}


inline std::unique_ptr<RadiativeTransfer> selectRadiativeTransfer(
  const ModuleSpec& spec,
  const size_t nb_grid_points,
  SpectralGrid* spectral_grid)
{
  const auto module_id = static_cast<rt_modules::id>(resolveModuleType(
    spec.type, rt_modules::description, rt_modules::description_short, "radiative_transfer"));

  ParamReader reader(spec, "radiative_transfer");
  std::unique_ptr<RadiativeTransfer> module;

  switch (module_id)
  {
    case rt_modules::disort :
      module = std::make_unique<DiscreteOrdinates>(
        spectral_grid, static_cast<int>(reader.getInt("nb_streams", 4)), nb_grid_points);
      break;

    case rt_modules::adding_doubling :
      module = std::make_unique<AddingDoubling>(
        spectral_grid, static_cast<int>(reader.getInt("nb_streams", 2)), nb_grid_points);
      break;
  }

  reader.finish();
  return module;
}


}
#endif
