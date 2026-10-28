/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _select_stellar_spectrum_h
#define _select_stellar_spectrum_h

#include <memory>
#include <string>
#include <vector>

#include "stellar_spectrum.h"
#include "../config/module_params.h"


namespace baerli {


// Host-star spectra, selected by a module spec {type, named parameters}. The spectrum's shape
// comes from the module; its normalisation (the instellation at the planet) and the incidence
// geometry are properties of the planet and are passed separately.
//
//   blackbody    temperature   stellar effective temperature in K (required)
//   tabulated    file          path to the spectrum file (required)
namespace stellar_modules{
  enum id {blackbody, tabulated};
  const std::vector<std::string> description {"blackbody", "tabulated"};
}


inline std::unique_ptr<StellarSpectrum> selectStellarSpectrum(
  const ModuleSpec& spec,
  const double instellation)
{
  const auto module_id = static_cast<stellar_modules::id>(resolveModuleType(
    spec.type, stellar_modules::description, {}, "stellar_spectrum"));

  ParamReader reader(spec, "stellar_spectrum");
  std::unique_ptr<StellarSpectrum> module;

  switch (module_id)
  {
    case stellar_modules::blackbody :
      module = std::make_unique<BlackbodyStar>(reader.requireDouble("temperature"), instellation);
      break;

    case stellar_modules::tabulated :
      module = std::make_unique<TabulatedStar>(reader.requireString("file"), instellation);
      break;
  }

  reader.finish();
  return module;
}


}
#endif
