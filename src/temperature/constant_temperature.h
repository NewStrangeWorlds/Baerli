/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _constant_temperature_h
#define _constant_temperature_h

#include "temperature.h"
#include "../atmosphere/atmosphere.h"

#include <vector>


namespace baerli {


class ConstantTemperature : public TemperatureProfile{
  public:
    ConstantTemperature() {nb_parameters = 1;}
    virtual ~ConstantTemperature() {}
    virtual void calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      Atmosphere& atmosphere) {
        atmosphere.temperature.assign(atmosphere.pressure.size(), parameters[0]);}
};


}
#endif
