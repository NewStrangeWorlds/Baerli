/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _temperature_profile_h
#define _temperature_profile_h

#include <vector>
#include <cstddef>

#include "../atmosphere/atmosphere.h"

namespace baerli {


class TemperatureProfile{
  public:
    virtual ~TemperatureProfile() {}
    virtual void calcProfile(
      const std::vector<double>& parameters,
      const double surface_gravity,
      Atmosphere& atmosphere) = 0;
    size_t nbParameters() {return nb_parameters;}
  protected:
    size_t nb_parameters {};
};


}
#endif
