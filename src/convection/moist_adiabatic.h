/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _moist_adiabatic_h
#define _moist_adiabatic_h

#include <iostream>

#include "convection.h"
#include "../atmosphere/atmosphere.h"


namespace baerli {


class MoistAdiabaticAdjustment : public Convection {
  public:
    MoistAdiabaticAdjustment(size_t max_sweeps_ = 10, double min_pressure_ = 1e-4)
      : max_sweeps(max_sweeps_), min_pressure(min_pressure_)
    {
      std::cout << "\n- Convection: moist adiabatic adjustment"
                << " (min pressure: " << min_pressure << " bar)\n\n";
    }
    virtual ~MoistAdiabaticAdjustment() {}
    void adjust(Atmosphere& atmosphere) const override;
    double convectiveGradient(
      const std::vector<double>& number_densities,
      double temperature,
      double pressure) const override;
  private:
    size_t max_sweeps;
    double min_pressure;
};


}

#endif
