/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _atmosphere_h
#define _atmosphere_h

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>



namespace baerli {


// The condensed phase of the equilibrium chemistry (FastChem with condensation), per level:
// densities[level][k] in formula units per cm^3 for condensate symbols[k] (FastChem's symbols, e.g.
// "MgSiO3(s,l)") with formula weights[k] in g/mol. In rainout mode a level's density is the amount
// that condenses AT that level from what arrives from below -- the condensation increment the cloud
// modules need. Empty without condensation.
struct CondensatePhase {
  std::vector<std::string> symbols;
  std::vector<double> weights;
  std::vector<std::vector<double>> densities;

  bool empty() const { return symbols.empty(); }
  int index(const std::string& symbol) const
  {
    for (size_t k = 0; k < symbols.size(); ++k) if (symbols[k] == symbol) return static_cast<int>(k);
    return -1;
  }
};


class Atmosphere {
  public:
    Atmosphere (
      const size_t nb_grid_points_,
      const std::vector<double>& atmos_boundaries);
    ~Atmosphere() {}

    const size_t nb_grid_points = 0;

    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> altitude;
    std::vector<double> scale_height;
    std::vector<double> mass_density;
    std::vector<double> mean_molecular_weight;
    std::vector< std::vector<double> > number_densities;
    std::vector<int> convective;
    std::vector<double> kzz;          // eddy diffusion coefficient [cm^2/s], set by the object
    CondensatePhase condensates;      // from the chemistry (empty without condensation)

    void calcAtmosphereStructure(
      const double surface_gravity,
      const double bottom_radius,
      const bool use_variable_gravity);

  private:
    void createPressureGrid(const std::vector<double>& domain_boundaries);
    void calcMassDensity();
    void calcAltitude(const double g);
    void calcAltitudeVariableGravity(
      const double g, const double bottom_radius);
    void calcScaleHeight(const double g);
};


}


#endif
