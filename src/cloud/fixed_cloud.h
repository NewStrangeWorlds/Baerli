/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#ifndef _fixed_cloud_h
#define _fixed_cloud_h

#include <string>
#include <vector>

#include "cloud_model.h"


namespace ngam {


// Prescribed cloud optical properties: a table of layer optical depth, single-scattering albedo and
// asymmetry parameter over (pressure, wavelength), read from a file or handed over from Python
// (Object::setCloudOpticalProperties -- the hook for an external cloud/microphysics code).
//
// Table convention: pressure LEVELS p_k [bar] in any order; the row (p_k, lambda) holds the optical
// depth of the input layer between p_k and the next level towards lower pressure (the optical depth
// in the lowest-pressure row is therefore ignored), plus that layer's w and g. The table is mapped
// onto the model grid CONSERVATIVELY: the input optical depth is integrated to a cumulative
// tau(p) from the bottom, interpolated linearly in ln p onto the model levels and differenced onto
// the model layers, so the column optical depth is the same on any model grid. w and g are
// interpolated in ln p at the layer's geometric-mean pressure. Everything is linear in wavelength
// and held constant outside the table's wavelength range (a single-wavelength table is grey).
//
// File format: whitespace-separated columns  pressure[bar]  wavelength[micron]  tau  w  g, '#'
// comments, one row per (pressure, wavelength) pair; the rows must cover the full grid of their
// distinct pressures and wavelengths.
class FixedCloud : public CloudModel {
  public:
    FixedCloud(SpectralGrid* spectral_grid_, const size_t nb_grid_points, const std::string& file);
    virtual ~FixedCloud() {}

    void update(const Atmosphere& atmosphere) override;

    // pressure[np] [bar], wavelength[nw] [micron], optical_depth/single_scattering/asymmetry
    // [np][nw]; see the table convention above
    void setTable(
      const std::vector<double>& pressure,
      const std::vector<double>& wavelength,
      const std::vector<std::vector<double>>& optical_depth,
      const std::vector<std::vector<double>>& single_scattering,
      const std::vector<std::vector<double>>& asymmetry);

    bool hasFile() const { return !file_name.empty(); }
    bool hasTable() const { return !table_pressure.empty(); }

  private:
    std::string file_name;

    // the table, sorted by decreasing pressure (bottom -> top like the model grid) and increasing
    // wavelength
    std::vector<double> table_pressure;
    std::vector<double> table_wavelength;
    std::vector<std::vector<double>> table_tau, table_ssa, table_g;   // [np][nw]

    void readFile(const std::string& file);
};


}
#endif
