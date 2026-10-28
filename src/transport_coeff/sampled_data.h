/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef SAMPLED_DATA_H
#define SAMPLED_DATA_H


#include <vector>
#include <string>
#include <cmath>


namespace baerli{


class CrossSectionFile {
  public:
    CrossSectionFile(const std::string filename_in, const bool log_data)
        : filename(filename_in)
        , is_data_log(log_data) 
        {}
    void loadFile();
    void unloadData();
    
    const std::string filename = "";

    const bool is_data_log = false;
    bool is_loaded = false;

    std::vector<double> cross_sections;
};



class SampledData{
  public:
    SampledData(
      const double temperature_data, 
      const double pressure_data, 
      const std::string file_name, 
      const bool log_data)
        : pressure(pressure_data)
        , log_pressure(std::log10(pressure_data))
        , temperature(temperature_data)
        , data_file(file_name, log_data)
        {}
    ~SampledData();
    void sampleCrossSections(
      const std::vector<size_t>& sampling_list, const double species_mass,
      const std::vector<double>* native_wavenumbers = nullptr,
      const double band_width = 0.0,
      const size_t nb_bands = 0);
    void deleteSampledData();

    const double pressure = 0.0;
    const double log_pressure = 0.0;
    const double temperature = 0.0;

    bool is_sampled = false;

    std::vector<double> cross_sections;
    //per-band integral int_band sigma dnu over the FULL native grid (band-closure correction),
    //same unit convention and log10 storage as cross_sections; computed during sampling, the
    //only moment the complete file is in memory
    std::vector<double> band_integrals;
    //per-band PEAK sigma over the native grid (escape-probability weighting of the closure:
    //the core optical depth to space decides whether the missing cores are thermalised);
    //same unit convention and log10 storage
    std::vector<double> band_peaks;
  private:
    CrossSectionFile data_file;
};


}

#endif
