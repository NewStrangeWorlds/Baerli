/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _fixed_chemistry_h
#define _fixed_chemistry_h

#include "chemistry.h"


#include <vector>
#include <string>


namespace baerli {


class FixedChemistry : public Chemistry{
  public:
    FixedChemistry(const std::string& file_path)
    {
      nb_parameters = 0;
      readFile(file_path);
    }
    virtual ~FixedChemistry() {}
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  protected:
    void readFile(const std::string& file_path);

    std::vector<std::vector<double>> fixed_mixing_ratios;
    std::vector<double> pressure;
};


}
#endif 
