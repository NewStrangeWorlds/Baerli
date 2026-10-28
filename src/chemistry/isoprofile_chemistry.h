/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _isoprofile_chemistry_h
#define _isoprofile_chemistry_h

#include "chemistry.h"


#include <vector>
#include <string>


namespace baerli {


class IsoprofileChemistry : public Chemistry{
  public:
    IsoprofileChemistry(const std::vector<std::string>& chemical_species);
    virtual ~IsoprofileChemistry() {}
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
  protected:
    bool sodium_free_parameter = false;
};


}
#endif 
