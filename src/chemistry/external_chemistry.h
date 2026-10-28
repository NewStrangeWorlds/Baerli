/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _external_chemistry_h
#define _external_chemistry_h

#include <string>
#include <vector>

#include "chemistry.h"
#include "chem_species.h"


namespace baerli {


// Composition supplied from outside the model, e.g. by a chemical kinetics code driven from
// Python (the neoVULCAN coupling). The module holds volume mixing ratio profiles on the model's own
// pressure grid for a subset of species and OVERWRITES those species in the composition produced by
// the preceding modules; everything else (species the external code does not know: TiO, VO, Na,
// K, ...) is left as it is. List it last: chemistry=[("equilibrium", {...}), ("external", {})],
// optionally with ("quench", {...}) in between as the first-pass guess for the kinetics.
//
// Until setMixingRatios() has been called the module is a no-op, so the same configuration serves
// the very first pass (equilibrium or quench composition) and every later one.
//
// Budget: the external mixing ratios are normalised to sum to one over the external code's own
// species, while the equilibrium solver allotted those same species a slightly smaller total (the
// rest sits in trace species outside the network). The supplied profiles are therefore rescaled
// per level so that the SUM over the supplied species equals what the preceding modules gave that
// same set -- the total stays conserved and the untouched species keep their share. The mean
// molecular weight is updated differentially (sum over changed species of df * m), for the same
// reason as in the quench module.
class ExternalChemistry : public Chemistry{
  public:
    ExternalChemistry();
    virtual ~ExternalChemistry() {}

    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);

    // symbols: species symbols as in chem_species.h (unknown ones are ignored, reported once);
    // mixing_ratios[level][k] for symbol k, on the model's pressure grid (index 0 = bottom).
    void setMixingRatios(
      const std::vector<std::string>& symbols,
      const std::vector<std::vector<double>>& mixing_ratios);

    // forget the external composition: the module becomes a no-op again
    void clear() { profiles.clear(); }

    bool isSet() const { return !profiles.empty(); }

  private:
    struct Profile { chemical_species_id id; std::vector<double> mixing_ratio; };
    std::vector<Profile> profiles;
    std::vector<std::string> ignored_symbols;   // reported once
    bool reported_ignored = false;
};


}
#endif
