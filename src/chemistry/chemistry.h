/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _chemistry_h
#define _chemistry_h

#include <vector>


#include "chem_species.h"
#include "../atmosphere/atmosphere.h"


namespace baerli {


class Chemistry{
  public:
    virtual ~Chemistry() {}
    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight) = 0;
    size_t nbParameters() {return nb_parameters;}

    // Called by the owning object after construction. Only modules that need the gravity (e.g.
    // the quench module's mixing timescale) override it.
    virtual void setSurfaceGravity(const double) {}

    // The model's eddy diffusion profile [cm^2/s, one value per level], passed by the owning object
    // before every outer-iteration chemistry evaluation. Only modules that need it override this.
    virtual void setKzz(const std::vector<double>&) {}

    // lagged = true while the temperature corrector probes TRIAL profiles: a module whose result
    // depends on T through anything the corrector's Jacobian does not model (the quench module's
    // quench points, for instance) must then reapply its last committed state instead of
    // re-deriving it from the trial profile. Modules that are cheap and smooth ignore this.
    virtual void setLagged(const bool) {}

    // The condensed phase this module computed alongside the gas (equilibrium chemistry with
    // condensation), copied onto the atmosphere by the owning object after every evaluation.
    // nullptr for modules without a condensed phase.
    virtual const CondensatePhase* condensates() const { return nullptr; }

    // A disabled module is skipped by the owning object (its composition passes through). Used for
    // a chemistry continuation: converge with equilibrium first, then switch the quench module on
    // and re-solve warm -- from a cold analytic start the quench composition can be so far from the
    // converged one (every family quenched at the grid bottom, a CH4/NH3 column up to the top) that
    // the damped Newton's linearisation of the deep level is useless (measured: warm Jupiter,
    // Guillot start, +1700 K at the bottom in one call, then a permanent crawl).
    bool enabled = true;

    // the module's own runtime parameters (e.g. metallicity and C/O for the equilibrium module,
    // the mixing ratios for the isoprofile module), set from its module spec at construction;
    // the objects pass these to calcChemicalComposition
    std::vector<double> parameters;
  protected:
    size_t nb_parameters {};
    std::vector<chemical_species_id> species;

    void meanMolecularWeight(
      const std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight) {
        for (size_t i=0; i<number_densities.size(); ++i)
        {
          double mu = 0;

          for (auto & j : constants::species_data)
            mu += number_densities[i][j.id]/number_densities[i][_TOTAL] * j.molecular_weight;

          mean_molecular_weight[i] = mu;
        }
      };

    bool checkMixingRatios(
      const std::vector<std::vector<double>>& number_densities) {
        for (size_t i=0; i<number_densities.size(); ++i)
        {
          double sum = 0;

          for (auto & j : constants::species_data)
            if (j.id != _TOTAL)
              sum += number_densities[i][j.id]/number_densities[i][_TOTAL];

          if (sum > 1+1e-10)
            return false;
        }

        return true;
      };
};


}
#endif 
