/*
* This file is part of the BeAR code (https://github.com/newstrangeworlds/BeAR).
* Copyright (C) 2024 Daniel Kitzmann
*
* BeAR is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* BeAR is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* BeAR directory under <LICENSE>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#ifndef _fastchem_chemistry_h
#define _fastchem_chemistry_h


#include "chemistry.h"

#include "../../_deps/fastchem-src/fastchem_src/fastchem.h"
#include "chem_species.h"

#include <memory>
#include <string>
#include <vector>



namespace ngam {


// FastChem equilibrium chemistry, optionally with condensation (doc/cloud_module_design.md
// Sec. 5.2). The FastChem parameter file has to list the condensate data file next to the gas-phase
// one for the condensation modes ("fastchem_data/logK.dat fastchem_data/logK_condensates.dat").
//
//   none          gas phase only (the default)
//   equilibrium   equilibrium condensation at every level independently: the condensate stays where
//                 it forms (the f_sed -> 0 limit of a cloud)
//   rainout       FastChem's rainout: bottom-up, the elements condensed at a level are removed from
//                 the gas above it; a level's condensate density is the amount condensing AT that
//                 level -- the increment the Ackerman & Marley cloud module integrates
//
// Operator split for the temperature corrector: the rainout march is sequential over the levels
// (no parallelism inside FastChem), so it runs only on the COMMITTED calls (once per outer
// iteration). While the corrector probes trial profiles (setLagged(true)) the gas phase is solved
// with the elements depleted by the committed degrees of condensation delta[level][element],
// epsilon (1 - delta), level by level in parallel with one FastChem instance per OpenMP thread.
// Exact at convergence (delta no longer changes), invisible to the Jacobians (they are built with
// the opacity frozen anyway), and smoother for the damped Newton than an inside rainout would be.
// The same split is applied to the equilibrium mode for uniformity.
class FastChemChemistry : public Chemistry{
  public:
    enum class Condensation {none, equilibrium, rainout};

    FastChemChemistry(
      const std::string& fastchen_parameter_file,
      const Condensation condensation_ = Condensation::none,
      const bool accelerate_ = true);
    virtual ~FastChemChemistry() {}

    virtual bool calcChemicalComposition(
      const std::vector<double>& parameters,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);

    void setLagged(const bool lagged_) override { lagged = lagged_; }
    const CondensatePhase* condensates() const override
    { return condensation == Condensation::none ? nullptr : &condensed_phase; }

    // the committed degrees of condensation [level][FastChem element] (empty before the first
    // committed call with condensation)
    const std::vector<std::vector<double>>& degreesOfCondensation() const { return cond_degree; }
    std::vector<std::string> elementSymbols();

  private:
    fastchem::FastChem<long double> fastchem;
    const Condensation condensation;
    const bool accelerate;           // Aitken acceleration of the lagged depletion
    bool lagged = false;

    std::vector<double> reference_element_abundances;
    std::vector<size_t> fastchem_species_indices;

    CondensatePhase condensed_phase;
    std::vector<std::vector<double>> cond_degree;        // committed delta[level][element]
    // the last raw rainout results, for the Aitken (Steffensen) acceleration of the lagged
    // depletion: the outer iteration maps delta -> T -> delta with a contraction factor measured at
    // 0.96 for the alkali depletion at the top of a 1300 K brown dwarf (60 slow iterations); every
    // third committed call extrapolates ln(1 - delta) geometrically from the last three values
    std::vector<std::vector<std::vector<double>>> degree_history;
    std::vector<std::vector<double>> accelerateDegrees(const std::vector<std::vector<double>>& raw);
    // one instance per OpenMP thread for the lagged per-level calls (each initialised from the
    // parameter file, see the constructor)
    std::vector<std::unique_ptr<fastchem::FastChem<long double>>> thread_copies;

    std::vector<double> scaledElementAbundances(const std::vector<double>& parameters) const;
    unsigned int runFull(
      const std::vector<double>& element_abundances,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      const bool with_condensation,
      fastchem::FastChemOutput& output);
    unsigned int runDepleted(
      const std::vector<double>& element_abundances,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      fastchem::FastChemOutput& output);
    void storeOutput(
      const fastchem::FastChemOutput& output,
      const std::vector<double>& temperature,
      const std::vector<double>& pressure,
      std::vector<std::vector<double>>& number_densities,
      std::vector<double>& mean_molecular_weight);
};


}
#endif
