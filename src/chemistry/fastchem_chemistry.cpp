/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#include "fastchem_chemistry.h"

#include "chem_species.h"
#include "../additional/exceptions.h"
#include "../additional/physical_const.h"
#include "../../_deps/fastchem-src/fastchem_src/fastchem.h"


#include <algorithm>
#include <chrono>
#include <memory>
#include <vector>
#include <omp.h>
#include <math.h>
#include <iostream>


namespace baerli {


FastChemChemistry::FastChemChemistry(
  const std::string& fastchen_parameter_file,
  const Condensation condensation_,
  const bool accelerate_)
  : fastchem(fastchen_parameter_file, 1)
  , condensation(condensation_)
  , accelerate(accelerate_)
{
  std::cout << "- Chemistry model: " << "equilibrium/FastChem" << "\n";
  std::cout << "  - Parameter file: " << fastchen_parameter_file << "\n";


  reference_element_abundances = fastchem.getElementAbundances();


  fastchem_species_indices.assign(constants::species_data.size(), fastchem::FASTCHEM_UNKNOWN_SPECIES);

  for (size_t i=0; i<constants::species_data.size(); ++i)
    fastchem_species_indices[i] = fastchem.getGasSpeciesIndex(constants::species_data[i].fastchem_symbol);


  //check if C, O, and H are present in FastChem
  if (fastchem_species_indices[_H] == fastchem::FASTCHEM_UNKNOWN_SPECIES
      || fastchem_species_indices[_O] == fastchem::FASTCHEM_UNKNOWN_SPECIES
      || fastchem_species_indices[_C] == fastchem::FASTCHEM_UNKNOWN_SPECIES)
  {
    std::string error_message = "Critical elements (H, C, or O) not found in FastChem\n";
    throw InvalidInput(std::string ("FastChemChemistry::FastChemChemistry"), error_message);
  }


  if (condensation != Condensation::none)
  {
    const unsigned int nb_condensates = fastchem.getCondSpeciesNumber();

    if (nb_condensates == 0)
      throw InvalidInput("FastChemChemistry",
        "condensation requested but the FastChem parameter file lists no condensate data; add the\n"
        "condensate file to the species data line, e.g.\n"
        "  fastchem_data/logK.dat fastchem_data/logK_condensates.dat\n"
        "(see fastchem_parameters_condensates.dat)\n");

    //formula weights from the stoichiometry and the element weights (FastChem's
    //getCondSpeciesWeight returns 0 for condensates)
    for (unsigned int k = 0; k < nb_condensates; ++k)
    {
      const std::vector<int> stoichiometry = fastchem.getCondSpeciesStoichiometry(k);
      double weight = 0.0;
      for (size_t e = 0; e < stoichiometry.size(); ++e)
        weight += stoichiometry[e] * fastchem.getElementWeight(e);

      if (!(weight > 0))
        throw InvalidInput("FastChemChemistry",
          "could not determine the formula weight of condensate " + fastchem.getCondSpeciesSymbol(k) + "\n");

      condensed_phase.symbols.push_back(fastchem.getCondSpeciesSymbol(k));
      condensed_phase.weights.push_back(weight);
    }

    //one instance per thread for the lagged, per-level gas-phase calls. Built from the parameter
    //file rather than copied: FastChem's copy constructor does not carry the is_initialised flag,
    //so a copy rejects every public calcDensities() call (FastChem's own thread copies only use
    //the private per-point routine and never notice).
    for (int t = 0; t < omp_get_max_threads(); ++t)
      thread_copies.push_back(
        std::make_unique<fastchem::FastChem<long double>>(fastchen_parameter_file, 0));

    std::cout << "  - Condensation: "
              << (condensation == Condensation::rainout ? "rainout" : "equilibrium")
              << " (" << nb_condensates << " condensates), lagged in the temperature corrector\n";
  }

  std::cout << "\n";

  nb_parameters = 2;
}



std::vector<std::string> FastChemChemistry::elementSymbols()
{
  std::vector<std::string> symbols;
  for (unsigned int e = 0; e < fastchem.getElementNumber(); ++e)
    symbols.push_back(fastchem.getElementSymbol(e));
  return symbols;
}



std::vector<double> FastChemChemistry::scaledElementAbundances(
  const std::vector<double>& parameters) const
{
  const double metallicity_factor = parameters[0];
  const double co_ratio = parameters[1];

  std::vector<double> element_abundances = reference_element_abundances;

  for (size_t i=0; i<element_abundances.size(); ++i)
    if (i != fastchem_species_indices[_H] && (i != fastchem_species_indices[_He] && fastchem_species_indices[_He] != fastchem::FASTCHEM_UNKNOWN_SPECIES) )
      element_abundances[i] *= metallicity_factor;

  element_abundances[fastchem_species_indices[_O]] = element_abundances[fastchem_species_indices[_C]] / co_ratio;

  return element_abundances;
}



//one FastChem call over the whole column, with or without condensation
unsigned int FastChemChemistry::runFull(
  const std::vector<double>& element_abundances,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  const bool with_condensation,
  fastchem::FastChemOutput& output)
{
  fastchem.setElementAbundances(element_abundances);

  fastchem::FastChemInput input;
  input.temperature = temperature;
  input.pressure = pressure;

  for (auto & i : input.temperature)
    if (i < 100) i = 100;

  if (with_condensation)
  {
    input.rainout_condensation = (condensation == Condensation::rainout);
    input.equilibrium_condensation = (condensation == Condensation::equilibrium);
  }

  return fastchem.calcDensities(input, output);
}



//gas phase only, level by level in parallel, the elements depleted by the committed degrees of
//condensation: epsilon (1 - delta[level])
unsigned int FastChemChemistry::runDepleted(
  const std::vector<double>& element_abundances,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  fastchem::FastChemOutput& output)
{
  const size_t nb_levels = temperature.size();
  const size_t nb_elements = element_abundances.size();
  const int nb_threads = static_cast<int>(thread_copies.size());

  output.number_densities.assign(nb_levels, std::vector<double>(fastchem.getGasSpeciesNumber(), 0.0));
  output.mean_molecular_weight.assign(nb_levels, 0.0);
  output.total_element_density.assign(nb_levels, 0.0);
  output.fastchem_flag.assign(nb_levels, fastchem::FASTCHEM_SUCCESS);

  //FastChem's calcDensities sets the OpenMP thread count for its own (here one-point) loop,
  //which would stick to the calling threads; restore it afterwards
  const int saved_threads = omp_get_max_threads();

  #pragma omp parallel for schedule(dynamic, 1) num_threads(nb_threads)
  for (size_t i = 0; i < nb_levels; ++i)
  {
    fastchem::FastChem<long double>& fc = *thread_copies[omp_get_thread_num()];

    std::vector<double> depleted(nb_elements);
    for (size_t e = 0; e < nb_elements; ++e)
      depleted[e] = element_abundances[e] * (1.0 - std::min(cond_degree[i][e], 1.0 - 1e-15));

    fc.setElementAbundances(depleted);

    fastchem::FastChemInput input;
    input.temperature = {std::max(temperature[i], 100.0)};
    input.pressure = {pressure[i]};

    fastchem::FastChemOutput out;
    output.fastchem_flag[i] = fc.calcDensities(input, out);

    if (out.number_densities.empty())
    {
      //FastChem returned before computing anything (initialisation / busy / input error)
      #pragma omp critical
      std::cout << "  WARNING: FastChem per-level call failed with flag " << output.fastchem_flag[i]
                << " at level " << i << "\n";
      if (output.fastchem_flag[i] == fastchem::FASTCHEM_SUCCESS)
        output.fastchem_flag[i] = fastchem::FASTCHEM_NO_CONVERGENCE;
      continue;
    }

    output.number_densities[i] = out.number_densities[0];
    output.mean_molecular_weight[i] = out.mean_molecular_weight[0];
    output.total_element_density[i] = out.total_element_density[0];
  }

  omp_set_num_threads(saved_threads);

  unsigned int status = fastchem::FASTCHEM_SUCCESS;
  for (const auto f : output.fastchem_flag)
    if (f != fastchem::FASTCHEM_SUCCESS) status = f;

  return status;
}



//Steffensen cycle on the committed degrees of condensation: two raw values are taken as they are,
//the third is replaced by the Aitken extrapolation of ln(1 - delta) (the gas remainder, which is
//what the opacity sees) wherever the sequence contracts geometrically and slowly (ratio in
//(0.5, 0.95)) with a small step (|d ln| < 0.1): the slow tail of the lag, not the early wild phase.
std::vector<std::vector<double>> FastChemChemistry::accelerateDegrees(
  const std::vector<std::vector<double>>& raw)
{
  if (!accelerate) return raw;

  if (degree_history.size() < 2)
  {
    degree_history.push_back(raw);
    return raw;
  }

  const auto& d0 = degree_history[0];
  const auto& d1 = degree_history[1];
  std::vector<std::vector<double>> accelerated = raw;
  size_t nb_accelerated = 0;

  auto logRemainder = [](double delta) { return std::log(std::max(1.0 - delta, 1e-300)); };

  for (size_t i = 0; i < raw.size() && i < d0.size() && i < d1.size(); ++i)
    for (size_t e = 0; e < raw[i].size(); ++e)
    {
      const double x0 = logRemainder(d0[i][e]), x1 = logRemainder(d1[i][e]), x2 = logRemainder(raw[i][e]);
      const double s0 = x1 - x0, s1 = x2 - x1;

      if (std::abs(s0) < 1e-12 || std::abs(s1) < 1e-12 || std::abs(s1) > 0.1) continue;
      const double ratio = s1 / s0;
      if (ratio <= 0.5 || ratio >= 0.95) continue;

      const double x_acc = std::min(0.0, x2 + s1 * ratio / (1.0 - ratio));
      accelerated[i][e] = 1.0 - std::exp(x_acc);
      ++nb_accelerated;
    }

  degree_history.clear();   // restart the cycle

  if (nb_accelerated > 0)
    std::cout << "  FastChem rainout: Aitken-accelerated the lagged depletion of " << nb_accelerated
              << " (level, element) pairs\n";

  return accelerated;
}



void FastChemChemistry::storeOutput(
  const fastchem::FastChemOutput& output,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  mean_molecular_weight = output.mean_molecular_weight;

  for (size_t j=0; j<temperature.size(); ++j)
    for (size_t i=0; i<constants::species_data.size(); ++i)
    {
      if (fastchem_species_indices[i] != fastchem::FASTCHEM_UNKNOWN_SPECIES)
        number_densities[j][i] = output.number_densities[j][fastchem_species_indices[i]];

      if (i == _TOTAL)
        number_densities[j][_TOTAL] = pressure[j] * 1.e6 / constants::boltzmann_k / temperature[j];
    }
}



bool FastChemChemistry::calcChemicalComposition(
  const std::vector<double>& parameters,
  const std::vector<double>& temperature,
  const std::vector<double>& pressure,
  std::vector<std::vector<double>>& number_densities,
  std::vector<double>& mean_molecular_weight)
{
  const std::vector<double> element_abundances = scaledElementAbundances(parameters);

  fastchem::FastChemOutput output;
  unsigned int status = fastchem::FASTCHEM_SUCCESS;

  if (condensation == Condensation::none)
  {
    status = runFull(element_abundances, temperature, pressure, false, output);
  }
  else if (lagged && cond_degree.size() == temperature.size())
  {
    //trial evaluation: gas phase with the committed depletion
    status = runDepleted(element_abundances, temperature, pressure, output);
  }
  else
  {
    //committed evaluation: the full condensation calculation
    const auto start = std::chrono::steady_clock::now();
    status = runFull(element_abundances, temperature, pressure, true, output);
    const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    if (status == fastchem::FASTCHEM_SUCCESS)
    {
      cond_degree = accelerateDegrees(output.element_cond_degree);
      condensed_phase.densities = output.number_densities_cond;

      //The gas phase is taken from the same depleted gas-only call the Newton's trial evaluations
      //use, so that the committed and the trial compositions agree by construction. The rainout's
      //own gas phase differs from it by FastChem's element-conservation tolerance, which for a
      //99 %-condensed element is a percent-level error on the gas remainder (measured: 2 % in the
      //Na/K depletion zone of a 1300 K brown dwarf, enough to floor the flux residual at 1e-3).
      status = runDepleted(element_abundances, temperature, pressure, output);

      static bool reported = false;
      if (!reported)
      {
        std::cout << "  FastChem " << (condensation == Condensation::rainout ? "rainout" : "equilibrium")
                  << " condensation over " << temperature.size() << " levels: "
                  << seconds << " s\n";
        reported = true;
      }
    }
    else
    {
      //a failed point stops FastChem's rainout march; keep the previous condensation state and
      //fall back to the depleted gas phase (or the plain gas phase if there is none yet)
      std::cout << "  WARNING: FastChem condensation failed (flag " << status
                << "); keeping the previous condensation state for this evaluation\n";
      if (cond_degree.size() == temperature.size())
        status = runDepleted(element_abundances, temperature, pressure, output);
      else
        status = runFull(element_abundances, temperature, pressure, false, output);
    }
  }

  if (status == fastchem::FASTCHEM_INITIALIZATION_FAILED)
  {
    std::string error_message = "FastChem initialisation failed!\n";
    throw InvalidInput(std::string ("FastChemChemistry::calcChemicalComposition"), error_message);
  }

  const bool neglect_model = (status != fastchem::FASTCHEM_SUCCESS);

  storeOutput(output, temperature, pressure, number_densities, mean_molecular_weight);

  return neglect_model;
}



}
