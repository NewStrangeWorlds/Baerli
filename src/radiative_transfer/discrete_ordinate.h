/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _discrete_ordinate_h
#define _discrete_ordinate_h


#include <vector>
#include <iostream>
#include <cmath>
#include <functional>

#include "radiative_transfer.h"
#include "../atmosphere/atmosphere.h"
#include "../spectral_grid/spectral_grid.h"


#include "../../_deps/disortpp-src/src/DisortFluxConfig.hpp"
#include "../../_deps/disortpp-src/src/FluxResult.hpp"
#include "../../_deps/disortpp-src/src/FluxSolver.hpp"


namespace baerli {


class DiscreteOrdinates : public RadiativeTransfer{
  public:
    DiscreteOrdinates(
      SpectralGrid* spectral_grid_ptr,
      const size_t nb_streams,
      const size_t nb_grid_points);
    virtual ~DiscreteOrdinates() {}
    
    virtual void calculate(
      const Atmosphere& atmosphere,
      const OpacityCalculation& opacity,
      RadiativeTransferOutput& output,
      const RadiativeBoundaryConditions& bc = RadiativeBoundaryConditions{}) override;
  private:
    size_t nb_streams = 0;
    size_t nb_grid_points = 0;
    size_t nb_layers = 0;

    int nb_threads = 0;

    //std::vector<disort_state> ds;
    //std::vector<disort_output> out;
    std::vector<disortpp::DisortFluxConfig> configs;
    std::vector<std::function<disortpp::FluxResult(const disortpp::DisortFluxConfig&)>> solvers;

    // scratch for the in-loop temperature-Jacobian accumulation (only used when
    // output.compute_jacobian is set): per-point trapezoidal weights and per-thread
    // accumulators reduced after the parallel spectral loop.
    std::vector<double> quad_weights;
    std::vector< std::vector<double> >                  heating_value_thread;       // [thread][i]
    std::vector< std::vector< std::vector<double> > >   jac_net_heating_thread;     // [thread][i][j]
    std::vector< std::vector< std::vector<double> > >   jac_net_flux_thread;        // [thread][i][j]
    std::vector< std::vector< std::vector<double> > >   jac_meanint_kappa_thread;   // [thread][i][j]

    void calculate(
      const OpacityCalculation& opacity,
      const std::vector<double>& vertical_grid,
      const double surface_albedo,
      const bool has_surface,
      const double incident_radiation,
      const double zenith_angle,
      const size_t nu_index,
      const double max_temperature,
      RadiativeTransferOutput& output);
    void calcTotalTransportCoeff(
      const OpacityCalculation& opacity,
      const std::vector<double>& vertical_grid, 
      const size_t nu_index,
      std::vector<double>& optical_depth,
      std::vector<double>& single_scattering_albedo,
      std::vector<double>& asymmetry_parameter);
    void setTemperatureStructure(
      const std::vector<double>& temperature_structure,
      const double& surface_temperature);
    void setDISORTParam(
      const int thread_id,
      const double wavenumber_input,
      const std::vector<double>& optical_depth,
      const std::vector<double>& single_scattering_albedo,
      const std::vector<double>& asymmetry_parameter,
      const double incident_radiation,
      const double zenith_angle,
      const double surface_albedo,
      const bool has_surface);
    void initDISORT();
};


}
#endif


