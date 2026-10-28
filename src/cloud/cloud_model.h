/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#ifndef _cloud_model_h
#define _cloud_model_h

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "../atmosphere/atmosphere.h"
#include "../spectral_grid/spectral_grid.h"


namespace ngam {


// Per-layer optical properties of ONE cloud module on the model's spectral grid, in the layout the
// radiative-transfer backends consume (opacity_calc.h): [spectral point][layer], with
// nb_layers = nb_grid_points - 1 and layer j between levels j and j+1 (0 = bottom).
//   optical_depth      extinction optical depth of the layer
//   single_scattering  scattering / extinction
//   asymmetry          Henyey-Greenstein asymmetry parameter of the scattered part
struct CloudOpticalProperties {
  std::vector<std::vector<double>> optical_depth;
  std::vector<std::vector<double>> single_scattering;
  std::vector<std::vector<double>> asymmetry;

  void assignZero(const size_t nb_spectral_points, const size_t nb_layers)
  {
    optical_depth.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));
    single_scattering.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));
    asymmetry.assign(nb_spectral_points, std::vector<double>(nb_layers, 0.0));
  }

  bool empty() const { return optical_depth.empty(); }
};


// A cloud module: anything that produces per-layer (tau, w, g) on the spectral grid -- a tabulated
// deck, a Mie slab, an Ackerman & Marley cloud, ... The owning object keeps a LIST of them (like the
// chemistry modules), sums their optical properties with the scattering-weighted mixing rule
// (combineCloudOptics below) and hands the result to the radiative transfer.
//
// Operator split (doc/cloud_module_design.md Sec. 3.4): update() runs ONCE per outer iteration from
// the committed atmosphere. The temperature corrector's trial evaluations never call it; they
// reapply the committed optics unchanged, so the corrector sees the cloud as frozen opacity --
// exactly like the gas opacity in its Jacobians. A module that relaxes or ramps its own structure
// between calls (the Ackerman & Marley module) does so inside update().
class CloudModel {
  public:
    CloudModel(SpectralGrid* spectral_grid_, const size_t nb_grid_points)
      : spectral_grid(spectral_grid_), nb_layers(nb_grid_points - 1)
    {
      committed.assignZero(spectral_grid->nbSpectralPoints(), nb_layers);
    }
    virtual ~CloudModel() {}

    // Recompute the cloud from the committed atmosphere (pressure, temperature, structure,
    // composition) and store it in `committed`.
    virtual void update(const Atmosphere& atmosphere) = 0;

    const CloudOpticalProperties& optics() const { return committed; }

    // The model's eddy diffusion profile [cm^2/s per level] and gravity, handed over by the owning
    // object before update(). Only modules that need them override these.
    virtual void setKzz(const std::vector<double>&) {}
    virtual void setSurfaceGravity(const double) {}

    // Diagnostic: the committed optical depth per layer at the spectral point closest to the given
    // wavelength [micron].
    std::vector<double> opticalDepth(const double wavelength) const
    {
      const std::vector<double>& wl = spectral_grid->wavelength_list;
      size_t best = 0;
      for (size_t k = 1; k < wl.size(); ++k)
        if (std::abs(wl[k] - wavelength) < std::abs(wl[best] - wavelength)) best = k;
      return committed.optical_depth[best];
    }

    // True for modules that read the chemistry's condensed phase (atmosphere.condensates); the
    // owning object checks at construction that a chemistry module provides it.
    virtual bool requiresCondensates() const { return false; }

    // Set by the owning object for the duration of compute(): the current profile is (near) a
    // converged state, so a module's start-up ramp is skipped.
    virtual void setWarmStart(const bool) {}

    // The temperature corrector's convergence residual of the previous outer iteration (< 0 =
    // unknown), handed over before update(); a module gates its start-up ramp on it.
    virtual void setResidual(const double) {}

    // Per-layer diagnostics (name -> nb_layers values) for the output file, e.g. the condensate
    // mass mixing ratio, particle radius and column density of a microphysical module.
    virtual std::map<std::string, std::vector<double>> diagnostics() const { return {}; }

    // Set by update(): the committed structure changed in this update (false when a module held
    // its previous structure). The owning object restarts its stall detection on a change, since
    // the problem the corrector solves has changed.
    bool structure_changed = false;

    // A disabled module is skipped by the owning object (no optical properties, no update).
    bool enabled = true;

    // the selector keyword this module was created from (for diagnostics and output)
    std::string type;

  protected:
    SpectralGrid* spectral_grid = nullptr;
    size_t nb_layers = 0;
    CloudOpticalProperties committed;
};


// Sum the optical properties of several cloud modules into one set of per-layer arrays:
//   tau = sum tau_k,   w = sum tau_k w_k / tau,   g = sum tau_k w_k g_k / sum tau_k w_k
// i.e. w weighted by extinction and g by scattering optical depth (the only consistent rule for
// a mix of absorbing and scattering decks). The output arrays must already have the
// [spectral point][layer] shape; they are overwritten.
inline void combineCloudOptics(
  const std::vector<const CloudOpticalProperties*>& parts,
  std::vector<std::vector<double>>& optical_depth,
  std::vector<std::vector<double>>& single_scattering,
  std::vector<std::vector<double>>& asymmetry)
{
  const size_t nb_spectral_points = optical_depth.size();

  #pragma omp parallel for
  for (size_t nu = 0; nu < nb_spectral_points; ++nu)
  {
    const size_t nb_layers = optical_depth[nu].size();

    for (size_t j = 0; j < nb_layers; ++j)
    {
      double tau = 0, tau_sca = 0, tau_sca_g = 0;

      for (const auto* p : parts)
      {
        const double t = p->optical_depth[nu][j];
        const double s = t * p->single_scattering[nu][j];
        tau += t;
        tau_sca += s;
        tau_sca_g += s * p->asymmetry[nu][j];
      }

      optical_depth[nu][j] = tau;
      single_scattering[nu][j] = tau > 0 ? std::min(tau_sca / tau, 1.0) : 0.0;
      asymmetry[nu][j] = tau_sca > 0 ? tau_sca_g / tau_sca : 0.0;
    }
  }
}


}
#endif
