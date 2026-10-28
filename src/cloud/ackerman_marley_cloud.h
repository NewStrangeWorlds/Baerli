/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _ackerman_marley_cloud_h
#define _ackerman_marley_cloud_h

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "cloud_model.h"
#include "cloud_optics.h"
#include "condensate_data.h"


namespace baerli {


// Ackerman & Marley (2001) eddy-sedimentation cloud of ONE condensate, in the form of virga's
// layer solver (Rooney et al. 2022) but fed by the equilibrium chemistry's condensation instead of
// vapour-pressure fits (doc/cloud_module_design.md Sec. 5).
//
// Condensation increments: the `equilibrium` chemistry with condensation = rainout stores at every
// level the condensate density that forms AT that level from the vapour arriving from below
// (atmosphere.condensates). Above the cloud base the vapour is at saturation in both the A&M
// picture and the rainout, so the only thing f_sed changes is the suspended condensate, which obeys
// the layer recursion (virga's closed form, Delta q_cond assigned to the layer below the level)
//   q_c,j = q_c,j-1 exp(-f_sed Delta z_j / L_j) + Delta q_cond,j
// with the mixing length L = H (or max(0.1, nabla/nabla_ad) H). The gas-phase depletion is the
// chemistry's job, so this module never touches the composition.
//
// Particle size per layer (A&M eqs. 13-18): the radius r_w whose fall speed equals the convective
// velocity scale w* = Kzz/L (Brent on fallVelocity), the local exponent alpha of v ~ (r/r_w)^alpha
// (least squares over the Mie table's radius grid, as virga), then the lognormal moments
//   r_g = r_w f_sed^(1/alpha) exp(-(alpha+6) l),  r_eff = r_g exp(5 l),
//   N   = 3 rho q_c Delta z / (4 pi rho_p r_g^3) exp(-9 l),   l = ln^2(sigma) / 2.
//
// Damping between outer iterations (the Kzz pattern, plus what the first cloudy runs taught): the
// structure is recomputed only while the corrector's residual of the previous iteration is below
// `gate` (a cloud that re-forms on every iteration from a Kzz that is itself still settling keeps
// the Newton from ever converging), the per-layer (N, r_g) are relaxed in log space towards the
// new values with their change per update capped at `max_change` dex, a layer the cloud enters
// starts at `seed` of its target and a layer it leaves decays by the cap. The structure is
// committed only when it moved by more than `tolerance` dex anywhere (dead band, so the cloud
// freezes once settled and the Newton can finish). The ONSET of the cloud opacity in the radiative
// transfer (nothing during the corrector's easy start, then an optical-depth-capped ramp) is the
// owning object's job for every cloud module (GenericObject::updateCloudScale).
class AckermanMarleyCloud : public CloudModel {
  public:
    struct Parameters {
      std::string condensate;          // name in condensate_data.h
      double fsed = 1.0;
      bool fsed_exponential = false;   // Rooney et al. (2022): f_sed(z) = eps + (f_sed - eps) exp(-(z - z_base)/(b H))
      double b = 1.0;                  // in scale heights
      double eps = 1e-2;
      double sigma = 2.0;              // lognormal geometric standard deviation
      double kzz_min = 1e5;            // cm^2/s, floor on the eddy diffusion (virga's kz_min)
      bool lapse_ratio_mixing_length = false;   // L = max(0.1, nabla/nabla_ad) H instead of H
      double relax = 0.5;              // log-space relaxation weight of the structure between updates (1 = none)
      double max_change = 0.3;         // cap [dex] on the change of N and r_g per update (0 = none)
      double tolerance = 0.05;         // dead band [dex] on N and r_g
      double gate = 1e-2;              // structure updates only while the corrector's residual is below this
      double seed = 1e-2;              // a layer the cloud enters starts at this fraction of its target
      double bulk_density = 0.0;       // g/cm^3, 0 = from the condensate table
      double molecule_diameter = 2.827e-8;   // cm, gas (H2) for the viscosity / mean free path
      double lj_temperature = 59.7;          // K
      double mie_resolution = 200.0;
      size_t mie_radii = 120;
    };

    AckermanMarleyCloud(
      SpectralGrid* spectral_grid_, const size_t nb_grid_points, const Parameters& parameters);
    virtual ~AckermanMarleyCloud() {}

    void update(const Atmosphere& atmosphere) override;
    void setKzz(const std::vector<double>& kzz_) override { kzz = kzz_; }
    void setSurfaceGravity(const double g) override { gravity = g; }
    void setWarmStart(const bool w) override { warm_start = w; }
    void setResidual(const double r) override { residual = r; }
    bool requiresCondensates() const override { return true; }
    std::map<std::string, std::vector<double>> diagnostics() const override;

    const std::string& fastchemSymbol() const { return fastchem_symbol; }

  private:
    Parameters par;
    std::string fastchem_symbol;
    double bulk_density = 0.0;
    double gravity = 0.0;
    bool warm_start = false;
    double residual = -1.0;          // the corrector's residual of the previous outer iteration
    size_t nb_updates = 0;

    std::unique_ptr<MieTable> mie_table;
    std::vector<double> kzz;

    // per layer: the physical (unrelaxed) A&M solution of the last update ...
    std::vector<double> q_c, r_eff, fsed_layer, alpha_layer, r_w_layer;
    // ... and the committed structure the optics are built from: the condensate MASS column per
    // layer [g/cm^2] (the conserved quantity, relaxed and capped between updates) and the
    // geometric mean radius (always the physical value where the layer has cloud -- relaxing the
    // radius lagged it 20x behind the physical value for dozens of iterations and made the deck
    // hundreds of times too thick); the column number density follows from the two
    std::vector<double> mass_column, r_g;
    std::vector<double> columnDensity() const;

    void computeOptics();
};


}
#endif
