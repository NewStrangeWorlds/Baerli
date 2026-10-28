/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _cloud_optics_h
#define _cloud_optics_h

#include <string>
#include <vector>

#include "../spectral_grid/spectral_grid.h"


namespace baerli {


// Mie efficiencies of one condensate, tabulated with LX-MIE (Kitzmann & Heng 2018) on a radius grid
// and a COARSE wavelength grid, and the machinery to turn a particle population into per-layer
// optical properties on the model's spectral grid.
//
// Why a coarse wavelength grid: Q_ext, Q_sca and g of a size DISTRIBUTION are smooth in wavelength
// (the single-size ripple averages out), so tabulating them at a resolving power of ~200 and
// interpolating linearly in ln(lambda) onto the opacity-sampling grid (10^4-10^5 points) costs
// nothing in accuracy and saves a factor ~100 in table size and per-iteration work.
//
// Radius grid: virga's convention -- nb_radii log-spaced radii between r_min and r_max (virga: 60,
// here 120 by default for a finer ln r interpolation), each bin averaged over `sub_bins` radii
// spread across the bin to smooth the single-size resonances (virga's calc_mie). A population is
// handed over as quadrature points (r_j, N_j) [cm, cm^-2]; the efficiencies are interpolated as
// ln Q in ln r between the bin centres (exact for the Rayleigh power laws), so a monodisperse
// population is a single point and a lognormal one a set of points from lognormalQuadrature().
// Accuracy against a direct Mie call for a single size: ~1 % (120 radii), apart from the
// resonance ripple of weak absorbers that the sub-bin average deliberately smooths (~2 %).
class MieTable {
  public:
    MieTable(
      const std::string& refractive_index_file,
      SpectralGrid* spectral_grid,
      const double mie_resolution = 200.0,
      const size_t nb_radii = 120,
      const double r_min = 1e-8,
      const double r_max = 5.4239131e-2,
      const unsigned int sub_bins = 6);

    // Column optical properties of a population on the coarse grid:
    //   tau_ext(l) = sum_j N_j pi r_j^2 Q_ext(r_j, l),  tau_sca likewise,
    //   tau_sca_g(l) = sum_j N_j pi r_j^2 Q_sca g  (the scattering-weighted asymmetry numerator)
    void integrate(
      const std::vector<double>& radius,
      const std::vector<double>& number_density,
      std::vector<double>& tau_ext,
      std::vector<double>& tau_sca,
      std::vector<double>& tau_sca_g) const;

    // a coarse-grid quantity on the model's spectral grid (linear in ln lambda)
    std::vector<double> toSpectralGrid(const std::vector<double>& coarse) const;

    // value of a coarse-grid quantity at an arbitrary wavelength [micron] (linear in ln lambda)
    double atWavelength(const std::vector<double>& coarse, const double wavelength) const;

    const std::vector<double>& wavelengths() const { return coarse_wavelength; }
    const std::vector<double>& radii() const { return radius_grid; }
    const std::string& file() const { return file_name; }

    // efficiencies at a bin centre
    double qExt(const size_t i_radius, const size_t i_wavelength) const { return q_ext[i_radius][i_wavelength]; }
    double qSca(const size_t i_radius, const size_t i_wavelength) const { return q_sca[i_radius][i_wavelength]; }
    double asymmetry(const size_t i_radius, const size_t i_wavelength) const { return g[i_radius][i_wavelength]; }

  private:
    std::string file_name;
    std::vector<double> coarse_wavelength;            // [micron], increasing
    std::vector<double> radius_grid, log_radius;      // bin centres [cm], ln r
    std::vector<std::vector<double>> q_ext, q_sca, g; // [radius][wavelength]

    // interpolation of the coarse grid onto the spectral grid: lower index and weight per point
    std::vector<size_t> grid_index;
    std::vector<double> grid_weight;

    void readRefractiveIndex(
      const std::string& file, std::vector<double>& wavelength,
      std::vector<double>& n, std::vector<double>& k) const;
    void buildTable(
      const std::vector<double>& wavelength_data, const std::vector<double>& n_data,
      const std::vector<double>& k_data, const unsigned int sub_bins);
    void radiusWeights(const double radius, size_t& lower, double& weight) const;
};


// Quadrature points (radius [cm], column number density [cm^-2]) of a lognormal size distribution
// with geometric mean radius r_g, geometric standard deviation sigma (> 1; sigma = 1 is a single
// point) and total column number density N: nb_points points equally spaced in ln r over
// +-4 ln(sigma), trapezoid weights renormalised to N.
void lognormalQuadrature(
  const double r_g, const double sigma, const double column_density,
  std::vector<double>& radius, std::vector<double>& number_density,
  const size_t nb_points = 41);


// Fraction of a column that falls into each model LAYER when the column is distributed uniformly
// in ln p between p_bottom and p_top (clipped to the model domain): nb_grid_points-1 values that
// sum to 1 (to the clipped fraction if the deck extends beyond the domain).
std::vector<double> layerFractionsUniformInLnP(
  const std::vector<double>& pressure, const double p_bottom, const double p_top);


}
#endif
