/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _quadrature_h
#define _quadrature_h

#include <vector>
#include <stdexcept>

namespace baerli{
namespace aux{



inline double quadratureTrapezoidal(const std::vector<double> &x, const std::vector<double> &y)
{

  if (x.size() != y.size())
    throw std::logic_error("trapezoidal quadrature: x and y must be the same size!\n");

  if (x.size() == 1)
    return y[0];


  double sum = 0.0;

  for (size_t i = 1; i < x.size(); ++i)
      sum += (x[i] - x[i-1]) * (y[i] + y[i-1]);


  return sum * 0.5;
}


// Per-point weights w_i of the trapezoidal rule, i.e. sum_i w_i y_i equals
// quadratureTrapezoidal(x, y). Used to accumulate spectral integrals (and their
// temperature Jacobians) point-by-point inside the monochromatic radiative-transfer
// loop, consistently with quadratureTrapezoidal used elsewhere.
inline std::vector<double> trapezoidalWeights(const std::vector<double>& x)
{
  const size_t n = x.size();
  std::vector<double> w(n, 0.0);

  if (n == 1) { w[0] = 1.0; return w; }

  w[0]   = 0.5 * (x[1] - x[0]);
  w[n-1] = 0.5 * (x[n-1] - x[n-2]);

  for (size_t i = 1; i + 1 < n; ++i)
    w[i] = 0.5 * (x[i+1] - x[i-1]);

  return w;
}


inline double quadratureTrapezoidal(
  const std::vector<double> &x,
  const std::vector<double> &y,
  const size_t idx_start,
  const size_t idx_end)
{

  if (x.size() != y.size())
    throw std::logic_error("trapezoidal quadrature: x and y must be the same size!\n");

  if (x.size() == 1)
    return y[0];


  double sum = 0.0;

  for (size_t i = idx_start+1; i < idx_end+1; ++i)
      sum += (x[i] - x[i-1]) * (y[i] + y[i-1]);


  return sum * 0.5;
}



}
}



#endif

