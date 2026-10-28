/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#ifndef _fall_velocity_h
#define _fall_velocity_h

#include <functional>


namespace ngam {


// Gas properties a falling particle sees (all cgs).
struct FallMedium {
  double temperature;        // K
  double pressure;           // dyn/cm^2
  double mean_molecular_weight;   // amu
  double gravity;            // cm/s^2
  double mean_free_path;     // cm
  double viscosity;          // dyn s/cm^2
};


// Mean free path and dynamic viscosity of the gas from kinetic theory with a Lennard-Jones
// correction (Ackerman & Marley 2001 eq. B2, Rosner 2000), as in virga:
//   mfp = 1 / (sqrt(2) n pi d^2),  eta = 5/16 sqrt(pi k T m) / (pi d^2) / (1.22 (T / eps_k)^-0.16)
// d: molecular diameter [cm], eps_k: Lennard-Jones well depth / k [K] (H2: 2.827e-8 cm, 59.7 K).
FallMedium fallMedium(
  const double temperature, const double pressure_dyn, const double mean_molecular_weight,
  const double gravity, const double molecule_diameter, const double lj_temperature);


// Terminal fall speed [cm/s] of a sphere of radius r [cm] and bulk density rho_p [g/cm^3]: Stokes
// flow with the Cunningham slip correction 1 + 1.26 Kn for Re < 1, the Pruppacher & Klett (1978)
// polynomial fit in ln(C_d Re^2) for 1 < Re <= 1e3, and a constant drag coefficient 0.45 above --
// virga's vfall (A. Ackerman's original, "vfall_legacy").
double fallVelocity(const double radius, const double bulk_density, const FallMedium& medium);


// Brent's root finder on [a, b] (f(a) and f(b) of opposite sign); relative tolerance on x.
double brentRoot(
  const std::function<double(double)>& f, double a, double b,
  const double rel_tol = 1e-10, const int max_iter = 200);


}
#endif
