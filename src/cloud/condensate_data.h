/*
* This file is part of the ngam code.
* Copyright (C) 2026 Daniel Kitzmann
*
* ngam is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#ifndef _condensate_data_h
#define _condensate_data_h

#include <string>
#include <vector>


namespace ngam {


// Bulk properties of the condensates the cloud modules know by name: the refractive-index file in
// data/refractive_indices/ (the Kitzmann & Heng 2018 compilation, references in each file's header)
// and the bulk density [g/cm^3] (virga's gas_properties where available, otherwise handbook values).
// A cloud module's `condensate` parameter is either one of these names or a path to a refractive-
// index file of the same format (then the density has to be given explicitly where it is needed).
struct CondensateData {
  std::string name;
  std::string refractive_index_file;   // file name in data/refractive_indices/
  double bulk_density;                 // g/cm^3
};


namespace constants {

const std::vector<CondensateData> condensate_data {
  {"Al2O3",         "Al2O3.dat",                      3.987},
  {"C",             "C.dat",                          2.26},    // graphite
  {"CaTiO3",        "CaTiO3.dat",                     3.987},
  {"Cr",            "Cr.dat",                         7.15},
  {"Fe",            "Fe.dat",                         7.875},
  {"FeO",           "FeO.dat",                        5.74},
  {"FeS",           "FeS.dat",                        4.84},
  {"Fe2O3",         "Fe2O3.dat",                      5.24},
  {"Fe2SiO4",       "Fe2SiO4.dat",                    4.39},
  {"H2O_ice",       "H2O_s.dat",                      0.93},
  {"KCl",           "KCl.dat",                        1.99},
  {"Mg2SiO4",       "Mg2SiO4_amorph_sol-gel.dat",     3.214},
  {"MgAl2O4",       "MgAl2O4.dat",                    3.64},
  {"MgFeSiO4",      "MgFeSiO4_amorph_glass.dat",      3.8},
  {"MgO",           "MgO.dat",                        3.58},
  {"MgSiO3",        "MgSiO3_amorph_glass.dat",        3.192},
  {"MgSiO3_sol-gel","MgSiO3_amorph_sol-gel.dat",      3.192},
  {"Mg08Fe02SiO3",  "Mg08Fe02SiO3_amorph_glass.dat",  3.3},
  {"Mg05Fe05SiO3",  "Mg05Fe05SiO3_amorph_glass.dat",  3.5},
  {"Mg04Fe06SiO3",  "Mg04Fe06SiO3_amorph_glass.dat",  3.6},
  {"Mg08Fe12SiO4",  "Mg08Fe12SiO4_amorph_glass.dat",  3.9},
  {"MnS",           "MnS.dat",                        4.0},
  {"Na2S",          "Na2S.dat",                       1.856},
  {"NaCl",          "NaCl.dat",                       2.165},
  {"SiC",           "SiC.dat",                        3.21},
  {"SiO",           "SiO.dat",                        2.13},
  {"SiO2",          "SiO2_amorph.dat",                2.648},
  {"SiO2_alpha",    "SiO2_alpha.dat",                 2.648},
  {"TiC",           "TiC.dat",                        4.93},
  {"TiO2",          "TiO2_anatase.dat",               4.25},
  {"tholin",        "Titan_tholin.dat",               1.0},
  {"ZnS",           "ZnS.dat",                        4.04},
};

}


// the condensate entry for a name, or nullptr if the name is not in the table
inline const CondensateData* findCondensate(const std::string& name)
{
  for (const auto& c : constants::condensate_data)
    if (c.name == name) return &c;
  return nullptr;
}


inline std::string condensateNames()
{
  std::string s;
  for (const auto& c : constants::condensate_data) s += (s.empty() ? "" : ", ") + c.name;
  return s;
}


// Resolve a `condensate` parameter to a refractive-index file: a table name maps to the shipped
// data file, anything else is taken as a path. Defined in cloud_optics.cpp (needs the data
// directory compiled into the library).
std::string refractiveIndexFile(const std::string& condensate);


}
#endif
