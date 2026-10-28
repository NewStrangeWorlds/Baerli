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
// data/refractive_indices/ (the Kitzmann & Heng 2018 compilation, references in each file's header),
// the bulk density [g/cm^3] (virga's gas_properties where available, otherwise handbook values) and
// the symbol of the condensate in FastChem's condensate data (fastchem_data/logK_condensates.dat),
// through which the Ackerman & Marley module reads the condensation increments.
// A cloud module's `condensate` parameter is either one of these names or a path to a refractive-
// index file of the same format (then the density has to be given explicitly where it is needed).
struct CondensateData {
  std::string name;
  std::string refractive_index_file;   // file name in data/refractive_indices/
  double bulk_density;                 // g/cm^3
  std::string fastchem_symbol;         // the condensate in FastChem's data ("" = not available)
};


namespace constants {

const std::vector<CondensateData> condensate_data {
  {"Al2O3",         "Al2O3.dat",                      3.987, "Al2O3(s,l)"},
  {"C",             "C.dat",                          2.26, "C(s)"},    // graphite
  {"CaTiO3",        "CaTiO3.dat",                     3.987, "CaTiO3(s)"},
  {"Cr",            "Cr.dat",                         7.15, "Cr(s,l)"},
  {"Fe",            "Fe.dat",                         7.875, "Fe(s,l)"},
  {"FeO",           "FeO.dat",                        5.74, "FeO(s,l)"},
  {"FeS",           "FeS.dat",                        4.84, "FeS(s,l)"},
  {"Fe2O3",         "Fe2O3.dat",                      5.24, "Fe2O3(s)"},
  {"Fe2SiO4",       "Fe2SiO4.dat",                    4.39, "Fe2SiO4(s)"},
  {"H2O_ice",       "H2O_s.dat",                      0.93, "H2O(s,l)"},
  {"KCl",           "KCl.dat",                        1.99, "KCl(s,l)"},
  {"Mg2SiO4",       "Mg2SiO4_amorph_sol-gel.dat",     3.214, "Mg2SiO4(s,l)"},
  {"MgAl2O4",       "MgAl2O4.dat",                    3.64, "MgAl2O4(s,l)"},
  {"MgFeSiO4",      "MgFeSiO4_amorph_glass.dat",      3.8, ""},
  {"MgO",           "MgO.dat",                        3.58, "MgO(s,l)"},
  {"MgSiO3",        "MgSiO3_amorph_glass.dat",        3.192, "MgSiO3(s,l)"},
  {"MgSiO3_sol-gel","MgSiO3_amorph_sol-gel.dat",      3.192, "MgSiO3(s,l)"},
  {"Mg08Fe02SiO3",  "Mg08Fe02SiO3_amorph_glass.dat",  3.3, ""},
  {"Mg05Fe05SiO3",  "Mg05Fe05SiO3_amorph_glass.dat",  3.5, ""},
  {"Mg04Fe06SiO3",  "Mg04Fe06SiO3_amorph_glass.dat",  3.6, ""},
  {"Mg08Fe12SiO4",  "Mg08Fe12SiO4_amorph_glass.dat",  3.9, ""},
  {"MnS",           "MnS.dat",                        4.0, "MnS(s)"},
  {"Na2S",          "Na2S.dat",                       1.856, "Na2S(s,l)"},
  {"NaCl",          "NaCl.dat",                       2.165, "NaCl(s,l)"},
  {"SiC",           "SiC.dat",                        3.21, "SiC(s)"},
  {"SiO",           "SiO.dat",                        2.13, "SiO(s)"},
  {"SiO2",          "SiO2_amorph.dat",                2.648, "SiO2(s,l)"},
  {"SiO2_alpha",    "SiO2_alpha.dat",                 2.648, "SiO2(s,l)"},
  {"TiC",           "TiC.dat",                        4.93, "TiC(s,l)"},
  {"TiO2",          "TiO2_anatase.dat",               4.25, "TiO2(s,l)"},
  {"tholin",        "Titan_tholin.dat",               1.0, ""},
  {"ZnS",           "ZnS.dat",                        4.04, ""},
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
