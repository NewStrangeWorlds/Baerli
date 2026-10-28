/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#include "sampled_data.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <assert.h>
#include <cmath>


namespace baerli{

void CrossSectionFile::loadFile()
{
  std::fstream file;


  file.open(filename.c_str(), std::ios::binary | std::ios::in | std::ios::ate);

  if (file.fail()) std::cout << "cross section file " << filename << " not found! :(((( \n";
  
  assert (!file.fail( ));


  int nb_data_points = file.tellg()/sizeof(float);
  file.seekg(0, std::ios::beg);

  
  cross_sections.resize(nb_data_points);


  for (int i=0; i<nb_data_points; ++i)
  {
    float x;

    file.read((char *) &x, sizeof x);

    cross_sections[i] = x;

    if (is_data_log) cross_sections[i] = std::pow(10.0, cross_sections[i]);
  }


  file.close();

  is_loaded = true;
}


void CrossSectionFile::unloadData()
{
  std::vector<double>().swap (cross_sections);

  is_loaded = false;
}


}