/*
 * This file is part of Bärli, a radiative-convective atmosphere model.
 * Copyright (C) 2024-2026 Daniel Kitzmann
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 * See the LICENSE file in the repository root for the full license text.
 */

#ifndef _EXCEPTIONS_H
#define	_EXCEPTIONS_H


#include <iostream>
#include <exception>
#include <stdexcept>


namespace baerli{

class FileNotFound : public std::runtime_error {
  public:
    FileNotFound(
      const std::string where, 
      const std::string what) : std::runtime_error("Critical Error - Aborting!") {
        std::cout << "File " << what << " in " << where << " not found!\n";
      }
};


class InvalidInput : public std::runtime_error {
  public:
    InvalidInput(
      const std::string where,
      const std::string what) : std::runtime_error("Invalid input in " + where + ": " + what) {
        std::cout << "Invalid input in " << where << ":\n" << what << "\n";
      }
};


}

#endif
 
