/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_SOLUTION_HPP_INCLUDED
#define BURNMAN_CORE_SOLUTION_HPP_INCLUDED

#include "burnman/core/mineral.hpp"

/**
 * Base class for solutions.

  TODO
  
 */ 
class Solution : public Mineral {

 public:

  virtual ~Solution() = default;

 protected:
  ;
 private:
  ;

};

#endif // BURNMAN_CORE_SOLUTION_HPP_INCLUDED