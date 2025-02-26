/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_CONSTANTS_HPP_INCLUDED
#define BURNMAN_UTILS_CONSTANTS_HPP_INCLUDED

#include <limits>

/**
 * Constants used in burnman
 */
namespace constants {

  /**
   * Physical constants
   * Mostly CODATA 2022 values
   */
  namespace physics {
    /**
     * R in [J/mol/K]
     */
    constexpr double gas_constant = 8.314462618;
    /**
     * k_B in [J/K]
     */
    constexpr double boltzmann = 1.380649e-23;
    /**
     * N_A [1/mol]
     */
    constexpr double avogadro = 6.02214076e23;
    /**
     * Reduced planck constant (Dirac) in Js.
     */
    constexpr double dirac = 1.054571817e-34;
    /**
     * Newtonian constant of gravitation, G in [m^3/kg/s^2].
     */
     constexpr double gravitation = 6.67430e-11;
  }

  /**
   * Machine epsilon, etc.
   */
  namespace precision {
    constexpr double double_eps = std::numeric_limits<double>::epsilon();
    constexpr double rel_tolerance_eps = 4.0 * double_eps;
    constexpr double abs_tolerance = 1.0e-12;
  }
}

#endif // BURNMAN_UTILS_CONSTANTS_HPP_INCLUDED
