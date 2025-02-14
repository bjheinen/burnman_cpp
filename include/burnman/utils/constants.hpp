/*
  TODO: Copyright Notice!
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
  }
}

#endif // BURNMAN_UTILS_CONSTANTS_HPP_INCLUDED
