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
#include <string>
#include <vector>

/**
 * Constants used in burnman
 */
namespace constants {

  /**
   * Physical constants
   * Using CODATA 2022 values
   * n.b. To match the python implementation / scipy > 1.15 we use
   *      more digits for R and hbar than the CODATA2022 database.
   *      R is defined as k_B * N_A, and hbar as h / (2pi), so we
   *      can calculate them to full double precision.
   */
  namespace physics {
    /**
     * R in [J/mol/K]
     * CODATA2022: 8.314 462 618 ... J mol-1 K-1
     */
    constexpr double gas_constant = 8.31446261815324;
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
     * CODATA2022: 1.054 571 817 ...  x 10-34 J s
     */
    constexpr double dirac = 1.0545718176461565e-34;
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
    constexpr double inverseish_eps = 1.0e-5;
    constexpr double logish_eps = 1.0e-7;
  }

  /**
   * Elements etc.
   */
  namespace chemistry {

    /**
     * IUPAC_element_order provides a list of all the elements.
     * Element order is based loosely on electronegativity,
     * following the scheme suggested by IUPAC, except that H
     * comes after the Group 16 elements, not before them.
    */
    inline const std::vector<std::string> IUPAC_element_order = {
      "v", "Og", "Rn", "Xe", "Kr", "Ar", "Ne", "He", // Group 18
      "Fr", "Cs", "Rb", "K", "Na", "Li", // Group 1 (not H)
      "Ra", "Ba", "Sr", "Ca", "Mg", "Be", // Group 2
      "Lr", "No", "Md", "Fm", "Es", "Cf", "Bk", "Cm",
      "Am", "Pu", "Np", "U", "Pa", "Th", "Ac", // Actinides
      "Lu", "Yb", "Tm", "Er", "Ho", "Dy", "Tb", "Gd", "Eu",
      "Sm", "Pm", "Nd", "Pr", "Ce", "La", // Lanthanides
      "Y", "Sc", // Group 3
      "Rf", "Hf", "Zr", "Ti", // Group 4
      "Db", "Ta", "Nb", "V", // Group 5
      "Sg", "W", "Mo", "Cr", // Group 6
      "Bh", "Re", "Tc", "Mn", // Group 7
      "Hs", "Os", "Ru", "Fe", // Group 8
      "Mt", "Ir", "Rh", "Co", // Group 9
      "Ds", "Pt", "Pd", "Ni", // Group 10
      "Rg", "Au", "Ag", "Cu", // Group 11
      "Cn", "Hg", "Cd", "Zn", // Group 12
      "Nh", "Tl", "In", "Ga", "Al", "B", // Group 13
      "Fl", "Pb", "Sn", "Ge", "Si", "C", // Group 14
      "Mc", "Bi", "Sb", "As", "P", "N", // Group 15
      "Lv", "Po", "Te", "Se", "S", "O", // Group 16
      "H", // Hydrogen
      "Ts", "At", "I", "Br", "Cl", "F" // Group 17
    };

  }

}

#endif // BURNMAN_UTILS_CONSTANTS_HPP_INCLUDED
