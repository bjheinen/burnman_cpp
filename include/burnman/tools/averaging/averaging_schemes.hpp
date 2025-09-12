/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_TOOLS_AVERAGING_AVERAGING_SCHEMES_HPP_INCLUDED
#define BURNMAN_TOOLS_AVERAGING_AVERAGING_SCHEMES_HPP_INCLUDED

// Convenience header to get all averaging schemes
#include "burnman/tools/averaging/averaging_base.hpp"
#include "burnman/tools/averaging/averaging_utils.hpp"
#include "burnman/tools/averaging/voigt.hpp"
#include "burnman/tools/averaging/reuss.hpp"
#include "burnman/tools/averaging/vrh.hpp"
#include "burnman/tools/averaging/hs_lower.hpp"
#include "burnman/tools/averaging/hs_upper.hpp"
#include "burnman/tools/averaging/hs.hpp"

// TODO: add to tools group / module in docs
namespace burnman {
  /**
   * @namespace burnman::averaging
   * @brief Tools and classes for averaging composite properties.
   *
   * This namespace collects classes that provide functions for
   * averaging elastic moduli and thermodynamic properties.
   * Implemented averaging schemes are:
   * - Voigt
   * - Reuss
   * - Voigt-Reuss-Hill
   * - Hashin-Shtrikman
   */
  namespace averaging {
    // blank for documentation only
  } // namespace averaging
} // namespace burnman

#endif // BURNMAN_TOOLS_AVERAGING_SCHEMES_HPP_INCLUDED
