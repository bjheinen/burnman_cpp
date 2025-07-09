/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_TYPES_HPP_INCLUDED
#define BURNMAN_UTILS_TYPES_HPP_INCLUDED

#include <string>
#include <unordered_map>

// Type alias for formula map
using FormulaMap = std::unordered_map<std::string, int>;

/**
 * Enum used to define EOS Type
 */
enum class EOSType {
  Auto, // Used to set EOS from params
  Custom, // Used when user passes custom EOSType
  Vinet,
  BM3,
  BM2,
  MGD2,
  MGD3,
  SLB2,
  SLB3,
  SLB3Conductive
};

#endif // BURNMAN_UTILS_TYPES_HPP_INCLUDED
