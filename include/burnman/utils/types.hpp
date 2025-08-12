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
#include <vector>

/**
 * Custom type alias for chemical formulae.
 */
using FormulaMap = std::unordered_map<std::string, double>;

// Overloaded operators to add and weight formulae
inline FormulaMap operator+(const FormulaMap& a, const FormulaMap& b) {
  FormulaMap result = a;
  for (const auto& [elem, count] : b) {
    result[elem] += count;
  }
  return result;
}

inline FormulaMap& operator+=(FormulaMap& a, const FormulaMap& b) {
  for (const auto& [elem, count] : b) {
    a[elem] += count;
  }
  return a;
}

inline FormulaMap operator*(const FormulaMap& a, double scalar) {
  FormulaMap result;
  for (const auto& [elem, count] : a) {
    result[elem] = scalar * count;
  }
  return result;
}

inline FormulaMap operator*(double scalar, const FormulaMap& a) {
  return a * scalar;
}

inline FormulaMap& operator*=(FormulaMap& a, double scalar) {
  for (auto& [_, count] : a) {
    count *= scalar;
  }
  return a;
}

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

/**
 * Enum used to define averaging scheme type
 */
enum class AveragingType {
  Voigt,
  Reuss,
  VRH,
  HashinShtrikmanLower,
  HashinShtrikmanUpper,
  HashinShtrikman
};

/**
 * Enum used to define fraction type
 */
enum class FractionType {
  Molar,
  Mass,
  Volume
};

/**
 * Forward declarations
 */
class Mineral;
class Solution;
class Assemblage;

/**
 * Type aliases for constructing solid solution models.
 */
using MineralFormulaPair = std::pair<Mineral, std::string>;
using PairedEndmemberList = std::vector<MineralFormulaPair>;

/**
 * Type aliases for phases when constructing Assemblage
 */
using PhaseVariant = std::variant<
  Mineral,
  Solution,
  Assemblage
>;
using PhaseVariantVector = std::vector<PhaseVariant>;

#endif // BURNMAN_UTILS_TYPES_HPP_INCLUDED
