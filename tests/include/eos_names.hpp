/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef TESTS_EOS_ENUM_STRING_HPP
#define TESTS_EOS_ENUM_STRING_HPP

#include <string>
#include "burnman/utils/types/simple_types.hpp"

inline constexpr const char* eos_string(burnman::types::EOSType eos_type) {
  switch (eos_type) {
    case burnman::types::EOSType::BM2:             return "BM2";
    case burnman::types::EOSType::BM3:             return "BM3";
    case burnman::types::EOSType::MT:              return "MT";
    case burnman::types::EOSType::Vinet:           return "Vinet";
    case burnman::types::EOSType::MGD2:            return "MGD2";
    case burnman::types::EOSType::MGD3:            return "MGD3";
    case burnman::types::EOSType::SLB2:            return "SLB2";
    case burnman::types::EOSType::SLB3:            return "SLB3";
    case burnman::types::EOSType::SLB3Conductive:  return "SLB3Conductive";
    case burnman::types::EOSType::HPTMT:           return "HPTMT";
    default:                                       return "Unknown";
  }
}

#endif // TESTS_EOS_ENUM_STRING_HPP
