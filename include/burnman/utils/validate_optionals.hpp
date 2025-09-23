/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_VALIDATE_OPTIONALS_INCLUDED
#define BURNMAN_UTILS_VALIDATE_OPTIONALS_INCLUDED

#include <optional>
#include <stdexcept>
#include <sstream>
#include <string>
#include "burnman/utils/string_utils.hpp"
#include "burnman/utils/warnings.hpp"

namespace burnman {
namespace utils {

/**
 * @brief Checks that an std::optional is set, else throws.
 *
 * @param opt The std::optional to check.
 * @param name Name of the parameter (for error message).
 * @param msg Optional custom message prefix (default: "Missing required parameter: ").
 * @throw std::invalid_argument if opt is not set.
 */
template <typename T>
void require_set(
  const std::optional<T>& opt,
  const std::string& name,
  const std::string& msg = "Missing required parameter: "
) {
  if (!opt.has_value()) {
    throw std::invalid_argument(
      "Missing required parameter: " + name
    );
  }
}

/**
 * @brief Sets a default value if an std::optional is not set.
 *
 * @param opt The std::optional to check/set.
 * @param value The default value to set if opt is not set.
 * @param verbose If true, prints a message when setting the default (default: false).
 * @param name Name of the parameter (for warning message).
 */
template <typename T, typename U>
void fallback_to_default(
  std::optional<T>& opt,
  const U& value,
  bool verbose = false,
  const std::string& name = "parameter"
) {
  if (!opt.has_value()) {
    opt = static_cast<T>(value);
    if (verbose) {
      utils::warn("Setting default value for " + name + ": " + utils::to_string(value));
    }
  }
}

template <typename T>
bool check_in_range(
  const std::optional<T>& opt,
  const T& min_value,
  const T& max_value,
  const std::string& name = "parameter",
  bool strict = false
) {
  const T& val = *opt;
  if (val < min_value || val > max_value) {
    std::string message =
      "Unusual value for " + name + ": "
      + std::to_string(val) + " is out of range "
      + "{ " + std::to_string(min_value) + " to " + std::to_string(max_value) + "}";
    if (strict) {
      throw std::runtime_error(message);
    } else {
      utils::warn(message);
    }
    return false;
  } else {
    return true;
  }
}

} // namespace utils
} // namespace burnman

#endif // BURNMAN_UTILS_VALIDATE_OPTIONALS_INCLUDED
