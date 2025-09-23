/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_STRING_UTILS_INCLUDED
#define BURNMAN_UTILS_STRING_UTILS_INCLUDED

#include <cstddef>
#include <string>
#include <regex>
#include <stdexcept>
#include <vector>

namespace burnman {
namespace utils {

  /**
   * @brief Extracts any leading number from a string.
   *
   * Gets numeric prefix
   *   (e.g.,"2SiO2" --> 2, "0.45SiO2" --> "0.45",
   *   "1/2Fe2O3" --> "1/2", "" or "Al2O3" --> "")
   * Works with e-2 but not E-2!
   *
   * @param s The string.
   * @return Prefix string (empty string if no prefix)
   */
  inline std::string extract_numeric_prefix(const std::string& s) {
    std::regex element_re("[A-Z][^A-Z]*");
    std::sregex_token_iterator it(s.begin(), s.end(), element_re, -1);
    std::sregex_token_iterator end;
    return (it != end) ? std::string(*it) : "";
  }

  /**
   * @brief Extended std::stod that works for fractions.
   *
   * e.g. "3" --> 3, "1e-2" --> 0.01, "1/2" --> 0.5
   *
   * @param s The string to convert.
   * @return Numerical value
   */
  inline double stod(const std::string& s) {
    auto slash_pos = s.find('/');
    // If no fraction, just use std::stod
    if (slash_pos == std::string::npos) {
      return std::stod(s);
    }
    // Else, split on position of / and use std::pos on both
    std::string numerator_str = s.substr(0, slash_pos);
    std::string denominator_str = s.substr(slash_pos + 1);
    double numerator = std::stod(numerator_str);
    double denominator = std::stod(denominator_str);
    if (denominator == 0.0) {
      throw std::invalid_argument("Division by zero in: " + s);
    }
    return numerator / denominator;
  }

  /**
   * @brief Convert value to string using ostringstream.
   *
   * Works for types that implement << operator.
   * Safe for use on strings/char* (unlike std::to_string).
   *
   * @param val The value to convert.
   * @return The string representation of the value.
   */
  template <typename T>
  inline std::string to_string(const T& val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
  }

  /**
   * @brief Join vector of strings (a la Python .join())
   *
   * @param string_list std::vector of string elements
   * @param delimiter optional delimiter to join strings
   *
   * @return Concatenated string
   */
  inline std::string join(
    const std::vector<std::string>& string_list,
    const std::string& delimiter = ""
  ) {
    std::string result;
    for (std::size_t i = 0; i < string_list.size(); ++i) {
      result += string_list[i];
      if (i + 1 < string_list.size()) {
        result += delimiter;
      }
    }
    return result;
  }

} // namespace utils
} // namespace burnman

#endif // BURNMAN_UTILS_STRING_UTILS_INCLUDED
