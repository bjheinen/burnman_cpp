/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_WARNINGS_HPP_INCLUDED
#define BURNMAN_UTILS_WARNINGS_HPP_INCLUDED

#include <iostream>
#include <string>

namespace burnman {
  namespace utils {
    /**
     * @brief Flag to suppress warnings across burnman.
     */
    inline bool supress_warnings = false;

    /**
     * @brief Prints a warning message.
     *
     * @param message The warning message to print.
     */
    inline void warn(const std::string& message) {
      // Simple implementation: print to std::cerr
      // Could switch to a logging framework if needed
      if (!supress_warnings) {
        std::cerr << "Warning: " << message << std::endl;
      }
    }
  } // namespace utils
} // namespace burnman

#endif // BURNMAN_UTILS_WARNINGS_HPP_INCLUDED
