/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_EXCEPTIONS_HPP_INCLUDED
#define BURNMAN_UTILS_EXCEPTIONS_HPP_INCLUDED

#include <stdexcept>
#include <string>

/**
 * Custom exception derived from std::logic_error
 * Usage: throw NotImplementedError(class_name, func_name)
 *        class_name from e.g. typeid(*this).name();
 *        func_name from e.g. __func__ identifier
 */
class NotImplementedError : public std::logic_error {
 public:
  explicit NotImplementedError(const std::string& class_name,
                               const std::string& func_name,
                               const std::string& message = "Function not implemented!")
    : std::logic_error("[" + class_name + "::" + func_name + "] " + message) {}
};

#endif  // BURNMAN_UTILS_EXCEPTIONS_HPP_INCLUDED