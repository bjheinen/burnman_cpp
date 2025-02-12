/*
  TODO: Copyright Notice!
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