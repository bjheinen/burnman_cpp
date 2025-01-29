/*
  TODO: Copyright Notice!
*/
#include <typeinfo>
#include "../include/burnman/classes/material.hpp"
#include "../include/burnman/util/exceptions.hpp"

[[noreturn]] void Material::throw_not_implemented_error(const std::string& method) const {
  throw NotImplementedError(get_class_name(), method);
}

std::string Material::get_class_name() const {
  return typeid(*this).name();
}

