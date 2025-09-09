/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <catch2/catch_test_macros.hpp>
#include "burnman/utils/exceptions.hpp"
#include <stdexcept>

TEST_CASE("NotImplementedError message", "[utils][exceptions]") {
  const std::string class_name = "AClass";
  const std::string func_name = "some_function";
  const std::string message = "Hello!";
  try {
    throw NotImplementedError(class_name, func_name, message);
  } catch (const NotImplementedError& e) {
    REQUIRE(std::string(e.what()) ==
            "[AClass::some_function] Hello!");
  } catch (...) {
    FAIL("NotImplementedError not thrown?");
  }
}

TEST_CASE("NotImplementedError default message", "[utils][exceptions]") {
  try {
    throw NotImplementedError("Class", "func");
  } catch (const std::logic_error& e) {
    REQUIRE(std::string(e.what()) ==
            "[Class::func] Function not implemented!");
  }
}

TEST_CASE("NotImplementedError derived from std::logic_error", "[utils][exceptions]") {
  NotImplementedError err("Class", "func", "msg");
  std::logic_error* base_ptr = &err;
  REQUIRE(base_ptr != nullptr);
}
