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
#include "burnman/core/material.hpp"
#include "burnman/core/equation_of_state.hpp"
#include "burnman/utils/eos.hpp"
#include "burnman/utils/exceptions.hpp"
#include <string>
#include <memory>

TEST_CASE("Get/Set state", "[material][core]") {
  double test_P = 24e9;
  double test_T = 2000.0;
  Material test_material;
  // Ensure P/T set correctly
  test_material.set_state(test_P, test_T);
  CHECK(test_material.get_pressure() == test_P);
  CHECK(test_material.get_temperature() == test_T);
  // Ensure reset doesn't wipe P & T
  test_material.reset();
  CHECK(test_material.get_pressure() == test_P);
  CHECK(test_material.get_temperature() == test_T);
}

TEST_CASE("Get/Set name", "[material][core]") {
  Material test_material;
  std::string test_name = "Test Name";
  std::string new_test_name = "Another Name";
  // Ensure default class name
  CHECK_FALSE(test_material.get_name().empty());
  CHECK_FALSE(test_material.get_name() == test_name);
  // Esure name can be set
  test_material.set_name(test_name);
  CHECK(test_material.get_name() == test_name);
  // Make sure it updates
  test_material.set_name(new_test_name);
  CHECK(test_material.get_name() == new_test_name);
  // Ensure reset() does not remove name
  test_material.reset();
  CHECK(test_material.get_name() == new_test_name);
}

TEST_CASE("Ensure default errors", "[material][core]") {
  Material test_material;
  CHECK_THROWS_AS(test_material.set_method(EOSType::Auto), NotImplementedError);
  CHECK_THROWS_AS(test_material.set_method(std::make_shared<EquationOfState>()), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_internal_energy(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_gibbs(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_helmholtz(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_mass(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_volume(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_volume_unmodified(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_density(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_entropy(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_enthalpy(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_isothermal_bulk_modulus_reuss(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_isentropic_bulk_modulus_reuss(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_isothermal_compressibility_reuss(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_isentropic_compressibility_reuss(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_shear_modulus(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_p_wave_velocity(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_bulk_sound_velocity(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_shear_wave_velocity(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_grueneisen_parameter(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_thermal_expansivity(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_heat_capacity_v(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_molar_heat_capacity_p(), NotImplementedError);
  CHECK_THROWS_AS(test_material.get_isentropic_thermal_gradient(), NotImplementedError);
}
