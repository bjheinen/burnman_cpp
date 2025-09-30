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
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "burnman/core/assemblage.hpp"
#include "burnman/utils/types/simple_types.hpp"
#include "solution_fixtures.hpp"

using namespace burnman;

TEST_CASE_METHOD(PyroliteAssemblageFixture, "Assemblage benchmarks", "[core][assemblage][!benchmark]") {
  // Setup assemblage object
  double P = 50.e9;
  double T = 2000.0;
  assemblage.set_method(types::EOSType::Auto);
  assemblage.set_state(P, T);
  // Run benchmark for clear_properties (runs reset_cache internally)
  BENCHMARK("clear_computed_properties") {
    assemblage.clear_computed_properties();
    return 0;
  };

  // Benchmark individual functions
  BENCHMARK("set_state") {
    assemblage.clear_computed_properties();
    return assemblage.set_state(P, T);
  };

  BENCHMARK("get_formula") {
    assemblage.clear_computed_properties();
    return assemblage.get_formula();
  };

  BENCHMARK("get_endmember_formulae") {
    assemblage.clear_computed_properties();
    return assemblage.get_endmember_formulae();
  };

  BENCHMARK("get_n_endmembers") {
    assemblage.clear_computed_properties();
    return assemblage.get_n_endmembers();
  };

  BENCHMARK("get_endmembers_per_phase") {
    assemblage.clear_computed_properties();
    return assemblage.get_endmembers_per_phase();
  };

  BENCHMARK("get_n_elements") {
    assemblage.clear_computed_properties();
    return assemblage.get_n_elements();
  };

  BENCHMARK("get_n_reactions") {
    assemblage.clear_computed_properties();
    return assemblage.get_n_reactions();
  };

  BENCHMARK("get_elements") {
    assemblage.clear_computed_properties();
    return assemblage.get_elements();
  };

  BENCHMARK("get_independent_element_indices") {
    assemblage.clear_computed_properties();
    return assemblage.get_independent_element_indices();
  };

  BENCHMARK("get_dependent_element_indices") {
    assemblage.clear_computed_properties();
    return assemblage.get_dependent_element_indices();
  };

  BENCHMARK("get_stoichiometric_matrix") {
    assemblage.clear_computed_properties();
    return assemblage.get_stoichiometric_matrix();
  };

  BENCHMARK("get_compositional_basis") {
    assemblage.clear_computed_properties();
    return assemblage.get_compositional_basis();
  };

  BENCHMARK("get_compositional_null_basis") {
    assemblage.clear_computed_properties();
    return assemblage.get_compositional_null_basis();
  };

  BENCHMARK("get_reaction_basis") {
    assemblage.clear_computed_properties();
    return assemblage.get_reaction_basis();
  };

  BENCHMARK("get_molar_internal_energy") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_internal_energy();
  };

  BENCHMARK("get_molar_gibbs") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_gibbs();
  };

  BENCHMARK("get_molar_helmholtz") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_helmholtz();
  };

  BENCHMARK("get_molar_mass") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_mass();
  };

  BENCHMARK("get_molar_volume") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_volume();
  };

  BENCHMARK("get_density") {
    assemblage.clear_computed_properties();
    return assemblage.get_density();
  };

  BENCHMARK("get_molar_entropy") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_entropy();
  };

  BENCHMARK("get_molar_enthalpy") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_enthalpy();
  };

  BENCHMARK("get_isothermal_bulk_modulus_reuss") {
    assemblage.clear_computed_properties();
    return assemblage.get_isothermal_bulk_modulus_reuss();
  };

  BENCHMARK("get_isentropic_bulk_modulus_reuss") {
    assemblage.clear_computed_properties();
    return assemblage.get_isentropic_bulk_modulus_reuss();
  };

  BENCHMARK("get_isothermal_compressibility_reuss") {
    assemblage.clear_computed_properties();
    return assemblage.get_isothermal_compressibility_reuss();
  };

  BENCHMARK("get_isentropic_compressibility_reuss") {
    assemblage.clear_computed_properties();
    return assemblage.get_isentropic_compressibility_reuss();
  };

  BENCHMARK("get_shear_modulus") {
    assemblage.clear_computed_properties();
    return assemblage.get_shear_modulus();
  };

  BENCHMARK("get_p_wave_velocity") {
    assemblage.clear_computed_properties();
    return assemblage.get_p_wave_velocity();
  };

  BENCHMARK("get_bulk_sound_velocity") {
    assemblage.clear_computed_properties();
    return assemblage.get_bulk_sound_velocity();
  };

  BENCHMARK("get_shear_wave_velocity") {
    assemblage.clear_computed_properties();
    return assemblage.get_shear_wave_velocity();
  };

  BENCHMARK("get_grueneisen_parameter") {
    assemblage.clear_computed_properties();
    return assemblage.get_grueneisen_parameter();
  };

  BENCHMARK("get_thermal_expansivity") {
    assemblage.clear_computed_properties();
    return assemblage.get_thermal_expansivity();
  };

  BENCHMARK("get_molar_heat_capacity_v") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_heat_capacity_v();
  };

  BENCHMARK("get_molar_heat_capacity_p") {
    assemblage.clear_computed_properties();
    return assemblage.get_molar_heat_capacity_p();
  };

}

TEST_CASE_METHOD(PyroliteAssemblageFixture, "Assemblage benchmarks - soft reset", "[core][assemblage][!benchmark]") {
  // Setup assemblage object
  double P = 50.e9;
  double T = 2000.0;
  assemblage.set_method(types::EOSType::Auto);
  assemblage.set_state(P, T);
  // Run benchmark for clear_properties (runs reset_cache internally)
  BENCHMARK("reset_cache") {
    assemblage.reset_cache();
    return 0;
  };

  // Benchmark individual functions
  BENCHMARK("set_state") {
    assemblage.reset_cache();
    return assemblage.set_state(P, T);
  };

  BENCHMARK("get_formula") {
    assemblage.reset_cache();
    return assemblage.get_formula();
  };

  BENCHMARK("get_endmember_formulae") {
    assemblage.reset_cache();
    return assemblage.get_endmember_formulae();
  };

  BENCHMARK("get_n_endmembers") {
    assemblage.reset_cache();
    return assemblage.get_n_endmembers();
  };

  BENCHMARK("get_endmembers_per_phase") {
    assemblage.reset_cache();
    return assemblage.get_endmembers_per_phase();
  };

  BENCHMARK("get_n_elements") {
    assemblage.reset_cache();
    return assemblage.get_n_elements();
  };

  BENCHMARK("get_n_reactions") {
    assemblage.reset_cache();
    return assemblage.get_n_reactions();
  };

  BENCHMARK("get_elements") {
    assemblage.reset_cache();
    return assemblage.get_elements();
  };

  BENCHMARK("get_independent_element_indices") {
    assemblage.reset_cache();
    return assemblage.get_independent_element_indices();
  };

  BENCHMARK("get_dependent_element_indices") {
    assemblage.reset_cache();
    return assemblage.get_dependent_element_indices();
  };

  BENCHMARK("get_stoichiometric_matrix") {
    assemblage.reset_cache();
    return assemblage.get_stoichiometric_matrix();
  };

  BENCHMARK("get_compositional_basis") {
    assemblage.reset_cache();
    return assemblage.get_compositional_basis();
  };

  BENCHMARK("get_compositional_null_basis") {
    assemblage.reset_cache();
    return assemblage.get_compositional_null_basis();
  };

  BENCHMARK("get_reaction_basis") {
    assemblage.reset_cache();
    return assemblage.get_reaction_basis();
  };

  BENCHMARK("get_molar_internal_energy") {
    assemblage.reset_cache();
    return assemblage.get_molar_internal_energy();
  };

  BENCHMARK("get_molar_gibbs") {
    assemblage.reset_cache();
    return assemblage.get_molar_gibbs();
  };

  BENCHMARK("get_molar_helmholtz") {
    assemblage.reset_cache();
    return assemblage.get_molar_helmholtz();
  };

  BENCHMARK("get_molar_mass") {
    assemblage.reset_cache();
    return assemblage.get_molar_mass();
  };

  BENCHMARK("get_molar_volume") {
    assemblage.reset_cache();
    return assemblage.get_molar_volume();
  };

  BENCHMARK("get_density") {
    assemblage.reset_cache();
    return assemblage.get_density();
  };

  BENCHMARK("get_molar_entropy") {
    assemblage.reset_cache();
    return assemblage.get_molar_entropy();
  };

  BENCHMARK("get_molar_enthalpy") {
    assemblage.reset_cache();
    return assemblage.get_molar_enthalpy();
  };

  BENCHMARK("get_isothermal_bulk_modulus_reuss") {
    assemblage.reset_cache();
    return assemblage.get_isothermal_bulk_modulus_reuss();
  };

  BENCHMARK("get_isentropic_bulk_modulus_reuss") {
    assemblage.reset_cache();
    return assemblage.get_isentropic_bulk_modulus_reuss();
  };

  BENCHMARK("get_isothermal_compressibility_reuss") {
    assemblage.reset_cache();
    return assemblage.get_isothermal_compressibility_reuss();
  };

  BENCHMARK("get_isentropic_compressibility_reuss") {
    assemblage.reset_cache();
    return assemblage.get_isentropic_compressibility_reuss();
  };

  BENCHMARK("get_shear_modulus") {
    assemblage.reset_cache();
    return assemblage.get_shear_modulus();
  };

  BENCHMARK("get_p_wave_velocity") {
    assemblage.reset_cache();
    return assemblage.get_p_wave_velocity();
  };

  BENCHMARK("get_bulk_sound_velocity") {
    assemblage.reset_cache();
    return assemblage.get_bulk_sound_velocity();
  };

  BENCHMARK("get_shear_wave_velocity") {
    assemblage.reset_cache();
    return assemblage.get_shear_wave_velocity();
  };

  BENCHMARK("get_grueneisen_parameter") {
    assemblage.reset_cache();
    return assemblage.get_grueneisen_parameter();
  };

  BENCHMARK("get_thermal_expansivity") {
    assemblage.reset_cache();
    return assemblage.get_thermal_expansivity();
  };

  BENCHMARK("get_molar_heat_capacity_v") {
    assemblage.reset_cache();
    return assemblage.get_molar_heat_capacity_v();
  };

  BENCHMARK("get_molar_heat_capacity_p") {
    assemblage.reset_cache();
    return assemblage.get_molar_heat_capacity_p();
  };

}
