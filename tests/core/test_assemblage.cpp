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
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "burnman/core/assemblage.hpp"
#include <memory>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "tolerances.hpp"
#include "solution_fixtures.hpp"

using namespace Catch::Matchers;
using namespace burnman;

struct MultiMineralFixture {
  CaPerovskiteFixture capv_fix;
  StishoviteFixture stish_fix;
};

TEST_CASE_METHOD(MultiMineralFixture, "Test assemblage creation", "[core][assemblage]") {
  // Make assemblage
  Assemblage rock;

  // Test name(s)
  REQUIRE_NOTHROW(rock.get_name());
  REQUIRE(rock.get_name() != "Test assemblage");
  REQUIRE_NOTHROW(rock.set_name("Test assemblage"));
  REQUIRE(rock.get_name() == "Test assemblage");

  // Test set averaging scheme
  REQUIRE_NOTHROW(rock.set_averaging_scheme(types::AveragingType::VRH));

  // Can't set fractions without phases
  REQUIRE_THROWS(rock.set_fractions((Eigen::ArrayXd(2) << 0.6, 0.4).finished()));

  // Can add phases
  auto capv_ptr = std::make_shared<Mineral>(capv_fix.ca_perovskite);
  auto stish_ptr = std::make_shared<Mineral>(stish_fix.stishovite);
  REQUIRE_NOTHROW(rock.add_phases({capv_ptr, stish_ptr}));
  // Check get phase
  REQUIRE(rock.get_phase(0) == capv_ptr);
  REQUIRE(rock.get_phase(1) == stish_ptr);
  // And set molar fractions
  REQUIRE_NOTHROW(rock.set_fractions((Eigen::ArrayXd(2) << 0.6, 0.4).finished()));
  REQUIRE(rock.get_n_endmembers() == 2);

  // Check computed formula
  types::FormulaMap expected_formula = {
    {"Si", 1.0},
    {"O", 2.6},
    {"Ca", 0.6}
  };
  types::FormulaMap computed_formula = rock.get_formula();
  REQUIRE(computed_formula.size() == computed_formula.size());
  for (const auto& [key, val] : computed_formula) {
    REQUIRE(expected_formula.find(key) != expected_formula.end());
    REQUIRE_THAT(val, WithinAbs(expected_formula.at(key), tol_abs) || WithinRel(expected_formula.at(key), tol_rel));
  }

  // Check state can be set
  REQUIRE_NOTHROW(rock.set_state(50.e9, 2000.0));
  // Confirm P, T set
  REQUIRE(rock.get_pressure() == 50.e9);
  REQUIRE(rock.get_temperature() == 2000.0);
  REQUIRE(rock.get_phase(0)->get_pressure() == 50.e9);
  REQUIRE(rock.get_phase(0)->get_temperature() == 2000.0);
  REQUIRE(rock.get_phase(1)->get_pressure() == 50.e9);
  REQUIRE(rock.get_phase(1)->get_temperature() == 2000.0);

  // Check get_names (and internal set_names from params) is working
  std::vector<std::string> expected_names = {"Ca-perovskite", "Stishovite"};
  REQUIRE(rock.get_endmember_names() == expected_names);

  // Check endmember_formulae
  std::vector<types::FormulaMap> expected_formulae = {
    {
      {"Ca", 1.0},
      {"Si", 1.0},
      {"O", 3.0}
    },
    {
      {"Si", 1.0},
      {"O", 2.0}
    }
  };
  REQUIRE(rock.get_endmember_formulae() == expected_formulae);

  // Test element list
  std::vector<std::string> expected_elements = {"Ca", "Si", "O"};
  REQUIRE(rock.get_elements() == expected_elements);

  // Test setting method for all phases
  REQUIRE_NOTHROW(rock.set_method(types::EOSType::Auto));

}

TEST_CASE_METHOD(BdgFperAssemblageFixture, "Test multi-solution assemblage", "[core][assemblage]") {
  REQUIRE(assemblage.get_name() == "Bdg + Fper assemblage");
  REQUIRE(assemblage.get_n_endmembers() == 5);
  REQUIRE(assemblage.get_n_elements() == 5);
  REQUIRE(assemblage.get_endmembers_per_phase() == std::vector<int>{3, 2});
  std::vector<std::string> expected_names = {
    "MgSiO3 perovskite in Bridgmanite",
    "FeSiO3 perovskite in Bridgmanite",
    "AlAlO3 perovskite in Bridgmanite",
    "Periclase in Ferro-periclase",
    "Wuestite in Ferro-periclase"
  };
  REQUIRE(assemblage.get_endmember_names() == expected_names);
  REQUIRE_NOTHROW(assemblage.set_state(50.e9, 2000.0));
  // Confirm P, T set
  REQUIRE(assemblage.get_pressure() == 50.e9);
  REQUIRE(assemblage.get_temperature() == 2000.0);
  REQUIRE(assemblage.get_phase(0)->get_pressure() == 50.e9);
  REQUIRE(assemblage.get_phase(0)->get_temperature() == 2000.0);
  REQUIRE(assemblage.get_phase(1)->get_pressure() == 50.e9);
  REQUIRE(assemblage.get_phase(1)->get_temperature() == 2000.0);
}

TEST_CASE_METHOD(PyroliteAssemblageFixture, "Test mixed phase type assemblage", "[core][assemblage]") {
  REQUIRE(assemblage.get_name() == "Bdg + Fper + CaPv assemblage");
  REQUIRE(assemblage.get_n_endmembers() == 6);
  REQUIRE(assemblage.get_endmembers_per_phase() == std::vector<int>{3, 2, 1});
  REQUIRE(assemblage.get_n_elements() == 6);
  REQUIRE_NOTHROW(assemblage.set_state(50.e9, 2000.0));
  std::vector<types::FormulaMap> expected_formulae = {
    {
      {"Mg", 1.0},
      {"Si", 1.0},
      {"O", 3.0}
    },
    {
      {"Fe", 1.0},
      {"Si", 1.0},
      {"O", 3.0}
    },
    {
      {"Al", 2.0},
      {"O", 3.0}
    },
    {
      {"Mg", 1.0},
      {"O", 1.0}
    },
    {
      {"Fe", 1.0},
      {"O", 1.0}
    },
    {
      {"Ca", 1.0},
      {"Si", 1.0},
      {"O", 3.0}
    }
  };
  REQUIRE(assemblage.get_endmember_formulae() == expected_formulae);
  types::FormulaMap expected_formula = {
    {"Mg", 0.796},
    {"Fe", 0.069},
    {"Ca", 0.1},
    {"Al", 0.07},
    {"Si", 0.765},
    {"O", 2.6}
  };
  types::FormulaMap computed_formula = assemblage.get_formula();
  REQUIRE(computed_formula.size() == computed_formula.size());
  for (const auto& [key, val] : computed_formula) {
    REQUIRE(expected_formula.find(key) != expected_formula.end());
    REQUIRE_THAT(val, WithinAbs(expected_formula.at(key), tol_abs) || WithinRel(expected_formula.at(key), tol_rel));
  }
}

TEST_CASE_METHOD(NestedAssemblageFixture, "Test nested assemblage", "[core][assemblage]") {
  REQUIRE(assemblage.get_name() == "Nested assemblage");
  REQUIRE(assemblage.get_n_endmembers() == 7);
  REQUIRE_NOTHROW(assemblage.set_state(50.e9, 2000.0));
  // Check access to nested assemblage phases
  REQUIRE(assemblage.get_pressure() == 50.e9);
  REQUIRE(assemblage.get_temperature() == 2000.0);
  REQUIRE(assemblage.get_phase(0)->get_pressure() == 50.e9);
  REQUIRE(assemblage.get_phase(1)->get_pressure() == 50.e9);
  // Check get_phase works with cast
  REQUIRE(assemblage.get_phase<Assemblage>(0)->get_phase(0)->get_pressure() == 50.e9);
}

TEST_CASE_METHOD(PyroliteAssemblageFixture, "Assemblage reference values", "[core][assemblage]") {
  assemblage.set_state(50.e9, 2000.0);
  assemblage.set_method(types::EOSType::Auto);
  // Assemblage functions
  int ref_n_endmembers = 6;
  int ref_n_elements = 6;
  int ref_n_reactions = 1;
  std::vector<std::string> ref_elements = {"Ca", "Mg", "Fe", "Al", "Si", "O"};
  std::vector<Eigen::Index> ref_independent_element_indices = {0,1,2,3,4};
  std::vector<Eigen::Index> ref_dependent_element_indices = {5};
  Eigen::MatrixXd ref_stoichiometric_matrix(6,6);
  ref_stoichiometric_matrix <<
    0, 1, 0, 0, 1, 3,
    0, 0, 1, 0, 1, 3,
    0, 0, 0, 2, 0, 3,
    0, 1, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 1,
    1, 0, 0, 0, 1, 3;
  Eigen::MatrixXd ref_compositional_basis(5,6);
  ref_compositional_basis <<
    1., 0., 0., 0., 0., 0.,
    0., 1., 0., 0., 0., 0.,
    0., 0., 1., 0., 0., 0.,
    0., 0., 0., 1., 0., 0.,
    0., 0., 0., 0., 0., 1.;
  Eigen::MatrixXd ref_compositional_null_basis(1,6);
  ref_compositional_null_basis <<
    -1.0, -1.0, -1.0, -1.5, -2.0, 1.0;
  Eigen::MatrixXd ref_reaction_basis(1,6);
  ref_reaction_basis <<
    1.0, -1.0, 0, -1.0, 1.0, 0;
  // Material functions
  double ref_molar_internal_energy = -974191.8211164112;
  double ref_molar_gibbs = 999561205.8828137;
  double ref_molar_helmholtz = -1419441.7062346758;
  double ref_molar_mass = 0.09218043750000002;
  double ref_molar_volume = 0.020019612951780966;
  double ref_density = 4.604506476824746;
  double ref_molar_entropy = 222.62494255913217;
  double ref_molar_enthalpy = 1000006455.7679318;
  double ref_isothermal_bulk_modulus_reuss = 1824185011747062.2;
  double ref_isentropic_bulk_modulus_reuss = 1898805807921876.2;
  double ref_isothermal_compressibility_reuss = 5.481900100923853e-16;
  double ref_isentropic_compressibility_reuss = 5.266467986499563e-16;
  double ref_shear_modulus = 173600057512.71997;
  double ref_p_wave_velocity = 20308376.231700655;
  double ref_bulk_sound_velocity = 20307138.536943752;
  double ref_shear_wave_velocity = 194170.56228630507;
  double ref_grueneisen_parameter = 6648.456796536628;
  double ref_thermal_expansivity = 1.9707992960705696e-08;
  double ref_molar_heat_capacity_v = 108.2545371299868;
  double ref_molar_heat_capacity_p = 114.40558099311642;

  CHECK(assemblage.get_n_endmembers() == ref_n_endmembers);
  CHECK(assemblage.get_n_elements() == ref_n_elements);
  CHECK(assemblage.get_n_reactions() == ref_n_reactions);
  CHECK(assemblage.get_elements() == ref_elements);
  CHECK(assemblage.get_independent_element_indices() == ref_independent_element_indices);
  CHECK(assemblage.get_dependent_element_indices() == ref_dependent_element_indices);
  CHECK(assemblage.get_stoichiometric_matrix().isApprox(ref_stoichiometric_matrix, tol_rel));
  CHECK(assemblage.get_compositional_basis().isApprox(ref_compositional_basis, tol_rel));
  CHECK(assemblage.get_compositional_null_basis().isApprox(ref_compositional_null_basis, tol_rel));
  CHECK(assemblage.get_reaction_basis().isApprox(ref_reaction_basis));
  CHECK_THAT(assemblage.get_molar_internal_energy(),
    WithinRel(ref_molar_internal_energy, tol_rel) || WithinAbs(ref_molar_internal_energy, tol_abs));
  CHECK_THAT(assemblage.get_molar_gibbs(),
    WithinRel(ref_molar_gibbs, tol_rel) || WithinAbs(ref_molar_gibbs, tol_abs));
  CHECK_THAT(assemblage.get_molar_helmholtz(),
    WithinRel(ref_molar_helmholtz, tol_rel) || WithinAbs(ref_molar_helmholtz, tol_abs));
  CHECK_THAT(assemblage.get_molar_mass(),
    WithinRel(ref_molar_mass, tol_rel) || WithinAbs(ref_molar_mass, tol_abs));
  CHECK_THAT(assemblage.get_molar_volume(),
    WithinRel(ref_molar_volume, tol_rel) || WithinAbs(ref_molar_volume, tol_abs));
  CHECK_THAT(assemblage.get_density(),
    WithinRel(ref_density, tol_rel) || WithinAbs(ref_density, tol_abs));
  CHECK_THAT(assemblage.get_molar_entropy(),
    WithinRel(ref_molar_entropy, tol_rel) || WithinAbs(ref_molar_entropy, tol_abs));
  CHECK_THAT(assemblage.get_molar_enthalpy(),
    WithinRel(ref_molar_enthalpy, tol_rel) || WithinAbs(ref_molar_enthalpy, tol_abs));
  CHECK_THAT(assemblage.get_isothermal_bulk_modulus_reuss(),
    WithinRel(ref_isothermal_bulk_modulus_reuss, tol_rel) || WithinAbs(ref_isothermal_bulk_modulus_reuss, tol_abs));
  CHECK_THAT(assemblage.get_isentropic_bulk_modulus_reuss(),
    WithinRel(ref_isentropic_bulk_modulus_reuss, tol_rel) || WithinAbs(ref_isentropic_bulk_modulus_reuss, tol_abs));
  CHECK_THAT(assemblage.get_isothermal_compressibility_reuss(),
    WithinRel(ref_isothermal_compressibility_reuss, tol_rel) || WithinAbs(ref_isothermal_compressibility_reuss, tol_abs));
  CHECK_THAT(assemblage.get_isentropic_compressibility_reuss(),
    WithinRel(ref_isentropic_compressibility_reuss, tol_rel) || WithinAbs(ref_isentropic_compressibility_reuss, tol_abs));
  CHECK_THAT(assemblage.get_shear_modulus(),
    WithinRel(ref_shear_modulus, tol_rel) || WithinAbs(ref_shear_modulus, tol_abs));
  CHECK_THAT(assemblage.get_p_wave_velocity(),
    WithinRel(ref_p_wave_velocity, tol_rel) || WithinAbs(ref_p_wave_velocity, tol_abs));
  CHECK_THAT(assemblage.get_bulk_sound_velocity(),
    WithinRel(ref_bulk_sound_velocity, tol_rel) || WithinAbs(ref_bulk_sound_velocity, tol_abs));
  CHECK_THAT(assemblage.get_shear_wave_velocity(),
    WithinRel(ref_shear_wave_velocity, tol_rel) || WithinAbs(ref_shear_wave_velocity, tol_abs));
  CHECK_THAT(assemblage.get_grueneisen_parameter(),
    WithinRel(ref_grueneisen_parameter, tol_rel) || WithinAbs(ref_grueneisen_parameter, tol_abs));
  CHECK_THAT(assemblage.get_thermal_expansivity(),
    WithinRel(ref_thermal_expansivity, tol_rel) || WithinAbs(ref_thermal_expansivity, tol_abs));
  CHECK_THAT(assemblage.get_molar_heat_capacity_v(),
    WithinRel(ref_molar_heat_capacity_v, tol_rel) || WithinAbs(ref_molar_heat_capacity_v, tol_abs));
  CHECK_THAT(assemblage.get_molar_heat_capacity_p(),
    WithinRel(ref_molar_heat_capacity_p, tol_rel) || WithinAbs(ref_molar_heat_capacity_p, tol_abs));
}


// TODO:
// reset_cache() (check over)

// add_phase (value)
// add_phases (value - variadic template)

// add_phase (pointer)
// add_phases (vector of ptrs)
// add_phases(init list of ptrs)

// set_fractions
// set_fractions using volume fractions
// set_averaging_scheme
// set_averaging_scheme(custom scheme with unique ptr)

// set_state
// set_method

// get_phase --> with the cast
// get_phase (as Material)

// get_endmembers_per_phase --> check set correctly

// // From Material
// set_name
// get_name
// get_formula
// get_pressure
// get_temperature

// // From composite material
// get_n_endmembers
// get_n_elements
// get_n_reactions
// get_elements
// get_endmember_names
// get_endmember_formulae
// get_independent_element_indices
// get_dependent_element_indices
// get_stoichiometric_matrix
// get_compositional_basis
// get_compositional_null_basis
// get_reaction_basis

// // From assemblage
// get_molar_internal_energy
// get_molar_gibbs
// get_molar_helmholtz
// get_molar_mass
// get_molar_volume
// get_density
// get_molar_entropy
// get_molar_enthalpy
// get_isothermal_bulk_modulus_reuss
// get_isentropic_bulk_modulus_reuss
// get_isothermal_compressibility_reuss
// get_isentropic_compressibility_reuss
// get_shear_modulus
// get_p_wave_velocity
// get_shear_wave_velocity
// get_bulk_sound_velocity
// get_grueneisen_parameter
// get_thermal_expansivity
// get_molar_heat_capacity_v
// get_molar_heat_capacity_p

// get_isentropic_therma_gradient
