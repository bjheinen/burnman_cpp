/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <memory>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "tolerances.hpp"
#include "solution_fixtures.hpp"
#include "burnman/core/assemblage.hpp"
#include "burnman/utils/types.hpp"

using namespace Catch::Matchers;

// Mineral, Mineral --> for making
// TODO:
// Do Solution, Solution next for checking interface a bit more
// Do Mineral, Solution just to check no errors, brief
// Do Mineral, Solution, Solution just to check no errors too, brief
// Do Mineral, Assemblage just for error check too
// Do Mineral, Solution, Solution again for the ref tests

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
  REQUIRE_NOTHROW(rock.set_averaging_scheme(AveragingType::VRH));

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
  FormulaMap expected_formula = {
    {"Si", 1.0},
    {"O", 2.6},
    {"Ca", 0.6}
  };
  FormulaMap computed_formula = rock.get_formula();
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
  std::vector<FormulaMap> expected_formulae = {
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
  REQUIRE_NOTHROW(rock.set_method(EOSType::Auto));

}

// TODO:
// reset() (check over)

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
