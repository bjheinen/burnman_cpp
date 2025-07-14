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
// #include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "tolerances.hpp"
#include "burnman/core/solution_model.hpp"
#include "burnman/core/mineral.hpp"

using namespace Catch::Matchers;

// TEST_CASE_METHOD(AveragingSchemeFixture, "Averaging::reuss_fn", "[core][averaging_schemes]") {
//   // Quick check of Reuss average
//   REQUIRE(Averaging::reuss_fn(v_a, X_a) == 1.6);
//   // Check reuss_fn normalises if needed
//   REQUIRE(Averaging::reuss_fn(v_a*34.235, X_a) == 1.6);
//   CHECK_THAT(Averaging::reuss_fn(v_b, X_b),
//     WithinRel(ref_reuss, tol_rel) ||
//     WithinAbs(ref_reuss, tol_abs));
// }

struct OlivineFixture {
  Mineral forsterite;
  Mineral fayalite;

  OlivineFixture() {
    // Set fo parameters
    forsterite.params.name = "Forsterite";
    forsterite.params.formula = FormulaMap{
      {"Mg", 2.0},
      {"Si", 1.0},
      {"O", 4.0}
    };
    forsterite.params.napfu = 7;
    forsterite.params.molar_mass = 0.1406931;
    forsterite.params.V_0 = 4.366e-05;
    forsterite.params.K_0 = 1.285e11;
    forsterite.params.Kprime_0 = 3.84;
    forsterite.params.Kdprime_0 = -3e-11;
    forsterite.params.H_0 = -2172590.0;
    forsterite.params.S_0 = 95.1;
    forsterite.params.Cp = CpParams{233.3, 0.001494, -603800.0, -1869.7};
    forsterite.params.a_0 = 2.85e-05;
    forsterite.params.equation_of_state = EOSType::SLB3; // TODO! Change!
    // set_method on mineral?.....

    // Fayalite
    fayalite.params.name = "Fayalite";
    fayalite.params.formula = FormulaMap{
      {"Fe", 2.0},
      {"Si", 1.0},
      {"O", 4.0}
    };
    fayalite.params.napfu = 7;
    fayalite.params.molar_mass = 0.2037731;
    fayalite.params.V_0 = 4.631e-05;
    fayalite.params.K_0 = 1.256e11;
    fayalite.params.Kprime_0 = 4.68;
    fayalite.params.Kdprime_0 = -3.7e-11;
    fayalite.params.H_0 = -1477720.0;
    fayalite.params.S_0 = 151.0;
    fayalite.params.Cp = CpParams{201.1, 0.01733, -1960600.0, -900.9};
    fayalite.params.a_0 = 2.82e-05;
    fayalite.params.equation_of_state = EOSType::SLB3; // TODO! Change!

  }
};

// Tests
// Fixture --> dummy solutions, up to 3 components
// IdealSolution to test main functionality
//   process_solution_chemistry
//   constructor for IdealSolution 
//   --> check can get endmembers, formulas, and solution_chemistry properties
//   check Cp_excess, alpha_V_excess, VoverKT_excess return 0
//   check endmember configurational entropies
// AsymmetricReg... 
//   Constructor
//   --> as IdealSolution + alphas available,
//   --> with W_e, W_v, W_s (square + correct)
// SymmetricRe...
//   Constructor --> alphas are 1, all else still available
