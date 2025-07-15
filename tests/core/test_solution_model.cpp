/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "tolerances.hpp"
#include "burnman/core/solution_model.hpp"
#include "burnman/core/mineral.hpp"

using namespace Catch::Matchers;

struct OlivineFixture {
  Mineral forsterite;
  Mineral fayalite;
  PairedEndmemberList olivine_endmembers;

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
    // Params not needed here... use for Solution test

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
    fayalite.params.equation_of_state = EOSType::SLB3; // TODO!

    // Paired list for solution model setup
    olivine_endmembers = {
      {forsterite, "[Mg]2SiO4"},
      {fayalite, "[Fe]2SiO4"}
    };
  }
};

struct MultiSiteFixture {
  Mineral min;
  PairedEndmemberList em;
  MultiSiteFixture() {
    em = {
      {min, "[Mg]3[Al]2Si3O12"},
      {min, "[Fe]3[Al]2Si3O12"},
      {min, "[Mg]3[Mg1/2Si1/2]2Si3O12"}
    };
  }
};

TEST_CASE_METHOD(OlivineFixture, "IdealSolution", "[core][solution_model]") {
  // Create IdealSolution model
  IdealSolution ol_ss(olivine_endmembers);
  SECTION("Derived parameters") {
    // Check endmembers (Mineral objects) can be accessed in order
    REQUIRE(ol_ss.endmembers[0].get_name() == "Forsterite");
    REQUIRE(ol_ss.endmembers[1].get_name() == "Fayalite");
    // ints
    CHECK(ol_ss.n_endmembers == 2);
    CHECK(ol_ss.n_occupancies == 2);
    CHECK(ol_ss.n_sites == 1);
    // strings
    CHECK(ol_ss.empty_formula == "[]2SiO4");
    CHECK(ol_ss.general_formula == "[Mg,Fe]2SiO4");
    CHECK(ol_ss.formulas == std::vector<std::string>{"[Mg]2SiO4", "[Fe]2SiO4"});
    CHECK(ol_ss.site_names == std::vector<std::string>{"Mg_A", "Fe_A"});
    CHECK(ol_ss.sites == std::vector<std::vector<std::string>>{{"Mg", "Fe"}});
    CHECK(ol_ss.solution_formulae ==
      std::vector<std::map<std::string, double>>{ {{"Mg", 2.0}, {"Fe", 2.0}} });
    // Eigen - can do CHECK(m1.isApprox(m2, rel_tolerance));
    CHECK(ol_ss.site_multiplicities.isApprox(
      (Eigen::ArrayXXd(2, 2) << 2.0, 2.0,
                                2.0, 2.0).finished(), tol_rel));
    CHECK(ol_ss.endmember_occupancies.isApprox(
      (Eigen::ArrayXXd(2, 2) << 1, 0,
                                0, 1).finished(), tol_rel));
    CHECK(ol_ss.endmember_n_occupancies.isApprox(
      (Eigen::ArrayXXd(2, 2) << 2, 0,
                                0, 2).finished(), tol_rel));
  }
  SECTION("Default returns") {
    CHECK(ol_ss.compute_Cp_excess() == 0);
    CHECK(ol_ss.compute_alphaV_excess() == 0);
    CHECK(ol_ss.compute_VoverKT_excess() == 0);
    // Dummy args
    double a = 0;
    double b = 0;
    Eigen::ArrayXd c = Eigen::ArrayXd::Zero(2);
    CHECK(ol_ss.compute_excess_partial_volumes(a, b, c).isZero());
    CHECK(ol_ss.compute_volume_hessian(a, b, c).isZero());
    CHECK(ol_ss.compute_activity_coefficients(a, b, c).isOnes());
  }
}

TEST_CASE_METHOD(OlivineFixture, "Single ss", "[core][solution_model]") {
  PairedEndmemberList fo_repeat = {
    {forsterite, "[Mg]2SiO4"},
    {forsterite, "[Mg]2SiO4"}
  };
  IdealSolution fo_ss(fo_repeat);
  // Ensure different objects
  REQUIRE(&fo_ss.endmembers[0] != &fo_ss.endmembers[1]);
  // ints
  CHECK(fo_ss.n_endmembers == 2);
  CHECK(fo_ss.n_occupancies == 1);
  CHECK(fo_ss.n_sites == 1);
  // strings
  CHECK(fo_ss.empty_formula == "[]2SiO4");
  CHECK(fo_ss.general_formula == "[Mg]2SiO4");
  CHECK(fo_ss.formulas == std::vector<std::string>{"[Mg]2SiO4", "[Mg]2SiO4"});
  CHECK(fo_ss.site_names == std::vector<std::string>{"Mg_A"});
  CHECK(fo_ss.sites == std::vector<std::vector<std::string>>{{"Mg"}});
  CHECK(fo_ss.solution_formulae ==
    std::vector<std::map<std::string, double>>{ {{"Mg", 2.0}, {"Mg", 2.0}} });
  CHECK(fo_ss.site_multiplicities.isApprox(
    (Eigen::ArrayXXd(2, 1) << 2.0,
                              2.0).finished(), tol_rel));
  CHECK(fo_ss.endmember_occupancies.isApprox(
    (Eigen::ArrayXXd(2, 1) << 1,
                              1).finished(), tol_rel));
  CHECK(fo_ss.endmember_n_occupancies.isApprox(
    (Eigen::ArrayXXd(2, 1) << 2,
                              2).finished(), tol_rel));
}

TEST_CASE_METHOD(MultiSiteFixture, "Check multi-site solutions", "[core][solution_model]") {
  IdealSolution sol(em);
  // ints
  CHECK(sol.n_endmembers == 3);
  CHECK(sol.n_occupancies == 5);
  CHECK(sol.n_sites == 2);
  // strings
  CHECK(sol.empty_formula == "[]3[]2Si3O12");
  CHECK(sol.general_formula == "[Mg,Fe]3[Al,Mg,Si]2Si3O12");
  CHECK(sol.formulas == std::vector<std::string>{"[Mg]3[Al]2Si3O12", "[Fe]3[Al]2Si3O12", "[Mg]3[Mg1/2Si1/2]2Si3O12"});
  CHECK(sol.site_names == std::vector<std::string>{"Mg_A", "Fe_A", "Al_B", "Mg_B", "Si_B"});
  CHECK(sol.sites == std::vector<std::vector<std::string>>{{"Mg", "Fe"}, {"Al", "Mg", "Si"}});
  CHECK(sol.solution_formulae ==
    std::vector<std::map<std::string, double>>{ {{"Mg", 3.0}, {"Al", 2.0}},
                                                {{"Fe", 3.0}, {"Al", 2.0}},
                                                {{"Mg", 4.0}, {"Si", 1.0}} });
  CHECK(sol.site_multiplicities.isApprox(
    (Eigen::ArrayXXd(3, 5) << 3, 3, 2, 2, 2,
                              3, 3, 2, 2, 2,
                              3, 3, 2, 2, 2).finished(), tol_rel));
  CHECK(sol.endmember_occupancies.isApprox(
    (Eigen::ArrayXXd(3, 2) << 1, 0, 1, 0, 0,
                              0, 1, 1, 0, 0,
                              1, 0, 0, 0.5, 0.5).finished(), tol_rel));
  CHECK(sol.endmember_n_occupancies.isApprox(
    (Eigen::ArrayXXd(3, 2) << 3, 0, 2, 0, 0,
                              0, 3, 2, 0, 0,
                              3, 0, 0, 1, 1).finished(), tol_rel));
}

TEST_CASE_METHOD(MultiSiteFixture, "Regular Solutions", "[core][solution_model]") {
  std::vector<double> alphas = {1.0, 2.0, 2.0};
  std::vector<double> a_ones = {1.0, 1.0, 1.0};
  std::vector<std::vector<double>> interactions = {
    {10.0e3, 5.0e3},
    {-10.0e3}
  };
  Eigen::MatrixXd ref_W(3,3);
  ref_W <<  0,     10000,   5000,
            0,         0, -10000,
            0,         0,      0;
  Eigen::MatrixXd ref_Wa(3,3);
  ref_Wa << 0, 10000*(2/3), 5000*(2/3),
            0,           0,       5000,
            0,           0,          0;

  SECTION("Symmetric Solution - alphas") {
    SymmetricRegularSolution sym_sol(em, interactions);
    AsymmetricRegularSolution asym_sol(em, a_ones, interactions);
    REQUIRE((sym_sol.alphas == asym_sol.alphas).all());
    REQUIRE(sym_sol.alphas.isOnes());
  }
  SECTION("Symmetric Solution - interactions") {
    SymmetricRegularSolution e_sol(em, interactions);
    SymmetricRegularSolution evs_sol(em, interactions, interactions, interactions);
    REQUIRE(e_sol.W_e.rows() == e_sol.W_e.cols());
    REQUIRE(e_sol.W_v.rows() == e_sol.W_v.cols());
    REQUIRE(e_sol.W_s.rows() == e_sol.W_s.cols());
    CHECK(e_sol.W_v.isZero());
    CHECK(e_sol.W_s.isZero());
    CHECK(e_sol.W_e.isApprox(ref_W, tol_rel));
    CHECK(evs_sol.W_e.isApprox(ref_W, tol_rel));
    CHECK(evs_sol.W_e.isApprox(evs_sol.W_v, tol_rel));
    CHECK(evs_sol.W_e.isApprox(evs_sol.W_s, tol_rel));
  }
  SECTION("Asymmetric Solution - interactions") {
    AsymmetricRegularSolution sol(em, alphas, interactions);
    REQUIRE((sol.alphas == (Eigen::ArrayXd(3) << 1, 2, 2).finished()).all());
    CHECK(sol.W_e.isApprox(ref_Wa, tol_rel));
  }
}
