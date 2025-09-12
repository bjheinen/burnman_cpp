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
#include "burnman/core/solution.hpp"
#include <memory>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/core/solution_model.hpp"
#include "tolerances.hpp"
#include "solution_fixtures.hpp"

using namespace Catch::Matchers;
using namespace burnman;

TEST_CASE_METHOD(BridgmaniteFixture, "Test interface", "[core][solution]") {
  // Make solution
  Solution bdg;
  // Test name(s)
  REQUIRE_NOTHROW(bdg.get_name());
  REQUIRE(bdg.get_name() != "Test solution");
  REQUIRE_NOTHROW(bdg.set_name("Test solution"));
  REQUIRE(bdg.get_name() == "Test solution");

  // Check composition can't be set if no solution model
  REQUIRE_THROWS(bdg.set_composition(molar_fractions));
  // Set solution_model
  REQUIRE_NOTHROW(bdg.set_solution_model(bdg_solution_model));

  // Check get_names (and internal set_names from params) is working
  std::vector<std::string> expected_names =
    {"MgSiO3 perovskite", "FeSiO3 perovskite", "AlAlO3 perovskite"};
  REQUIRE(bdg.get_endmember_names() == expected_names);
  // Check endmember_formulae
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
    }
  };
  REQUIRE(bdg.get_endmember_formulae() == expected_formulae);
  // Test element list
  std::vector<std::string> expected_elements = {"Mg", "Fe", "Al", "Si", "O"};
  REQUIRE(bdg.get_elements() == expected_elements);

  // Test setting method for all endmembers
  REQUIRE_NOTHROW(bdg.set_method(types::EOSType::Auto));

  // Test composition setting
  // Molar fractions of endmembers (Mg, Fe, Al)
  // Check bad molar fractions array
  REQUIRE_THROWS(bdg.set_composition((Eigen::ArrayXd(2) << 1, 1).finished()));
  REQUIRE_THROWS(bdg.set_composition((Eigen::ArrayXd(3) << 0.1, 0.1, 0.1).finished()));
  // Check set molar_fractions
  REQUIRE_NOTHROW(bdg.set_composition(molar_fractions));
  // Test formula
  types::FormulaMap expected_solution_formula = {
    {"Mg", 0.88},
    {"Fe", 0.07},
    {"Al", 0.1},
    {"Si", 0.95},
    {"O", 3.0}
  };
  REQUIRE(bdg.get_formula() == expected_solution_formula);
  // Test state
  REQUIRE_NOTHROW(bdg.set_state(P, T));
  REQUIRE_THAT(bdg.get_pressure(),
    WithinRel(P, tol_rel) || WithinAbs(P, tol_abs));
  REQUIRE_THAT(bdg.get_temperature(),
    WithinRel(T, tol_rel) || WithinAbs(T, tol_abs));
  // Test reset runs
  REQUIRE_NOTHROW(bdg.reset());
}

TEST_CASE_METHOD(BridgmaniteFixture, "Test reference values", "[core][solution]") {
  // Reference data
  // Solution functions
  double ref_excess_gibbs = -10757.973872013106;
  double ref_excess_volume = 0.0;
  double ref_excess_entropy = 5.378986936006553;
  double ref_excess_enthalpy = 0.0;
  Eigen::ArrayXd ref_activities(3);
  ref_activities << 0.836, 0.0665, 0.0025;
  Eigen::ArrayXd ref_activity_coefficients(3);
  ref_activity_coefficients << 1., 1., 1.;
  Eigen::ArrayXd ref_excess_partial_gibbs(3);
  ref_excess_partial_gibbs << -2978.683935037304, -45073.58869554721, -99631.61600983949;
  Eigen::ArrayXd ref_excess_partial_volumes(3);
  ref_excess_partial_volumes << 0., 0., 0.;
  Eigen::ArrayXd ref_excess_partial_entropies(3);
  ref_excess_partial_entropies << 1.48934196751865, 22.5367943477736, 49.81580800491975;
  Eigen::ArrayXd ref_partial_gibbs(3);
  ref_partial_gibbs << -740818.9553017676, -429520.9819551775, -988122.8420793302;
  Eigen::ArrayXd ref_partial_entropies(3);
  ref_partial_entropies << 250.55024892868906, 277.1966162401331, 301.75523166155415;
  Eigen::MatrixXd ref_gibbs_hessian(3,3);
  ref_gibbs_hessian <<
      3142.787305426345, -15753.718644921928, -33257.85047261296,
    -15753.718644921928,   221802.3561594563, -33257.85047261296,
    -33257.85047261296,  -33257.85047261296,  631899.1589796463;
  Eigen::MatrixXd ref_entropy_hessian(3,3);
  ref_entropy_hessian <<
    -1.57139365271317,    7.87685932246096,  16.62892523630648,
    7.87685932246096, -110.90117807972815,  16.62892523630648,
    16.62892523630648,   16.62892523630648, -315.9495794898231;
  Eigen::MatrixXd ref_volume_hessian = Eigen::MatrixXd::Zero(3, 3); // DO ISZERO INSTEAD
  // Composite Material functions
  int ref_n_endmembers = 3;
  int ref_n_elements = 5;
  int ref_n_reactions = 0;
  std::vector<Eigen::Index> ref_independent_element_indices = {0, 1, 2};
  std::vector<Eigen::Index> ref_dependent_element_indices = {3, 4};
  Eigen::MatrixXd ref_stoichiometric_matrix(3,5);
  ref_stoichiometric_matrix <<
    1, 0, 0, 1, 3,
    0, 1, 0, 1, 3,
    0, 0, 2, 0, 3;
  Eigen::MatrixXd ref_compositional_basis(3,3);
  ref_compositional_basis <<
    1., 0., 0.,
    0., 1., 0.,
    0., 0., 1.;
  Eigen::MatrixXd ref_compositional_null_basis(2,5);
  ref_compositional_null_basis <<
    -1, -1,    0, 1, 0,
    -3, -3, -1.5, 0, 1;
  Eigen::MatrixXd ref_reaction_basis(0,3);
  //ref_reaction_basis << []; // ? Check what to do here....
  // Material functions
  double ref_molar_internal_energy = -1118021.6956653593;
  double ref_molar_gibbs = -731393.2915063845;
  double ref_molar_helmholtz = -1627973.183219626;
  double ref_molar_mass = 0.102675125;
  double ref_molar_volume = 2.2414497292831042e-05;
  double ref_density = 4580.746275886329;
  double ref_molar_entropy = 254.97574377713337;
  double ref_molar_enthalpy = -221441.80395211757;
  double ref_isothermal_bulk_modulus_reuss = 359295734163.2916;
  double ref_isentropic_bulk_modulus_reuss = 381900523901.2901;
  double ref_isothermal_compressibility_reuss = 2.783222579385045e-12;
  double ref_isentropic_compressibility_reuss = 2.618482922685045e-12;
  double ref_shear_modulus = 193760526087.90207;
  double ref_p_wave_velocity = 11822.408460993425;
  double ref_bulk_sound_velocity = 9130.761700467749;
  double ref_shear_wave_velocity = 6503.76040770877;
  double ref_grueneisen_parameter = 1.4347327542714985;
  double ref_thermal_expansivity = 2.1925393775858323e-05;
  double ref_molar_heat_capacity_v = 123.07148837536111;
  double ref_molar_heat_capacity_p = 130.81442783426158;
  double ref_isentropic_thermal_gradient = 7.513646431353596e-09;

  // Make solution
  Solution bdg;
  bdg.set_solution_model(bdg_solution_model);
  bdg.set_method(types::EOSType::Auto);
  bdg.set_composition(molar_fractions);
  bdg.set_state(P, T);

  CHECK_THAT(bdg.get_excess_gibbs(),
    WithinRel(ref_excess_gibbs, tol_rel) || WithinAbs(ref_excess_gibbs, tol_abs));
  CHECK_THAT(bdg.get_excess_volume(),
    WithinRel(ref_excess_volume, tol_rel) || WithinAbs(ref_excess_volume, tol_abs));
  CHECK_THAT(bdg.get_excess_entropy(),
    WithinRel(ref_excess_entropy, tol_rel) || WithinAbs(ref_excess_entropy, tol_abs));
  CHECK_THAT(bdg.get_excess_enthalpy(),
    WithinRel(ref_excess_enthalpy, tol_rel) || WithinAbs(ref_excess_enthalpy, tol_abs));
  CHECK_THAT(bdg.get_molar_internal_energy(),
    WithinRel(ref_molar_internal_energy, tol_rel) || WithinAbs(ref_molar_internal_energy, tol_abs));
  CHECK_THAT(bdg.get_molar_gibbs(),
    WithinRel(ref_molar_gibbs, tol_rel) || WithinAbs(ref_molar_gibbs, tol_abs));
  CHECK_THAT(bdg.get_molar_helmholtz(),
    WithinRel(ref_molar_helmholtz, tol_rel) || WithinAbs(ref_molar_helmholtz, tol_abs));
  CHECK_THAT(bdg.get_molar_mass(),
    WithinRel(ref_molar_mass, tol_rel) || WithinAbs(ref_molar_mass, tol_abs));
  CHECK_THAT(bdg.get_molar_volume(),
    WithinRel(ref_molar_volume, tol_rel) || WithinAbs(ref_molar_volume, tol_abs));
  CHECK_THAT(bdg.get_density(),
    WithinRel(ref_density, tol_rel) || WithinAbs(ref_density, tol_abs));
  CHECK_THAT(bdg.get_molar_entropy(),
    WithinRel(ref_molar_entropy, tol_rel) || WithinAbs(ref_molar_entropy, tol_abs));
  CHECK_THAT(bdg.get_molar_enthalpy(),
    WithinRel(ref_molar_enthalpy, tol_rel) || WithinAbs(ref_molar_enthalpy, tol_abs));
  CHECK_THAT(bdg.get_isothermal_bulk_modulus_reuss(),
    WithinRel(ref_isothermal_bulk_modulus_reuss, tol_rel) || WithinAbs(ref_isothermal_bulk_modulus_reuss, tol_abs));
  CHECK_THAT(bdg.get_isentropic_bulk_modulus_reuss(),
    WithinRel(ref_isentropic_bulk_modulus_reuss, tol_rel) || WithinAbs(ref_isentropic_bulk_modulus_reuss, tol_abs));
  CHECK_THAT(bdg.get_isothermal_compressibility_reuss(),
    WithinRel(ref_isothermal_compressibility_reuss, tol_rel) || WithinAbs(ref_isothermal_compressibility_reuss, tol_abs));
  CHECK_THAT(bdg.get_isentropic_compressibility_reuss(),
    WithinRel(ref_isentropic_compressibility_reuss, tol_rel) || WithinAbs(ref_isentropic_compressibility_reuss, tol_abs));
  CHECK_THAT(bdg.get_shear_modulus(),
    WithinRel(ref_shear_modulus, tol_rel) || WithinAbs(ref_shear_modulus, tol_abs));
  CHECK_THAT(bdg.get_p_wave_velocity(),
    WithinRel(ref_p_wave_velocity, tol_rel) || WithinAbs(ref_p_wave_velocity, tol_abs));
  CHECK_THAT(bdg.get_bulk_sound_velocity(),
    WithinRel(ref_bulk_sound_velocity, tol_rel) || WithinAbs(ref_bulk_sound_velocity, tol_abs));
  CHECK_THAT(bdg.get_shear_wave_velocity(),
    WithinRel(ref_shear_wave_velocity, tol_rel) || WithinAbs(ref_shear_wave_velocity, tol_abs));
  CHECK_THAT(bdg.get_grueneisen_parameter(),
    WithinRel(ref_grueneisen_parameter, tol_rel) || WithinAbs(ref_grueneisen_parameter, tol_abs));
  CHECK_THAT(bdg.get_thermal_expansivity(),
    WithinRel(ref_thermal_expansivity, tol_rel) || WithinAbs(ref_thermal_expansivity, tol_abs));
  CHECK_THAT(bdg.get_molar_heat_capacity_v(),
    WithinRel(ref_molar_heat_capacity_v, tol_rel) || WithinAbs(ref_molar_heat_capacity_v, tol_abs));
  CHECK_THAT(bdg.get_molar_heat_capacity_p(),
    WithinRel(ref_molar_heat_capacity_p, tol_rel) || WithinAbs(ref_molar_heat_capacity_p, tol_abs));
  CHECK_THAT(bdg.get_isentropic_thermal_gradient(),
    WithinRel(ref_isentropic_thermal_gradient, tol_rel) || WithinAbs(ref_isentropic_thermal_gradient, tol_abs));
  CHECK(bdg.get_n_endmembers() == ref_n_endmembers);
  CHECK(bdg.get_n_elements() == ref_n_elements);
  CHECK(bdg.get_n_reactions() == ref_n_reactions);
  CHECK(bdg.get_independent_element_indices() == ref_independent_element_indices);
  CHECK(bdg.get_dependent_element_indices() == ref_dependent_element_indices);
  CHECK(bdg.get_excess_partial_volumes().isZero());
  CHECK(bdg.get_volume_hessian().isZero());
  CHECK(bdg.get_activity_coefficients().isOnes());
  CHECK(bdg.get_compositional_basis().isIdentity());
  CHECK(bdg.get_activities().isApprox(ref_activities, tol_rel));
  CHECK(bdg.get_excess_partial_gibbs().isApprox(ref_excess_partial_gibbs, tol_rel));
  CHECK(bdg.get_excess_partial_entropies().isApprox(ref_excess_partial_entropies, tol_rel));
  CHECK(bdg.get_partial_gibbs().isApprox(ref_partial_gibbs, tol_rel));
  CHECK(bdg.get_partial_entropies().isApprox(ref_partial_entropies, tol_rel));
  CHECK(bdg.get_gibbs_hessian().isApprox(ref_gibbs_hessian, tol_rel));
  CHECK(bdg.get_entropy_hessian().isApprox(ref_entropy_hessian, tol_rel));
  CHECK(bdg.get_stoichiometric_matrix().isApprox(ref_stoichiometric_matrix, tol_rel));
  CHECK(bdg.get_compositional_basis().isApprox(ref_compositional_basis, tol_rel));
  CHECK(bdg.get_compositional_null_basis().isApprox(ref_compositional_null_basis, tol_rel));
  Eigen::MatrixXd reaction_basis = bdg.get_reaction_basis();
  CHECK(reaction_basis.size() == 0);
  CHECK(reaction_basis.rows() == ref_reaction_basis.rows());
  CHECK(reaction_basis.cols() == ref_reaction_basis.cols());
}

TEST_CASE_METHOD(BridgmaniteFixture, "Test reactions", "[core][solution]") {
  // Make solution model & solution
  types::PairedEndmemberList em = {
    {mg_si_perovskite, "[Mg][Si]O3"},
    {fe_si_perovskite, "[Fe][Si]O3"},
    {mg_si_perovskite, "[Mg][Si]O3"},
    {fe_si_perovskite, "[Fe][Si]O3"}
  };
  std::shared_ptr<IdealSolution> sol_model;
  sol_model = std::make_shared<IdealSolution>(em);
  Solution sol;
  sol.set_solution_model(sol_model);
  REQUIRE(sol.get_n_endmembers() == 4);
  REQUIRE(sol.get_n_reactions() == 2);

  Eigen::MatrixXd expected_stoich_mat(4,4);
  expected_stoich_mat <<
    1, 0, 1, 3,
    0, 1, 1, 3,
    1, 0, 1, 3,
    0, 1, 1, 3;
  Eigen::MatrixXd expected_comp_basis(2, 4);
  expected_comp_basis <<
    1, 0, 0, 0,
    0, 1, 0, 0;
  Eigen::MatrixXd expected_comp_null_basis(2, 4);
  expected_comp_null_basis <<
    -1, -1, 1, 0,
    -3, -3, 0, 1;
  Eigen::MatrixXd expected_reac_basis(2, 4);
  expected_reac_basis <<
    -1, 0, 1, 0,
    0, -1, 0, 1;
  REQUIRE((sol.get_stoichiometric_matrix().array() == expected_stoich_mat.array()).all());
  REQUIRE((sol.get_compositional_basis().array() == expected_comp_basis.array()).all());
  REQUIRE((sol.get_compositional_null_basis().array() == expected_comp_null_basis.array()).all());
  REQUIRE((sol.get_reaction_basis().array() == expected_reac_basis.array()).all());
}
