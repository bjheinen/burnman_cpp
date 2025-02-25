/*
  TODO: Copyright Notice!
*/
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "tolerances.hpp"
#include "burnman/eos/vinet.hpp"
#include "burnman/utils/eos.hpp"
#include <cmath>
#include <vector>

using namespace Catch::Matchers;



TEST_CASE("Check reference volume", "[vinet][eos]") {
  // Set up test params
  MineralParams params;
  params.equation_of_state = EOSType::Vinet;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  // Set up test P, T, V
  auto P = GENERATE(0.0, 10.0, 25.e9);
  auto T = GENERATE(300.0, 2000.0);
  double V = *params.V_0;
  Vinet vinet;
  CHECK_THAT(vinet.compute_isothermal_bulk_modulus_reuss(P, T, V, params),
    WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
  CHECK_THAT(vinet.compute_isentropic_bulk_modulus_reuss(P, T, V, params),
    WithinRel(*params.K_0, tol_rel) || WithinAbs(*params.K_0, tol_abs));
  CHECK_THAT(vinet.compute_molar_internal_energy(P, T, V, params),
    WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
  CHECK_THAT(vinet.compute_gibbs_free_energy(P, T, V, params),
    WithinRel(P*V, tol_rel) || WithinAbs(P*V, tol_abs));
  CHECK_THAT(vinet.compute_volume(0.0, T, params),
    WithinRel(V, tol_rel) || WithinAbs(V, tol_abs));
}

TEST_CASE("Check hard-coded returns", "[vinet][eos]") {
  // Set up test P, T, V
  auto P = GENERATE(0.0, 25.e9);
  auto T = GENERATE(300.0, 2000.0);
  auto V = GENERATE(5.0e-6, 10.0e-6);
  // Set up test params
  MineralParams params;
  params.equation_of_state = EOSType::Vinet;
  params.T_0 = 300.0;
  params.P_0 = 0.0;
  params.V_0 = 11.24e-6;
  params.K_0 = 161.0e9;
  params.Kprime_0 = 3.8;
  params.molar_mass = 0.0403;
  params.napfu = 2;
  // Unused parameters
  params.G_0 = 131.0e9;
  params.Gprime_0 = 2.1;
  params.grueneisen_0 = 1.5;
  Vinet vinet;
  CHECK(vinet.compute_grueneisen_parameter(P, T, V, params) == 0);
  CHECK(vinet.compute_shear_modulus(P, T, V, params) == 0);
  CHECK(vinet.compute_thermal_expansivity(P, T, V, params) == 0);
  CHECK(vinet.compute_entropy(P, T, V, params) == 0);
  CHECK(vinet.compute_molar_heat_capacity_v(P, T, V, params) == 1.0e99);
  CHECK(vinet.compute_molar_heat_capacity_v(P, T, V, params) == 1.0e99);
}