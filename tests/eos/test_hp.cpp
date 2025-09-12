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
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "burnman/eos/hp.hpp"
#include <cmath>
#include <map>
#include <tuple>
#include "burnman/utils/types/mineral_params.hpp"
#include "tolerances.hpp"

using namespace Catch::Matchers;
using namespace burnman;

TEST_CASE("Test validate parameters", "[eos][hp]") {
  types::MineralParams params;
  params.molar_mass = 0.1003887;
  params.napfu = 5;
  params.V_0 = 2.445e-05;
  params.K_0 = 251000e6;
  params.Kprime_0 = 4.14;
  params.Kdprime_0 = -1.6e-11;
  params.Cp = CpParams{149.3, 0.002918, -2983000.0, -799.1};
  params.a_0 = 1.87e-05;
  params.H_0 = -1443030.0;
  params.S_0 = 62.6;
  eos::HP_TMT hp;
  SECTION("Missing V_0") {
    params.V_0.reset();
    REQUIRE_THROWS(hp.validate_parameters(params));
  }
  SECTION("Missing K_0") {
    params.K_0.reset();
    REQUIRE_THROWS(hp.validate_parameters(params));
  }
  SECTION("Missing Kprime_0") {
    params.Kprime_0.reset();
    REQUIRE_THROWS(hp.validate_parameters(params));
  }
  SECTION("Missing Kdprime_0") {
    params.Kdprime_0.reset();
    REQUIRE_THROWS(hp.validate_parameters(params));
  }
  SECTION("Missing molar mass") {
    params.molar_mass.reset();
    REQUIRE_THROWS(hp.validate_parameters(params));
  }
  SECTION("Missing napfu") {
    params.napfu.reset();
    REQUIRE_THROWS(hp.validate_parameters(params));
  }
  SECTION("Missing Cp params") {
    params.Cp.reset();
    REQUIRE_THROWS(hp.validate_parameters(params));
  }
  SECTION("Missing a_0") {
    params.a_0.reset();
    REQUIRE_THROWS(hp.validate_parameters(params));
  }
  SECTION("Optional Parameters") {
    params.H_0.reset();
    params.S_0.reset();
    REQUIRE_FALSE(params.T_0.has_value());
    REQUIRE_FALSE(params.P_0.has_value());
    REQUIRE_FALSE(params.E_0.has_value());
    REQUIRE_FALSE(params.H_0.has_value());
    REQUIRE_FALSE(params.F_0.has_value());
    REQUIRE_FALSE(params.G_0.has_value());
    REQUIRE_FALSE(params.Gprime_0.has_value());
    REQUIRE_FALSE(params.T_einstein.has_value());
    REQUIRE_NOTHROW(hp.validate_parameters(params));
    REQUIRE(*params.T_0 == 298.15);
    REQUIRE(*params.P_0 == 1.0e5);
    REQUIRE(*params.E_0 == 0);
    REQUIRE(std::isnan(*params.H_0));
    REQUIRE(std::isnan(*params.S_0));
    REQUIRE(std::isnan(*params.G_0));
    REQUIRE(std::isnan(*params.Gprime_0));
    REQUIRE(std::isnan(*params.T_einstein));
  }
  SECTION("Optional T_einstein") {
    params.S_0 = 53147.8;
    hp.validate_parameters(params);
    CHECK_THAT(*params.T_einstein,
      WithinRel(1.0, tol_rel) ||
      WithinAbs(1.0, tol_abs));
  }
}

TEST_CASE("Check reference conditions", "[eos][hp]") {
  // Set up test params
  types::MineralParams params;
  params.molar_mass = 0.1003887;
  params.napfu = 5;
  params.V_0 = 2.445e-05;
  params.K_0 = 251000e6;
  params.Kprime_0 = 4.14;
  params.Kdprime_0 = -1.6e-11;
  params.Cp = CpParams{149.3, 0.002918, -2983000.0, -799.1};
  params.a_0 = 1.87e-05;
  params.H_0 = -1443030.0;
  params.S_0 = 62.6;
  eos::HP_TMT hp;
  REQUIRE_NOTHROW(hp.validate_parameters(params));
  SECTION("V-T dependent functions") {
    double T = *params.T_0;
    double V = *params.V_0;
    CHECK_THAT(hp.compute_pressure(T, V, params),
      WithinRel(*params.P_0, tol_rel) || WithinAbs(*params.P_0, tol_abs));
  }
  SECTION("P-T dependent functions") {
    auto P = *params.P_0;
    auto T = *params.T_0;
    double V = GENERATE(2.445e-05, 1.0e-5, 2.0);
    CHECK_THAT(hp.compute_volume(P, T, params),
      WithinRel(*params.V_0, tol_rel) ||
      WithinAbs(*params.V_0, tol_abs));
    CHECK_THAT(hp.compute_entropy(P, T, V, params),
      WithinRel(*params.S_0, tol_rel) ||
      WithinAbs(*params.S_0, tol_abs));
    CHECK_THAT(hp.compute_enthalpy(P, T, V, params),
      WithinRel(*params.H_0, tol_rel) ||
      WithinAbs(*params.H_0, tol_abs));
    double G_0 = *params.H_0 - (*params.T_0 * *params.S_0);
    CHECK_THAT(hp.compute_gibbs_free_energy(P, T, V, params),
      WithinRel(G_0, tol_rel) ||
      WithinAbs(G_0, tol_abs));
    CHECK_THAT(hp.compute_thermal_expansivity(P, T, V, params),
      WithinRel(*params.a_0, tol_rel) ||
      WithinAbs(*params.a_0, tol_abs));
    CHECK_THAT(hp.compute_isothermal_bulk_modulus_reuss(P, T, V, params),
      WithinRel(*params.K_0, tol_rel) ||
      WithinAbs(*params.K_0, tol_abs));
  }
  SECTION("P-T-V dependent functions") {
    double P = *params.P_0;
    double T = *params.T_0;
    double V = *params.V_0;
    double F_0 = *params.H_0
      - (*params.T_0 * *params.S_0)
      - (*params.P_0 * *params.V_0);
    CHECK_THAT(hp.compute_helmholtz_free_energy(P, T, V, params),
      WithinRel(F_0, tol_rel) ||
      WithinAbs(F_0, tol_abs));
  }
  SECTION("Zero returns") {
    CHECK(hp.compute_shear_modulus(1.0, 1.0, 1.0, params) == 0);
  }
  SECTION("Zero T known functions") {
    auto P = GENERATE(0.0, 10.0, 25.e9);
    double V = GENERATE(2.445e-05, 1.0e-5, 2.0);
    double T = 0.0;
    CHECK_THAT(hp.compute_thermal_expansivity(P, T, V, params),
      WithinRel(0.0, tol_rel) || WithinAbs(0.0, tol_abs));
  }
}

TEST_CASE("Test volume", "[eos][hp]") {
  types::MineralParams params;
  params.molar_mass = 0.1003887;
  params.napfu = 5;
  params.V_0 = 2.445e-05;
  params.K_0 = 251000e6;
  params.Kprime_0 = 4.14;
  params.Kdprime_0 = -1.6e-11;
  params.Cp = CpParams{149.3, 0.002918, -2983000.0, -799.1};
  params.a_0 = 1.87e-05;
  params.H_0 = -1443030.0;
  params.S_0 = 62.6;
  eos::HP_TMT hp;
  REQUIRE_NOTHROW(hp.validate_parameters(params));
  double T_a = 800.0;
  double T_b = 2000.0;
  double V_a = 0.9 * (*params.V_0);
  double V_b = 0.5 * (*params.V_0);
  double P_aa = hp.compute_pressure(T_a, V_a, params);
  double P_ab = hp.compute_pressure(T_a, V_b, params);
  double P_ba = hp.compute_pressure(T_b, V_a, params);
  double P_bb = hp.compute_pressure(T_b, V_b, params);
  CHECK_THAT(hp.compute_volume(P_aa, T_a, params),
    WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
  CHECK_THAT(hp.compute_volume(P_ab, T_a, params),
    WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
  CHECK_THAT(hp.compute_volume(P_ba, T_b, params),
    WithinRel(V_a, tol_rel) || WithinAbs(V_a, tol_abs));
  CHECK_THAT(hp.compute_volume(P_bb, T_b, params),
    WithinRel(V_b, tol_rel) || WithinAbs(V_b, tol_abs));
}

TEST_CASE("HP python reference values", "[eos][hp]") {
  types::MineralParams params;
  params.molar_mass = 0.1003887;
  params.napfu = 5;
  params.V_0 = 2.445e-05;
  params.K_0 = 251000e6;
  params.Kprime_0 = 4.14;
  params.Kdprime_0 = -1.6e-11;
  params.Cp = CpParams{149.3, 0.002918, -2983000.0, -799.1};
  params.a_0 = 1.87e-05;
  params.H_0 = -1443030.0;
  params.S_0 = 62.6;
  eos::HP_TMT hp;
  REQUIRE_NOTHROW(hp.validate_parameters(params));

  SECTION("T-V dependent functions") {
    double T1 = 300, T2 = 800, T3 = 2500;
    double x1 = 0.99, x2 = 0.80, x3 = 0.40;
    // Reference data
    std::map<std::tuple<double, double>, double> ref_P = {
      {{T1, x1}, 2584473797.106938},
      {{T1, x2}, 88540858732.01479},
      {{T1, x3}, 1157572417438.9905},
      {{T2, x1}, 5396882398.7405405},
      {{T2, x2}, 91353267333.64839},
      {{T2, x3}, 1160384826040.6243},
      {{T3, x1}, 15888372858.290174},
      {{T3, x2}, 101844757793.19803},
      {{T3, x3}, 1170876316500.1738}
    };
    // Generate combos
    auto T = GENERATE(300.0, 800.0, 2500.0);
    auto x = GENERATE(0.99, 0.8, 0.4);
    CAPTURE(T);
    CAPTURE(x);
    double V = *params.V_0 * x;
    auto key = std::make_tuple(T, x);
    CHECK_THAT(hp.compute_pressure(T, V, params),
      WithinRel(ref_P[key], tol_rel) ||
      WithinAbs(ref_P[key], tol_abs));
  }
  SECTION("P-T dependent functions") {
    double V = 0.8 * (*params.V_0);
    double P1 = 1.0e9, P2 = 25.0e9, P3 = 125.0e9;
    double T1 = 300, T2 = 800, T3 = 2500;
    // Reference data
    std::map<std::tuple<double, double>, double> ref_V = {
      {{P1, T1}, 2.4354413813599477e-05},
      {{P1, T2}, 2.4630795305961834e-05},
      {{P1, T3}, 2.5831159689724342e-05},
      {{P2, T1}, 2.2486262470184822e-05},
      {{P2, T2}, 2.2670479187011876e-05},
      {{P2, T3}, 2.3429630586858282e-05},
      {{P3, T1}, 1.846281636193151e-05},
      {{P3, T2}, 1.8538591264063576e-05},
      {{P3, T3}, 1.883306402824874e-05}
    };
    std::map<std::tuple<double, double>, double> ref_S = {
      {{P1, T1}, 62.582620959889574},
      {{P1, T2}, 160.31082953959867},
      {{P1, T3}, 308.54961025974063},
      {{P2, T1}, 53.78394368426048},
      {{P2, T2}, 148.54406131658308},
      {{P2, T3}, 293.59688570391285},
      {{P3, T1}, 34.834192119861},
      {{P3, T2}, 123.74246577654709},
      {{P3, T3}, 264.9771227489722}
    };
    std::map<std::tuple<double, double>, double> ref_G = {
      {{P1, T1}, -1437410.3765602056},
      {{P1, T2}, -1495244.2221249093},
      {{P1, T3}, -1913961.3910165555},
      {{P2, T1}, -876799.7503743},
      {{P2, T2}, -929252.0500285919},
      {{P2, T3}, -1325243.6074070022},
      {{P3, T1}, 1141586.3457654482},
      {{P3, T2}, 1100600.7102598934},
      {{P3, T3}, 750305.9791734531}
    };
    std::map<std::tuple<double, double>, double> ref_alpha = {
      {{P1, T1}, 1.846299219677111e-05},
      {{P1, T2}, 2.4657610677261276e-05},
      {{P1, T3}, 3.133433379492222e-05},
      {{P2, T1}, 1.3452770701588973e-05},
      {{P2, T2}, 1.7688838234736712e-05},
      {{P2, T3}, 2.086821179738008e-05},
      {{P3, T1}, 6.815622480750097e-06},
      {{P3, T2}, 8.797504962203125e-06},
      {{P3, T3}, 9.588778354759562e-06}
    };
    std::map<std::tuple<double, double>, double> ref_Cp = {
      {{P1, T1}, 70.63144262748978},
      {{P1, T2}, 118.60787420371608},
      {{P1, T3}, 139.77953742761602},
      {{P2, T1}, 65.57197594679991},
      {{P2, T2}, 116.66064253797676},
      {{P2, T3}, 134.66581629115453},
      {{P3, T1}, 54.87355506059076},
      {{P3, T2}, 113.50212206177662},
      {{P3, T3}, 129.6261485045782}
    };
    std::map<std::tuple<double, double>, double> ref_KT = {
      {{P1, T1}, 255095763379.74985},
      {{P1, T2}, 243433361881.57788},
      {{P1, T3}, 198706472520.83252},
      {{P2, T1}, 350101194258.324},
      {{P2, T2}, 339337438868.3897},
      {{P2, T3}, 298364565092.2997},
      {{P3, T1}, 691034619656.831},
      {{P3, T2}, 682294024149.0497},
      {{P3, T3}, 649335578195.814}
    };
    std::map<std::tuple<double, double>, double> ref_H = {
      {{P1, T1}, -1418635.5902722387},
      {{P1, T2}, -1366995.5584932303},
      {{P1, T3}, -1142587.365367204},
      {{P2, T1}, -860664.5672690219},
      {{P2, T2}, -810416.8009753254},
      {{P2, T3}, -591251.3931472201},
      {{P3, T1}, 1152036.6034014064},
      {{P3, T2}, 1199594.682881131},
      {{P3, T3}, 1412748.7860458835}
    };
    // Generate combos
    auto P = GENERATE(1.0e9, 25.0e9, 125.0e9);
    auto T = GENERATE(300, 800, 2500);
    CAPTURE(P);
    CAPTURE(T);
    auto key = std::make_tuple(P, T);
    CHECK_THAT(hp.compute_volume(
        P, T, params
      ),
      WithinRel(ref_V[key], tol_rel) ||
      WithinAbs(ref_V[key], tol_abs));
    CHECK_THAT(hp.compute_entropy(
        P, T, V, params
      ),
      WithinRel(ref_S[key], tol_rel) ||
      WithinAbs(ref_S[key], tol_abs));
    CHECK_THAT(hp.compute_gibbs_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_G[key], tol_rel) ||
      WithinAbs(ref_G[key], tol_abs));
    CHECK_THAT(hp.compute_thermal_expansivity(
        P, T, V, params
      ),
      WithinRel(ref_alpha[key], tol_rel) ||
      WithinAbs(ref_alpha[key], tol_abs));
    CHECK_THAT(hp.compute_molar_heat_capacity_p(
        P, T, V, params
      ),
      WithinRel(ref_Cp[key], tol_rel) ||
      WithinAbs(ref_Cp[key], tol_abs));
    CHECK_THAT(hp.compute_isothermal_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KT[key], tol_rel) ||
      WithinAbs(ref_KT[key], tol_abs));
    CHECK_THAT(hp.compute_enthalpy(
        P, T, V, params
      ),
      WithinRel(ref_H[key], tol_rel) ||
      WithinAbs(ref_H[key], tol_abs));
  }
  SECTION("P-T-V dependent functions") {
    double P1 = 1.e9, P2 = 54.e9;
    double T1 = 800, T2 = 2500;
    double x1 = 0.99, x2 = 0.4;
    std::map<std::tuple<double, double, double>, double> ref_gamma = {
        {{P1, T1, x1}, 1.2553212527890947},
        {{P1, T1, x2}, 0.49982435028360145},
        {{P1, T2, x1}, 1.1776785174232474},
        {{P1, T2, x2}, 0.45103200390765347},
        {{P2, T1, x1}, 1.2774046649808242},
        {{P2, T1, x2}, 0.5119319566866284},
        {{P2, T2, x1}, 1.1929603871089238},
        {{P2, T2, x2}, 0.469295659521815}
    };
    std::map<std::tuple<double, double, double>, double> ref_Cv = {
        {{P1, T1, x1}, 115.74180861898594},
        {{P1, T1, x2}, 117.44986790685542},
        {{P1, T2, x1}, 127.97342236627493},
        {{P1, T2, x2}, 135.00938992808426},
        {{P2, T1, x1}, 113.74089681897487},
        {{P2, T1, x2}, 114.67208317564322},
        {{P2, T2, x1}, 126.33407777028971},
        {{P2, T2, x2}, 129.75520751174346}
    };
    std::map<std::tuple<double, double, double>, double> ref_KS = {
        {{P1, T1, x1}, 249461399537.01764},
        {{P1, T1, x2}, 245833512438.9918},
        {{P1, T2, x1}, 217038024765.32724},
        {{P1, T2, x2}, 205727163330.12305},
        {{P2, T1, x1}, 452659755086.5626},
        {{P2, T1, x2}, 448983964288.3453},
        {{P2, T2, x1}, 427233348513.5859},
        {{P2, T2, x2}, 415968901073.136}
    };
    std::map<std::tuple<double, double, double>, double> ref_F = {
        {{P1, T1, x1}, -1519449.7221249093},
        {{P1, T1, x2}, -1505024.2221249093},
        {{P1, T2, x1}, -1938166.8910165555},
        {{P1, T2, x2}, -1923741.3910165555},
        {{P2, T1, x1}, -1603757.2595452142},
        {{P2, T1, x2}, -824780.2595452144},
        {{P2, T2, x1}, -1981570.1899917815},
        {{P2, T2, x2}, -1202593.1899917815}
    };
    // Generate combos
    auto P = GENERATE(1.e9, 54.e9);
    auto T = GENERATE(800.0, 2500.0);
    auto x = GENERATE(0.99, 0.4);
    CAPTURE(P);
    CAPTURE(T);
    CAPTURE(x);
    double V = *params.V_0 * x;
    auto key = std::make_tuple(P, T, x);
    CHECK_THAT(hp.compute_grueneisen_parameter(
        P, T, V, params
      ),
      WithinRel(ref_gamma[key], tol_rel) ||
      WithinAbs(ref_gamma[key], tol_abs));
    CHECK_THAT(hp.compute_molar_heat_capacity_v(
        P, T, V, params
      ),
      WithinRel(ref_Cv[key], tol_rel) ||
      WithinAbs(ref_Cv[key], tol_abs));
    CHECK_THAT(hp.compute_isentropic_bulk_modulus_reuss(
        P, T, V, params
      ),
      WithinRel(ref_KS[key], tol_rel) ||
      WithinAbs(ref_KS[key], tol_abs));
    CHECK_THAT(hp.compute_helmholtz_free_energy(
        P, T, V, params
      ),
      WithinRel(ref_F[key], tol_rel) ||
      WithinAbs(ref_F[key], tol_abs));
  }
}
