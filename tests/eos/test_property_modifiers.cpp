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
#include "tolerances.hpp"
#include "burnman/eos/property_modifiers.hpp"
#include "burnman/utils/eos.hpp"
using namespace Catch::Matchers;

TEST_CASE("Test Landau", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;
  SECTION("Landau A") {
    double G_ref = 1137.8563285203;
    double dGdT_ref = 0.44067528178863924;
    double dGdP_ref = 9.594794189827874e-09;
    double d2GdT2_ref = -0.0015112986893857335;
    double d2GdP2_ref = -1.35546123054467e-19;
    double d2GdPdT_ref = 1.4885723933047452e-11;
    excesses::LandauParams params = {800.0, 1.0e-7, 5.0};
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
  SECTION("Landau B") {
    double G_ref = 1019.371598207319;
    double dGdT_ref = 0.44133754667682723;
    double dGdP_ref = 8.987151749036388e-09;
    double d2GdT2_ref = -0.0012366676242447335;
    double d2GdP2_ref = -1.090365427626495e-19;
    double d2GdPdT_ref = 1.2587509746776752e-11;
    excesses::LandauParams params = {1200.0, 1.0e-7, 5.0};
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
}

TEST_CASE("Test Landau HP", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;
  SECTION("Landau HP A") {
    double G_ref = -2534.1849669400353;
    double dGdT_ref = 3.5398427428570116;
    double dGdP_ref = -7.079685485714023e-08;
    double d2GdT2_ref = -0.0020833344907417047;
    double d2GdP2_ref = -8.333337962966819e-19;
    double d2GdPdT_ref = 4.1666689814834094e-11;
    excesses::LandauHPParams params = {
      298.15, 1.0e5, 800.0, 1.0e-7, 5.0
    };
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
  SECTION("Landau HP B") {
    double G_ref = -1696.3556185260459;
    double dGdT_ref = 2.4354537839959924;
    double dGdP_ref = -4.870907567991984e-08;
    double d2GdT2_ref = -0.0015386443366256072;
    double d2GdP2_ref = -6.154577346502429e-19;
    double d2GdPdT_ref = 3.077288673251214e-11;
    excesses::LandauHPParams params = {
      298.15, 1.0e5, 1200.0, 1.0e-7, 5.0
    };
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
}

TEST_CASE("Test Landau SLB", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;

  SECTION("Landau SLB A") {
    double G_ref = -1333.333333333334;
    double dGdT_ref = 2.5;
    double dGdP_ref = -5e-08;
    double d2GdT2_ref = -0.0005952380952380953;
    double d2GdP2_ref = -8.333333333333332e-19;
    double d2GdPdT_ref = 1.1904761904761904e-11;
    excesses::LandauSLB2022Params params = {800.0, 1.0e-7, 5.0};
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
  SECTION("Landau SLB B") {
    double G_ref = -929.380272332839;
    double dGdT_ref = 1.7700320038632995;
    double dGdP_ref = -3.5400640077266006e-08;
    double d2GdT2_ref = -0.0005769913639656223;
    double d2GdP2_ref = -6.154574548966636e-19;
    double d2GdPdT_ref = 1.1539827279312444e-11;
    excesses::LandauSLB2022Params params = {1200.0, 1.0e-7, 5.0};
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
}

TEST_CASE("Test Linear", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;
  double G_ref = 6200.0;
  double dGdT_ref = -5.0;
  double dGdP_ref = 1e-07;
  double d2GdT2_ref = 0.0;
  double d2GdP2_ref = 0.0;
  double d2GdPdT_ref = 0.0;
  excesses::LinearParams params = {1.0e-7, 5.0, 1200.0};
  excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
  CHECK_THAT(calc_excess.G,
    WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdT,
    WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdP,
    WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdT2,
    WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdP2,
    WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdPdT,
    WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
}

TEST_CASE("Test Bragg-Williams", "[prop_mod][eos]") {
  excesses::BraggWilliamsParams params = {
    1, 0.8, 1000.0, 1.0e-7, 1000.0, 1.0e-7};
  SECTION("BW High P") {
    double P = 1.0e11;
    double T = 1000.0;
    double G_ref = -551.347145218419;
    double dGdT_ref = -2.5551956240190066;
    double dGdP_ref = 1.8216804392068298e-08;
    double d2GdT2_ref = -0.007757316802781132;
    double d2GdP2_ref = -1.3642420526593924e-18;
    double d2GdPdT_ref = 7.051539796520956e-11;
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
  SECTION("BW Low P") {
    // No root for lower P. Make sure we are catching Brent
    // exception in order_gibbs.
    double P = 5.e10;
    double T = 2000.0;
    double G_ref = -12442.06822926073;
    double dGdT_ref = -9.221034114634676;
    double dGdP_ref = 1.0000000111176633e-07;
    double d2GdT2_ref = 0.0;
    double d2GdP2_ref = 0.0;
    double d2GdPdT_ref = 0.0;
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
}

TEST_CASE("Test Magnetic Chs", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;
  SECTION("Magnetic Chs A") {
    double G_ref = -14069.155333136323;
    double dGdT_ref = 19.194307974919088;
    double dGdP_ref = -2.2610537234163914e-07;
    double d2GdT2_ref = -0.0068214881170452;
    double d2GdP2_ref = -9.82760289180444e-19;
    double d2GdPdT_ref = 9.425343695892169e-11;
    excesses::MagneticChsParams params = {
      0.4, 800.0, 1.0e-8, 2.2, 1.0e-10};
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
  SECTION("Magnetic Chs B") {
    double G_ref = -21582.56374006573;
    double dGdT_ref = 20.218806497890633;
    double dGdP_ref = -2.5337465107094005e-07;
    double d2GdT2_ref = -0.003710226214794699;
    double d2GdP2_ref = -7.123458198725806e-19;
    double d2GdPdT_ref = 7.62289687062366e-11;
    excesses::MagneticChsParams params = {
      0.4, 1200.0, 1.0e-8, 2.2, 1.0e-10};
    excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
    CHECK_THAT(calc_excess.G,
      WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdT,
      WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
    CHECK_THAT(calc_excess.dGdP,
      WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdT2,
      WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdP2,
      WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
    CHECK_THAT(calc_excess.d2GdPdT,
      WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
  }
}

TEST_CASE("Test Debye", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;
  double G_ref = -565.3149805314973;
  double dGdT_ref = -1.186112668544687;
  double dGdP_ref = 0.0;
  double d2GdT2_ref = -0.0009315448135567599;
  double d2GdP2_ref = 0.0;
  double d2GdPdT_ref = 0.0;
  excesses::DebyeParams params = {1.0, 1200.0};
  excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
  CHECK_THAT(calc_excess.G,
    WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdT,
    WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdP,
    WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdT2,
    WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdP2,
    WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdPdT,
    WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
}

TEST_CASE("Test Debye Delta", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;
  double G_ref = -620.7976880131898;
  double dGdT_ref = -0.9315448135567599;
  double dGdP_ref = 0.0;
  double d2GdT2_ref = -0.0001301242749913705;
  double d2GdP2_ref = 0.0;
  double d2GdPdT_ref = 0.0;
  excesses::DebyeDeltaParams params = {1.0, 1200.0};
  excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
  CHECK_THAT(calc_excess.G,
    WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdT,
    WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdP,
    WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdT2,
    WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdP2,
    WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdPdT,
    WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
}

TEST_CASE("Test Einstein", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;
  double G_ref = -358.3824178604341;
  double dGdT_ref = -0.8755977306924339;
  double dGdP_ref = 0.0;
  double d2GdT2_ref = -0.0008881700552263031;
  double d2GdP2_ref = 0.0;
  double d2GdPdT_ref = 0.0;
  excesses::EinsteinParams params = {1.0, 1200.0};
  excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
  CHECK_THAT(calc_excess.G,
    WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdT,
    WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdP,
    WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdT2,
    WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdP2,
    WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdPdT,
    WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
}

TEST_CASE("Test Einstein Delta", "[prop_mod][eos]") {
  double P = 1.0e11;
  double T = 1000.0;
  double G_ref = -517.2153128319998;
  double dGdT_ref = -0.8881700552263031;
  double dGdP_ref = 0.0;
  double d2GdT2_ref = -0.00020821426174273138;
  double d2GdP2_ref = 0.0;
  double d2GdPdT_ref = 0.0;
  excesses::EinsteinDeltaParams params = {1.0, 1200.0};
  excesses::Excesses calc_excess = excesses::compute_excesses(P, T, params);
  CHECK_THAT(calc_excess.G,
    WithinRel(G_ref, tol_rel) || WithinAbs(G_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdT,
    WithinRel(dGdT_ref, tol_rel) || WithinAbs(dGdT_ref, tol_abs));
  CHECK_THAT(calc_excess.dGdP,
    WithinRel(dGdP_ref, tol_rel) || WithinAbs(dGdP_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdT2,
    WithinRel(d2GdT2_ref, tol_rel) || WithinAbs(d2GdT2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdP2,
    WithinRel(d2GdP2_ref, tol_rel) || WithinAbs(d2GdP2_ref, tol_abs));
  CHECK_THAT(calc_excess.d2GdPdT,
    WithinRel(d2GdPdT_ref, tol_rel) || WithinAbs(d2GdPdT_ref, tol_abs));
}
