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
#include "burnman/utils/chemistry_utils.hpp"
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"

using namespace burnman;

TEST_CASE("sort_element_list_to_IUPAC_order (sort)", "[utils][chemistry_utils]") {
  std::unordered_set<std::string> input = {
    "H", "Si", "O", "He", "Na", "Ni", "K",
    "Li", "Fe", "Ti", "Mn", "Al", "C", "Pt"};
  std::vector<std::string> expected = {
    "He", "K", "Na", "Li", "Ti", "Mn", "Fe",
    "Pt", "Ni", "Al", "Si", "C", "O", "H"};
  std::vector<std::string> result = utils::sort_element_list_to_IUPAC_order(input);
  REQUIRE(result.size() == expected.size());
  REQUIRE(result == expected);
}

TEST_CASE("sort_element_list_to_IUPAC_order (empty input)", "[utils][chemistry]") {
  std::unordered_set<std::string> input;
  std::vector<std::string> result = utils::sort_element_list_to_IUPAC_order(input);
  REQUIRE(result.empty());
}

TEST_CASE("sum_formulae (weights)", "[utils][chemistry_utils]") {
  types::FormulaMap f1 = {{"Si", 1}, {"O", 2}};
  types::FormulaMap f2 = {{"Mg", 1}, {"O", 1}};
  types::FormulaMap f3 = {{"Fe", 2}, {"O", 3}};
  std::vector<types::FormulaMap> formulae = {f1, f2, f3};
  Eigen::ArrayXd weights(3);
  weights << 0.5, 0.25, 0.25;
  types::FormulaMap result = utils::sum_formulae(formulae, weights);
  CHECK(result["Si"] == 0.5);
  CHECK(result["O"] == 2.0); // 0.5 * 2 + 1.0 * 1
  CHECK(result["Mg"] == 0.25);
  CHECK(result["Fe"] == 0.5);
}

TEST_CASE("sum_formulae (errors)", "[utils][chemistry_utils]") {
  types::FormulaMap f1 = {{"Si", 1}, {"O", 2}};
  std::vector<types::FormulaMap> formulae = {f1};
  Eigen::ArrayXd weights(2);
  weights << 1.0, 1.0;
  CHECK_THROWS_AS(utils::sum_formulae(formulae, weights), std::invalid_argument);
}

TEST_CASE("sum_formulae (no weights)", "[utils][chemistry_utils]") {
  types::FormulaMap f1 = {{"Si", 1}, {"O", 2}};
  types::FormulaMap f2 = {{"Mg", 1}, {"O", 1}};
  std::vector<types::FormulaMap> formulae = {f1, f2};
  Eigen::ArrayXd weights(2);
  weights << 1.0, 1.0;
  types::FormulaMap result = utils::sum_formulae(formulae);
  types::FormulaMap result_w = utils::sum_formulae(formulae, weights);
  CHECK(result == result_w);
}
