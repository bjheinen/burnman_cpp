/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/core/composite_material.hpp"
#include <algorithm>
#include <utility>
#include <unordered_set>
#include "burnman/utils/chemistry_utils.hpp"
#include "burnman/utils/matrix_utils.hpp"

namespace burnman {

void CompositeMaterial::reset() {
  // Reset caches Material properties
  Material::reset();
  // Reset cached CompositeMaterial properties
  partial_gibbs.reset();
}

int CompositeMaterial::get_n_endmembers() const {
  if (!n_endmembers.has_value()) {
    n_endmembers = compute_n_endmembers();
  }
  return *n_endmembers;
}

int CompositeMaterial::get_n_elements() const {
  if (!n_elements.has_value()) {
    n_elements = compute_n_elements();
  }
  return *n_elements;
}

int CompositeMaterial::get_n_reactions() const {
  if (!n_reactions.has_value()) {
    n_reactions = compute_n_reactions();
  }
  return *n_reactions;
}

const std::vector<std::string>& CompositeMaterial::get_elements() const {
  if (!elements.has_value()) {
    elements = compute_elements();
  }
  return *elements;
}

const std::vector<Eigen::Index>& CompositeMaterial::get_independent_element_indices() const {
  if (!independent_element_indices.has_value()) {
    independent_element_indices = compute_independent_element_indices();
  }
  return *independent_element_indices;
}

const std::vector<Eigen::Index>& CompositeMaterial::get_dependent_element_indices() const {
  if (!dependent_element_indices.has_value()) {
    dependent_element_indices = compute_dependent_element_indices();
  }
  return *dependent_element_indices;
}

const Eigen::MatrixXd& CompositeMaterial::get_stoichiometric_matrix() const {
  if (!stoichiometric_matrix.has_value()) {
    stoichiometric_matrix = compute_stoichiometric_matrix();
  }
  return *stoichiometric_matrix;
}

const Eigen::MatrixXd& CompositeMaterial::get_reduced_stoichiometric_matrix() const {
  if (!reduced_stoichiometric_matrix.has_value()) {
    reduced_stoichiometric_matrix = compute_reduced_stoichiometric_matrix();
  }
  return *reduced_stoichiometric_matrix;
}

const Eigen::MatrixXd& CompositeMaterial::get_compositional_basis() const {
  if (!compositional_basis.has_value()) {
    compositional_basis = compute_compositional_basis();
  }
  return *compositional_basis;
}

const Eigen::MatrixXd& CompositeMaterial::get_compositional_null_basis() const {
  if (!compositional_null_basis.has_value()) {
    compositional_null_basis = compute_compositional_null_basis();
  }
  return *compositional_null_basis;
}

const Eigen::MatrixXd& CompositeMaterial::get_reaction_basis() const {
  if (!reaction_basis.has_value()) {
    reaction_basis = compute_reaction_basis();
  }
  return *reaction_basis;
}

const std::vector<std::string>& CompositeMaterial::get_endmember_names() const {
  if (!endmember_names.has_value()) {
    setup_endmember_names();
  }
  return *endmember_names;
}

const std::vector<FormulaMap>& CompositeMaterial::get_endmember_formulae() const {
  if (!endmember_formulae.has_value()) {
    setup_endmember_formulae();
  }
  return *endmember_formulae;
}

const Eigen::ArrayXd& CompositeMaterial::get_partial_gibbs() const {
  if (!partial_gibbs.has_value()) {
    partial_gibbs = compute_partial_gibbs();
  }
  return *partial_gibbs;
}

// Setters
void CompositeMaterial::set_endmember_names(std::vector<std::string> names) const {
  endmember_names = std::move(names);
}
void CompositeMaterial::set_endmember_formulae(std::vector<FormulaMap> formulae) const {
  endmember_formulae = std::move(formulae);
}

// Compute functions common to all composite materials.

int CompositeMaterial::compute_n_elements() const {
  return static_cast<int>(get_elements().size());
}

int CompositeMaterial::compute_n_reactions() const {
  return static_cast<int>(get_reaction_basis().rows());
}

std::vector<std::string> CompositeMaterial::compute_elements() const {
  // Grab set of all elements in formulae
  std::unordered_set<std::string> all_elements;
  for (const auto& embr : get_endmember_formulae()) {
    for (const auto& [element, n] : embr) {
      all_elements.insert(element);
    }
  }
  // Sort elements into IUPAC order
  return utils::sort_element_list_to_IUPAC_order(all_elements);
}

std::vector<Eigen::Index> CompositeMaterial::compute_independent_element_indices() const {
  return utils::get_independent_row_indices(get_stoichiometric_matrix());
}

std::vector<Eigen::Index> CompositeMaterial::compute_dependent_element_indices() const {
  const std::vector<std::string>& elems = get_elements();
  const std::vector<Eigen::Index>& indep = get_independent_element_indices();
  std::vector<Eigen::Index> dep;
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(elems.size()); ++i) {
    if (std::find(
      indep.begin(), indep.end(), i) == indep.end()
    ) {
      dep.push_back(i);
    }
  }
  return dep;
}

Eigen::MatrixXd CompositeMaterial::compute_stoichiometric_matrix() const {
  int n_elem = get_n_elements();
  int n_embr = get_n_endmembers();
  const std::vector<std::string> elems = get_elements();
  const std::vector<FormulaMap>& embr_formulae = get_endmember_formulae();
  Eigen::MatrixXd stoich_mat;
  stoich_mat.resize(n_embr, n_elem);
  stoich_mat.setZero();
  for (int i = 0; i < n_embr; ++i) {
    const auto& formula_i = embr_formulae[i];
    for (int j = 0; j < n_elem; ++j) {
      const std::string& element = elems[j];
      auto it = formula_i.find(element);
      if (it != formula_i.end()) {
        stoich_mat(i, j) = it->second;
      }
    }
  }
  return stoich_mat;
}

Eigen::MatrixXd CompositeMaterial::compute_reduced_stoichiometric_matrix() const {
  return (get_stoichiometric_matrix())(Eigen::all, get_independent_element_indices());
}

Eigen::MatrixXd CompositeMaterial::compute_compositional_basis() const {
  const Eigen::MatrixXd& reac_basis_mat = get_reaction_basis();
  return utils::complete_basis(reac_basis_mat)
    (Eigen::seq(get_n_reactions(), Eigen::last), Eigen::all);
}

Eigen::MatrixXd CompositeMaterial::compute_compositional_null_basis() const {
  const Eigen::MatrixXd& stoich_mat = get_stoichiometric_matrix();
  Eigen::MatrixXd nullspace = utils::nullspace(stoich_mat);
  // Maybe consider threshold?
  // Check matrix:
  // const std::vector<Eigen::Index>& dep = get_dependent_element_indices();
  // Eigen::MatrixXd M(nullspace.rows(), static_cast<Eigen::Index>(dep.size()));
  // M = nullspace(Eigen::all, dep);
  // assert(M.cols() == M.rows());
  // assert(M.isIdentity(constants::precision::abs_tolerance));
  return nullspace;
}

Eigen::MatrixXd CompositeMaterial::compute_reaction_basis() const {
  const Eigen::MatrixXd& stoich_mat = get_stoichiometric_matrix();
  Eigen::MatrixXd nullspace = utils::nullspace(stoich_mat.transpose());
  if (nullspace.rows() == 0) {
    return Eigen::MatrixXd(0, get_n_endmembers());
  } else {
    return nullspace;
  }
}

} // namespace burnman
