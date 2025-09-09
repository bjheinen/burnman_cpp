/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_COMPOSITE_MATERIAL_HPP_INCLUDED
#define BURNMAN_CORE_COMPOSITE_MATERIAL_HPP_INCLUDED

#include <optional>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "burnman/core/material.hpp"
#include "burnman/utils/types/simple_types.hpp"

/**
 * @class CompositeMaterial
 * @brief Abtract base class for composite materials (solutions, assemblages).
 *
 * This class is an intermediate abstraction layer for composite materials.
 * It defines the interface and implementation for several properties
 * calculated from stoichiometry, e.g. stoichiometric matrices,
 * reaction basis, number of endmembers, element lists, etc.
 *
 */
class CompositeMaterial : public Material {

 public:

  virtual ~CompositeMaterial() = default;

  // Utility functions?
  // map_to_array etc.
  // get_partial_volumes

  // Public getter functions for stored composite material properties

  /**
   * @brief Number of endmembers in the material.
   */
  int get_n_endmembers() const;

  /**
   * @brief Number of different chemical elements in the material.
   */
  int get_n_elements() const;

  /**
   * @brief Number of possible reactions from reaction basis.
   */
  int get_n_reactions() const;

  /**
   * @brief Vector of chemical elements in the material in IUPAC order.
   */
  const std::vector<std::string>& get_elements() const;

  /**
   * @brief Vector of names of all endmembers in the material.
   */
  const std::vector<std::string>& get_endmember_names() const;

  /**
   * @brief Vector of chemical formula of all endmembers in the material.
   *
   * Each formula is a FormulaMap, where
   * FormulaMap = std::unordered_map<std::string, double>>;
   */
  const std::vector<FormulaMap>& get_endmember_formulae() const;

  /**
   * @brief The independent set of element indices.
   *
   * If the amounts of these elements are known the amounts of other
   * elements can be inferred by:
   *   -compositional_null_basis[independent_element_indices].dot(element_amounts)
   */
  const std::vector<Eigen::Index>& get_independent_element_indices() const;

  /**
   * @brief The element indices not in the independent list.
   */
  const std::vector<Eigen::Index>& get_dependent_element_indices() const;

  /**
   * @brief The matrix describing material stoichiometry.
   *
   * Each element M[i,j] corresponds to the number of
   * atoms of element j in endmember i.
   */
  const Eigen::MatrixXd& get_stoichiometric_matrix() const;

  /**
   * @brief The compositional basis of the material.
   */
  const Eigen::MatrixXd& get_compositional_basis() const;

  /**
   * @brief The compositional null basis of the material.
   *
   * The matrix where N[b] = 0 for all bulk compositions that
   * can be produced with a linear sum of the endmembers in the material.
   */
  const Eigen::MatrixXd& get_compositional_null_basis() const;

  /**
   * @brief The reaction basis of the material.
   *
   * Each element M[i,j] corresponds to the number of moles
   * of endmember j involved in reaction i.
   */
  const Eigen::MatrixXd& get_reaction_basis() const;

 protected:

  // Pure virtual funtions that must be implemented in derived classes

  virtual int compute_n_endmembers() const = 0;
  // Void functions should use protected setters!
  virtual void setup_endmember_formulae() const = 0;
  virtual void setup_endmember_names() const = 0;

  // Protected setters
  void set_endmember_names(std::vector<std::string> names) const;
  void set_endmember_formulae(std::vector<FormulaMap> formulae) const; 

 private:

  // Cached CompositeMaterial properties
  // Note: Not cleared by reset()!
  mutable std::optional<int> n_endmembers;
  mutable std::optional<int> n_elements;
  mutable std::optional<int> n_reactions;
  mutable std::optional<std::vector<Eigen::Index>> independent_element_indices;
  mutable std::optional<std::vector<Eigen::Index>> dependent_element_indices;
  mutable std::optional<std::vector<std::string>> endmember_names;
  mutable std::optional<std::vector<std::string>> elements;
  mutable std::optional<std::vector<FormulaMap>> endmember_formulae;
  mutable std::optional<Eigen::MatrixXd> stoichiometric_matrix;
  // stoichiometric_array?
  mutable std::optional<Eigen::MatrixXd> compositional_basis;
  mutable std::optional<Eigen::MatrixXd> compositional_null_basis;
  mutable std::optional<Eigen::MatrixXd> reaction_basis;

  // Compute functions common to all composite materials.
  int compute_n_elements() const;
  int compute_n_reactions() const;
  std::vector<std::string> compute_elements() const;
  std::vector<Eigen::Index> compute_independent_element_indices() const;
  std::vector<Eigen::Index> compute_dependent_element_indices() const;
  Eigen::MatrixXd compute_stoichiometric_matrix() const;
  Eigen::MatrixXd compute_compositional_basis() const;
  Eigen::MatrixXd compute_compositional_null_basis() const;
  Eigen::MatrixXd compute_reaction_basis() const;

  // TODO: combined helpers to set all properties
  // e.g. (from old Solution form)
  // void Solution::setup_solution_properties() {
  //   setup_endmember_names();
  //   setup_endmember_formulae();
  //   setup_elements();
  //   setup_stoichiometric_matrix();
  //   setup_independent_element_indices();
  //   setup_dependent_element_indices();
  //   setup_reaction_basis();
  //   setup_n_reactions();
  //   setup_compositional_basis();
  //   setup_compositional_null_basis();
  // }

};

#endif // BURNMAN_CORE_COMPOSITE_MATERIAL_HPP_INCLUDED
