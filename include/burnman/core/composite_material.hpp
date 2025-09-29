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
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/core/material.hpp"

namespace burnman {

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

  // Override of reset_cache to include additional solution properties
  void reset_cache() override;

  /**
   * @brief Computes and stores composite material properties.
   *
   * Convenience function to call getters/setup functions for
   * all stored material properties. This is useful to call
   * before cloning a material so computations only have to be
   * run once. Properties only need to be recomputed if changing
   * the endmember phases or chemistry.
   */
  void setup_composite_material_properties() const;

  /**
   * @brief Clears all stored composite material properties.
   *
   * Convenience function to clear all computed properties.
   * This includes material properties based on endmember
   * phases/chemistry (element list, stoichiometric matrix etc.)
   * as well as cached properties for the material P-T state.
   *
   * This function should only be required if changing phases
   * in the composite material. Cached properties can be cleared
   * with `reset_cache()'.
   *
   * @note Derived classes which store other computed properties
   * should extend this function.
   */
  void clear_computed_properties() override;

  // Public getter functions for stored composite material properties
  /**
   * @brief Number of endmembers in the material.
   */
  Eigen::Index get_n_endmembers() const;

  /**
   * @brief Number of different chemical elements in the material.
   */
  Eigen::Index get_n_elements() const;

  /**
   * @brief Number of possible reactions from reaction basis.
   */
  Eigen::Index get_n_reactions() const;

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
  const std::vector<types::FormulaMap>& get_endmember_formulae() const;

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
   * @brief Stoichiometric matrix restricted to independent elements.
   */
  const Eigen::MatrixXd& get_reduced_stoichiometric_matrix() const;

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

  /**
   * @brief Retrieves the endmember partial molar gibbs free energy.
   *
   * Uses a cached value if available, or calls
   * `CompositeMaterial::compute_partial_gibbs()` and caches the result.
   *
   * @note Use `CompositeMaterial::reset_cache()` to clear cached values.
   *
   * @return Partial molar gibbs free energies in [J/mol].
   */
  const Eigen::ArrayXd& get_partial_gibbs() const;

 protected:

  // Pure virtual funtions that must be implemented in derived classes
  virtual Eigen::Index compute_n_endmembers() const = 0;
  // Void functions should use protected setters!
  virtual void setup_endmember_formulae() const = 0;
  virtual void setup_endmember_names() const = 0;

  /**
   * @brief Computes endmember partial molar gibbs free energies.
   *
   * @return Partial molar gibbs free energies in [J/mol].
   */
  virtual Eigen::ArrayXd compute_partial_gibbs() const = 0;

  // Protected setters
  void set_endmember_names(std::vector<std::string> names) const;
  void set_endmember_formulae(std::vector<types::FormulaMap> formulae) const;

 private:

  // Cached CompositeMaterial properties
  // Note: Not cleared by reset_cache()! Use clear_computed_properties() instead.
  // Use setup_composite_material_properties() to initialise all values.
  mutable std::optional<Eigen::Index> n_endmembers;
  mutable std::optional<Eigen::Index> n_elements;
  mutable std::optional<Eigen::Index> n_reactions;
  mutable std::optional<std::vector<Eigen::Index>> independent_element_indices;
  mutable std::optional<std::vector<Eigen::Index>> dependent_element_indices;
  mutable std::optional<std::vector<std::string>> endmember_names;
  mutable std::optional<std::vector<std::string>> elements;
  mutable std::optional<std::vector<types::FormulaMap>> endmember_formulae;
  mutable std::optional<Eigen::MatrixXd> stoichiometric_matrix;
  mutable std::optional<Eigen::MatrixXd> reduced_stoichiometric_matrix;
  mutable std::optional<Eigen::MatrixXd> compositional_basis;
  mutable std::optional<Eigen::MatrixXd> compositional_null_basis;
  mutable std::optional<Eigen::MatrixXd> reaction_basis;

  // Cached properties cleared by reset_cache
  mutable std::optional<Eigen::ArrayXd> partial_gibbs;

  // Compute functions common to all composite materials.
  Eigen::Index compute_n_elements() const;
  Eigen::Index compute_n_reactions() const;
  std::vector<std::string> compute_elements() const;
  std::vector<Eigen::Index> compute_independent_element_indices() const;
  std::vector<Eigen::Index> compute_dependent_element_indices() const;
  Eigen::MatrixXd compute_stoichiometric_matrix() const;
  Eigen::MatrixXd compute_reduced_stoichiometric_matrix() const;
  Eigen::MatrixXd compute_compositional_basis() const;
  Eigen::MatrixXd compute_compositional_null_basis() const;
  Eigen::MatrixXd compute_reaction_basis() const;

};

} // namespace burnman

#endif // BURNMAN_CORE_COMPOSITE_MATERIAL_HPP_INCLUDED
