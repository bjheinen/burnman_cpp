/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_SOLUTION_HPP_INCLUDED
#define BURNMAN_CORE_SOLUTION_HPP_INCLUDED

#include <optional>
#include <memory>
#include <Eigen/Dense>
#include "burnman/core/material.hpp"
#include "burnman/core/mineral.hpp"
#include "burnman/core/solution_model.hpp"

/**
 * @class Solution
 * @brief Base class for solid solutions.
 *
 * Base class for all solutions. 
 *
 * Site occupancies, endmember activities and the constant
 * pressure and temperature dependencies of the excess properties
 * require set_composition().
 * Solution states equire set_state().
 *
 * Uses an instance of `SolutionModel' (or derived class) to calculate
 * interaction terms between endmembers.
 *
 * All solution parameters and properties are in SI units.
 * e.g. Interaction parameters in J/mol, with T & P derivatives in J/K/mol
 * and m^3/mol.
 *
 */
class Solution : public Material {

 public:

  virtual ~Solution() = default;

  // Override of reset to include additional solution properties
  void reset() override;

  // Utility functions
  /**
   * @brief Utility function to map endmember properties to Array.
   *
   * Usage:
   *  map_endmembers_to_array(&Mineral::get_property);
   *    where, get_property() is the function needed
   *  alternatively, with a lambda function, e.g.
   * map_endmembers_to_array([](const Mineral& m) {
   *  return m.get_property() + 10.0;
   * });
   *
   * @param func Pointer to public getter function in `Mineral'
   * @return properties Array of endmember properties.
   */
  template <typename Func>
  Eigen::ArrayXd map_endmembers_to_array(Func&& func) const {
    const auto& em_ref = solution_model->endmembers;
    Eigen::ArrayXd mapped_properties(solution_model->n_endmembers);
    std::transform(
      em_ref.begin(), em_ref.end(), mapped_properties.data(),
      [&func](const auto& em) {
        return (em.*func)();
    });
    return mapped_properties;
  }

  // Public getters for extra Solution functions
  /**
   * @brief Retrieves molar excess gibbs free energy the of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_excess_gibbs()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Excess gibbs free energy in [J/mol].
   */
  double get_excess_gibbs() const;

  /**
   * @brief Retrieves the excess molar volume of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_excess_volume()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Excess volume in [m^3/mol].
   */
  double get_excess_volume() const;

  /**
   * @brief Retrieves the excess molar entropy of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_excess_entropy()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Excess entropy in [J/K/mol].
   */
  double get_excess_entropy() const;

  /**
   * @brief Retrieves the excess molar enthalpy of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_excess_enthalpy()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Excess enthalpy in [J/mol].
   */
  double get_excess_enthalpy() const;

  /**
   * @brief Retrieves the endmember activities.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_activities()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Endmember activities [unitless].
   */
  Eigen::ArrayXd get_activities() const;

  /**
   * @brief Retrieves endmember activity coefficients of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_activity_coefficients()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return (gamma = activity / ideal activity) [unitless].
   */
  Eigen::ArrayXd get_activity_coefficients() const;

  /**
   * @brief Retrieves the excess partial molar gibbs free energies of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_excess_partial_gibbs()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Excess partial gibbs in [J/mol].
   */
  Eigen::ArrayXd get_excess_partial_gibbs() const;

  /**
   * @brief Retrieves the excess partial volumes of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_excess_partial_volumes()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Excess partial volumes in [m^3/mol].
   */
  Eigen::ArrayXd get_excess_partial_volumes() const;

  /**
   * @brief Retrieves the excess partial entropies of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_excess_partial_entropies()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Excess partial entropies in [J/K].
   */
  Eigen::ArrayXd get_excess_partial_entropies() const;

  /**
   * @brief Retrieves the endmember partial molar gibbs free energy.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_partial_gibbs()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Partial molar gibbs free energies in [J/mol].
   */
  Eigen::ArrayXd get_partial_gibbs() const;

  /**
   * @brief Retrieves the endmember partial volumes of the solid solution.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_partial_volumes()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Partial volumes in [m^3].
   */
  Eigen::ArrayXd get_partial_volumes() const;

  /**
   * @brief Retrieves the endmember partial entropies.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_partial_entropies()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Partial entropies in [J/K].
   */
  Eigen::ArrayXd get_partial_entropies() const;

  /**
   * @brief Retrieves the second compositional derivative of the Gibbs free energy.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_gibbs_hessian()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Second derivative of Gibbs in [J].
   */
  Eigen::MatrixXd get_gibbs_hessian() const;

  /**
   * @brief Retrieves the second compositional derivative of the entropy.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_entropy_hessian()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Second derivative of entropy [J/K].
   */
  Eigen::MatrixXd get_entropy_hessian() const;

  /**
   * @brief Retrieves the second compositional derivative of the volume.
   *
   * Uses a cached value if available, or calls
   * `Solution::compute_volume_hessian()` and caches the result.
   *
   * @note Use `Solution::reset()` to clear cached values.
   *
   * @return Second derivatives of volume in [m^3].
   */
  Eigen::MatrixXd get_volume_hessian() const;

  // Public getter functions for stored solution properties
  int get_n_reactions() const { return n_reactions; }
  const Eigen::MatrixXd& get_stoichiometric_matrix() const { return stoichiometric_matrix; }
  const Eigen::MatrixXd& get_compositional_basis() const { return compositional_basis; }
  const Eigen::MatrixXd& get_compositional_null_basis() const { return compositional_null_basis; }
  const Eigen::MatrixXd& get_reaction_basis() const { return reaction_basis; }
  const std::vector<int>& get_independent_element_indices() const { return independent_element_indices }
  const std::vector<int>& get_dependent_element_indices() const { return dependent_element_indices; }
  const std::vector<std::string>& get_endmember_names() const { return endmember_names; }
  const std::vector<std::string>& get_elements() const { return elements; }
  const std::vector<FormulaMap>& get_endmember_formulae() const { return endmember_formulae; }

 protected:

  // Overrides of defaults from Material
  double compute_molar_internal_energy() const override;
  double compute_molar_gibbs() const override;
  double compute_molar_helmholtz() const override;
  double compute_molar_mass() const override;
  double compute_molar_volume() const override;
  double compute_density() const override;
  double compute_molar_entropy() const override;
  double compute_molar_enthalpy() const override;
  double compute_isothermal_bulk_modulus_reuss() const override;
  double compute_isentropic_bulk_modulus_reuss() const override;
  double compute_isothermal_compressibility_reuss() const override;
  double compute_isentropic_compressibility_reuss() const override;
  double compute_shear_modulus() const override;
  double compute_p_wave_velocity() const override;
  double compute_bulk_sound_velocity() const override;
  double compute_shear_wave_velocity() const override;
  double compute_grueneisen_parameter() const override;
  double compute_thermal_expansivity() const override;
  double compute_molar_heat_capacity_v() const override;
  double compute_molar_heat_capacity_p() const override;

  // New functions
  /**
   * @brief Computes molar excess gibbs free energy the of the solid solution.
   *
   * @return Excess gibbs free energy in [J/mol].
   */
  double compute_excess_gibbs() const;

  /**
   * @brief Computes excess molar volume of the solid solution.
   *
   * @return Excess volume in [m^3/mol].
   */
  double compute_excess_volume() const;

  /**
   * @brief Computes excess molar entropy of the solid solution.
   *
   * @return Excess entropy in [J/K/mol].
   */
  double compute_excess_entropy() const;

  /**
   * @brief Computes excess molar enthalpy of the solid solution.
   *
   * @return Excess enthalpy in [J/mol].
   */
  double compute_excess_enthalpy() const;

  /**
   * @brief Computes endmember activities.
   *
   * @return Endmember activities [unitless].
   */
  Eigen::ArrayXd compute_activities() const;

  /**
   * @brief Computes endmember activity coefficients of the solid solution.
   *
   * @return Endmember activity coefficients [unitless].
   */
  Eigen::ArrayXd compute_activity_coefficients() const;

  /**
   * @brief Computes excess partial molar gibbs free energies.
   *
   * @return Excess partial gibbs in [J/mol].
   */
  Eigen::ArrayXd compute_excess_partial_gibbs() const;

  /**
   * @brief Computes excess partial volumes of the solid solution.
   *
   * @return Excess partial volumes in [m^3/mol].
   */
  Eigen::ArrayXd compute_excess_partial_volumes() const;

  /**
   * @brief Computes excess partial entropies of the solid solution.
   *
   * @return Excess partial entropies in [J/K].
   */
  Eigen::ArrayXd compute_excess_partial_entropies() const;

  /**
   * @brief Computes endmember partial molar gibbs free energies.
   *
   * @return Partial molar gibbs free energies in [J/mol].
   */
  Eigen::ArrayXd compute_partial_gibbs() const;

  /**
   * @brief Computes endmember partial volumes of the solid solution.
   *
   * @return Partial volumes in [m^3].
   */
  Eigen::ArrayXd compute_partial_volumes() const;

  /**
   * @brief Computes endmember partial entropies.
   *
   * @return Partial entropies in [J/K].
   */
  Eigen::ArrayXd compute_partial_entropies() const;

  /**
   * @brief Hessian of gibbs free energy.
   *
   * @return Second compositional derivative [J].
   */
  Eigen::MatrixXd compute_gibbs_hessian() const;

  /**
   * @brief Hessian of entropy.
   *
   * @return Second compositional derivative [J/K].
   */
  Eigen::MatrixXd compute_entropy_hessian() const;

  /**
   * @brief Hessian of volumes.
   *
   * @return Second compositional derivative [m^3].
   */
  Eigen::MatrixXd compute_volume_hessian() const;

 private:
  // Cached properties
  mutable std::optional<double> excess_gibbs;
  mutable std::optional<double> excess_volume;
  mutable std::optional<double> excess_entropy;
  mutable std::optional<double> excess_enthalpy;
  // 1D Eigen arrays (could be made Vectors)
  mutable std::optional<Eigen::ArrayXd> activities;
  mutable std::optional<Eigen::ArrayXd> activity_coefficients;
  mutable std::optional<Eigen::ArrayXd> excess_partial_gibbs;
  mutable std::optional<Eigen::ArrayXd> excess_partial_volumes;
  mutable std::optional<Eigen::ArrayXd> excess_partial_entropies;
  mutable std::optional<Eigen::ArrayXd> partial_gibbs;
  mutable std::optional<Eigen::ArrayXd> partial_volumes;
  mutable std::optional<Eigen::ArrayXd> partial_entropies;
  // 2D Matrices
  mutable std::optional<Eigen::MatrixXd> gibbs_hessian;
  mutable std::optional<Eigen::MatrixXd> entropy_hessian;
  mutable std::optional<Eigen::MatrixXd> volume_hessian;
  // site_occupancies - Map?

  // Solution properties set on initialisation
  Eigen::MatrixXd stoichiometric_matrix;
  Eigen::MatrixXd compositional_basis;
  Eigen::MatrixXd compositional_null_basis;
  Eigen::MatrixXd reaction_basis;
  int n_reactions;
  std::vector<int> independent_element_indices;
  std::vector<int> dependent_element_indices;
  std::vector<std::string> endmember_names;
  std::vector<std::string> elements;
  std::vector<FormulaMap> endmember_formulae;

  // Molar fractions - not cached so not deleted by reset()
  // If Solution is subclassed will need a public getter function
  Eigen::ArrayXd molar_fractions;

  // Shared pointer to solution model class
  std::shared_ptr<SolutionModel> solution_model;

  // Compute functions for intialised properties
  /**
   * @brief Stores human readable endmember names
   */
  void setup_endmember_names();

  /**
   * @brief Stores a vector of endmember formulae.
   *
   * Each formula is a FormulaMap, where
   * FormulaMap = std::unordered_map<std::string, int>;
   */
  void setup_endmember_formulae();

  /**
   * @brief Stores vector of elementes sorted in IUPAC order.
   */
  void setup_elements();

  /**
   * @brief Initialises matrix detailing stoichiometry.
   *
   * Each element M[i,j] corresponds to to the number of
   * atoms of element j in endmember i.
   */
  void setup_stoichiometric_matrix();

  /**
   * @brief Stores the independent set of element indices.
   *
   * If the amounts of these elements are known the amounts of other
   * elements can be inferred by:
   *   -compositional_null_basis[independent_element_indices].dot(element_amounts)
   */
  void setup_independent_element_indices();

  /**
   * @brief Stores the element indices not in the independent list.
   */
  void setup_dependent_element_indices();

  /**
   * @brief Stores the reaction basis matrix.
   *
   * Each element M[i,j] corresponds to the number of moles
   * of endmember j involved in reaction i.
   */
  void setup_reaction_basis();

  /**
   * @bried Stores the number of possible reactions from reaction basis.
   */
  void setup_n_reactions();

  /**
   * @brief Stores the compositional basis matrix.
   */
  void setup_compositional_basis();

  /**
   * @brief Stores the compositional null basis.
   *
   * Stores the matrix where N[b] = 0 for all bulk compositions that
   * can be produced with a linear sum of the endmembers in the solution.
   */
  void setup_compositional_null_basis();

  /**
   * @brief Helper function to setup stored solution properties.
   */
  void setup_solution_properties();

};

#endif // BURNMAN_CORE_SOLUTION_HPP_INCLUDED
