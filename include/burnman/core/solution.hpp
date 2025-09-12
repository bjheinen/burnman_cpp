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

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/core/mineral.hpp"
#include "burnman/core/solution_model.hpp"
#include "burnman/core/composite_material.hpp"

namespace burnman {

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
 * Uses an instance of `solution_models::SolutionModel' (or derived class) to calculate
 * interaction terms between endmembers.
 *
 * All solution parameters and properties are in SI units.
 * e.g. Interaction parameters in J/mol, with T & P derivatives in J/K/mol
 * and m^3/mol.
 *
 * For composition and stoichiometric dependent properties
 * see the abstract class `CompositeMaterial'.
 *
 */
class Solution : public CompositeMaterial {

 public:

  virtual ~Solution() = default;

  // Override of reset to include additional solution properties
  void reset() override;

  // Utility functions
  // TODO: std::invoke to deal with lambda functions?
  /**
   * @brief Utility function to map endmember properties to Eigen::Array.
   *
   * Usage:
   *  map_endmembers_to_array(&Mineral::get_property);
   *    where, get_property() is the function needed
   *  alternatively, with a lambda function, e.g.
   * map_endmembers_to_array([](const Mineral& m) {
   *  return m.get_property() + 10.0;
   * });
   *
   * @note Properties expected as type double.
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

  /**
   * @brief Utility function to map endmember properties to std::vector.
   *
   * Usage: see `Solution::map_endmembers_to_array'.
   *
   * @param func Pointer to public getter function in `Mineral'
   * @return vector of endmember properties (of type T)
   */
  template <typename T, typename Func>
  std::vector<T> map_endmembers_to_vector(Func&& func) const {
    const auto& em_ref = solution_model->endmembers;
    std::vector<T> mapped_properties;
    mapped_properties.reserve(static_cast<std::size_t>(solution_model->n_endmembers));
    std::transform(
      em_ref.begin(), em_ref.end(), std::back_inserter(mapped_properties),
      [&func](const auto& em) {
        return (em.*func)();
      });
    return mapped_properties;
  }

  // Override public methods
  void set_state(double new_pressure, double new_temperature) override;
  void set_method(types::EOSType new_method) override;
  void set_method(std::shared_ptr<EquationOfState> new_method) override;

  // Public setters for Solution
  /**
   * @brief Sets the solution model for the material.
   *
   * @param model Shared pointer to the solution model.
   */
  void set_solution_model(std::shared_ptr<solution_models::SolutionModel> model);

  /**
   * @brief Sets the molar amounts of each endmember.
   *
   * @param composition_vector Molar fractions of endmembers.
   * @throws RuntimeError if solution model not set.
   * @throws RuntimeError if sum(composition_vector) != 1.
   * @throws RuntimeError if length(composition_vector) != n_endmembers.
   */
  void set_composition(const Eigen::ArrayXd& composition_vector);

  // Public getters for extra Solution functions
  /**
   * @brief Retrieves molar fractions array.
   */
  Eigen::ArrayXd get_molar_fractions() const;

  /**
   * @brief Returns fractional site occupancy matrix from solution model.
   */
  Eigen::ArrayXXd get_endmember_occupancies() const;

  /**
   * @brief Returns site occupancy matrix (absolute number of atoms) from solution model.
   */
  Eigen::ArrayXXd get_endmember_n_occupancies() const;

  /**
   * @brief Returns site names from solution model.
   */
  std::vector<std::string> get_site_names() const;

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

  // get_partial_gibbs() provided by CompositeMaterial

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
  double compute_isentropic_thermal_gradient() const override;
  types::FormulaMap compute_formula() const override;

  // Overrides of pure virtual functions from CompositeMaterial
  int compute_n_endmembers() const override;
  void setup_endmember_names() const override;
  void setup_endmember_formulae() const override;
  Eigen::ArrayXd compute_partial_gibbs() const override;

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

  // Molar fractions - not cached so not deleted by reset()
  Eigen::ArrayXd molar_fractions;

  // Shared pointer to solution model class
  std::shared_ptr<solution_models::SolutionModel> solution_model;

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
  // partial_gibbs provided by CompositeMaterial
  mutable std::optional<Eigen::ArrayXd> partial_volumes;
  mutable std::optional<Eigen::ArrayXd> partial_entropies;
  // 2D Matrices
  mutable std::optional<Eigen::MatrixXd> gibbs_hessian;
  mutable std::optional<Eigen::MatrixXd> entropy_hessian;
  mutable std::optional<Eigen::MatrixXd> volume_hessian;
  // site_occupancies - Map?

};

} // namespace burnman

#endif // BURNMAN_CORE_SOLUTION_HPP_INCLUDED
